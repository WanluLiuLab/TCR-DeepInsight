from dataclasses import dataclass
from functools import partial
import torch
from torch import nn
from typing import List, Optional, Tuple, Union
from transformers.activations import ACT2FN
from transformers.modeling_outputs import (
    BaseModelOutputWithPast, 
    CausalLMOutputWithPast
)
from transformers.modeling_utils import PreTrainedModel
from transformers.utils import logging
from transformers.utils.model_parallel_utils import assert_device_map, get_device_map
from transformers.modeling_outputs import BaseModelOutput

from scatlasvae.utils._utilities import exists, absent


from ._config import TCRGenConfig
from ._modeling_utils import ModelOutput
from ._primitives import Linear, RotaryEmbedding, ChunkLayer, permute_final_dims

import importlib 
import einops

fa_is_installed = importlib.util.find_spec("flash_attn") is not None
if fa_is_installed:
    try:
        from flash_attn import flash_attn_qkvpacked_func, flash_attn_func  
    except:
        fa_is_installed = False

logger = logging.get_logger(__name__)

class TCRGenAttentionBase(nn.Module):
    def __init__(self, config):
        super().__init__()
        max_positions = config.max_position_embeddings
        self.register_buffer(
            "bias",
            torch.tril(torch.ones(
                max_positions, max_positions, dtype=torch.bool
            )).view(
                1, 1, max_positions, max_positions
            )
        )
        self.register_buffer(
            "mask_bias",
            torch.tensor(config.neg_inf)
        )

        self.dropout_attention = nn.Dropout(config.dropout_rate_attention)
        self.dropout_hidden = nn.Dropout(config.dropout_rate_hidden)

        self.hidden_dim = config.hidden_dim
        self.num_attention_heads = config.num_attention_heads
        self.head_dim = self.hidden_dim // config.num_attention_heads
        if self.head_dim * self.num_attention_heads != self.hidden_dim:
            raise ValueError(
                f"hidden_dim ({self.hidden_dim}) must be divisible by num_attention_heads ({self.num_attention_heads})"
            )
        
        self.scale_attn_weights = torch.sqrt(
            torch.tensor(self.head_dim, dtype=torch.float32)
        )
        self.qkv_proj = Linear(
            self.hidden_dim, 
            self.hidden_dim * 3, 
            bias=False,
            init='normal'
        )
        self.out_proj = Linear(
            self.hidden_dim, 
            self.hidden_dim, 
            bias=False,
            init='normal'
        )
        self.rotary_dim = config.rotary_dim if exists(config.rotary_dim) else None

class TCRGenAttention(TCRGenAttentionBase):
    def __init__(self, config):
        super().__init__(config)
    
    def _attn(self, 
        q: torch.Tensor, 
        k: torch.Tensor, 
        v: torch.Tensor, 
        attention_mask: Optional[torch.Tensor], 
        head_mask: Optional[torch.Tensor], 
    ):
        '''
        Main attention function
        
        :param q: query tensor. Shape: [batch_size, query_length, hidden_dim, num_attention_heads]
        :param k: key tensor. Shape: [batch_size, key_length, hidden_dim, num_attention_heads]
        :param v: value tensor. Shape: [batch_size, value_length, hidden_dim, num_attention_heads]
        :param attention_mask: attention mask. 
        :param head_mask: head mask.
        '''
        q_l = q.size(-2)
        k_l = k.size(-2)
        # [1, 1, query_length, key_length]
        causal_mask = self.bias[:, :, k_l - q_l : k_l, :k_l]

        # fp32
        q = q.to(torch.float32)
        k = k.to(torch.float32)

        a = torch.matmul(q, k.transpose(-1, -2))

        a = a / self.scale_attn_weights
        a = torch.where(causal_mask, a, self.mask_bias.to(a.dtype))

        if exists(attention_mask):
            a = a + attention_mask

        a = nn.Softmax(dim=-1)(a)
        a = a.to(v.dtype)
        a = self.dropout_attention(a)

        if exists(head_mask):
            a = a * head_mask

        # return output and attention weights
        return torch.matmul(a, v), a

    def _flash_attn(self,
        q: torch.Tensor, 
        k: torch.Tensor, 
        v: torch.Tensor, 
        attention_mask: Optional[torch.Tensor], 
        head_mask: Optional[torch.Tensor], 
    ):  
        raise NotImplementedError("flash_attn is not yet implemented")
        o, softmax_lse, S_dmask = flash_attn_func(
            q,
            k,
            v,
            dropout_p=0.0, 
            softmax_scale=None, 
            causal=True,
            return_attn_probs=True
        )
        return o, softmax_lse

    def _split_heads(self, x, n_head, head_dim, mp_num):
        reshaped = x.reshape(x.shape[:-1] + (n_head//mp_num, head_dim))
        reshaped = reshaped.reshape(x.shape[:-2] + (-1, ) + reshaped.shape[-1:])
        return reshaped
    
    def _merge_heads(self, x, num_attention_heads, head_dim):
        if len(x.shape) == 4:
            x = x.permute(0, 2, 1, 3).contiguous()
        elif len(x.shape) == 5:
            x = x.permute(0, 1, 3, 2, 4).contiguous()
        else:
            raise ValueError(f"Unexpected input with shape {x.shape}")
        return x.view(x.shape[:-2] + (num_attention_heads * head_dim, ))

    def forward(
        self,
        hidden_states: torch.Tensor,
        layer_past: Optional[Tuple[torch.Tensor]] = None,
        attention_mask: Optional[torch.Tensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        use_cache: Optional[bool] = None,
        output_attentions: Optional[bool] = None
    ):
        MP_NUM = 8
        # [batch_size, seq_length, 3 * hidden_dim]
        qkv = self.qkv_proj(hidden_states)
        # [batch_size, seq_length, num_attention_heads, hidden_dim * 3 / num_attention_heads]
        qkv_split = qkv.reshape(qkv.shape[:-1] + (MP_NUM, -1))
        # [batch_size, num_attention_heads, seq_length, num_attention_heads, head_dim]
        q,k,v = torch.split(
            qkv_split, 
            self.head_dim * self.num_attention_heads // MP_NUM,
            dim=-1
        )
        # [batch_size, num_attention_heads, seq_length, head_dim]
        q = self._split_heads(q, self.num_attention_heads, self.head_dim, MP_NUM)
        k = self._split_heads(k, self.num_attention_heads, self.head_dim, MP_NUM)
        v = self._split_heads(v, self.num_attention_heads, self.head_dim, MP_NUM)
        # [batch_size, seq_length, num_attention_heads, head_dim]

        v = v.permute(0, 2, 1, 3)

        seq_len = k.shape[1]
        offset = 0

        if exists(layer_past):
            offset = layer_past[0].shape[-2]
            seq_len += offset

        if exists(self.rotary_dim):
            k_rot = k[:, :, :, :self.rotary_dim ]
            k_pass = k[:, :, :, self.rotary_dim:]
            q_rot = q[:, :, :, :self.rotary_dim ]
            q_pass = q[:, :, :, self.rotary_dim:]
            sincos = RotaryEmbedding.fixed_position_embedding(
                k_rot, seq_dim=1, seq_len=seq_len
            )
            k_rot = RotaryEmbedding.apply_rotary_pos_emb(
                k_rot, sincos, offset=offset
            )
            q_rot = RotaryEmbedding.apply_rotary_pos_emb(
                q_rot, sincos, offset=offset
            )
            k = torch.cat([k_rot, k_pass], dim=-1)
            q = torch.cat([q_rot, q_pass], dim=-1)
        else: 
            sincos = RotaryEmbedding.fixed_position_embedding(
                k, seq_dim=1, seq_len=seq_len
            )
            k = RotaryEmbedding.apply_rotary_pos_emb(
                k, sincos, offset=offset
            )
            q = RotaryEmbedding.apply_rotary_pos_emb(
                q, sincos, offset=offset
            )

        k = k.permute(0, 2, 1, 3)
        q = q.permute(0, 2, 1, 3)

        if exists(layer_past):
            past_k, past_v = layer_past[0], layer_past[1]
            k = torch.cat([past_k, k], dim=-2)
            v = torch.cat([past_v, v], dim=-2)
        
        if use_cache:
            present = (k, v)
        else:
            present = None

        attn_out, attn_weights = self._attn(
            q, k, v, attention_mask, head_mask
        )

        attn_out = self._merge_heads(
            attn_out, self.num_attention_heads, self.head_dim
        )

        attn_out = self.out_proj(attn_out)
        attn_out = self.dropout_hidden(attn_out)

        outputs = ModelOutput({
            "attention_output": attn_out,
            "present": present,
        })

        if output_attentions:
            outputs["attention_weights"] = attn_weights

        return outputs
    
class TCRGenRowAttention(TCRGenAttentionBase):
    def __init__(self, config):
        super().__init__()
        
        self.register_buffer(
            "bias",
            torch.tril(torch.ones(
                self.max_positions, self.max_positions, dtype=torch.bool
            )).view(
                1, 1, self.max_positions, self.max_positions
            )
        )
        self.register_buffer(
            "mask_bias",
            torch.tensor(config.neg_inf)
        )

class TCRGenColumnAttention(TCRGenAttentionBase):
    def __init__(self, config):
        super().__init__()

class TCRGenTriangularAttention(TCRGenAttentionBase):
    def __init__(self, config, starting: bool = True):
        super().__init__()

        self.starting = starting
        self.ln = nn.LayerNorm(config.hidden_dim, eps=config.layer_norm_epsilon)
        self.linear = Linear(
            config.hidden_dim, 
            config.num_attention_heads, 
            bias=False, 
            init='normal'
        )
        self.mha = TCRGenAttention(config)

    @torch.jit.ignore
    def _chunk(self,
        hidden_states: torch.Tensor,
        biases: List[torch.Tensor],
        chunk_size: int,
    ):
        mha_inputs = {
            "hidden_states": hidden_states,
            "biases": biases,
        }

        return ChunkLayer.chunk_layer(
            partial(self.mha),
            mha_inputs,
            chunk_size=chunk_size,
            no_batch_dims=len(hidden_states.shape[:-2]),
            _out=hidden_states,
        )
        

    def forward(self,
        hidden_states: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        chunk_size: Optional[int] = None,
    ):
        """
        :param hidden_states: [batch_size, seq_length, seq_length, hidden_dim]
        """
        if absent(attention_mask):
            attention_mask = hidden_states.new_ones(
                hidden_states.shape[:-1] 
            )
        if not self.starting:
            hidden_states = hidden_states.transpose(-2, -3)
            attention_mask = attention_mask.transpose(-1, -2)

        hidden_states = self.ln(hidden_states)

        # [batch_size, seq_length, 1, 1, seq_length]
        mask_bias = (self.config.neg_inf * (1 - attention_mask))[..., :, None, None, :]
        # [batch_size, num_heads, seq_length, seq_length]
        triangular_bias = permute_final_dims(self.linear(hidden_states), (2,0,1))
        # [batch_size, 1, num_heads, seq_length, seq_length]
        triangular_bias = triangular_bias.unsquueze(-4)

        biases = [mask_bias, triangular_bias]

class TCRGenTransition(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.hidden_dim = config.hidden_dim
        self.intermediate_size = config.intermediate_size
        self.fc_in = Linear(
            self.hidden_dim, 
            self.intermediate_size,
            bias=True,
            init='normal'
        )
        self.fc_out = Linear(
            self.intermediate_size,
            self.hidden_dim,
            bias=True,
            init='normal'
        )
        self.act = ACT2FN[config.activation_function]
        self.dropout_hidden = nn.Dropout(config.dropout_rate_hidden)

    def forward(self, hidden_states):
        hidden_states = self.fc_in(hidden_states)
        hidden_states = self.act(hidden_states)
        hidden_states = self.fc_out(hidden_states)
        hidden_states = self.dropout_hidden(hidden_states)
        return hidden_states
    
class TCRGenBlock(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.attn = TCRGenAttention(config)
        self.ln = nn.LayerNorm(config.hidden_dim, eps=config.layer_norm_epsilon)
        self.transition = TCRGenTransition(config)

    def forward(
        self,
        hidden_states: torch.Tensor,
        layer_past: Optional[Tuple[torch.Tensor]] = None,
        attention_mask: Optional[torch.Tensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        use_cache: Optional[bool] = None,
        output_attentions: Optional[bool] = None
    ):
        residual = hidden_states
        hidden_states = self.ln(hidden_states)
        attn_outputs = self.attn(
            hidden_states, 
            layer_past=layer_past, 
            attention_mask=attention_mask, 
            head_mask=head_mask, 
            use_cache=use_cache, 
            output_attentions=output_attentions
        )
        attn_output = attn_outputs["attention_output"]
        ffd_output = self.transition(attn_output)
        hidden_states = attn_output + residual + ffd_output
        outputs = ModelOutput({
            "hidden_states": hidden_states,
            "present": attn_outputs["present"] if use_cache else None,
            "attention_output": attn_outputs["attention_output"] if output_attentions else None,
        })
        return outputs 


class TCRGenPreTrainedModel(PreTrainedModel):
    config_class = TCRGenConfig
    base_model_prefix = "transformer"
    is_parallelizable = False

    def __init__(self, *inputs, **kwargs):
        super().__init__(*inputs, **kwargs)

    def _init_weights(self, module):
        if isinstance(module, nn.Embedding):
            module.weight.data.normal_(mean=0.0, std=self.config.initializer_range)
            if module.padding_idx is not None:
                module.weight.data[module.padding_idx].zero_()
        elif isinstance(module, nn.LayerNorm):
            module.bias.data.zero_()
            module.weight.data.fill_(1.0)
        
class TCRGenModel(TCRGenPreTrainedModel):
    def __init__(self, config):
        super().__init__(config)

        self.hidden_dim = config.hidden_dim
        self.vocab_size = config.vocab_size
        self.word_to_embedding = nn.Embedding(self.vocab_size, self.hidden_dim)
        self.dropout_embedding = nn.Dropout(config.dropout_rate_embedding)

        self.blocks = nn.ModuleList([
            TCRGenBlock(config) for _ in range(config.num_blocks)
        ])

        self.ln_f = nn.LayerNorm(config.hidden_dim, eps=config.layer_norm_epsilon)

        self.rotary_dim = min(
            config.rotary_dim, 
            config.context_size // config.num_attention_heads
        )
        self.init_weights()
        self.model_parallel = False

    def parallelize(self, device_map=None):
        
        self.device_map = get_device_map(len(self.blocks), range(torch.cuda.device_count())) if absent(device_map) else device_map

        assert_device_map(device_map, len(self.blocks))
        self.blocks = nn.parallel.replicate(
            self.blocks, 
            device_ids=device_map
        )
        self.model_parallel = True
        self.first_device = 'cpu' if 'cpu' in self.device_map.keys() else 'cuda:' + str(min(self.device_map.keys()))
        self.last_device = 'cuda:' + str(max(self.device_map.keys()))
        self.wte = self.wte.to(self.first_device)
        # load blocks to device
        for k,v in self.device_map.items():
            for block in v:
                cuda_device = 'cuda:' + str(k)
                self.blocks[block] = self.blocks[block].to(cuda_device)

        self.ln_f = self.ln_f.to(self.last_device)

    def deparellelize(self):
        self.blocks = nn.parallel.scatter(
            self.blocks, 
            device_ids=list(self.device_map.keys())
        )
        self.device_map = None 
        self.model_parallel = False
        self.first_device = 'cpu'
        self.last_device = 'cpu'
        self.wte = self.wte.to('cpu')
        self.ln_f = self.ln_f.to('cpu')
        for block in self.blocks:
            block = block.to('cpu')
        torch.cuda.empty_cache()

    @property
    def input_embeddings(self):
        return self.word_to_embedding

    @input_embeddings.setter
    def input_embeddings(self, value):
        self.word_to_embedding = value

    def forward(
        self,
        input_ids: Optional[torch.Tensor],
        past_key_values: Optional[Tuple[torch.Tensor]] = None,
        attention_mask: Optional[torch.Tensor] = None,
        token_type_ids: Optional[torch.Tensor] = None,
        position_ids: Optional[torch.Tensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        inputs_embeddings: Optional[torch.Tensor] = None,
        use_cache: Optional[bool] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ):
        output_attentions = output_attentions if exists(output_attentions) else self.config.output_attentions
        output_hidden_states = output_hidden_states if exists(output_hidden_states) else self.config.output_hidden_states
        use_cache = use_cache if exists(use_cache) else self.config.use_cache

        if exists(input_ids) and exists(inputs_embeddings):
            raise ValueError("You cannot specify both input_ids and inputs_embeddings at the same time")
        elif exists(input_ids):
            input_shape = input_ids.shape
            input_ids = input_ids.view(-1, input_shape[-1])
            batch_size = input_ids.shape[0]
        elif exists(inputs_embeddings):
            input_shape = inputs_embeddings.shape[:-1]
            batch_size = inputs_embeddings.shape[0]
        else: 
            raise ValueError("You have to specify either input_ids or inputs_embeddings")
        
        device = input_ids.device if exists(input_ids) else inputs_embeddings.device

        if exists(token_type_ids):
            token_type_ids = token_type_ids.view(-1, input_shape[-1])

        if exists(position_ids):
            position_ids = position_ids.view(-1, input_shape[-1])

        if absent(past_key_values):
            past_length = 0
            past_key_values = tuple([None] * len(self.blocks))
        else: 
            past_length = past_key_values[0][0].shape[-2]

        if absent(position_ids):
            position_ids = torch.arange(
                past_length,
                input_shape[-1] + past_length,
                dtype=torch.long,
                device=device
            )
            position_ids = position_ids.unsqueeze(0).view(-1, input_shape[-1])
        
        if exists(attention_mask):
            if batch_size == 0:
                raise ValueError("batch_size == 0")
            attention_mask = attention_mask.view(batch_size, -1)
            # broadcast to [batch_size, num_heads, from_seq_length, to_seq_length]
            attention_mask = attention_mask[:, None, None, :]
            attention_mask = attention_mask.to(dtype=self.dtype)
            attention_mask = (1.0 - attention_mask) * -10000.0

        head_mask = self.get_head_mask(head_mask, self.config.num_blocks)

        if absent(inputs_embeddings):
            inputs_embeddings = self.word_to_embedding(input_ids)

        hidden_states = inputs_embeddings

        if exists(token_type_ids):
            token_type_embeddings = self.token_type_embeddings(token_type_ids)
            hidden_states = hidden_states + token_type_embeddings

        hidden_states = self.dropout_embedding(hidden_states)

        output_shape = input_shape + (hidden_states.shape[-1], )

        presents = () if use_cache else None 
        all_self_attentions = () if output_attentions else None
        all_hidden_states = () if output_hidden_states else None

        for i, (block, layer_past) in enumerate(zip(
            self.blocks,
            past_key_values
        )):
            # Model parallel
            if self.model_parallel:
                torch.cuda.set_device(hidden_states.device)
                if exists(layer_past):
                    layer_past = tuple([p.to(hidden_states.device) for p in layer_past])
                if exists(attention_mask):
                    attention_mask = attention_mask.to(hidden_states.device)
                if isinstance(head_mask, torch.Tensor):
                    head_mask = head_mask.to(hidden_states.device)

            if output_hidden_states:
                all_hidden_states = all_hidden_states + (hidden_states, )
            
            if getattr(self.config, "gradient_checkpointing", False) and self.training:
                if use_cache:
                    logger.warning(
                        "`use_cache=True` is incompatible with `config.gradient_checkpointing=True`. Setting "
                        "`use_cache=False`..."
                    )
                    use_cache = False
                def create_custom_forward(module):
                    def custom_forward(*inputs):
                        return module(
                            *inputs,
                            use_cache,
                            output_attentions
                        )
                    return custom_forward
                outputs = torch.utils.checkpoint.checkpoint(
                    create_custom_forward(block),
                    hidden_states,
                    None,
                    attention_mask,
                    head_mask[i],
                )
            else:
                outputs = block(
                    hidden_states,
                    layer_past=layer_past,
                    attention_mask=attention_mask,
                    head_mask=head_mask[i],
                    use_cache=use_cache,
                    output_attentions=output_attentions
                )

            hidden_states = outputs["hidden_states"]
            if use_cache:
                presents = presents + (outputs["present"], )

            if output_attentions:
                all_self_attentions = all_self_attentions + (outputs["attention_output"], )

            if self.model_parallel:
                pass 


        hidden_states = self.ln_f(hidden_states)
        hidden_states = hidden_states.view(*output_shape)
        if output_hidden_states:
            all_hidden_states = all_hidden_states + (hidden_states, )
        if not return_dict:
            return tuple(v for v in [hidden_states, presents, all_hidden_states, all_self_attentions] if exists(v))
        
        return BaseModelOutputWithPast(
            last_hidden_state=hidden_states,
            past_key_values=presents,
            hidden_states=all_hidden_states,
            attentions=all_self_attentions
        )
                
@dataclass
class ConstrastiveCausalLMOutputWithPast(BaseModelOutput):
    loss: Optional[torch.FloatTensor] = None
    contrastive_loss: Optional[torch.FloatTensor] = None
    logits: torch.FloatTensor = None
    past_key_values: Optional[Tuple[Tuple[torch.FloatTensor]]] = None
    hidden_states: Optional[Tuple[torch.FloatTensor]] = None
    attentions: Optional[Tuple[torch.FloatTensor]] = None
    cross_attentions: Optional[Tuple[torch.FloatTensor]] = None

class TCRGenForCausalLM(TCRGenPreTrainedModel):
    def __init__(self, config):
        super().__init__(config)
        self.tcrgen = TCRGenModel(config)
        self.lm_head = Linear(
            config.hidden_dim, 
            config.vocab_size, 
            bias=False,
            init='final'
        )
    
    def parallize(self, device_map=None):
        self.device_map = (
            get_device_map(len(self.tcrgen.blocks), range(torch.cuda.device_count()))
            if absent(device_map) else device_map
        )
        assert_device_map(device_map, len(self.tcrgen.blocks))
        self.tcrgen.parallelize(self.device_map)
        self.lm_head = self.lm_head.to(self.tcrgen.first_device)
        self.model_parallel = True

    def deparallelize(self):
        self.tcrgen.deparallelize()
        self.lm_head = self.lm_head.to('cpu')
        self.tcrgen = self.tcrgen.to('cpu')
        self.model_parallel = False
        self.device_map = None
        torch.cuda.empty_cache()

    def prepare_inputs_for_generation(self, input_ids, past_key_values=None, **kwargs):
        token_type_ids = kwargs.get("token_type_ids", None)
        # only last token for inputs_ids if past is defined in kwargs
        if exists(past_key_values):
            input_ids = input_ids[:, -1].unsqueeze(-1)
            if token_type_ids is not None:
                token_type_ids = token_type_ids[:, -1].unsqueeze(-1)

        attention_mask = kwargs.get("attention_mask", None)
        position_ids = kwargs.get("position_ids", None)

        if attention_mask is not None and position_ids is None:
            # create position_ids on the fly for batch generation
            position_ids = attention_mask.long().cumsum(-1) - 1
            position_ids.masked_fill_(attention_mask == 0, 1)
            if past_key_values:
                position_ids = position_ids[:, -1].unsqueeze(-1)
        else:
            position_ids = None
        return {
            "input_ids": input_ids,
            "past_key_values": past_key_values,
            "use_cache": kwargs.get("use_cache"),
            "position_ids": position_ids,
            "attention_mask": attention_mask,
            "token_type_ids": token_type_ids,
        }
    
    def forward(
        self,
        input_ids: Optional[torch.Tensor],
        past_key_values: Optional[Tuple[torch.Tensor]] = None,
        attention_mask: Optional[torch.Tensor] = None,
        token_type_ids: Optional[torch.Tensor] = None,
        position_ids: Optional[torch.Tensor] = None,
        head_mask: Optional[torch.Tensor] = None,
        inputs_embeddings: Optional[torch.Tensor] = None,
        labels: Optional[torch.Tensor] = None,
        label_weights: Optional[torch.Tensor] = None,
        contrastive_labels: Optional[torch.Tensor] = None,
        contrastive_mask: Optional[torch.Tensor] = None,
        binder: Optional[torch.Tensor] = None,
        use_cache: Optional[bool] = None,
        output_attentions: Optional[bool] = None,
        output_hidden_states: Optional[bool] = None,
        return_dict: Optional[bool] = None,
    ):
        return_dict = return_dict if return_dict is not None else self.config.use_return_dict

        outputs = self.tcrgen(
            input_ids,
            past_key_values=past_key_values,
            attention_mask=attention_mask,
            token_type_ids=token_type_ids,
            position_ids=position_ids,
            head_mask=head_mask,
            inputs_embeddings=inputs_embeddings,
            use_cache=use_cache,
            output_attentions=output_attentions,
            output_hidden_states=output_hidden_states,
            return_dict=return_dict
        )

        hidden_states = outputs[0]

        lm_logits = self.lm_head(hidden_states).to(torch.float32)

        loss = None 
        contrastive_loss = torch.tensor(0.0, device=lm_logits.device)

        if exists(labels):
            shift_logits = lm_logits[..., :-1, :].contiguous()
            shift_labels = labels[..., 1:].contiguous()
            loss_fct = nn.CrossEntropyLoss(reduction='none')
            loss = loss_fct(
                shift_logits.view(-1, shift_logits.size(-1)), 
                shift_labels.view(-1)
            )

            if exists(label_weights):
                loss = (loss * label_weights.repeat_interleave(
                    shift_logits.shape[-2]
                )).mean()

            if exists(binder):
                binder = binder.cpu().to(torch.float32).apply_(lambda x: x if x == 1 else -0.01).to(lm_logits.device)
                loss = (loss * binder.repeat_interleave(
                    shift_logits.shape[-2]
                )).mean() + 10
            loss = loss.to(hidden_states.dtype)

        if exists(contrastive_labels):
            shift_logits = lm_logits[..., :-1, :].contiguous()
            loss_fct = nn.CrossEntropyLoss(reduction='none')
            contrastive_mask = contrastive_mask[:,:shift_logits.shape[-2]].reshape(-1)
            for i in range(contrastive_labels.shape[1]):
                shift_contrastive_label = contrastive_labels[:,i,:shift_logits.shape[-2]+1][..., 1:].contiguous()
                contrastive_loss += (loss_fct(
                    shift_logits.view(-1, shift_logits.size(-1)), 
                    shift_contrastive_label[:,:shift_logits.shape[-2]].view(-1)
                ) * contrastive_mask).mean()
            contrastive_loss = contrastive_loss / shift_labels.shape[0]
            if exists(binder):
                loss = loss - contrastive_loss*0.5
            else:
                loss = loss + (10-contrastive_loss*0.5)

        if not return_dict:
            output = (lm_logits,) + outputs[1:]
            return ((loss,) + output) if exists(loss) else output
        
        return ConstrastiveCausalLMOutputWithPast(
            loss=loss,
            contrastive_loss=contrastive_loss,
            logits=lm_logits,
            past_key_values=outputs.past_key_values,
            hidden_states=outputs.hidden_states,
            attentions=outputs.attentions
        )
    
    @staticmethod 
    def _reorder_cache(
        past: Tuple[Tuple[torch.Tensor]],
        beam_idx: torch.Tensor
    ):
        return tuple(
            tuple(
                past_state.index_select(0, beam_idx.to(past_state.device))
                for past_state in layer_past
            )
            for layer_past in past
        )

    def generate(self, *args, **kwargs):
        return super().generate(*args, **kwargs)
        
    