from transformers.configuration_utils import PretrainedConfig
from transformers.utils import logging

logger = logging.get_logger(__name__)

class TCRGenConfig(PretrainedConfig):
    model_type = "tcrgen"

    def __init__(
        self,
        vocab_size=50400,
        context_size=512,
        max_position_embeddings=2048,
        hidden_dim=1024,
        rotary_dim=64,
        intermediate_size=4096,
        num_attention_heads=8,
        num_blocks=32,
        dropout_rate_attention=0.0,
        dropout_rate_embedding=0.0,
        dropout_rate_hidden=0.0,
        activation_function="gelu_new",
        neg_inf=-1e9,
        layer_norm_epsilon=1e-5,
        initializer_range=0.02,
        use_cache=True,
        use_flash_attention=False,
        gradient_checkpointing=False,
        bos_token_id=50256,
        eos_token_id=50256,
        **kwargs
    ):
        super().__init__(
            bos_token_id=bos_token_id, 
            eos_token_id=eos_token_id, 
            **kwargs
        )
        self.vocab_size = vocab_size
        self.context_size = context_size
        self.max_position_embeddings = max_position_embeddings
        self.hidden_dim = hidden_dim
        self.rotary_dim = rotary_dim
        self.intermediate_size = intermediate_size
        self.num_attention_heads = num_attention_heads
        self.num_blocks = num_blocks
        self.dropout_rate_attention = dropout_rate_attention
        self.dropout_rate_embedding = dropout_rate_embedding
        self.dropout_rate_hidden = dropout_rate_hidden
        self.activation_function = activation_function
        self.neg_inf = neg_inf 
        self.layer_norm_epsilon = layer_norm_epsilon
        self.initializer_range = initializer_range
        self.use_cache = use_cache
        self.use_flash_attention = use_flash_attention
        self.gradient_checkpointing = gradient_checkpointing

        self.bos_token_id = bos_token_id
        self.eos_token_id = eos_token_id

    @property
    def num_hidden_layers(self):
        return self.num_blocks