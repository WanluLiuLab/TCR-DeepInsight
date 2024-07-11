import torch
import torch.nn as nn
from typing import Optional, List, Tuple
import math
import einops
from functools import partialmethod, partial

from scatlasvae.model._primitives import Linear

from .._model_utils import permute_final_dims, softmax_no_cast, flatten_final_dims, add


class Attention(nn.Module):
    def __init__(
        self,
        c_q: int,
        c_k: int,
        c_v: int,
        c_hidden: int, 
        num_heads: int,
        gating: bool = True
    ):
        super(Attention, self).__init__()
        
        self.c_q = c_q 
        self.c_k = c_k
        self.c_v = c_v 
        self.c_hidden = c_hidden
        self.num_heads = num_heads
        self.gating = gating

        self.linear_q = Linear(
            self.c_q, self.c_hidden * self.num_heads, bias=False, init="glorot"
        )
        
        self.linear_k = Linear(
            self.c_k, self.c_hidden * self.num_heads, bias=False, init="glorot"
        )
        
        self.linear_v = Linear(
            self.c_v, self.c_hidden * self.num_heads, bias=False, init="glorot"
        )
        
        self.linear_o = Linear(
            self.c_hidden * self.num_heads, self.c_q, init="final"
        )
        
        self.linear_g = None
        if self.gating:
            self.linear_g = Linear(
                self.c_q, self.c_hidden * self.num_heads, init="gating"
            )

        self.sigmoid = nn.Sigmoid()


    @staticmethod
    def _attention(
        query: torch.Tensor,
        key: torch.Tensor, 
        value: torch.Tensor,
        biases: List[torch.Tensor],
        output_attentions: bool = False
    ) -> torch.Tensor:
        # [*, Head, C_hidden, K]
        # equivalent to permute_final_dims(key, (1, 0))
        key = permute_final_dims(key, (1, 0))

        # [*, H, Q, K]
        a = torch.matmul(query, key)

        for b in biases:
            a += b

        a = softmax_no_cast(a, -1)

        # [*, H, Q, C_hidden]
        o = torch.matmul(a, value)
        
        if output_attentions:
            return o, a
        return o

    def _prep_qkv(
        self, 
        q_x: torch.Tensor,
        kv_x: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        # [*, Q/K/V, H * C_hidden]
        q = self.linear_q(q_x)
        k = self.linear_k(kv_x)
        v = self.linear_v(kv_x)

        # [*, Q/K, H, C_hidden]
        q = q.view(q.shape[:-1] + (self.num_heads, -1))
        k = k.view(k.shape[:-1] + (self.num_heads, -1))
        v = v.view(v.shape[:-1] + (self.num_heads, -1))

        # [*, H, Q/K, C_hidden]
        q = q.transpose(-2, -3)
        k = k.transpose(-2, -3)
        v = v.transpose(-2, -3)

        q /= math.sqrt(self.c_hidden)

        return q, k, v
    
    
    def _wrap_up(self,
        o: torch.Tensor, 
        q_x: torch.Tensor
    ) -> torch.Tensor:
        if(self.linear_g is not None):
            g = self.sigmoid(self.linear_g(q_x))
        
            # [*, Q, H, C_hidden]
            g = g.view(g.shape[:-1] + (self.num_heads, -1))
            o = o * g

        # [*, Q, H * C_hidden]
        o = flatten_final_dims(o, 2)

        # [*, Q, C_q]
        o = self.linear_o(o)

        return o

    def forward(
        self,
        q_x: torch.Tensor,
        kv_x: torch.Tensor,
        biases: Optional[List[torch.Tensor]] = None,
        output_attentions: bool = False
    ) -> torch.Tensor:
        """
        Args:
            q_x:
                [*, Q, C_q] query data
            kv_x:
                [*, K, C_k] key data
            biases:
                List of biases that broadcast to [*, Head, Q, K]
        """
        if biases is None:
            biases = []
        q, k, v = self._prep_qkv(q_x, kv_x)
        o = Attention._attention(q, k, v, biases, output_attentions=output_attentions)
        if output_attentions:
            o, a = o
        o = o.transpose(-2, -3)
        o = self._wrap_up(o, q_x)
        if output_attentions:
            return o, a
        return o


class JointEmbedder(nn.Module):
    """
    Embeds single-representation from single representations
    """

    def __init__(
        self, 
        tf_dim: int,
        c_m: int,
        c_z: int,
        relpos_k: int = 128
    ):
        """
        Args:
            tf_dim:
                Final dimension of the target features
            c_j:
                Joint embedding dimension
        """
        super(JointEmbedder, self).__init__()

        self.tf_dim = tf_dim
        self.c_m = c_m
        self.c_z = c_z


        self.linear_tf_z_i = Linear(tf_dim, c_z)
        self.linear_tf_z_j = Linear(tf_dim, c_z)

        # RPE stuff
        self.relpos_k = relpos_k
        self.no_bins = 2 * relpos_k + 1
        self.linear_relpos = Linear(self.no_bins, c_z)

        self.layer_norm_z = nn.LayerNorm(c_z)

    def relpos(self, ri: torch.Tensor, rj: torch.Tensor):
        d = ri[..., None] - rj[..., None, :]
        boundaries = torch.arange(
            start=-self.relpos_k, end=self.relpos_k + 1, device=d.device
        ) 
        reshaped_bins = boundaries.view(((1,) * len(d.shape)) + (len(boundaries),))
        d = d[..., None] - reshaped_bins
        d = torch.abs(d)
        d = torch.argmin(d, dim=-1)
        d = nn.functional.one_hot(d, num_classes=len(boundaries)).float()
        d = d.to(ri.dtype)
        return self.linear_relpos(d)
    
    def forward(
        self, 
        tf_1: torch.Tensor,
        tf_2: torch.Tensor,
        mask_1: torch.Tensor,
        mask_2: torch.Tensor,
        inplace_safe: bool = False,
    ): 
        """
        Args:
            tf_1:
                "single" features of shape [*, N_res_1, tf_dim]
            tf_2:
                "single" features of shape [*, N_res_2, tf_dim]

            mask_1:
                attention mask of shape [*, N_res_1, 1]
            mask_2:
                attention mask of shape [*, N_res_2, 1]
        Returns:
            joint_emb:
                [*, N_res_1, N_res_1, C_j] joint embedding
        """

        # [*, N_res, c_j]
        tf_emb_i = self.linear_tf_z_i(tf_1)
        tf_emb_j = self.linear_tf_z_j(tf_2)

        # [*, N_res_1, N_res_2]
        joint_mask = torch.matmul(
            mask_1.unsqueeze(-1), 
            mask_2.unsqueeze(-2)
        )

        # [*, N_res_1, N_res_2, c_j]
        joint_emb = self.relpos(
            torch.arange(tf_emb_i.shape[1]).to(tf_emb_i.device).type(tf_emb_i.dtype) + \
            tf_emb_j.shape[0],
            torch.arange(tf_emb_j.shape[1]).to(tf_emb_j.device).type(tf_emb_j.dtype)
        )

        joint_emb = add(
            joint_emb,
            tf_emb_i[..., None, :],
            inplace=inplace_safe
        )

        joint_emb = add(
            joint_emb,
            tf_emb_j[..., None, :, :],
            inplace=inplace_safe
        )

        return joint_emb, joint_mask 

class OuterProductMean(nn.Module):
    """
    Implements Algorithm 10.
    """

    def __init__(self, c_m, c_z, c_hidden, eps=1e-3):
        """
        Args:
            c_m:
                MSA embedding channel dimension
            c_z:
                Pair embedding channel dimension
            c_hidden:
                Hidden channel dimension
        """
        super(OuterProductMean, self).__init__()

        self.c_m = c_m
        self.c_z = c_z
        self.c_hidden = c_hidden
        self.eps = eps

        self.layer_norm_1 = nn.LayerNorm(c_m)
        self.layer_norm_2 = nn.LayerNorm(c_m)

        self.linear_1 = Linear(c_m, c_hidden)
        self.linear_2 = Linear(c_m, c_hidden)
    
        self.linear_out = Linear(c_hidden**2, c_z, init="final")

    def _opm(self, a, b):
        # [*, N_res, N_res, C, C]
        outer = torch.einsum("...bac,...dae->...bdce", a, b)

        # [*, N_res, N_res, C * C]
        outer = outer.reshape(outer.shape[:-2] + (-1,))

        # [*, N_res, N_res, C_z]
        outer = self.linear_out(outer)

        return outer

    def _forward(
        self,
        m_1: torch.Tensor,
        m_2: torch.Tensor,
        mask_1: Optional[torch.Tensor] = None,
        mask_2: Optional[torch.Tensor] = None,
        inplace_safe: bool = False,
    ) -> torch.Tensor:
        """
        Args:
            m:
                [*, N_seq, N_res, C_m] MSA embedding
            mask:
                [*, N_seq, N_res] MSA mask
        Returns:
            [*, N_res, N_res, C_z] pair embedding update
        """
        if mask_1 is None:
            mask_1 = m_1.new_ones(m_1.shape[:-1])
        if mask_2 is None:
            mask_2 = m_2.new_ones(m_2.shape[:-1])

        # [*, N_seq, N_res, C_m]
        ln_1 = self.layer_norm_1(m_1)
        ln_2 = self.layer_norm_2(m_2)

        # [*, N_seq, N_res, C]
        mask_1 = mask_1.unsqueeze(-1)
        mask_2 = mask_2.unsqueeze(-1)

        a = self.linear_1(ln_1)
        a = a * mask_1

        b = self.linear_2(ln_2)
        b = b * mask_2

        del ln_1
        del ln_2

        a = a.transpose(-2, -3)
        b = b.transpose(-2, -3)

        outer = self._opm(a, b)

        # [*, N_res, N_res, 1]
        norm = torch.einsum("...abc,...adc->...bdc", mask_1, mask_2)
        norm = norm + self.eps

        # [*, N_res, N_res, C_z]
        if inplace_safe:
            outer /= norm
        else:
            outer = outer / norm

        return outer

    def forward(
        self,
        m_1: torch.Tensor,
        m_2: torch.Tensor,
        mask_1: Optional[torch.Tensor] = None,
        mask_2: Optional[torch.Tensor] = None,
        inplace_safe: bool = False,
    ) -> torch.Tensor:
        return self._forward(m_1, m_2, mask_1, mask_2, inplace_safe)

class TriangleAttention(nn.Module):
    def __init__(
        self, c_in, c_hidden, no_heads, starting=True, inf=1e9
    ):
        """
        Args:
            c_in:
                Input channel dimension
            c_hidden:
                Overall hidden channel dimension (not per-head)
            no_heads:
                Number of attention heads
        """
        super(TriangleAttention, self).__init__()

        self.c_in = c_in
        self.c_hidden = c_hidden
        self.no_heads = no_heads
        self.starting = starting
        self.inf = inf

        self.layer_norm = nn.LayerNorm(self.c_in)

        self.linear = Linear(c_in, self.no_heads, bias=False, init="normal")

        self.mha = Attention(
            self.c_in, self.c_in, self.c_in, self.c_hidden, self.no_heads
        )

    def forward(self, 
        x: torch.Tensor, 
        mask: Optional[torch.Tensor] = None,
        output_attentions: bool = False,
    ) -> torch.Tensor:
        """
        Args:
            x:
                [*, I, J, C_in] input tensor (e.g. the pair representation)
        Returns:
            [*, I, J, C_in] output tensor
        """ 
        if mask is None:
            # [*, I, J]
            mask = x.new_ones(
                x.shape[:-1],
            )

        if(not self.starting):
            x = x.transpose(-2, -3)
            mask = mask.transpose(-1, -2)

        # [*, I, J, C_in]
        x = self.layer_norm(x) 

        # [*, I, 1, 1, J]
        mask_bias = (self.inf * (mask - 1))[..., :, None, None, :]

        biases = [mask_bias]
        
        x = self.mha(
            q_x=x, 
            kv_x=x, 
            biases=biases, 
            output_attentions=output_attentions
        )
        if output_attentions:
            x, a = x

        if(not self.starting):
            x = x.transpose(-2, -3)

        if output_attentions:
            return x, a
        return x


TriangleAttentionStartingNode = TriangleAttention


class TriangleAttentionEndingNode(TriangleAttention):
    """
    Implements Algorithm 14.
    """
    __init__ = partialmethod(TriangleAttention.__init__, starting=False)



class TriangleMultiplicativeUpdate(nn.Module):
    """
    Implements Algorithms 11 and 12.
    """
    def __init__(self, c_z, c_hidden, _outgoing=True):
        """
        Args:
            c_z:
                Input channel dimension
            c:
                Hidden channel dimension
        """
        super(TriangleMultiplicativeUpdate, self).__init__()
        self.c_z = c_z
        self.c_hidden = c_hidden
        self._outgoing = _outgoing

        self.linear_a_p = Linear(self.c_z, self.c_hidden)
        self.linear_a_g = Linear(self.c_z, self.c_hidden, init="gating")
        self.linear_b_p = Linear(self.c_z, self.c_hidden)
        self.linear_b_g = Linear(self.c_z, self.c_hidden, init="gating")
        self.linear_g = Linear(self.c_z, self.c_z, init="gating")
        self.linear_z = Linear(self.c_hidden, self.c_z, init="final")

        self.layer_norm_in = nn.LayerNorm(self.c_z)
        self.layer_norm_out = nn.LayerNorm(self.c_hidden)

        self.sigmoid = nn.Sigmoid()

    def _combine_projections(self,
        a: torch.Tensor,
        b: torch.Tensor,
    ) -> torch.Tensor:
        if(self._outgoing):
            a = permute_final_dims(a, (2, 0, 1))
            b = permute_final_dims(b, (2, 1, 0))
        else:
            a = permute_final_dims(a, (2, 1, 0))
            b = permute_final_dims(b,  (2, 0, 1))

        p = torch.matmul(a, b)

        return permute_final_dims(p, (1, 2, 0))

    def _inference_forward(self,
        z: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        inplace_chunk_size: Optional[int] = None,
        with_add: bool = True,
    ):
        if mask is None:
            mask = z.new_ones(z.shape[:-1])

        mask = mask.unsqueeze(-1)
       
        def compute_projection_helper(pair, mask, a=True):
            if(a):
                linear_g = self.linear_a_g
                linear_p = self.linear_a_p
            else:
                linear_g = self.linear_b_g
                linear_p = self.linear_b_p
            
            pair = self.layer_norm_in(pair)
            p = linear_g(pair)
            p.sigmoid_()
            p *= linear_p(pair)
            p *= mask
            p = permute_final_dims(p, (2, 0, 1))
            return p

        def compute_projection(pair, mask, a=True, chunked=True): 
            need_transpose = self._outgoing ^ a
            if(not chunked):
                p = compute_projection_helper(pair, mask, a)
                if(need_transpose):
                    p = p.transpose(-1, -2)
            else:
                # This computation is chunked so as not to exceed our 2.5x 
                # budget with a large intermediate tensor
                linear_g = self.linear_a_g if a else self.linear_b_g
                c = linear_g.bias.shape[-1]
                out_shape = pair.shape[:-3] + (c,) + pair.shape[-3:-1]
                p = pair.new_zeros(out_shape)
                for i in range(0, pair.shape[-3], inplace_chunk_size):
                    pair_chunk = pair[..., i: i + inplace_chunk_size, :, :]
                    mask_chunk = mask[..., i: i + inplace_chunk_size, :, :]
                    pair_chunk = compute_projection_helper(
                        pair[..., i: i + inplace_chunk_size, :, :],
                        mask[..., i: i + inplace_chunk_size, :, :], 
                        a,
                    )
                    if(need_transpose):
                        pair_chunk = pair_chunk.transpose(-1, -2)
                        p[..., i: i + inplace_chunk_size] = pair_chunk
                    else:
                        p[..., i: i + inplace_chunk_size, :] = pair_chunk
                    
                    del pair_chunk

            return p

        # We start by fully manifesting a. In addition to the input, this
        # brings total memory consumption to 2x z (disregarding size of chunks)
        # [*, N, N, c]
        a = compute_projection(z, mask, True, chunked=True)

        b = compute_projection(z, mask, False, False)
        x = torch.matmul(a, b)
        x = self.layer_norm_out(x)
        x = self.linear_z(x)
        g = self.linear_g(z)
        g.sigmoid_()
        x *= g
        if(with_add):
            z += x
        else:
            z = x

        return z

    def forward(self, 
        z: torch.Tensor, 
        mask: Optional[torch.Tensor] = None,
        inplace_safe: bool = False,
        _add_with_inplace: bool = False,
    ) -> torch.Tensor:
        """
        Args:
            x:
                [*, N_res, N_res, C_z] input tensor
            mask:
                [*, N_res, N_res] input mask
        Returns:
            [*, N_res, N_res, C_z] output tensor
        """
        if(inplace_safe):
            x = self._inference_forward(
                z, 
                mask, 
                with_add=_add_with_inplace,
            )
            return x

        if mask is None:
            mask = z.new_ones(z.shape[:-1])

        mask = mask.unsqueeze(-1)
        
        z = self.layer_norm_in(z)
        a = mask
        a = a * self.sigmoid(self.linear_a_g(z)) 
        a = a * self.linear_a_p(z)

        b = mask
        b = b * self.sigmoid(self.linear_b_g(z))
        b = b * self.linear_b_p(z)

        x = self._combine_projections(a, b)
        
        del a, b
        x = self.layer_norm_out(x)
        x = self.linear_z(x)
        g = self.sigmoid(self.linear_g(z))
        x = x * g

        return x


class TriangleMultiplicationOutgoing(TriangleMultiplicativeUpdate):
    """
    Implements Algorithm 11.
    """
    __init__ = partialmethod(TriangleMultiplicativeUpdate.__init__, _outgoing=True)


class TriangleMultiplicationIncoming(TriangleMultiplicativeUpdate):
    """
    Implements Algorithm 12.
    """
    __init__ = partialmethod(TriangleMultiplicativeUpdate.__init__, _outgoing=False)

class PairTransition(nn.Module):
    """
    Implements Algorithm 15.
    """

    def __init__(self, c_z, n):
        """
        Args:
            c_z:
                Pair transition channel dimension
            n:
                Factor by which c_z is multiplied to obtain hidden channel
                dimension
        """
        super(PairTransition, self).__init__()

        self.c_z = c_z
        self.n = n

        self.layer_norm = nn.LayerNorm(self.c_z)
        self.linear_1 = Linear(self.c_z, self.n * self.c_z, init="relu")
        self.relu = nn.ReLU()
        self.linear_2 = Linear(self.n * self.c_z, c_z, init="final")

    def _transition(self, z, mask):
        # [*, N_res, N_res, C_z]
        z = self.layer_norm(z)
        
        # [*, N_res, N_res, C_hidden]
        z = self.linear_1(z)
        z = self.relu(z)

        # [*, N_res, N_res, C_z]
        z = self.linear_2(z) 
        z = z * mask

        return z
    

    def forward(self, 
        z: torch.Tensor, 
        mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Args:
            z:
                [*, N_res, N_res, C_z] pair embedding
        Returns:
            [*, N_res, N_res, C_z] pair embedding update
        """
        if mask is None:
            mask = z.new_ones(z.shape[:-1])

        # [*, N_res, N_res, 1]
        mask = mask.unsqueeze(-1)

        z = self._transition(z=z, mask=mask)

        return z