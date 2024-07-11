from functools import partial
import math
from typing import Dict, Optional, Callable, List, Tuple, Sequence, Literal, Mapping, Any
import numpy as np
import torch
import torch.nn as nn
from scipy.stats import truncnorm
import einops
from functools import partial

####################################################
#                    Tensors                       #
####################################################

def dict_map(fn, dic, leaf_type):
    new_dict = {}
    for k, v in dic.items():
        if type(v) is dict:
            new_dict[k] = dict_map(fn, v, leaf_type)
        else:
            new_dict[k] = tree_map(fn, v, leaf_type)

    return new_dict

def tree_map(fn, tree, leaf_type):
    if isinstance(tree, dict):
        return dict_map(fn, tree, leaf_type)
    elif isinstance(tree, list):
        return [tree_map(fn, x, leaf_type) for x in tree]
    elif isinstance(tree, tuple):
        return tuple([tree_map(fn, x, leaf_type) for x in tree])
    elif isinstance(tree, leaf_type):
        return fn(tree)
    else:
        print(type(tree))
        raise ValueError("Not supported")


tensor_tree_map = partial(tree_map, leaf_type=torch.Tensor)


def dict_multimap(fn, dicts):
    first = dicts[0]
    new_dict = {}
    for k, v in first.items():
        all_v = [d[k] for d in dicts]
        if type(v) is dict:
            new_dict[k] = dict_multimap(fn, all_v)
        else:
            new_dict[k] = fn(all_v)

    return new_dict

def permute_final_dims(tensor: torch.Tensor, inds: List[int]):
    zero_index = -1 * len(inds)
    first_inds = list(range(len(tensor.shape[:zero_index])))
    return tensor.permute(first_inds + [zero_index + i for i in inds])


####################################################
#                    Linear                        #
####################################################

def _prod(nums):
    out = 1
    for n in nums:
        out = out * n
    return out

def _calculate_fan(linear_weight_shape, fan="fan_in"):
    fan_out, fan_in = linear_weight_shape

    if fan == "fan_in":
        f = fan_in
    elif fan == "fan_out":
        f = fan_out
    elif fan == "fan_avg":
        f = (fan_in + fan_out) / 2
    else:
        raise ValueError("Invalid fan option")

    return f

def _trunc_normal_init_(weights, scale=1.0, fan="fan_in"):
    shape = weights.shape
    f = _calculate_fan(shape, fan)
    scale = scale / max(1, f)
    a = -2
    b = 2
    std = math.sqrt(scale) / truncnorm.std(a=a, b=b, loc=0, scale=1)
    size = _prod(shape)
    samples = truncnorm.rvs(a=a, b=b, loc=0, scale=std, size=size)
    samples = np.reshape(samples, shape)
    with torch.no_grad():
        weights.copy_(torch.tensor(samples, device=weights.device))
        
def lecun_normal_init_(weights):
    _trunc_normal_init_(weights, scale=1.0)


def he_normal_init_(weights):
    _trunc_normal_init_(weights, scale=2.0)


def glorot_uniform_init_(weights):
    nn.init.xavier_uniform_(weights, gain=1)


def final_init_(weights):
    with torch.no_grad():
        weights.fill_(0.0)


def gating_init_(weights):
    with torch.no_grad():
        weights.fill_(0.0)


def normal_init_(weights):
    torch.nn.init.kaiming_normal_(weights, nonlinearity="linear")
    
class Linear(nn.Linear):
    """
    A Linear layer with built-in nonstandard initializations. Called just
    like torch.nn.Linear.
    """
    
    def __init__(
        self,
        in_dim: int,
        out_dim: int,
        bias: bool = True,
        init: Literal["default","final","gating","glorot","normal","relu"] = "default",
        init_fn: Optional[Callable[[torch.Tensor, torch.Tensor], None]] = None,
    ):
        """
        Args:
            in_dim:
                The final dimension of inputs to the layer
            out_dim:
                The final dimension of layer outputs
            bias:
                Whether to learn an additive bias. True by default
            init:
                The initializer to use. Choose from:

                "default": LeCun fan-in truncated normal initialization
                "relu": He initialization w/ truncated normal distribution
                "glorot": Fan-average Glorot uniform initialization
                "gating": Weights=0, Bias=1
                "normal": Normal initialization with std=1/sqrt(fan_in)
                "final": Weights=0, Bias=0

                Overridden by init_fn if the latter is not None.
            init_fn:
                A custom initializer taking weight and bias as inputs.
                Overrides init if not None.
        """
        super(Linear, self).__init__(in_dim, out_dim, bias=bias)

        if bias:
            with torch.no_grad():
                self.bias.fill_(0)

        with torch.no_grad():
            if init_fn is not None:
                init_fn(self.weight, self.bias)
            elif init == "default":
                lecun_normal_init_(self.weight)
            elif init == "final":
                final_init_(self.weight)
            elif init == "gating":
                gating_init_(self.weight)
                if bias:
                    self.bias.fill_(1.0)
            elif init == "glorot":
                glorot_uniform_init_(self.weight)
            elif init == "normal":
                normal_init_(self.weight)
            elif init == "relu":
                he_normal_init_(self.weight)
            else:
                raise ValueError("Invalid init string.")

###################################################
#                  Embeddings                     #
###################################################

class RotaryEmbedding:
    @staticmethod
    def fixed_position_embedding(
        x: torch.Tensor,
        seq_dim: int,
        seq_len: Optional[int] = None,
    ):
        d = x.shape[-1]
        if seq_len is None:
            seq_len = x.shape[seq_dim]
        # create constant "pe" matrix with values dependant on    
        inv_freq = 1.0 / (10000 ** (torch.arange(0, d, 2.0) / d))
        # compute positional encodings for each position and feature
        sinusoid_inp = torch.einsum("i,j->ij", torch.arange(seq_len), inv_freq).to(
            device=x.device
        ).float()
        return torch.sin(sinusoid_inp), torch.cos(sinusoid_inp)
    
    @staticmethod
    def rotate_every_two(x: torch.Tensor):
        x1 = x[:, :, :, 0::2]
        x2 = x[:, :, :, 1::2]
        x = torch.stack((-x2, x1), dim=-1)
        return einops.rearrange(x, "... d j -> ... (d j)")
    
    @staticmethod
    def apply_rotary_pos_emb(
        x,
        sincos,
        offset=0
    ):
        sin, cos = map(lambda t: einops.repeat(
            t[offset:x.shape[1] + offset, :], 
            "n d -> () n () (d j)", j = 2
        ), sincos)
        return (x * cos) + (RotaryEmbedding.rotate_every_two(x) * sin)
    
    
####################################################
#                    Chunks                        #
####################################################

class ChunkLayer:
    def __init__(
        self,
        chunk_size: int,
        _low_mem: bool = False,
        _add_into_out: bool = False 
    ):
        self.chunk_size = chunk_size
        self._low_mem = _low_mem 
        self._add_into_out = _add_into_out

    @staticmethod
    @torch.jit.ignore
    def _flat_idx_to_idx(
        flat_idx: int,
        dims: Tuple[int],
    ) -> Tuple[int]:
        idx = []
        for d in reversed(dims):
            idx.append(flat_idx % d)
            flat_idx = flat_idx // d

        return tuple(reversed(idx))

    @staticmethod
    @torch.jit.ignore
    def _get_minimal_slice_set(
        start: Sequence[int],
        end: Sequence[int],
        dims: int,
        start_edges: Optional[Sequence[bool]] = None,
        end_edges: Optional[Sequence[bool]] = None,
    ) -> Sequence[Tuple[int]]:
        """ 
            Produces an ordered sequence of tensor slices that, when used in
            sequence on a tensor with shape dims, yields tensors that contain every
            leaf in the contiguous range [start, end]. Care is taken to yield a 
            short sequence of slices, and perhaps even the shortest possible (I'm 
            pretty sure it's the latter).
            
            end is INCLUSIVE. 
        """
        # start_edges and end_edges both indicate whether, starting from any given
        # dimension, the start/end index is at the top/bottom edge of the
        # corresponding tensor, modeled as a tree
        def reduce_edge_list(l):
            tally = 1
            for i in range(len(l)):
                reversed_idx = -1 * (i + 1)
                l[reversed_idx] *= tally
                tally = l[reversed_idx]

        if(start_edges is None):
            start_edges = [s == 0 for s in start]
            reduce_edge_list(start_edges)
        if(end_edges is None):
            end_edges = [e == (d - 1) for e,d in zip(end, dims)]
            reduce_edge_list(end_edges)        

        # Base cases. Either start/end are empty and we're done, or the final,
        # one-dimensional tensor can be simply sliced
        if(len(start) == 0):
            return [tuple()]
        elif(len(start) == 1):
            return [(slice(start[0], end[0] + 1),)]

        slices = []
        path = []
    
        # Dimensions common to start and end can be selected directly
        for s,e in zip(start, end):
            if(s == e):
                path.append(slice(s, s + 1))
            else:
                break

        path = tuple(path)
        divergence_idx = len(path)

        # start == end, and we're done
        if(divergence_idx == len(dims)):
            return [tuple(path)]

        def upper():
            sdi = start[divergence_idx]
            return [
                path + (slice(sdi, sdi + 1),) + s for s in 
                ChunkLayer._get_minimal_slice_set(
                    start[divergence_idx + 1:],
                    [d - 1 for d in dims[divergence_idx + 1:]],
                    dims[divergence_idx + 1:],
                    start_edges=start_edges[divergence_idx + 1:],
                    end_edges=[1 for _ in end_edges[divergence_idx + 1:]]
                )
            ]

        def lower():
            edi = end[divergence_idx]
            return [
                path + (slice(edi, edi + 1),) + s for s in 
                ChunkLayer._get_minimal_slice_set(
                    [0 for _ in start[divergence_idx + 1:]],
                    end[divergence_idx + 1:],
                    dims[divergence_idx + 1:],
                    start_edges=[1 for _ in start_edges[divergence_idx + 1:]],
                    end_edges=end_edges[divergence_idx + 1:],
                )
            ]

        # If both start and end are at the edges of the subtree rooted at
        # divergence_idx, we can just select the whole subtree at once
        if(start_edges[divergence_idx] and end_edges[divergence_idx]):
            slices.append(
                path + (slice(start[divergence_idx], end[divergence_idx] + 1),)
            )
        # If just start is at the edge, we can grab almost all of the subtree, 
        # treating only the ragged bottom edge as an edge case
        elif(start_edges[divergence_idx]):
            slices.append(
                path + (slice(start[divergence_idx], end[divergence_idx]),)
            )
            slices.extend(lower())
        # Analogous to the previous case, but the top is ragged this time
        elif(end_edges[divergence_idx]):
            slices.extend(upper())
            slices.append(
                path + (slice(start[divergence_idx] + 1, end[divergence_idx] + 1),)
            )
        # If both sides of the range are ragged, we need to handle both sides
        # separately. If there's contiguous meat in between them, we can index it
        # in one big chunk
        else:
            slices.extend(upper())
            middle_ground = end[divergence_idx] - start[divergence_idx]
            if(middle_ground > 1):
                slices.append(
                    path + (slice(start[divergence_idx] + 1, end[divergence_idx]),)
                )
            slices.extend(lower())

        return [tuple(s) for s in slices]


    @torch.jit.ignore
    def _chunk_slice(
        self,
        t: torch.Tensor,
        flat_start: int,
        flat_end: int,
        no_batch_dims: int,
    ) -> torch.Tensor:
        """
            Equivalent to
            
                t.reshape((-1,) + t.shape[no_batch_dims:])[flat_start:flat_end]

            but without the need for the initial reshape call, which can be 
            memory-intensive in certain situations. The only reshape operations
            in this function are performed on sub-tensors that scale with
            (flat_end - flat_start), the chunk size.
        """

        batch_dims = t.shape[:no_batch_dims]
        start_idx = list(ChunkLayer._flat_idx_to_idx(flat_start, batch_dims))
        # _get_minimal_slice_set is inclusive
        end_idx = list(ChunkLayer._flat_idx_to_idx(flat_end - 1, batch_dims))

        # Get an ordered list of slices to perform
        slices = ChunkLayer._get_minimal_slice_set(
            start_idx,
            end_idx,
            batch_dims,
        )

        sliced_tensors = [t[s] for s in slices]

        return torch.cat(
            [s.view((-1,) + t.shape[no_batch_dims:]) for s in sliced_tensors]
        )
    
    def _prep_inputs(self, t, no_batch_dims, orig_batch_dims):
        if not self._low_mem:
            if not sum(t.shape[:no_batch_dims]) == no_batch_dims:
                t = t.expand(orig_batch_dims + t.shape[no_batch_dims:])
            t = t.reshape(-1, *t.shape[no_batch_dims:])
        else:
            t = t.expand(orig_batch_dims + t.shape[no_batch_dims:])
        return t
    
    def _fetch_dims_from_tree(self, tree: Optional[Mapping[str, Any]]):
        shapes = []
        if isinstance(tree, dict):
            for v in tree.values():
                shapes.extend(self._fetch_dims_from_tree(v))
        elif isinstance(tree, list) or isinstance(tree, tuple):
            for t in tree: 
                shapes.extend(self._fetch_dims_from_tree(t))
        elif isinstance(tree, torch.Tensor):
            shapes.append(tree.shape)
        else:
            raise ValueError("Tree type not supported")
        return shapes
    
    @staticmethod
    def _fetch_dims(tree):
        shapes = []
        tree_type = type(tree)
        if tree_type is dict:
            for v in tree.values():
                shapes.extend(_fetch_dims(v))
        elif tree_type is list or tree_type is tuple:
            for t in tree:
                shapes.extend(_fetch_dims(t))
        elif tree_type is torch.Tensor:
            shapes.append(tree.shape)
        else:
            raise ValueError("Not supported")

        return shapes
    
    @staticmethod
    def chunk_layer(
        layer: Callable,
        inputs: Dict[str, Any],
        chunk_size: int,
        no_batch_dims: int,
        low_mem: bool = False,
        _out: Any = None,
        _add_into_out: bool = False,
    ) -> Any:
        """
        Implements the "chunking" procedure described in section 1.11.8.

        Layer outputs and inputs are assumed to be simple "pytrees,"
        consisting only of (arbitrarily nested) lists, tuples, and dicts with
        torch.Tensor leaves.

        Args:
            layer:
                The layer to be applied chunk-wise
            inputs:
                A (non-nested) dictionary of keyworded inputs. All leaves must
                be tensors and must share the same batch dimensions.
            chunk_size:
                The number of sub-batches per chunk. If multiple batch
                dimensions are specified, a "sub-batch" is defined as a single
                indexing of all batch dimensions simultaneously (s.t. the
                number of sub-batches is the product of the batch dimensions).
            no_batch_dims:
                How many of the initial dimensions of each input tensor can
                be considered batch dimensions.
            low_mem:
                Avoids flattening potentially large input tensors. Unnecessary
                in most cases, and is ever so slightly slower than the default
                setting.
        Returns:
            The reassembled output of the layer on the inputs.
        """
        if not (len(inputs) > 0):
            raise ValueError("Must provide at least one input")

        initial_dims = [shape[:no_batch_dims] for shape in ChunkLayer._fetch_dims(inputs)]
        orig_batch_dims = tuple([max(s) for s in zip(*initial_dims)])

        def _prep_inputs(t):
            if(not low_mem):
                if not sum(t.shape[:no_batch_dims]) == no_batch_dims:
                    t = t.expand(orig_batch_dims + t.shape[no_batch_dims:])
                t = t.reshape(-1, *t.shape[no_batch_dims:])
            else:
                t = t.expand(orig_batch_dims + t.shape[no_batch_dims:])
            return t

        prepped_inputs = tensor_tree_map(_prep_inputs, inputs)
        prepped_outputs = None
        if(_out is not None):
            reshape_fn = lambda t: t.view([-1] + list(t.shape[no_batch_dims:]))
            prepped_outputs = tensor_tree_map(reshape_fn, _out)

        flat_batch_dim = 1
        for d in orig_batch_dims:
            flat_batch_dim *= d

        no_chunks = flat_batch_dim // chunk_size + (
            flat_batch_dim % chunk_size != 0
        )

        i = 0
        out = prepped_outputs
        for _ in range(no_chunks):
            # Chunk the input
            if(not low_mem):
                select_chunk = (
                    lambda t: t[i : i + chunk_size] if t.shape[0] != 1 else t
                )
            else:
                select_chunk = (
                    partial(
                        ChunkLayer._chunk_slice, 
                        flat_start=i, 
                        flat_end=min(flat_batch_dim, i + chunk_size), 
                        no_batch_dims=len(orig_batch_dims)
                    )
                )

            chunks = tensor_tree_map(select_chunk, prepped_inputs)

            # Run the layer on the chunk
            output_chunk = layer(**chunks)

            # Allocate space for the output
            if out is None:
                allocate = lambda t: t.new_zeros((flat_batch_dim,) + t.shape[1:])
                out = tensor_tree_map(allocate, output_chunk)

            # Put the chunk in its pre-allocated space
            out_type = type(output_chunk)
            if out_type is dict:
                def assign(d1, d2):
                    for k, v in d1.items():
                        if type(v) is dict:
                            assign(v, d2[k])
                        else:
                            if(_add_into_out):
                                v[i: i + chunk_size] += d2[k]
                            else:
                                v[i: i + chunk_size] = d2[k]

                assign(out, output_chunk)
            elif out_type is tuple:
                for x1, x2 in zip(out, output_chunk):
                    if(_add_into_out):
                        x1[i: i + chunk_size] += x2
                    else:
                        x1[i : i + chunk_size] = x2
            elif out_type is torch.Tensor:
                if(_add_into_out):
                    out[i: i + chunk_size] += output_chunk
                else:
                    out[i: i + chunk_size] = output_chunk
            else:
                raise ValueError("Not supported")

            i += chunk_size

        reshape = lambda t: t.view(orig_batch_dims + t.shape[1:])
        out = tensor_tree_map(reshape, out)

        return out



    def __call__(
        self,
        layer: Callable,
        inputs: Optional[Mapping[str, Any]],
        no_batch_dims: int,
        _out: Any = None,
        _add_into_out: bool = False,
    ) -> Any: 
        """
        Layer outputs and inputs are assumed to be simple "pytrees,"
        consisting only of (arbitrarily nested) lists, tuples, and dicts with
        torch.Tensor leaves.

        Args:
            layer:
                The layer to be applied chunk-wise
            inputs:
                A (non-nested) dictionary of keyworded inputs. All leaves must
                be tensors and must share the same batch dimensions.
            chunk_size:
                The number of sub-batches per chunk. If multiple batch
                dimensions are specified, a "sub-batch" is defined as a single
                indexing of all batch dimensions simultaneously (s.t. the
                number of sub-batches is the product of the batch dimensions).
            no_batch_dims:
                How many of the initial dimensions of each input tensor can
                be considered batch dimensions.
            low_mem:
                Avoids flattening potentially large input tensors. Unnecessary
                in most cases, and is ever so slightly slower than the default
                setting.
        Returns:
            The reassembled output of the layer on the inputs.
        """
        if len(inputs) == 0:
            raise ValueError("Must provide at least one input")
        
        initial_dims = [shape[:no_batch_dims] for shape in self._fetch_dims_from_tree(inputs)]
        orig_batch_dims = tuple([max(s) for s in zip(*initial_dims)])
        prepared_inputs = tensor_tree_map(
            partial(
                self._prep_inputs, 
                no_batch_dims=no_batch_dims, 
                orig_batch_dims=orig_batch_dims
            ), 
            inputs
        )
        prepared_outputs = None
        if _out is not None:
            reshape_fn = lambda t: t.view([-1] + list(t.shape[no_batch_dims:]))
            prepared_outputs = tensor_tree_map(reshape_fn, _out)
        flat_batch_dim = 1
        for d in orig_batch_dims:
            flat_batch_dim *= d

        no_chunks = flat_batch_dim // self.chunk_size + (
            flat_batch_dim % self.chunk_size != 0
        )

        i = 0
        out = prepared_outputs
        for _ in range(no_chunks):
            if not self._low_mem:
                select_chunk = lambda t: t[i : i + self.chunk_size] if t.shape[0] != 1 else t
            else: 
                select_chunk = (
                    partial(
                        self._chunk_slice, 
                        flat_start=i, 
                        flat_end=min(flat_batch_dim, i + self.chunk_size), 
                        no_batch_dims=len(orig_batch_dims)
                    )
                )
            chunks = tensor_tree_map(select_chunk, prepared_inputs)
            # Run the layer on the chunk
            output_chunk = layer(**chunks)
            
            # Allocate space for the output
            if out is None:
                allocate = lambda t: t.new_zeros((flat_batch_dim,) + t.shape[1:])
                out = tensor_tree_map(allocate, output_chunk)
                
            # Put the chunk in its pre-allocated space
            out_type = type(output_chunk)
            if out_type is dict:
                def assign(d1, d2):
                    for k, v in d1.items():
                        if type(v) is dict:
                            assign(v, d2[k])
                        else:
                            if(_add_into_out):
                                v[i: i + self.chunk_size] += d2[k]
                            else:
                                v[i: i + self.chunk_size] = d2[k]

                assign(out, output_chunk)
            elif out_type is tuple:
                for x1, x2 in zip(out, output_chunk):
                    if(_add_into_out):
                        x1[i: i + self.chunk_size] += x2
                    else:
                        x1[i : i + self.chunk_size] = x2
            elif out_type is torch.Tensor:
                if(_add_into_out):
                    out[i: i + self.chunk_size] += output_chunk
                else:
                    out[i: i + self.chunk_size] = output_chunk
            else:
                raise ValueError("Not supported")

            i += self.chunk_size

        reshape = lambda t: t.view(orig_batch_dims + t.shape[1:])
        out = tensor_tree_map(reshape, out)

        return out

