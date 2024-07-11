import scanpy as sc 
import datasets
import torch 
import numpy as np
import pandas as pd
from typing import Any, List, Mapping
from collections import OrderedDict
from functools import partial

from scatlasvae.utils._decorators import typed

from ..utils._logger import get_tqdm, mt, mw
from ..utils._compat import Literal

class ModelOutput(OrderedDict):
    def __init__(self, __dict, **kwargs):
        super().__init__(**kwargs)
        self.__dict__ = self 
        for k,v in __dict.items():
            self[k] = v

    def __getitem__(self, __key: Any) -> Any:
        if isinstance(__key, int):
            return self.__dict__[list(self.__dict__.keys())[__key]]
        else:
            return super().__getitem__(__key)
        
class EarlyStopping():
    """
    Early stopping to stop the training when the loss does not improve after
    certain epochs.
    """
    def __init__(self, patience=5, min_delta=0):
        """
        :param patience: how many epochs to wait before stopping when loss is
               not improving
        :param min_delta: minimum difference between new loss and old loss for
               new loss to be considered as an improvement
        """
        self.patience = patience
        self.min_delta = min_delta
        self.counter = 0
        self.best_loss = None
        self.early_stop = False
    def __call__(self, val_loss):
        if self.best_loss == None:
            self.best_loss = val_loss
        elif self.best_loss - val_loss > self.min_delta:
            self.best_loss = val_loss
            # reset counter if validation loss improves
            self.counter = 0
        elif self.best_loss - val_loss < self.min_delta:
            self.counter += 1
            mt(f"INFO: Early stopping counter {self.counter} of {self.patience}")
            if self.counter >= self.patience:
                mt('INFO: Early stopping')
                self.early_stop = True

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


def add(m1, m2, inplace):
    # The first operation in a checkpoint can't be in-place, but it's
    # nice to have in-place addition during inference. Thus...
    if(not inplace):
        m1 = m1 + m2
    else:
        m1 += m2

    return m1

def permute_final_dims(tensor: torch.Tensor, inds: List[int]):
    zero_index = -1 * len(inds)
    first_inds = list(range(len(tensor.shape[:zero_index])))
    return tensor.permute(first_inds + [zero_index + i for i in inds])

def flatten_final_dims(t: torch.Tensor, no_dims: int):
    return t.reshape(t.shape[:-no_dims] + (-1,))

@torch.jit.ignore
def softmax_no_cast(t: torch.Tensor, dim: int = -1) -> torch.Tensor:
    """
        Softmax, but without automatic casting to fp32 when the input is of
        type bfloat16
    """
    d = t.dtype

    if(d is torch.bfloat16):
        with torch.cuda.amp.autocast(enabled=False):
            s = torch.nn.functional.softmax(t, dim=dim)
    else:
        s = torch.nn.functional.softmax(t, dim=dim)

    return s

def partial_load_state_dict(model, state_dict: Mapping[str, torch.Tensor]):
    """
    Partially load the state dict
    :param state_dict: Mapping[str, torch.Tensor]. State dict to load
    """
    original_state_dict = model.state_dict()
    ignored_keys = []
    for k,v in state_dict.items():
        if k not in original_state_dict.keys():
            mt(f"Warning: {k} not found in the model. Ignoring {k} in the provided state dict.")
            ignored_keys.append(k)
        elif v.shape != original_state_dict[k].shape:
            mw(f"Warning: shape of {k} does not match. \n" + \
                ' '*40 + "\tOriginal:" + f" {original_state_dict[k].shape},\n" + \
                ' '*40 + f"\tNew: {v.shape}")
            state_dict[k] = original_state_dict[k]
    for k,v in original_state_dict.items():
        if k not in state_dict.keys():
            mw(f"Warning: {k} not found in the provided state dict. " + \
                 f"Using {k} in the original state dict.")
            state_dict[k] = v
    for i in ignored_keys:
        state_dict.pop(i)
    model.load_state_dict(state_dict)