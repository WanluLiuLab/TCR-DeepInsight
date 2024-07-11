from typing import Any, List
from collections import OrderedDict
import torch
from functools import partial

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

