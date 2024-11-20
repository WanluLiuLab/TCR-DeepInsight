from abc import ABC, ABCMeta
from enum import Enum, EnumMeta, unique
from functools import wraps
from typing import Any, Callable

from ..utils._definitions import PrettyEnum, ModeEnum

@unique
class TCR_BERT_ENCODING(ModeEnum):
    VJCDR3 = 'vjcdr3'
    CDR123 = 'cdr123'

@unique
class TCR_BERT_POOLING(ModeEnum):
    MEAN = 'mean'
    SUM = 'sum'
    CLS = 'cls'
    MAX = 'max'
    TRA = 'tra'
    TRB = 'trb'
    WEIGHTED = 'weighted'
    POOL = 'pool'