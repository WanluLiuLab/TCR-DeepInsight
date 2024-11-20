from abc import ABC, ABCMeta
from enum import Enum, EnumMeta, unique
from functools import wraps
from typing import Any, Callable

from ..utils._definitions import PrettyEnum, ModeEnum

@unique
class FAISS_INDEX_BACKEND(ModeEnum):
    KMEANS = "kmeans"
    FLAT = "flat"
