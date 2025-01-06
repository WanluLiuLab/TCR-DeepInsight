from abc import ABC, ABCMeta
from enum import Enum, EnumMeta, unique
from functools import wraps
from typing import Any, Callable

from ..utils._definitions import PrettyEnum, ModeEnum

@unique
class FAISS_INDEX_BACKEND(ModeEnum):
    KMEANS = "kmeans"
    FLAT = "flat"

@unique 
class TDI_RESULT_FIELD(PrettyEnum):
    NUMBER_OF_UNIQUE_TCR = "number_of_unique_tcrs"
    NUMBER_OF_INDIVIDUAL = "number_of_individuals"
    NUMBER_OF_TCR = "number_of_tcrs"
    NUMBER_OF_CELL = "number_of_cells"
    CONVERGENCE = "tcr_convergence_score"
    DISEASE_ASSOCIATION = "disease_association_score"