from abc import ABC, ABCMeta
from enum import Enum, EnumMeta, unique
from functools import wraps
from typing import Any, Callable


class PrettyEnum(Enum):
    """Enum with a pretty :meth:`__str__` and :meth:`__repr__`."""

    @property
    def v(self) -> Any:
        """Alias for :attr`value`."""
        return self.value

    def __repr__(self) -> str:
        return f"{self.value!r}"

    def __str__(self) -> str:
        return f"{self.value!s}"
    
class ModeEnum(str, PrettyEnum, metaclass=EnumMeta):
    """Enum with a pretty :meth:`__str__` and :meth:`__repr__`."""


class FCDEF(Enum):
    ENCODER = 0x0
    DECODER = 0x1


@unique
class SPECIES(ModeEnum):
    HUMAN = 'human'
    MOUSE = 'mouse'
