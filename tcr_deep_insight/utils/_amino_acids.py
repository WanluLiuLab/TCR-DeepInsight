import json
from enum import Enum
from collections import OrderedDict
from pathlib import Path
import pandas as pd

_AMINO_ACIDS = ['R',
                'H',
                'K',
                'D',
                'E',
                'S',
                'T',
                'N',
                'Q',
                'C',
                'U',
                'G',
                'P',
                'A',
                'V',
                'I',
                'L',
                'M',
                'F',
                'Y',
                'W']

_AMINO_ACIDS_MOTIF3 = []
for i in _AMINO_ACIDS:
    for j in _AMINO_ACIDS:
        for n in _AMINO_ACIDS:
            _AMINO_ACIDS_MOTIF3.append(i+j+n)

_AMINO_ACIDS_ADDITIONALS = OrderedDict(
    PAD = '.',
    MASK = '#',
    UNK = '*',
    SEP = '|',
    CLS = '^',
    GAP = ':'
)

_AMINO_ACIDS_INDEX = dict(
    zip(
        _AMINO_ACIDS + list(_AMINO_ACIDS_ADDITIONALS.values()), 
        list(range(len(_AMINO_ACIDS) + len(_AMINO_ACIDS_ADDITIONALS)))
    )
)

_AMINO_ACIDS_INDEX_REVERSE = {v:k for k,v in _AMINO_ACIDS_INDEX.items()}

                      
_factor_names = ['Factor I', 'Factor II', 'Factor III', 'Factor IV', 'Factor V']

_AMINO_ACIDS_ONE_HOT = _AMINO_ACIDS_INDEX