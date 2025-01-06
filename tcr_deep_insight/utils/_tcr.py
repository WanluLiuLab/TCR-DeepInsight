from typing import Literal
from Bio import SeqIO
from pathlib import Path
import difflib
import numpy as np
import pickle

from ._amino_acids import _AMINO_ACIDS_INDEX, _AMINO_ACIDS_INDEX_REVERSE, seq_encode_aa, seq_decode_aa
from ._tcr_definitions import TRAV2INDEX, TRAJ2INDEX, TRBV2INDEX, TRBJ2INDEX

MODULE_PATH = Path(__file__).parent

class TCR:
    """
    TCR class to store TCR information and provide utility functions.

    :param cdr3a: CDR3 alpha sequence
    :param cdr3b: CDR3 beta sequence
    :param trav: TRAV gene
    :param trbv: TRBV gene
    :param traj: TRAJ gene
    :param trbj: TRBJ gene
    :param individual: Individual identifier
    :param species: Species identifier
    """
    def __init__(self, 
        cdr3a: str = None,
        cdr3b: str = None,
        trav: str = None,
        trbv: str = None,
        traj: str = None,
        trbj: str = None,
        individual: str = None,
        species: Literal['human','mouse'] = "human",
    ):
        self._species = species
        self._individual = individual
        self._cdr3a, self._cdr3b = None, None
        self._cdr3a_length, self._cdr3b_length = None, None
        self._trav, self._trbv, self._traj, self._trbj = None, None, None, None
        if cdr3a is not None:
            self._cdr3a = seq_encode_aa(cdr3a)
            self._cdr3a_length = len(cdr3a)
        if cdr3b is not None:
            self._cdr3b = seq_encode_aa(cdr3b)
            self._cdr3b_length = len(cdr3b)
        if trav is not None:
            self._trav = TRAV2INDEX[self._species][trav]
        if trbv is not None:
            self._trbv = TRBV2INDEX[self._species][trbv]
        if traj is not None:
            self._traj = TRAJ2INDEX[self._species][traj]
        if trbj is not None:
            self._trbj = TRBJ2INDEX[self._species][trbj]

    @property
    def species(self):
        return self._species

    @property
    def cdr3a(self):
        return seq_decode_aa(self._cdr3a, self._cdr3a_length)

    @property
    def cdr3b(self):
        return seq_decode_aa(self._cdr3b, self._cdr3b_length)

    @property
    def trav(self):
        return list(TRAV2INDEX[self._species].keys())[list(TRAV2INDEX[self._species].values()).index(self._trav)]

    @property
    def trbv(self):
        return list(TRBV2INDEX[self._species].keys())[list(TRBV2INDEX[self._species].values()).index(self._trbv)]

    @property
    def traj(self):
        return list(TRAJ2INDEX[self._species].keys())[list(TRAJ2INDEX[self._species].values()).index(self._traj)]

    @property
    def trbj(self):
        return list(TRBJ2INDEX[self._species].keys())[list(TRBJ2INDEX[self._species].values()).index(self._trbj)]

    @property
    def individual(self):
        return self._individual

    def __eq__(self, other):
        if self._individual is None and other._individual is None:
            return (
                self.cdr3a == other.cdr3a
                and self.cdr3b == other.cdr3b
                and self.trav == other.trav
                and self.trbv == other.trbv
                and self.traj == other.traj
                and self.trbj == other.trbj
            )
        else:
            return (
                self.cdr3a == other.cdr3a
                and self.cdr3b == other.cdr3b
                and self.trav == other.trav
                and self.trbv == other.trbv
                and self.traj == other.traj
                and self.trbj == other.trbj
                and self.individual == other.individual
            )

    def __hash__(self):
        return hash((self.cdr3a, self.cdr3b, self.trav, self.trbv, self.traj, self.trbj))

    def __str__(self):
        return f"TCR(species='{self._species}', cdr3a='{self.cdr3a}', cdr3b='{self.cdr3b}', trav='{self.trav}', trbv='{self.trbv}', traj='{self.traj}', trbj='{self.trbj}')"

    def to_string(self):
        return "{}={}={}={}={}={}={}={}".format(
            self.cdr3a,
            self.cdr3b,
            self.trav,
            self.traj,
            self.trbv,
            self.trbj,
            self.individual,
            self.species
        )

    def to_tcr_string(self):
        return "{}={}={}={}={}".format(
            self.cdr3a,
            self.cdr3b,
            self.trav,
            self.traj,
            self.trbv,
            self.trbj
        )

    @classmethod
    def from_string(cls, string: str):
        cdr3a, cdr3b, trav, traj, trbv, trbj, individual, *species = string.split("=")
        species = species[0] if species else "human"
        return cls(
            cdr3a=cdr3a,
            cdr3b=cdr3b,
            trav=trav,
            trbv=trbv,
            traj=traj,
            trbj=trbj,
            individual=individual,
            species=species,
        )

    def serialize(self):
        return pickle.dumps(self)

    @classmethod
    def deserialize(cls, string):
        return pickle.loads(string)
