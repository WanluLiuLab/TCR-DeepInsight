import pandas as pd
import torch
from torch.nn import functional as F
from typing import Union, Iterable, List, Tuple, Optional
from einops import rearrange, repeat
from sklearn.model_selection import train_test_split
import datasets
import tqdm
import numpy as np
from pathlib import Path
        
from tokenizers import Tokenizer, Encoding

from transformers import (
    PreTrainedTokenizerBase,
)

from scatlasvae.utils._tensor_utils import get_k_elements

from ...utils._tcr import (
    HumanTCRAnnotations,
    MouseTCRAnnotations
)
from ...utils._definitions import (
    TRAV2CDR1a,
    TRAV2CDR2a,
    TRBV2CDR1b,
    TRBV2CDR2b,
)
from ...utils._compat import Literal
from ...utils._decorators import deprecated
from ...utils._amino_acids import (
    _AMINO_ACIDS_ADDITIONALS,
    _AMINO_ACIDS_INDEX_REVERSE,
    _AMINO_ACIDS_INDEX
)

MODULE_PATH = Path(__file__).parent

class AminoAcidTokenizer(PreTrainedTokenizerBase):
    """Tokenizer for amino acids. The amino acid to token index follows the same layout as tcr-bert"""
    def __init__(self, 
        *,
        model_max_length: int,
        append_cls_token: bool = True,
        append_eos_token: bool = True,
        **kwargs
    ) -> None:
        # A special token representing an out-of-vocabulary token.
        if "pad_token" not in kwargs.keys() or not kwargs.get("pad_token", None):
            kwargs["pad_token"] = _AMINO_ACIDS_ADDITIONALS["PAD"] 
        # A special token representing the class of the input
        if "unk_token" not in kwargs.keys()  or not kwargs.get("unk_token", None):
            kwargs["unk_token"] = _AMINO_ACIDS_ADDITIONALS["UNK"]
        # A special token representing a masked token (used by masked-language modeling pretraining objectives, like BERT). 
        if "mask_token" not in kwargs.keys()  or not kwargs.get("mask_token", None):
            kwargs["mask_token"] = _AMINO_ACIDS_ADDITIONALS["PAD"] 
        if "sep_token" not in kwargs.keys()  or not kwargs.get("sep_token", None):
            kwargs["sep_token"] = _AMINO_ACIDS_ADDITIONALS["SEP"] 
        if "cls_token" not in kwargs.keys()  or not kwargs.get("cls_token", None):
            kwargs["cls_token"] = _AMINO_ACIDS_ADDITIONALS["CLS"] 

        kwargs["model_max_length"] = model_max_length
        super(AminoAcidTokenizer, self).__init__(**kwargs)
        self._vocab_size = len(_AMINO_ACIDS_INDEX)
        self.append_cls_token = append_cls_token
        self.append_eos_token = append_eos_token
        self._pad_token_id = _AMINO_ACIDS_INDEX[kwargs["pad_token"]]
        self._unk_token_id = _AMINO_ACIDS_INDEX[kwargs["unk_token"]]
        self._mask_token_id = _AMINO_ACIDS_INDEX[kwargs["mask_token"]]
        self._sep_token_id = _AMINO_ACIDS_INDEX[kwargs["sep_token"]]
        self._cls_token_id = _AMINO_ACIDS_INDEX[kwargs["cls_token"]]

    @property
    def is_fast(self) -> bool:
        return False

    @property
    def vocab_size(self) -> int:
        return self._vocab_size

    def _encode(self, aa: str, max_length:int = None) -> torch.Tensor:
        if self.append_cls_token and not aa.startswith(
            _AMINO_ACIDS_ADDITIONALS["CLS"]
        ):
            aa = _AMINO_ACIDS_ADDITIONALS["CLS"] + aa

        if self.append_eos_token and not aa.startswith(
            _AMINO_ACIDS_ADDITIONALS["SEP"]
        ):
            aa += _AMINO_ACIDS_ADDITIONALS["SEP"]

        max_length = max_length or self.model_max_length
        if len(aa) < max_length:
            aa = aa + self._pad_token * (max_length - len(aa))
        return torch.Tensor(list(map(lambda a: _AMINO_ACIDS_INDEX[a], aa)))
    
    def _unpad(self, s) -> str:
        return ''.join(list(filter(lambda x: x != self.pad_token, s)))

    def _decode(self, ids) -> str:
        return self._unpad(list(map(lambda t: _AMINO_ACIDS_INDEX_REVERSE[t], ids)))

    def convert_tokens_to_ids(self, sequence: Union[Iterable[str], str]) -> torch.Tensor:
        if isinstance(sequence, str):
            ids = rearrange(self._encode(sequence), '(n h) -> n h', n = 1).type(torch.LongTensor)
        elif isinstance(sequence, Iterable):
            ids = torch.vstack(list(map(lambda x: self._encode(x), sequence))).type(torch.LongTensor)
        mask = self.convert_ids_to_mask(ids)
        return {"indices": ids, "mask": mask}

    def convert_ids_to_tokens(
        self, 
        ids: Union[torch.Tensor, np.ndarray, Iterable[int]]
    ) -> Iterable[str]:
        ids = ids.detach().cpu().numpy().astype(np.int64)
        if len(ids.shape) == 1:
            return self._decode(ids)
        else:
            return list(map(self._decode, ids))
    
    def convert_ids_to_mask(self, ids: torch.Tensor) -> torch.Tensor:
        return (ids != self._pad_token_id) & (ids != self._cls_token_id) & (ids != self._sep_token_id)

    def to_dataset(
        self, 
        ids: Iterable[str],
        chains: Iterable[str], 
        
        split: bool = False
    ) -> datasets.DatasetDict:
        # sourcery skip: remove-unnecessary-else, swap-if-else-branches
        assert(len(ids) == len(chains))
        tokenized = self.convert_tokens_to_ids(chains)
        if split:
            train_idx, test_idx = train_test_split(list(range(len(ids))))
            return datasets.DatasetDict({
                "train": datasets.Dataset.from_dict({
                    "id": np.array(ids)[train_idx],
                    "chains": np.array(chains)[train_idx],
                    "input_ids": tokenized["indices"][train_idx],
                    "attention_mask": tokenized["mask"][train_idx],
                }),
                "test": datasets.Dataset.from_dict({
                    "id": np.array(ids)[test_idx],
                    "chains": np.array(chains)[test_idx],
                    "input_ids": tokenized["indices"][test_idx],
                    "attention_mask": tokenized["mask"][test_idx],
                })
            })
        else: 
            return datasets.DatasetDict({
                "train": datasets.Dataset.from_dict({
                    "id": np.array(ids),
                    "chains": np.array(chains),
                    "input_ids": tokenized["indices"],
                    "attention_mask": tokenized["mask"],
                }),
            })

class TRabTokenizerForVJCDR3(AminoAcidTokenizer):
    """Tokenizer for TRA and TRB sequence. Encode V,J genes into tokens, and CDR3 into amino acids"""
    def __init__(self,
        *,
        tra_max_length:int,
        trb_max_length:int,
        pad_token: Optional[str] = None,
        unk_token: Optional[str] = None,
        mask_token: Optional[str] = None,
        cls_token: Optional[str] = None,
        sep_token: Optional[str] = None,
        species: Literal['human', 'mouse'] = 'human',
        **kwargs
    ) -> None:
        # A special token representing an out-of-vocabulary token.
        kwargs["pad_token"] = pad_token or _AMINO_ACIDS_ADDITIONALS["PAD"]
        # A special token representing the class of the input
        kwargs["unk_token"] = unk_token or _AMINO_ACIDS_ADDITIONALS["UNK"]
        # A special token representing a masked token (used by masked-language modeling pretraining objectives, like BERT). 
        kwargs["mask_token"] = mask_token or _AMINO_ACIDS_ADDITIONALS["MASK"]

        kwargs["cls_token"] = cls_token or _AMINO_ACIDS_ADDITIONALS["CLS"]
        kwargs["sep_token"] = sep_token or _AMINO_ACIDS_ADDITIONALS["SEP"]

        super(TRabTokenizerForVJCDR3, self).__init__(model_max_length = tra_max_length + trb_max_length, **kwargs)
        self.tra_max_length = tra_max_length
        self.trb_max_length = trb_max_length
        if species == 'human':
            self.VJ_GENES2INDEX = HumanTCRAnnotations.VJ_GENES2INDEX
            self.VJ_GENES2INDEX_REVERSE = HumanTCRAnnotations.VJ_GENES2INDEX_REVERSE
        elif species == 'mouse':
            self.VJ_GENES2INDEX = MouseTCRAnnotations.VJ_GENES2INDEX
            self.VJ_GENES2INDEX_REVERSE = MouseTCRAnnotations.VJ_GENES2INDEX_REVERSE


    def _encode(self, aa: str, v_gene: str = None, j_gene: str = None, max_length: int = None) -> torch.Tensor:
        aa = list(aa)
        if self.append_cls_token:
            if aa[0] != _AMINO_ACIDS_ADDITIONALS["CLS"]:
                aa = [_AMINO_ACIDS_ADDITIONALS["CLS"]] + aa
            if aa[-1] != _AMINO_ACIDS_ADDITIONALS["SEP"]:
                aa += [_AMINO_ACIDS_ADDITIONALS["SEP"]]
        max_length = max_length or self.model_max_length
        if v_gene and j_gene:
            aa = [v_gene] + aa + [j_gene]
        if len(aa) < max_length:
            aa += list(self._pad_token * (max_length - len(aa)))
        return torch.Tensor(list(map(lambda a: _AMINO_ACIDS_INDEX.get(a) if a in _AMINO_ACIDS_INDEX.keys() else self.VJ_GENES2INDEX.get(a, 0) + len(_AMINO_ACIDS_INDEX), aa)))

    def convert_tokens_to_ids(self, sequence: Union[List[Tuple[str]], Tuple[str]], alpha_vj: Optional[Union[List[Tuple[str]], Tuple[str]]] = None, beta_vj: Optional[Union[List[Tuple[str]], Tuple[str]]] = None):
        # sourcery skip: none-compare, swap-if-else-branches
        if isinstance(sequence, list):
            if alpha_vj != None and beta_vj != None:
                ids = torch.hstack([
                    torch.vstack(
                        list(map(lambda x: self._encode(x[0], x[1][0], x[1][1], max_length=self.tra_max_length), 
                        zip(get_k_elements(sequence, 0), alpha_vj)))
                    ), 
                    torch.vstack(
                        list(map(lambda x: self._encode(x[0], x[1][0], x[1][1], max_length=self.trb_max_length), 
                        zip(get_k_elements(sequence, 1), beta_vj)))
                    )]).type(torch.LongTensor) 
                    
            else:
                ids = torch.hstack([
                    torch.vstack(list(map(lambda x: self._encode(x[0], max_length=self.tra_max_length), sequence))), 
                    torch.vstack(list(map(lambda x: self._encode(x[1], max_length=self.trb_max_length), sequence)))
                ]).type(torch.LongTensor)

        elif alpha_vj != None and beta_vj != None:
            ids = rearrange(torch.hstack([self._encode(sequence[0], alpha_vj[0], alpha_vj[1], max_length=self.tra_max_length), self._encode(sequence[1], beta_vj[0], beta_vj[1], max_length=self.tra_max_length)]), '(n h) -> n h', n=1).type(torch.LongTensor)

        else:
            ids = rearrange(torch.hstack([self._encode(sequence[0], max_length=self.tra_max_length), self._encode(sequence[1], max_length=self.tra_max_length)]), '(n h) -> n h', n=1).type(torch.LongTensor)

        mask = self.convert_ids_to_mask(ids)
        token_type_ids = torch.hstack([torch.zeros(ids.shape[0], self.tra_max_length), torch.ones(ids.shape[0], self.trb_max_length)]).type(torch.LongTensor)

        return {"indices": ids, "mask": mask, "token_type_ids": token_type_ids}

    def _decode(self, ids):
        return list(map(lambda t: _AMINO_ACIDS_INDEX_REVERSE[t] if t in _AMINO_ACIDS_INDEX_REVERSE.keys() else self.VJ_GENES2INDEX_REVERSE[t - len(_AMINO_ACIDS_INDEX)], ids))

    def _trab_decode(self, ids):
        dec = self._decode(ids)
        return (self._unpad(dec[:self.tra_max_length]), self._unpad(dec[self.tra_max_length:]))

    def convert_ids_to_tokens(self, ids: torch.Tensor):
        ids = ids.detach().cpu().numpy().astype(np.int64)
        if len(ids.shape) == 1:
            return self._trab_decode(ids)
        else:
            return list(map(lambda x: self._trab_decode(x), ids))

    def to_dataset(
        self,
        df: pd.DataFrame = None, 
        ids: Iterable[str] = None, 
        alpha_chains: Iterable[str] = None, 
        beta_chains: Iterable[str] = None, 
        alpha_v_genes: Iterable[str] = None,
        alpha_j_genes: Iterable[str] = None,
        beta_v_genes: Iterable[str] = None,
        beta_j_genes: Iterable[str] = None,
        pairing: Iterable[int] = None,
        split: bool = False
    ):
        if df is not None:
            ids = df['id']
            alpha_chains = df['CDR3a']
            beta_chains = df['CDR3b']
            if 'TRAV' in df.columns and 'TRAJ' in df.columns and 'TRBV' in df.columns and 'TRBJ' in df.columns:
                alpha_v_genes = df['TRAV']
                alpha_j_genes = df['TRAJ']
                beta_v_genes = df['TRBV']
                beta_j_genes = df['TRBJ']
            if 'pairing' in df.columns:
                pairing = df['pairing']
            return self._to_dataset(
                ids, 
                alpha_chains, 
                beta_chains, 
                alpha_v_genes, 
                alpha_j_genes, 
                beta_v_genes, 
                beta_j_genes, 
                pairing, 
                split
            )
        else:
            return self._to_dataset(
                ids, 
                alpha_chains, 
                beta_chains, 
                alpha_v_genes, 
                alpha_j_genes, 
                beta_v_genes, 
                beta_j_genes, 
                pairing, 
                split
            )
        
    def _to_dataset(
            self, 
            ids: Iterable[str],
            alpha_chains: Iterable[str], 
            beta_chains: Iterable[str], 
            alpha_v_genes: Iterable[str] = None,
            alpha_j_genes: Iterable[str] = None,
            beta_v_genes: Iterable[str] = None,
            beta_j_genes: Iterable[str] = None,
            pairing: Iterable[int] = None,
            split: bool = False
        ) -> datasets.DatasetDict:
        if not len(ids) == len(alpha_chains) == len(beta_chains):
            raise ValueError("Length of ids(%d), alpha_chains(%d) and beta_chains(%d) do not match" % (len(ids), len(alpha_chains), len(beta_chains)))
        alpha_chains = list(alpha_chains)
        beta_chains = list(beta_chains)
        if all(map(lambda x: x is not None, [alpha_v_genes, alpha_j_genes, beta_v_genes, beta_j_genes])):
            tokenized = self.convert_tokens_to_ids(
                list(zip(alpha_chains, beta_chains)), 
                list(zip(alpha_v_genes, alpha_j_genes)),
                list(zip(beta_v_genes, beta_j_genes)),
            )
        else:
            tokenized = self.convert_tokens_to_ids(list(zip(alpha_chains, beta_chains)))
        if split:
            train_idx, test_idx = train_test_split(list(range(len(ids))))
            return datasets.DatasetDict({
                "train": datasets.Dataset.from_dict({
                    "id": np.array(ids)[train_idx],
                    "alpha_chains": np.array(alpha_chains)[train_idx],
                    "beta_chains": np.array(beta_chains)[train_idx],
                    "input_ids": tokenized["indices"][train_idx],
                    "token_type_ids": tokenized["token_type_ids"][train_idx],
                    "attention_mask": tokenized["mask"][train_idx],
                    "pairing": np.array(pairing)[train_idx] if pairing is not None else np.ones(len(train_idx), dtype=np.uint8)
                }),
                "test": datasets.Dataset.from_dict({
                    "id": np.array(ids)[test_idx],
                    "alpha_chains": np.array(alpha_chains)[test_idx],
                    "beta_chains": np.array(beta_chains)[test_idx],
                    "input_ids": tokenized["indices"][test_idx],
                    "token_type_ids": tokenized["token_type_ids"][test_idx],
                    "attention_mask": tokenized["mask"][test_idx],
                    "pairing": np.array(pairing)[test_idx] if pairing is not None else np.ones(len(test_idx), dtype=np.uint8)
                })
            })
        else:
            return datasets.DatasetDict({
                "train": datasets.Dataset.from_dict({
                    "id": ids,
                    "alpha_chains": alpha_chains,
                    "beta_chains": beta_chains,
                    "input_ids": tokenized["indices"],
                    "token_type_ids": tokenized["token_type_ids"],
                    "attention_mask": tokenized["mask"],
                    "pairing": np.array(pairing) if pairing is not None else np.ones(len(alpha_chains), dtype=np.uint8)
                })
            })
    
    def __call__(self):
        raise NotImplementedError

class TCRabTokenizerForCDR123(AminoAcidTokenizer):
    @deprecated(ymd=(2024, 6, 1))
    def __init__(
        self,
        tra_max_length:int,
        trb_max_length:int,
        pad_token: Optional[str] = None,
        unk_token: Optional[str] = None,
        mask_token: Optional[str] = None,
        cls_token: Optional[str] = None,
        sep_token: Optional[str] = None,
        **kwargs
    ) -> None:
        # A special token representing an out-of-vocabulary token.
        kwargs["pad_token"] = pad_token or _AMINO_ACIDS_ADDITIONALS["PAD"]
        # A special token representing the class of the input
        kwargs["unk_token"] = unk_token or _AMINO_ACIDS_ADDITIONALS["UNK"]
        # A special token representing a masked token (used by masked-language modeling pretraining objectives, like BERT). 
        kwargs["mask_token"] = mask_token or _AMINO_ACIDS_ADDITIONALS["MASK"]

        kwargs["cls_token"] = cls_token or _AMINO_ACIDS_ADDITIONALS["CLS"]
        kwargs["sep_token"] = sep_token or _AMINO_ACIDS_ADDITIONALS["SEP"]

        super(TCRabTokenizerForCDR123, self).__init__(model_max_length = tra_max_length + trb_max_length, **kwargs)
        self.tra_max_length = tra_max_length
        self.trb_max_length = trb_max_length

    def _encode(self, aa: str, v_gene: str = None, j_gene: str = None, max_length: int = None) -> torch.Tensor:
        cdr1 = TRAV2CDR1a[v_gene] if v_gene in TRAV2CDR1a.keys() else TRBV2CDR1b[v_gene]
        cdr2 = TRAV2CDR2a[v_gene] if v_gene in TRAV2CDR2a else TRBV2CDR2b[v_gene] 
        aa = cdr1 + '.' + cdr2 + '.' + aa
        aa = list(aa)

        if self.append_cls_token:
            if aa[0] != _AMINO_ACIDS_ADDITIONALS["CLS"]:
                aa = [_AMINO_ACIDS_ADDITIONALS["CLS"]] + aa
            if aa[-1] != _AMINO_ACIDS_ADDITIONALS["SEP"]:
                aa += [_AMINO_ACIDS_ADDITIONALS["SEP"]]
        max_length = max_length or self.model_max_length
        if len(aa) < max_length:
            aa += list(self._pad_token * (max_length - len(aa)))
        return torch.Tensor(list(map(lambda a: _AMINO_ACIDS_INDEX.get(a), aa)))
    

    def convert_tokens_to_ids(self, sequence: Union[List[Tuple[str]], Tuple[str]], alpha_vj: Optional[Union[List[Tuple[str]], Tuple[str]]] = None, beta_vj: Optional[Union[List[Tuple[str]], Tuple[str]]] = None):
        # sourcery skip: none-compare, swap-if-else-branches
        if isinstance(sequence, list):
            if alpha_vj != None and beta_vj != None:
                ids = torch.hstack([
                    torch.vstack(
                        list(map(lambda x: self._encode(x[0], x[1][0], x[1][1], max_length=self.tra_max_length), 
                        zip(get_k_elements(sequence, 0), alpha_vj)))
                    ), 
                    torch.vstack(
                        list(map(lambda x: self._encode(x[0], x[1][0], x[1][1], max_length=self.trb_max_length), 
                        zip(get_k_elements(sequence, 1), beta_vj)))
                    )]).type(torch.LongTensor) 
                    
            else:
                ids = torch.hstack([
                    torch.vstack(list(map(lambda x: self._encode(x[0], max_length=self.tra_max_length), sequence))), 
                    torch.vstack(list(map(lambda x: self._encode(x[1], max_length=self.trb_max_length), sequence)))
                ]).type(torch.LongTensor)

        elif alpha_vj != None and beta_vj != None:
            ids = rearrange(torch.hstack([self._encode(sequence[0], alpha_vj[0], alpha_vj[1], max_length=self.tra_max_length), self._encode(sequence[1], beta_vj[0], beta_vj[1], max_length=self.tra_max_length)]), '(n h) -> n h', n=1).type(torch.LongTensor)

        else:
            ids = rearrange(torch.hstack([self._encode(sequence[0], max_length=self.tra_max_length), self._encode(sequence[1], max_length=self.tra_max_length)]), '(n h) -> n h', n=1).type(torch.LongTensor)

        mask = self.convert_ids_to_mask(ids)
        token_type_ids = torch.hstack([torch.zeros(ids.shape[0], self.tra_max_length), torch.ones(ids.shape[0], self.trb_max_length)]).type(torch.LongTensor)

        return {"indices": ids, "mask": mask, "token_type_ids": token_type_ids}


    def _decode(self, ids):
        return list(map(lambda t: _AMINO_ACIDS_INDEX_REVERSE[t] if t in _AMINO_ACIDS_INDEX_REVERSE.keys() else self.VJ_GENES2INDEX_REVERSE[t - len(_AMINO_ACIDS_INDEX)], ids))
    
    def convert_ids_to_tokens(self, ids) -> Iterable[str]:
        return super().convert_ids_to_tokens(ids)
    
    def convert_ids_to_mask(self, ids: torch.Tensor) -> torch.Tensor:
        return super().convert_ids_to_mask(ids)
    
    def to_dataset(
        self,
        df: pd.DataFrame = None, 
        ids: Iterable[str] = None, 
        alpha_chains: Iterable[str] = None, 
        beta_chains: Iterable[str] = None, 
        alpha_v_genes: Iterable[str] = None,
        alpha_j_genes: Iterable[str] = None,
        beta_v_genes: Iterable[str] = None,
        beta_j_genes: Iterable[str] = None,
        pairing: Iterable[int] = None,
        split: bool = False
    ):
        if df is not None:
            ids = df['id']
            alpha_chains = df['CDR3a']
            beta_chains = df['CDR3b']
            if 'TRAV' in df.columns and 'TRAJ' in df.columns and 'TRBV' in df.columns and 'TRBJ' in df.columns:
                alpha_v_genes = df['TRAV']
                alpha_j_genes = df['TRAJ']
                beta_v_genes = df['TRBV']
                beta_j_genes = df['TRBJ']
            if 'pairing' in df.columns:
                pairing = df['pairing']
            return self._to_dataset(
                ids, 
                alpha_chains, 
                beta_chains, 
                alpha_v_genes, 
                alpha_j_genes, 
                beta_v_genes, 
                beta_j_genes, 
                pairing, 
                split
            )
        else:
            return self._to_dataset(
                ids, 
                alpha_chains, 
                beta_chains, 
                alpha_v_genes, 
                alpha_j_genes, 
                beta_v_genes, 
                beta_j_genes, 
                pairing, 
                split
            )
        
    def _to_dataset(
        self, 
        ids: Iterable[str],
        alpha_chains: Iterable[str], 
        beta_chains: Iterable[str], 
        alpha_v_genes: Iterable[str] = None,
        alpha_j_genes: Iterable[str] = None,
        beta_v_genes: Iterable[str] = None,
        beta_j_genes: Iterable[str] = None,
        pairing: Iterable[int] = None,
        split: bool = False
    ) -> datasets.DatasetDict:
        if not len(ids) == len(alpha_chains) == len(beta_chains):
            raise ValueError("Length of ids(%d), alpha_chains(%d) and beta_chains(%d) do not match" % (len(ids), len(alpha_chains), len(beta_chains)))
        alpha_chains = list(alpha_chains)
        beta_chains = list(beta_chains)
        if all(map(lambda x: x is not None, [alpha_v_genes, alpha_j_genes, beta_v_genes, beta_j_genes])):
            tokenized = self.convert_tokens_to_ids(
                list(zip(alpha_chains, beta_chains)), 
                list(zip(alpha_v_genes, alpha_j_genes)),
                list(zip(beta_v_genes, beta_j_genes)),
            )
        else:
            tokenized = self.convert_tokens_to_ids(list(zip(alpha_chains, beta_chains)))

        if split:
            train_idx, test_idx = train_test_split(list(range(len(ids))))
            return datasets.DatasetDict({
                "train": datasets.Dataset.from_dict({
                    "id": np.array(ids)[train_idx],
                    "alpha_chains": np.array(alpha_chains)[train_idx],
                    "beta_chains": np.array(beta_chains)[train_idx],
                    "input_ids": tokenized["indices"][train_idx],
                    "token_type_ids": tokenized["token_type_ids"][train_idx],
                    "attention_mask": tokenized["mask"][train_idx],
                    "pairing": np.array(pairing)[train_idx] if pairing is not None else np.ones(len(train_idx), dtype=np.uint8)
                }),
                "test": datasets.Dataset.from_dict({
                    "id": np.array(ids)[test_idx],
                    "alpha_chains": np.array(alpha_chains)[test_idx],
                    "beta_chains": np.array(beta_chains)[test_idx],
                    "input_ids": tokenized["indices"][test_idx],
                    "token_type_ids": tokenized["token_type_ids"][test_idx],
                    "attention_mask": tokenized["mask"][test_idx],
                    "pairing": np.array(pairing)[test_idx] if pairing is not None else np.ones(len(test_idx), dtype=np.uint8)
                })
            })
        else:
            return datasets.DatasetDict({
                "train": datasets.Dataset.from_dict({
                    "id": ids,
                    "alpha_chains": alpha_chains,
                    "beta_chains": beta_chains,
                    "input_ids": tokenized["indices"],
                    "token_type_ids": tokenized["token_type_ids"],
                    "attention_mask": tokenized["mask"],
                    "pairing": np.array(pairing) if pairing is not None else np.ones(len(alpha_chains), dtype=np.uint8)
                })
            })

TRabTokenizerForPseudoSequence = Tokenizer

trab_tokenizer_for_pseudosequence = TRabTokenizerForPseudoSequence.from_str(
    open(
        MODULE_PATH / "tokenizer.json"
    ).read()
)
trab_tokenizer_for_pseudosequence.pad_token_id = 0
trab_tokenizer_for_pseudosequence.pad_token = "<|bos|>"
trab_tokenizer_for_pseudosequence.eos_token_id = 1
trab_tokenizer_for_pseudosequence.eos_token = "<|eos|>"
trab_tokenizer_for_pseudosequence.bos_token_id = 2
trab_tokenizer_for_pseudosequence.bos_token = "<|bos|>"
trab_tokenizer_for_pseudosequence.bos_token_id = 3

default_max_length = {
    'CDR3a': 36,
    'CDR3b': 36,
    'CDR1b': 8,
    'CDR2b': 8,
    'CDR1a': 8,
    'CDR2a': 8,
}

def tcr_pseudo_sequence_weight(
    cdr1a: float = 1.,
    cdr2a: float = 1.,
    cdr3a: float = 3.,
    cdr1b: float = 1.,
    cdr2b: float = 1.,
    cdr3b: float = 5.,
):
    weight = torch.zeros(110)
    weight[1:9] = cdr1a
    weight[10:18] = cdr2a
    weight[19:55] = cdr3a
    weight[56:64] = cdr1b
    weight[65:73] = cdr2b
    weight[74:110] = cdr3b
    return weight 

def tokenize_tcr_pseudo_sequence_to_fixed_length(
    pseudo_sequence, 
    max_length = default_max_length, 
    pad_token = '.',
    show_progress=False
):
    _f = False
    if isinstance(pseudo_sequence, str):
        pseudo_sequences = [pseudo_sequence]
        _f = True
    elif isinstance(pseudo_sequence, Iterable):
        pseudo_sequences = pseudo_sequence
    tokens, ids, attention_mask = [], [], []
    if show_progress:
        pbar = tqdm.tqdm(total=len(pseudo_sequences))
    
    for pseudo_sequence in pseudo_sequences:
        cdr1b, cdr2b, cdr3b, cdr1a, cdr2a, cdr3a = pseudo_sequence.split(':')
        cdr1a = cdr1a + pad_token * (max_length['CDR1a'] - len(cdr1a))
        cdr2a = cdr2a + pad_token * (max_length['CDR2a'] - len(cdr2a))
        cdr1b = cdr1b + pad_token * (max_length['CDR1b'] - len(cdr1b))
        cdr2b = cdr2b + pad_token * (max_length['CDR2b'] - len(cdr2b))
        cdr3a = cdr3a + pad_token * (max_length['CDR3a'] - len(cdr3a))
        cdr3b = cdr3b + pad_token * (max_length['CDR3b'] - len(cdr3b))
        tokens.append('^' + ':'.join([cdr1a, cdr2a, cdr3a, cdr1b, cdr2b, cdr3b]))
        encoding = trab_tokenizer_for_pseudosequence.encode(tokens[-1])
        ids.append(np.array(encoding.ids))
        attention_mask.append( ids[-1] >= 6 )
        if show_progress:
            pbar.update(1)
    if show_progress:
        pbar.close()
    if _f:
        return tokens[0], ids[0], attention_mask[0]
    else:
        return tokens, ids, attention_mask
    
def tokenize_to_fixed_length(
    sequence,
    max_length,
    pad_token = '.',
):
    _f = False
    if isinstance(sequence, str):
        sequences = [sequence]
        _f = True
    elif isinstance(sequence, Iterable):
        sequences = sequence

    tokens, ids, attention_mask = [], [], []

    
    for sequence in sequences:
        sequence = sequence + pad_token * (max_length - len(sequence))
        tokens.append('^' + sequence)
        encoding = trab_tokenizer_for_pseudosequence.encode(tokens[-1])
        ids.append(np.array(encoding.ids))
        attention_mask.append( ids[-1] >= 6 )

    if _f:
        return tokens[0], ids[0], attention_mask[0]
    else:
        return tokens, ids, attention_mask
    
