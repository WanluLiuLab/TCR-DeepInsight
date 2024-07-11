import random
from typing import Literal, Tuple
import torch
import numpy as np

from ...utils._amino_acids import (_AMINO_ACIDS, _AMINO_ACIDS_ADDITIONALS,
                                  _AMINO_ACIDS_INDEX_REVERSE,
                                  _AMINO_ACIDS_INDEX)

from ...utils._tcr import HumanTCRAnnotations, MouseTCRAnnotations


class AminoAcidsCollator:
    def __init__(
        self, 
        mask_token_id: int, 
        max_length: int, 
        mlm_probability: float = 0.2
    ) -> None:
        self.mask_token_id = mask_token_id
        self.max_length = max_length
        self.mlm_probability = mlm_probability

    def _rand_mask(self, input_ids, attention_mask) -> torch.Tensor:
        ids = input_ids.copy()
        mask = np.array(attention_mask)
        mask_idx = list(
            map(lambda z: z[0], list(filter(lambda x: x[1], enumerate(attention_mask))))
        )
        n_mask = int(len(mask_idx) * self.mlm_probability)
        mask_idx = np.random.choice(mask_idx, n_mask)
        ids[mask_idx] = self.mask_token_id
        mask[mask_idx] = False
        return torch.tensor(ids.astype(np.int64)), torch.tensor(mask.astype(np.bool))

    def __call__(self, input_ids, attention_mask) -> torch.Tensor:
        if len(input_ids.shape) == 1:
            return self._rand_mask(input_ids, attention_mask)
        else:
            ret = [
                self._rand_mask(i, a)
                for i, a in zip(input_ids, attention_mask)
            ]
            return torch.vstack(list(map(lambda x: x[0], ret))), torch.vstack(
                list(map(lambda x: x[1], ret)))


class TRabCollatorForVJCDR3:
    def __init__(
        self,
        tra_max_length: int,
        trb_max_length: int,
        mask_token_id: int = _AMINO_ACIDS_INDEX[
            _AMINO_ACIDS_ADDITIONALS["MASK"]],
        mlm_probability: float = 0.1,
        mask_trb_probability: float = 0.5,
        species: Literal["human", "mouse"] = "human",
    ) -> None:
        self.mask_token_id = mask_token_id
        self.tra_max_length = tra_max_length
        self.trb_max_length = trb_max_length
        self.mlm_probability = mlm_probability
        self.mask_trb_probability = mask_trb_probability
        if species == "human":
            self.VJ_GENES2INDEX = HumanTCRAnnotations.VJ_GENES2INDEX
        elif species == "mouse":
            self.VJ_GENES2INDEX = MouseTCRAnnotations.VJ_GENES2INDEX
        else:
            raise ValueError("Species must be either human or mouse")

    def _rand_mask(self, input_ids, attention_mask, mr_mask = None) -> Tuple[torch.Tensor]:
        mask_tr = random.random()
        if mask_tr < self.mask_trb_probability:
            ids = np.array(input_ids[self.tra_max_length:])
            mask = np.array(attention_mask[self.tra_max_length:])
            if mr_mask is not None:
                mr = np.array(mr_mask[self.tra_max_length:])
                if ids[0]-27 in self.VJ_GENES2INDEX.values():
                    mr[0] = True 
            else:
                mr = np.ones(mask.shape)
            mask_idx = list(
                map(
                    lambda z: z[0],
                    list(
                        filter(
                            lambda x:
                            ((_AMINO_ACIDS_INDEX_REVERSE.get(x[1][0], None) in
                              _AMINO_ACIDS) or x[1][0] > len(
                                  _AMINO_ACIDS_INDEX_REVERSE)) and x[1][1] and x[1][2],
                            enumerate(zip(ids, mask, mr))
                        )
                    )
                )
            )

            if not mask_idx:
                return torch.tensor(input_ids,
                                    dtype=torch.int64), torch.tensor(
                                        attention_mask, dtype=torch.bool)

            n_mask = int(len(mask_idx) * self.mlm_probability)
            _mask_idx = np.random.choice(mask_idx, n_mask, replace=False)
            
            # Mask V gene
            if ids[mask_idx[0]]-27 in self.VJ_GENES2INDEX.values(
            ) and mask_idx[0] not in _mask_idx:
                _mask_idx = np.hstack([mask_idx[0], _mask_idx])

            

            ids[_mask_idx] = self.mask_token_id
            mask[_mask_idx] = False

            return torch.tensor(
                np.hstack([input_ids[:self.tra_max_length], ids]),
                dtype=torch.int64), torch.tensor(np.hstack(
                    [attention_mask[:self.tra_max_length], mask]),
                                                 dtype=torch.bool)

        else:
            ids = np.array(input_ids[:self.tra_max_length])
            mask = np.array(attention_mask[:self.trb_max_length])
            if mr_mask is not None:
                mr = np.array(mr_mask[:self.trb_max_length])
                if ids[0]-27 in self.VJ_GENES2INDEX.values():
                    mr[0] = True 
            else: 
                mr = np.ones(mask.shape)
            mask_idx = list(
                map(
                    lambda z: z[0],
                    list(
                        filter(
                            lambda x:
                            ((_AMINO_ACIDS_INDEX_REVERSE.get(x[1][0], None) in
                              _AMINO_ACIDS) or x[1][0] > len(
                                  _AMINO_ACIDS_INDEX_REVERSE)) and x[1][1] and x[1][2],
                            enumerate(zip(ids, mask, mr))
                        )
                    )
                )
            )

            if not mask_idx:
                return torch.tensor(input_ids,
                                    dtype=torch.int64), torch.tensor(
                                        attention_mask, dtype=torch.bool)

            n_mask = int(len(mask_idx) * self.mlm_probability)
            _mask_idx = np.random.choice(mask_idx, n_mask, replace=False)

            # Mask V gene
            if ids[mask_idx[0]]-27 in self.VJ_GENES2INDEX.values(
            ) and mask_idx[0] not in _mask_idx:
                _mask_idx = np.hstack([mask_idx[0], _mask_idx])
            
            ids[_mask_idx] = self.mask_token_id
            mask[_mask_idx] = False
            return torch.tensor(
                np.hstack([ids, input_ids[self.tra_max_length:]]),
                dtype=torch.int64), torch.tensor(np.hstack(
                    [mask, attention_mask[self.tra_max_length:]]),
                                                 dtype=torch.bool)

    def __call__(self, input_ids, attention_mask, mr_mask) -> torch.Tensor:
        if len(input_ids.shape) == 1:
            return self._rand_mask(input_ids, attention_mask, mr_mask)
        else:
            ret = []
            for i, a, m in zip(input_ids, attention_mask, mr_mask):
                ret.append(self._rand_mask(i, a, m))
            return torch.vstack(list(map(lambda x: x[0], ret))), torch.vstack(
                list(map(lambda x: x[1], ret)))


class TRABMutator:
    def __init__(self,
                 tra_max_length: int,
                 trb_max_length: int,
                 max_mutation_aa: int = 2,
                 mutate_trb_probability: float = 1,
                 substitution_probability: float = 0.5,
                 insertion_probability: float = 0.2,
                 max_insertion_aa: int = 1,
                 deletion_probability: float = 0.2,
                 max_deletion_aa: int = 1,
                 is_full_length: bool = False) -> None:
        self.tra_max_length = tra_max_length
        self.trb_max_length = trb_max_length
        self.max_mutation_aa = max_mutation_aa
        self.mutate_trb_probability = mutate_trb_probability
        self.substitution_probability = substitution_probability
        self.insertion_probability = insertion_probability
        self.max_insertion_aa = max_deletion_aa
        self.deletion_probability = deletion_probability
        self.max_deletion_aa = max_deletion_aa
        self.is_full_length = False

    def _mutate(self, input_ids, attention_mask) -> Tuple[torch.Tensor]:
        mask_tr = random.random()
        if mask_tr < self.mutate_trb_probability:
            # mutate TRB
            ids = np.array(input_ids[self.tra_max_length:])
            mask = np.array(attention_mask[self.tra_max_length:])
            seq_len = max(
                map(lambda x: x[0], filter(lambda z: z[1], enumerate(mask))))

            mutation_idx = list(
                map(
                    lambda z: z[0],
                    list(
                        filter(
                            lambda x: (_AMINO_ACIDS_INDEX_REVERSE.get(
                                x[1][0], None) in _AMINO_ACIDS
                                       if x[1][0] in _AMINO_ACIDS_INDEX_REVERSE
                                       .keys() else False) and x[1][1],
                            enumerate(zip(ids, mask))
                        )
                    )
                )
            )
            if not self.is_full_length:
                # We focus on the middle region of the CDR3 region
                mutation_idx = list(
                    filter(lambda x: x > 3 and x < seq_len - 3, mutation_idx))

            if random.random() < self.substitution_probability:
                mutation_idx_substitution = np.random.choice(
                    mutation_idx, min(len(mutation_idx), self.max_mutation_aa))
                for i in mutation_idx_substitution:
                    ids[i] = _AMINO_ACIDS_INDEX[np.random.choice(
                        _AMINO_ACIDS, 1)[0]]
            else:
                if random.random() < self.insertion_probability:
                    mutation_idx_insertion = np.random.choice(
                        mutation_idx,
                        min(len(mutation_idx), self.max_insertion_aa))
                    for i in mutation_idx_insertion:
                        ids = np.insert(
                            ids, i, _AMINO_ACIDS_INDEX[np.random.choice(
                                _AMINO_ACIDS, 1)[0]])[:-1]

                if random.random() < self.deletion_probability:
                    mutation_idx_deletion = np.random.choice(
                        mutation_idx,
                        min(len(mutation_idx), self.max_deletion_aa))
                    for i in mutation_idx_deletion:
                        ids = np.hstack([
                            np.delete(ids, i),
                            np.array([
                                _AMINO_ACIDS_INDEX[
                                    _AMINO_ACIDS_ADDITIONALS["PAD"]]
                            ])
                        ])

            return torch.tensor(
                np.hstack([input_ids[:self.tra_max_length], ids]),
                dtype=torch.int64), torch.tensor(attention_mask,
                                                 dtype=torch.bool)
        else:
            # mutate TRA
            ids = np.array(input_ids[:self.tra_max_length])
            mask = np.array(attention_mask[:self.tra_max_length])
            seq_len = max(
                map(lambda x: x[0], filter(lambda z: z[1], enumerate(mask))))
            mutation_idx = list(
                map(
                    lambda z: z[0],
                    list(
                        filter(
                            lambda x: (_AMINO_ACIDS_INDEX_REVERSE.get(
                                x[1][0], None) in _AMINO_ACIDS
                                       if x[1][0] in _AMINO_ACIDS_INDEX_REVERSE
                                       .keys() else False) and x[1][1],
                            enumerate(zip(ids, mask))
                        )
                    )
                )
            )
            if not self.is_full_length:
                mutation_idx = list(
                    filter(lambda x: x > 2 and x < seq_len - 3, mutation_idx))
            mutation_idx = np.random.choice(
                mutation_idx, min(len(mutation_idx), self.max_mutation_aa))
            for i in mutation_idx:
                ids[i] = _AMINO_ACIDS_INDEX[np.random.choice(_AMINO_ACIDS,
                                                             1)[0]]
            return torch.tensor(
                np.hstack([ids, input_ids[self.tra_max_length:]]),
                dtype=torch.int64), torch.tensor(attention_mask,
                                                 dtype=torch.bool)

    def __call__(self, input_ids, attention_mask) -> torch.Tensor:
        if len(input_ids.shape) == 1:
            return self._mutate(input_ids, attention_mask)
        else:
            ret = []
            for i, a in zip(input_ids, attention_mask):
                ret.append(self._mutate(i, a))
            return torch.vstack(list(map(lambda x: x[0], ret))), torch.vstack(
                list(map(lambda x: x[1], ret)))
