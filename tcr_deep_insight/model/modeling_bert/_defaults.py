from typing import Literal, Union
import torch 
from . import TRabModelingBertForPseudoSequence, TRabModelingBertForVJCDR3
from ..tokenizers import (
    TRabTokenizerForVJCDR3,
    TRabTokenizerForPseudoSequence
)
from ._collator import (
    AminoAcidsCollator,
    TRabCollatorForVJCDR3,
)

def default_optimizer(model, *args, **kwargs):
    optimizer = torch.optim.AdamW(model.parameters(), *args, **kwargs), 
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.1, patience=5)
    return (optimizer, scheduler)

def default_collator(
    model: Union[TRabModelingBertForVJCDR3, TRabModelingBertForPseudoSequence],
    species:Literal['human','mouse'],
    tra_max_length:int=48,
    trb_max_length:int=48,
):
    if isinstance(model, TRabModelingBertForPseudoSequence):
        return AminoAcidsCollator(
            mask_token_id=4,
            max_length=110,
            mlm_probability=0.15,
        )
    elif isinstance(model, TRabModelingBertForVJCDR3):
        return TRabCollatorForVJCDR3(
            tra_max_length=tra_max_length, 
            trb_max_length=trb_max_length,
            species=species
        )
    else:
        raise ValueError("model must be an instance of TRabModelingBertForPseudoSequence or TRabModelingBertForVJCDR3")

