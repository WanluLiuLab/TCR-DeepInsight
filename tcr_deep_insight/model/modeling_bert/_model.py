# Pytorch
import torch
from torch import nn, optim
from torch.nn import functional as F
from torch.utils.data import DataLoader
from torch.utils.data.sampler import SubsetRandomSampler
from torch.distributions import kl_divergence as kld
from torch.optim.lr_scheduler import ReduceLROnPlateau


# Hugginface Transformers
from transformers import (
    BertConfig,
    PreTrainedModel,
    BertForMaskedLM
)

# Third Party Transformers
from .externals.bert_layers import (
    BertForMaskedLM as TritonBertForMaskedLM
)

# Third Party packages
from sklearn.neighbors import KNeighborsClassifier
from joblib import dump, load
import datasets
import numpy as np
from pathlib import Path
import umap
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from sklearn.utils import class_weight

# Built-in
import time
from collections import Counter
from itertools import chain
from copy import deepcopy
import json
from typing import Callable, Mapping, Union, Iterable, Tuple, Optional, Mapping
import os
import warnings

# Friendly dependencites
from scatlasvae.model._primitives import *
from scatlasvae.model._primitives import Linear, FCLayer, SAE
from scatlasvae.utils._loss import LossFunction
from scatlasvae.utils._parallelizer import Parallelizer

# Package
from ..tokenizers._tokenizer import TRabTokenizerForVJCDR3
from ._config import get_config
from ._layers import (
    OuterProductMean,
    TriangleAttentionStartingNode,
    TriangleAttentionEndingNode,
    PairTransition,
    JointEmbedder,
)
from .._model_utils import add, ModelOutput
from ...utils._logger import mt, Colors, get_tqdm
from ...utils._compat import Literal


MODULE_PATH = Path(__file__).parent
warnings.filterwarnings("ignore")


class ModuleBase(nn.Module):
    def to(self, device:str):
        super(ModuleBase, self).to(device)
        self.device=device
        return self
    

class TRabModelingBertForVJCDR3(ModuleBase):
    def __init__(self,
        bert_config: BertConfig,
        pooling: Literal["cls", "mean", "max", "pool", "trb", "tra", "weighted"] = "mean",
        pooling_cls_position: int = 0,
        pooling_weight = (0.1,0.9),
        labels_number: int = 1,
        device = "cuda",
    ) -> None:
        """
        TRABModel is a BERT model that takes in a abTCR sequence
        :param bert_config: BertConfig
        :param pooling: Pooling method, one of "cls", "mean", "max", "pool", "trb", "tra", "weighted"
        :param pooling_cls_position: Position of the cls token
        :param pooling_weight: Weight of the cls token
        :param hidden_layers: Hidden layers of the classifier
        :param labels_number: Number of labels

        :example:
            >>> from t_deep_insight as tdi
            >>> model = tdi.model.TCRModel(
            >>>    tdi.model.config.get_human_config(),
            >>>    labels_number=4
            >>> )
        """
        super(TRabModelingBertForVJCDR3, self).__init__()
        self.model = BertForMaskedLM(bert_config)
        self.pooler = nn.Sequential(
            nn.Linear(bert_config.hidden_size,bert_config.hidden_size),
            nn.Tanh()
        )

        self.pooling = pooling
        self.pooling_cls_position = pooling_cls_position
        self.pooling_weight = pooling_weight

        self.config = bert_config

        self.labels_number = labels_number

        self.classifier = nn.Sequential(
            nn.Dropout(0.1),
            Linear(bert_config.hidden_size, labels_number, init='final')
        )
        self.device = device
        self.to(device)

    def __repr__(self):
        return f'{Colors.GREEN}TRABModel{Colors.NC} object containing:\n' + \
            f'    bert_config: {self.config}\n' + \
            f'    pooling: {self.pooling}\n' + \
            f'    pooling_cls_position: {self.pooling_cls_position}\n' + \
            f'    pooling_weight: {self.pooling_weight}\n' + \
            f'    labels_number: {self.labels_number}\n'

    def forward(self,
        *,
        input_ids: torch.Tensor,
        attention_mask: torch.Tensor,
        labels: torch.Tensor,
        token_type_ids: torch.Tensor,
        output_hidden_states = True,
    ):
        '''
        Forward pass of the model

        :param input_ids: Input ids
        :param attention_mask: Attention mask
        :param labels: Labels
        :param token_type_ids: Token type ids

        :return: Output of the model
        '''

        output = self.model.forward(
                input_ids = input_ids,
                attention_mask = attention_mask,
                return_dict = True,
                token_type_ids = token_type_ids,
                labels = labels,
                output_hidden_states = output_hidden_states
        )
        hidden_states = None
        hidden_states_length = int(output.hidden_states[-1].shape[1]/2)
        if output_hidden_states:
            if self.pooling == "mean":
                hidden_states = output.hidden_states[-1][:,1:,:].mean(1)
            elif self.pooling == "max":
                hidden_states = output.hidden_states[-1][:,1:,:].max(1)
            elif self.pooling == "cls":
                hidden_states = output.hidden_states[-1][:,self.pooling_cls_position,:]
            elif self.pooling == 'pool':
                hidden_states = self.pooler(output.hidden_states[-1][:,self.pooling_cls_position,:])
            elif self.pooling == 'tra':
                if self.pooling_cls_position == 1:
                    hidden_states = output.hidden_states[-1][
                        :,
                        torch.hstack([
                            torch.tensor([self.pooling_cls_position]),
                            torch.arange(2,hidden_states_length),
                        ]),
                        :
                    ].mean(1) + output.hidden_states[-1][:, 0] + output.hidden_states[-1][:, hidden_states_length]
                else:
                    hidden_states = output.hidden_states[-1][
                        :,
                        torch.arange(1,hidden_states_length),
                        :
                    ].mean(1)
            elif self.pooling == 'trb':
                if self.pooling_cls_position == 1:
                    hidden_states = output.hidden_states[-1][
                        :,
                        torch.hstack([
                            torch.tensor([self.pooling_cls_position]),
                            torch.arange(hidden_states_length+2,hidden_states_length*2),
                        ]),
                        :
                    ].mean(1) + output.hidden_states[-1][:, 0] + output.hidden_states[-1][:, hidden_states_length]
                else:
                    hidden_states = output.hidden_states[-1][
                        :,
                        torch.arange(hidden_states_length+1,hidden_states_length*2),
                        :
                    ].mean(1)
            elif self.pooling == 'weighted':
                if self.pooling_cls_position == 1:
                    hidden_states = (output.hidden_states[-1][
                        :,
                        torch.hstack([
                            torch.tensor([self.pooling_cls_position]),
                            torch.arange(2,hidden_states_length)
                        ]),
                        :
                    ] * self.pooling_weight[0] + output.hidden_states[-1][
                        :,
                        torch.hstack([
                            torch.tensor([self.pooling_cls_position]),
                            torch.arange(hidden_states_length+2,hidden_states_length*2),
                        ]),
                        :
                    ] * self.pooling_weight[1]).mean(1) + output.hidden_states[-1][:, 0] + output.hidden_states[-1][:, hidden_states_length]
                else:
                    hidden_states = (output.hidden_states[-1][
                        :,
                        torch.arange(hidden_states_length+1,hidden_states_length*2),
                        :
                    ] * self.pooling_weight[0] + output.hidden_states[-1][
                        :,
                        torch.arange(1,hidden_states_length)
                        :
                    ] * self.pooling_weight[1]).mean(1)
            else:
                raise ValueError("Unrecognized pool strategy")

        prediction_out = self.classifier(hidden_states)

        return {
            "output": output,
            "hidden_states": hidden_states,
            "prediction_out": prediction_out,
        }

class TRabModelingBertForPseudoSequence(nn.Module):
    def __init__(self,
        bert_config: BertConfig,
        pooling: Union[Callable, Literal["cls", "mean", "max", "cdr3a", "cdr3b", "weighted"]] = "mean",
        pooling_weight: Optional[torch.Tensor] = None,
        labels_number: int = 1,
        use_triton: bool = False,
        device = "cuda",
    ) -> None:
        super(TRabModelingBertForPseudoSequence, self).__init__()
        if use_triton:
            self.model = TritonBertForMaskedLM(bert_config)
        else:
            self.model = BertForMaskedLM(bert_config)
        self.classifier = nn.Sequential(
            nn.Dropout(0.1),
            Linear(bert_config.hidden_size, labels_number, init='final')
        )
        self.pooling = pooling
        self.pooling_weight = pooling_weight.to(device) if pooling_weight is not None else None
        self.config = bert_config
        self.device = device
        self.labels_number = labels_number
        self.to(device)

    def forward(self,
        *,
        input_ids: torch.Tensor,
        attention_mask: torch.Tensor,
        labels: torch.Tensor = None,
        output_hidden_states = True,
    ):
        '''
        Forward pass of the model

        :param input_ids: Input ids
        :param attention_mask: Attention mask
        :param labels: Labels
        :param token_type_ids: Token type ids

        :return: Output of the model
        '''

        output = self.model.forward(
            input_ids=input_ids,
            attention_mask=attention_mask,
            return_dict=True,
            labels=labels,
            output_hidden_states=output_hidden_states,
        )

        hidden_states = None
        if output_hidden_states:
            if isinstance(self.pooling, str):
                if self.pooling == "mean":
                    hidden_states = output.hidden_states[-1][:, 1:, :].mean(1)
                elif self.pooling == "max":
                    hidden_states = output.hidden_states[-1][:,1:,:].max(1)
                elif self.pooling == "cls":
                    hidden_states = output.hidden_states[-1][:,0,:]
                elif self.pooling == 'cdr3b':
                    hidden_states = output.hidden_states[-1][:,74:,:]
                elif self.pooling == 'cdr3a':
                    hidden_states = output.hidden_states[-1][:,19:55,:]
                elif self.pooling == 'weighted':
                    hidden_states = (output.hidden_states[-1] * self.pooling_weight.to(input_ids.device)).mean(1)
            elif callable(self.pooling):
                hidden_states = self.pooling(output.hidden_states[-1])
            else:
                raise ValueError("Unrecognized pool strategy")

        prediction_out = self.classifier(hidden_states)

        return {
            "output": output,
            "hidden_states": hidden_states,
            "prediction_out": prediction_out,
        }

class TCRpMHCPairUpdateBlock(nn.Module):
    def __init__(self, bert_config) -> None:
        super(TCRpMHCPairUpdateBlock, self).__init__()

        self.opm = OuterProductMean(
            bert_config.hidden_size,
            bert_config.hidden_size,
            bert_config.hidden_size,
        )

        self.tri_att_start = TriangleAttentionStartingNode(
            bert_config.hidden_size,
            bert_config.hidden_size,
            bert_config.num_attention_heads
        )

        self.tri_att_end = TriangleAttentionEndingNode(
            bert_config.hidden_size,
            bert_config.hidden_size,
            bert_config.num_attention_heads
        )

        self.transition = PairTransition(
            bert_config.hidden_size, 1
        )

    def forward(self,
        z: torch.Tensor,
        z_mask: torch.Tensor,
        tcr_hidden_states: torch.Tensor,
        pmhc_hidden_states: torch.Tensor,
        tcr_attention_mask: torch.Tensor,
        pmhc_attention_mask: torch.Tensor,
        inplace_safe=False,
        output_attentions=False,
    ):
        z = z + self.opm(
            tcr_hidden_states, 
            pmhc_hidden_states, 
            tcr_attention_mask, 
            pmhc_attention_mask
        )
        if output_attentions:
            o, a1 = self.tri_att_start(z, z_mask, output_attentions=output_attentions)
            z = add(z, o, inplace_safe)
            o, a2 = self.tri_att_end(z, z_mask, output_attentions=output_attentions)
            z = add(z, o, inplace_safe)
        else:
            z = add(z, self.tri_att_start(z, z_mask, output_attentions=output_attentions), inplace_safe)
            z = add(z, self.tri_att_end(z, z_mask, output_attentions=output_attentions), inplace_safe)

        z = add(z, self.transition(z, z_mask), inplace_safe)
        if output_attentions:
            return z, a1, a2
        return z

class ContactDistogramHead(nn.Module):
    def __init__(self, c_z, no_bins=64):
        """
        Args:
            c_z:
                Input channel dimension
            no_bins:
                Number of distogram bins
        """
        super(ContactDistogramHead, self).__init__()

        self.c_z = c_z
        self.no_bins = no_bins

        self.linear = Linear(self.c_z, self.no_bins, init="final")

    def _forward(self, z):  # [*, N, N, C_z]
        """
        Args:
            z:
                [*, N_res, N_res, C_z] pair embedding
        Returns:
            [*, N, N, no_bins] distogram probability distribution
        """
        # [*, N, N, no_bins]
        logits = self.linear(z)
        return logits
    
    def forward(self, z): 
        return self._forward(z)


class TCRpMHCPairAttention(nn.Module):
    def __init__(self,
        bert_config: BertConfig,
        use_triton: bool = False,
    ) -> None:
        super(TCRpMHCPairAttention, self).__init__()
        if use_triton:
            self.model = TritonBertForMaskedLM(bert_config)
        else:
            self.model = BertForMaskedLM(bert_config)
        self.embedder = JointEmbedder(
            bert_config.hidden_size,
            bert_config.hidden_size,
            bert_config.hidden_size,
        )
        self.layers = nn.ModuleList([
            TCRpMHCPairUpdateBlock(bert_config)
            for _ in range(bert_config.num_hidden_layers)
        ])
        self.distogram_head = ContactDistogramHead(bert_config.hidden_size)
        self.contact_head = Linear(bert_config.hidden_size, 1, init="final")


    def forward(self,
        *,
        input_ids: torch.Tensor,
        attention_mask: torch.Tensor,
        tcr_hidden_states: torch.Tensor,
        tcr_attention_mask: torch.Tensor,
        output_hidden_states = True,
        output_triangular_attentions = True,
    ):
        output = self.model.forward(
            input_ids=input_ids,
            attention_mask=attention_mask,
            return_dict=True,
            output_hidden_states=output_hidden_states,
        )

        pmhc_hidden_states = output.hidden_states[-1]

        z, z_mask = self.embedder(
            tcr_hidden_states,
            pmhc_hidden_states,
            tcr_attention_mask,
            attention_mask
        )
        triangular_attentions = []
        for layer in self.layers:
            z = layer(
                z, 
                z_mask,
                tcr_hidden_states, 
                pmhc_hidden_states, 
                tcr_attention_mask, 
                attention_mask,
                output_attentions=output_triangular_attentions
            )
            if output_triangular_attentions:
                z, a1, a2 = z
                triangular_attentions.append((a1, a2))

        distogram_logits = self.distogram_head(z)
        contact_logits = self.contact_head(z[:,0,0,:])
        return ModelOutput(dict(
            z=z,
            z_mask=z_mask,
            distogram_logits=distogram_logits,
            contact_logits=contact_logits,
            triangular_attentions=triangular_attentions
        ))

def softmax_cross_entropy(logits, labels):
    loss = -1 * torch.sum(
        labels * torch.nn.functional.log_softmax(logits, dim=-1),
        dim=-1,
    )
    return loss

class TRabModelingBertForPseudoSequenceWithContactModule(TRabModelingBertForPseudoSequence):
    def __init__(
        self,
        bert_config: BertConfig,
        pooling: Literal["cls", "mean", "max"] = "mean",
        labels_number: int = 1,
        use_triton: bool = False,
        pretrained_checkpoint: Optional[str] = None,
        device = "cuda",
    ) -> None:
        super(TRabModelingBertForPseudoSequenceWithContactModule, self).__init__(
            bert_config,
            pooling,
            labels_number,
            use_triton
        )
        if pretrained_checkpoint is not None:
            self.load_state_dict(torch.load(pretrained_checkpoint))
        
        self.contact_module = TCRpMHCPairAttention(
            bert_config,
            use_triton
        )
        self.device = device
        self.to(device)
        
    def forward(self,
        *,
        tcr_input_ids: torch.Tensor,
        tcr_attention_mask: torch.Tensor,
        pmhc_input_ids: torch.Tensor,
        pmhc_attention_mask: torch.Tensor,
        tcr_pmhc_distogram: Optional[torch.Tensor] = None,
        tcr_pmhc_binding: Optional[torch.Tensor] = None,
        output_triangular_attentions: bool = False,
    ):
        tcr_output = self.model.forward(
            input_ids=tcr_input_ids,
            attention_mask=tcr_attention_mask,
            return_dict=True,
            output_hidden_states=True,
        )

        contact_output = self.contact_module(
            input_ids=pmhc_input_ids,
            attention_mask=pmhc_attention_mask,
            tcr_hidden_states=tcr_output.hidden_states[-1],
            tcr_attention_mask=tcr_attention_mask,
            output_triangular_attentions=output_triangular_attentions,
        )

        if tcr_pmhc_distogram is not None:
            contact_loss = (softmax_cross_entropy(
                contact_output.distogram_logits,
                tcr_pmhc_distogram
            ) * contact_output.z_mask).mean()

        if tcr_pmhc_binding is not None:
            binding_loss = torch.nn.functional.binary_cross_entropy_with_logits(
                contact_output.contact_logits,
                tcr_pmhc_binding
            )

        return {
            "tcr_output": tcr_output,
            "contact_output": contact_output,
            "contact_loss": contact_loss if tcr_pmhc_distogram is not None else None,
            "binding_loss": binding_loss if tcr_pmhc_binding is not None else None,
        }