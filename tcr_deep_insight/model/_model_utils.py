import scanpy as sc 
import datasets
import torch 
import tqdm 
import numpy as np
import pandas as pd

from ._model import TRABModelMixin
from ._tokenizer import TRABTokenizer
from ..utils._decorators import typed
from ..utils._compat import Literal

@typed({
    "adata": sc.AnnData,
    "tokenizer": TRABTokenizer
})
def tcr_adata_to_datasets(adata: sc.AnnData, tokenizer: TRABTokenizer) -> datasets.arrow_dataset.Dataset :
    """
    Convert adata to tcr datasets
    :param adata: AnnData
    :param tokenizer: tokenizer
    :return: tcr datasets
    """
    for i in ['TRAV', 'TRAJ', 'TRBV', 'TRBJ', 'CDR3a', 'CDR3b']:
        if i not in adata.obs.columns:
            raise ValueError(f"Column {i} not found in adata.obs.columns")
    tcr_dataset = tokenizer.to_dataset(
        ids=adata.obs.index,
        alpha_v_genes=list(adata.obs['TRAV']),
        alpha_j_genes=list(adata.obs['TRAJ']),
        beta_v_genes=list(adata.obs['TRBV']),
        beta_j_genes=list(adata.obs['TRBJ']),
        alpha_chains=list(adata.obs['CDR3a']),
        beta_chains=list(adata.obs['CDR3b']),
    )
    return tcr_dataset

@typed({
    "df": pd.DataFrame,
    "tokenizer": TRABTokenizer
})
def tcr_dataframe_to_datasets(
    df: pd.DataFrame,
    tokenizer: TRABTokenizer
) -> datasets.arrow_dataset.Dataset :
    """
    Convert dataframe to tcr datasets
    :param df: dataframe
    :param tokenizer: tokenizer
    :return: tcr datasets
    """
    for i in ['TRAV', 'TRAJ', 'TRBV', 'TRBJ', 'CDR3a', 'CDR3b']:
        if i not in df.columns:
            raise ValueError(f"Column {i} not found in dataframe columns")
    tcr_dataset = tokenizer.to_dataset(
        ids=df.index,
        alpha_v_genes=list(df['TRAV']),
        alpha_j_genes=list(df['TRAJ']),
        beta_v_genes=list(df['TRBV']),
        beta_j_genes=list(df['TRBJ']),
        alpha_chains=list(df['CDR3a']),
        beta_chains=list(df['CDR3b']),
    )
    return tcr_dataset


def to_embedding_tcr_only(
    model: TRABModelMixin, 
    eval_dataset: datasets.arrow_dataset.Dataset, 
    k: str, 
    device: str = 'cuda', 
    n_per_batch: int = 64, 
    progress: bool = False, 
    mask_tr: Literal['tra','trb','none'] = 'none'
):
    """
    Get embedding from model
    :param model: model
    :param eval_dataset: eval_dataset
    :param k: k
    :param device: device
    :param n_per_batch: n_per_batch
    :param progress: progress
    :param mask_tr: mask_tr
    :return: embedding
    """
    
    model.eval()
    all_embedding = []
    with torch.no_grad():
        if progress:
            for_range = tqdm.trange(0,len(eval_dataset),n_per_batch)
        else:
            for_range = range(0,len(eval_dataset),n_per_batch)
        for j in for_range:
           
            tcr_input_ids = torch.tensor(
                eval_dataset[j:j+n_per_batch]['input_ids'] if 'input_ids' in eval_dataset.features.keys() else  eval_dataset[j:j+n_per_batch]['tcr_input_ids']
            ).to(device)
            tcr_attention_mask = torch.tensor(
                eval_dataset[j:j+n_per_batch]['attention_mask'] if 'attention_mask' in eval_dataset.features.keys() else  eval_dataset[j:j+n_per_batch]['tcr_attention_mask']
            ).to(device)
            indices_length = int(tcr_attention_mask.shape[0]/2)
            if mask_tr == 'tra':
                # tcr_input_ids[:,2:indices_length] = _AMINO_ACIDS_INDEX[_AMINO_ACIDS_ADDITIONALS['MASK']]
                tcr_attention_mask[:,2:indices_length] = False
            elif mask_tr == 'trb':
                # tcr_input_ids[:,indices_length+2:indices_length*2] = _AMINO_ACIDS_INDEX[_AMINO_ACIDS_ADDITIONALS['MASK']]
                tcr_attention_mask[:,indices_length+2:indices_length*2] = False
            tcr_token_type_ids = torch.tensor(
                eval_dataset[j:j+n_per_batch]['token_type_ids'] if 'token_type_ids' in eval_dataset.features.keys() else  eval_dataset[j:j+n_per_batch]['tcr_token_type_ids']
            ).to(device)
      
            output = model(
                input_ids = tcr_input_ids,
                attention_mask = tcr_attention_mask,
                labels = tcr_input_ids,
                token_type_ids = tcr_token_type_ids,
            ) 
            all_embedding.append(output[k].detach().cpu().numpy())
    all_embedding = np.vstack(all_embedding)
    model.train()
    return all_embedding


def to_embedding_tcr_only_from_pandas_v2(
    model, 
    df, 
    tokenizer, 
    device,  
    n_per_batch=64, 
    mask_tr='none'
):
    all_embedding = []
    for i in tqdm.trange(0,len(df), n_per_batch):
        ds = tcr_dataframe_to_datasets(df.iloc[i:i+n_per_batch,:], tokenizer)['train']
        all_embedding.append(to_embedding_tcr_only(model, ds, 'hidden_states', device, mask_tr=mask_tr))
    return np.vstack(all_embedding)