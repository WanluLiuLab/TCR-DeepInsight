import torch
import scanpy as sc 
from sklearn.decomposition import PCA, KernelPCA

# Faiss
import faiss

from sklearn.neighbors import KNeighborsClassifier
from joblib import dump, load
import datasets

# Built-in
import time
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from sklearn.utils import class_weight
from collections import Counter
from itertools import chain
from copy import deepcopy
import json
from typing import Any, Callable, Mapping, Union, Iterable, Tuple, Optional, Mapping
import os
from pathlib import Path
import warnings
import tqdm
import umap
from transformers import (
    BertConfig,
    PreTrainedModel,
    BertForMaskedLM
)

from scatlasvae.model._primitives import *
from scatlasvae.utils._tensor_utils import one_hot, get_k_elements
from scatlasvae.utils._decorators import typed, deprecated
from scatlasvae.utils._loss import LossFunction
from scatlasvae.utils._parallelizer import Parallelizer

# TCR Deep Insight
from ._deep_insight_result import TDIResult
from ..utils._compat import Literal
from ..model import GEXModelingVAE
from ..utils._tcr_definitions import TCR, _get_hla_pseudo_sequence
from ..model.modeling_bert._model import (
    TRabModelingBertForVJCDR3,
    TRabModelingBertForPseudoSequence
)
from ..model.tokenizers._tokenizer import (
    TRabTokenizerForVJCDR3,
    TRabTokenizerForPseudoSequence,
    trab_tokenizer_for_pseudosequence,
    tokenize_tcr_pseudo_sequence_to_fixed_length,
)
from ..model._model_utils import partial_load_state_dict

from ..model.modeling_bert._config import get_config, get_human_config, get_mouse_config
from ..utils._compat import Literal
from ..utils._logger import mt, mw, get_tqdm
from ..utils._utilities import random_subset_by_key, euclidean

from ..utils._definitions import (
    TRA_DEFINITION_ORIG,
    TRAB_DEFINITION_ORIG,
    TRB_DEFINITION_ORIG,
    TRA_DEFINITION,
    TRAB_DEFINITION,
    TRB_DEFINITION,
)
from ..utils._logger import mt
from ..utils._utilities import default_aggrf, default_pure_criteria
from ..utils._utilities import nearest_neighbor_eucliean_distances

from ._constants import FAISS_INDEX_BACKEND 

MODULE_PATH = Path(__file__).parent
warnings.filterwarnings("ignore")


@typed({
    "df": pd.DataFrame,
})
def add_tcr_pseudosequence_to_dataframe(
    df: pd.DataFrame,
) -> None:
    """
    Add TCR pseudosequence to dataframe

    :param df: dataframe

    :note:
        This method modifies the `df` inplace and adds the following columns:
            - tcr_pseudosequence
    """
    for i in ["TRAV", "TRAJ", "TRBV", "TRBJ", "CDR3a", "CDR3b"]:
        if i not in df.columns:
            raise ValueError(f"Column {i} not found in adata.obs.columns")
    df['tcr_pseudosequence'] = [TCR(
        alpha_v_gene=trav,
        alpha_j_gene=traj,
        beta_v_gene=trbv,
        beta_j_gene=trbj,
        alpha_cdr3=cdr3a,
        beta_cdr3=cdr3b,
    ).pseudo_sequence() for trav, traj, trbv, trbj, cdr3a, cdr3b in zip(
        df["TRAV"],
        df["TRAJ"],
        df["TRBV"],
        df["TRBJ"],
        df["CDR3a"],
        df["CDR3b"],
    )]

@typed({
    "df": pd.DataFrame,
})
def add_pmhc_pseudosequence_to_dataframe(
    df: pd.DataFrame,
) -> None:
    """
    Add PMHC pseudosequence to dataframe

    :param df: dataframe

    :note:
        This method modifies the `df` inplace and adds the following columns:
            - hla_pseudosequence
            - pmhc_pseudosequence
    """
    hla_pseudo_sequence = _get_hla_pseudo_sequence()
    for i in ["HLA", "peptide"]:
        if i not in df.columns:
            raise ValueError(f"Column {i} not found in adata.obs.columns")
    df["hla_pseudosequence"] = [
        hla_pseudo_sequence[x + ":01"]["pseudosequence"] for x in df["HLA"]
    ]
    df["pmhc_pseudosequence"] = list(
        map(":".join, zip(df["hla_pseudosequence"], df["peptide"]))
    )

@typed({
    "adata": sc.AnnData,
})
def tcr_adata_to_datasets(
    adata: sc.AnnData,
    tokenizer: Optional[Union[TRabTokenizerForVJCDR3, TRabTokenizerForPseudoSequence]] = None,
    show_progress: bool = False,
) -> datasets.arrow_dataset.Dataset:
    """
    Convert adata to tcr datasets
    :param adata: AnnData
    :param tokenizer: tokenizer
    :return: tcr datasets
    """
    for i in ["TRAV", "TRAJ", "TRBV", "TRBJ", "CDR3a", "CDR3b"]:
        if i not in adata.obs.columns:
            raise ValueError(f"Column {i} not found in adata.obs.columns")

    if isinstance(tokenizer, TRabTokenizerForVJCDR3):
        tcr_dataset = tokenizer.to_dataset(
            ids=adata.obs.index,
            alpha_v_genes=list(adata.obs["TRAV"]),
            alpha_j_genes=list(adata.obs["TRAJ"]),
            beta_v_genes=list(adata.obs["TRBV"]),
            beta_j_genes=list(adata.obs["TRBJ"]),
            alpha_chains=list(adata.obs["CDR3a"]),
            beta_chains=list(adata.obs["CDR3b"]),
        )
    elif isinstance(tokenizer, TRabTokenizerForPseudoSequence) or tokenizer is None:
        if "tcr_pseudosequence" not in adata.obs.columns:
            raise ValueError("Column 'tcr_pseudosequence' not found in adata.obs.columns")
        tokenize_result = tokenize_tcr_pseudo_sequence_to_fixed_length(
            list(filter(lambda x: type(x) == str, adata.obs["tcr_pseudosequence"])),
            show_progress=show_progress,
        )
        tcr_dataset = datasets.Dataset.from_dict(
            {"input_ids": tokenize_result[1], "attention_mask": tokenize_result[2]}
        )
    else:
        raise ValueError("Invalid tokenizer")
    return tcr_dataset


def tcr_dataframe_to_datasets(
    df: pd.DataFrame,
    tokenizer: Union[TRabTokenizerForVJCDR3, TRabTokenizerForPseudoSequence],
    show_progress: bool = False,
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
        
    if isinstance(tokenizer, TRabTokenizerForVJCDR3):
        tcr_dataset = tokenizer.to_dataset(
            ids=df.index,
            alpha_v_genes=list(df['TRAV']),
            alpha_j_genes=list(df['TRAJ']),
            beta_v_genes=list(df['TRBV']),
            beta_j_genes=list(df['TRBJ']),
            alpha_chains=list(df['CDR3a']),
            beta_chains=list(df['CDR3b']),
        )["train"]
    elif isinstance(tokenizer, TRabTokenizerForPseudoSequence):
        if 'tcr_pseudosequence' not in df.columns:
            raise ValueError("Column 'tcr_pseudosequence' not found in dataframe columns")
        tokenize_result = tokenize_tcr_pseudo_sequence_to_fixed_length(
            list(filter(lambda x: type(x) == str, df['tcr_pseudosequence'])),
            show_progress=show_progress,
        )
        tcr_dataset = datasets.Dataset.from_dict(
            {
                'input_ids': tokenize_result[1], 
                'attention_mask': tokenize_result[2]
            }
        )
    else:
        raise ValueError("Invalid tokenizer")
    return tcr_dataset

def to_embedding_tcr_only(
    model: Union[TRabModelingBertForVJCDR3, TRabModelingBertForPseudoSequence], 
    tcr_dataset: datasets.arrow_dataset.Dataset, 
    k: str = 'hidden_states', 
    device: str = 'cuda', 
    n_per_batch: int = 64, 
    show_progress: bool = False, 
) -> np.ndarray:
    """
    Get embedding from model
    :param model: nn.Module. The TCR model.
    :param tcr_dataset: datasets.arrow_dataset.Dataset. evaluation datasets.
    :param k: str. 'hidden_states' or 'last_hidden_state'. 
    :param device: str. 'cuda' or 'cpu'. If 'cuda', use GPU. If 'cpu', use CPU.
    :param n_per_batch: int. Number of samples per batch.
    :param show_progress: bool. If True, show progress bar.

    :return: embedding
    """

    if not isinstance(model, (TRabModelingBertForVJCDR3, TRabModelingBertForPseudoSequence)):
        raise ValueError(
            "Invalid model class. Should be one of tdi.model.modeling_bert.TRabModelingBertForVJCDR3 " + \
            "or tdi.model.modeling_bert.TRabModelingBertForPseudoSequence"
        )

    model.eval()
    all_embedding = []
    with torch.no_grad():
        if show_progress:
            import tqdm
            for_range = tqdm.trange(0, len(tcr_dataset), n_per_batch)
        else:
            for_range = range(0, len(tcr_dataset), n_per_batch)
        for j in for_range:
            tcr_input_ids = torch.tensor(
                tcr_dataset[j : j + n_per_batch]["input_ids"]
                if "input_ids" in tcr_dataset.features.keys()
                else tcr_dataset[j : j + n_per_batch]["tcr_input_ids"]
            ).to(device)
            tcr_attention_mask = torch.tensor(
                tcr_dataset[j : j + n_per_batch]["attention_mask"]
                if "attention_mask" in tcr_dataset.features.keys()
                else tcr_dataset[j : j + n_per_batch]["tcr_attention_mask"]
            ).to(device)
            indices_length = int(tcr_attention_mask.shape[0] / 2)

            if "token_type_ids" in tcr_dataset.features.keys():
                tcr_token_type_ids = torch.tensor(
                    tcr_dataset[j : j + n_per_batch]["token_type_ids"]
                    if "token_type_ids" in tcr_dataset.features.keys()
                    else tcr_dataset[j : j + n_per_batch]["tcr_token_type_ids"]
                ).to(device)

            if isinstance(model, TRabModelingBertForVJCDR3):
                output = model(
                    input_ids=tcr_input_ids,
                    attention_mask=tcr_attention_mask,
                    labels=tcr_input_ids,
                    token_type_ids=tcr_token_type_ids,
                )
            elif isinstance(model, TRabModelingBertForPseudoSequence):
                output = model(
                    input_ids=tcr_input_ids,
                    attention_mask=tcr_attention_mask,
                    labels=tcr_input_ids,
                )
            else:
                raise ValueError("Invalid model")
            all_embedding.append(output[k].detach().cpu().numpy())
    all_embedding = np.vstack(all_embedding)
    return all_embedding


def to_embedding_tcr_only_from_pandas(
    model: Union[TRabModelingBertForVJCDR3, TRabModelingBertForPseudoSequence],
    df: pd.DataFrame,
    tokenizer: Union[TRabTokenizerForVJCDR3, TRabTokenizerForPseudoSequence],
    device,
    n_per_batch=64,
    pooling: Literal["cls", "mean", "max", "pool", "trb", "tra", "weighted"] = "mean",
):
    original_pooling = None
    if pooling != model.pooling:
        original_pooling = model.pooling
        model.pooling = pooling
    if isinstance(model, TRabModelingBertForPseudoSequence):
        if pooling in [ "pool", "trb", "tra", "weighted"]:
            raise ValueError("Pooling method not supported for pseudo sequence")
        
    all_embedding = []
    for i in tqdm.trange(0, len(df), n_per_batch):
        ds = tcr_dataframe_to_datasets(
            df.iloc[i : i + n_per_batch, :], tokenizer
        )
        all_embedding.append(
            to_embedding_tcr_only(model, ds, "hidden_states", device)
        )

    if original_pooling is not None:
        model.pooling = original_pooling
    return np.vstack(all_embedding)

def get_pretrained_tcr_embedding(
    tcr_adata: sc.AnnData, 
    bert_config: Mapping[str, Any],
    encoding: Literal['vjcdr3','cdr123'] = 'cdr123',
    pooling: Literal["cls", "mean", "max", "pool", "trb", "tra", "weighted"] = "mean",
    species: Literal['human','mouse'] = 'human',
    checkpoint_path=os.path.join(MODULE_PATH, '../data/pretrained_weights/bert_tcr_768.ckpt'),
    pca_path: Optional[str] = None, 
    use_pca:bool=True, 
    use_kernel_pca: bool = False,
    use_faiss_pca: bool = True,
    pca_n_components: int = 50,
    device='cuda:0',
    n_per_batch: int = 256,
):
    """
    Get TCR embedding from pretrained BERT model

    .. note::
        This method modifies the `tcr_adata` inplace. It adds the following fields to `tcr_adata.obsm`:
            - X_tcr: TCR embedding
            - X_tcr_pca: PCA of TCR embedding

    :param tcr_adata: AnnData object containing TCR data
    :param bert_config: BERT config
    :param encoding: Encoding type. Default: 'cdr123'.
    :param pooling: Pooling method. Default: 'mean'
    :param species: Species. Default: 'human'
    :param checkpoint_path: Path to pretrained BERT model. Default: None
    :param pca_path: Path to PCA model, if previously saved. Default: None
    :param use_pca: Whether to use PCA. Default: True
    :param use_kernel_pca: Whether to use Kernel PCA. High memory required for large dataset. Default: False
    :param use_faiss_pca: Whether to use Faiss PCA instead fo scikit-learn. Default: True
    :param pca_n_components: Number of PCA components. Default: 50
    :param device: Device for Faiss PCA. Default: 'cuda:0'
    :param n_per_batch: Number of samples per batch in getting TCR embeddings. Default: 256
    
    """
    mt("Building BERT model")

    if encoding == 'vjcdr3':
        tcr_tokenizer = TRabTokenizerForVJCDR3(
            tra_max_length=48, 
            trb_max_length=48,
            species=species
        )
        tcr_model = TRabModelingBertForVJCDR3(
            bert_config,
            pooling_cls_position=1,
            labels_number=1,
            pooling=pooling,
            device=device
        )
    elif encoding == 'cdr123':
        tcr_tokenizer = trab_tokenizer_for_pseudosequence
        tcr_model = TRabModelingBertForPseudoSequence(
            bert_config,
            pooling = pooling,
            labels_number=1,
            device=device
        )

    else:
        raise ValueError("Invalid encoding")
    
    mt("Loading BERT model checkpoints...")
    try:
        partial_load_state_dict(tcr_model, torch.load(checkpoint_path, map_location=device))
    except RuntimeError as e:
        mt("Failed to load the full pretrained BERT model. Please make sure the checkpoint path is correct.")

    
    mt("Computing TCR Embeddings...")
    all_embedding = to_embedding_tcr_only_from_pandas(
        tcr_model,
        tcr_adata.obs,
        tcr_tokenizer,
        device,
        n_per_batch=n_per_batch,
        pooling=pooling
    )

    if use_pca and os.path.exists(pca_path):
        mt("Loading PCA model...")
        pca = load(pca_path)
        if use_kernel_pca:
            assert(isinstance(pca, KernelPCA), "PCA model is not a KernelPCA model")
        elif use_faiss_pca:
            assert(isinstance(pca, faiss.PCAMatrix), "PCA model is not a Faiss PCA model")
        else:
            assert(isinstance(pca, PCA), "PCA model is not a PCA model")
            
        mt("Performing PCA...")
        if use_faiss_pca:
            all_embedding_pca = pca.apply(all_embedding)
        else: 
            all_embedding_pca = pca.transform(all_embedding)
    elif use_pca:
        if use_kernel_pca:
            pca = KernelPCA(
                n_components=pca_n_components, 
                kernel="rbf", 
                gamma=10, 
                fit_inverse_transform=True, 
                alpha=0.1
            )
            all_embedding_pca = np.array(pca.fit_transform(all_embedding))
        elif use_faiss_pca:
            mat = faiss.PCAMatrix (all_embedding.shape[1], pca_n_components)
            mat.train(all_embedding)
            all_embedding_pca = mat.apply(all_embedding)
        else:
            pca = PCA(n_components=pca_n_components).fit(all_embedding)
            all_embedding_pca = np.array(pca.transform(all_embedding))
        if pca_path is not None:
            mt("Saving PCA model...")
            dump(pca, pca_path)
    else:
        all_embedding_pca = all_embedding


    tcr_adata.obsm["X_tcr"] = all_embedding
    tcr_adata.obsm["X_tcr_pca"] = all_embedding_pca

def _prepare_tcr_embedding(
    tcr_adata: sc.AnnData,
    layer_norm: bool = True,
    use_gex: bool = True,
    _tcr_embedding_weight: float = 6.,
) -> Tuple[np.ndarray, int, int]:
    if layer_norm:
        # LayerNorm for TCR and GEX
        ln_1 = torch.nn.LayerNorm(tcr_adata.obsm["X_tcr_pca"].shape[1])
        X_tcr_pca = ln_1(torch.tensor(tcr_adata.obsm["X_tcr_pca"]).float()).detach().numpy()
        ln_2 = torch.nn.LayerNorm(tcr_adata.obsm["X_gex"].shape[1])
        X_gex = ln_2(torch.tensor(tcr_adata.obsm["X_gex"]).float()).detach().numpy()
        if not use_gex:
            all_tcr_gex_embedding = X_tcr_pca
        else:
            all_tcr_gex_embedding = np.hstack([
                X_tcr_pca,
                X_gex
            ])
    else:
        X_tcr_pca = tcr_adata.obsm["X_tcr_pca"]
        X_gex = tcr_adata.obsm["X_gex"]
        if not use_gex:
            all_tcr_gex_embedding = tcr_adata.obsm["X_tcr_pca"]
        else:
            all_tcr_gex_embedding = np.hstack([
                tcr_adata.obsm["X_tcr_pca"],
                _tcr_embedding_weight*tcr_adata.obsm["X_gex"]
            ])
    return all_tcr_gex_embedding, X_tcr_pca.shape[1], X_gex.shape[1]


def cluster_tcr(
    tcr_adata: sc.AnnData,
    label_key: str = None,  
    include_hla_keys: Optional[Iterable[str]] = None,
    gpu=0,
    pure_label: bool = True,
    pure_criteria: Callable = default_pure_criteria,
    layer_norm: bool = True,
    max_distance: float = 3,
    max_cluster_size: int = 40,
    use_gex: bool = True,
    filter_intersection_fraction: float = 0.7,
    k: int = -1
):
    """
    Cluster TCRs by joint TCR-GEX embedding.

    :param tcr_adata: AnnData object containing TCR data
    :param label_key: Key of the label to cluster
    :param gpu: GPU ID. Default: 0
    :param pure_label: Whether to use pure label. Default: True
    :param pure_criteria: Pure criteria. Default: default_pure_criteria
    :param layer_norm: Whether to use LayerNorm. Default: True
    :param max_distance: Maximum distance. Default: 25
    :param max_cluster_size: Maximum cluster size for dTCR clusters. Default: 40
    :param use_gex: Whether to use GEX embedding. Default: True
    :param filter_intersection_fraction: Filter intersection fraction. Default: 0.7
    :param k: Number of nearest neighbors for background. Default: -1, which means background neighbors equal to cluster size

    :return: sc.AnnData containing clustered TCRs.

    .. note::
        The `gpu` parameter indicates GPU to use for clustering. If `gpu` is 0, CPU is used.

    """
    all_tcr_gex_embedding, tcr_dim, gex_dim = _prepare_tcr_embedding(
        tcr_adata,
        layer_norm=layer_norm,
        use_gex=use_gex
    )
    result = _cluster_tcr_by_label_core(
        tcr_adata.obs, 
        all_tcr_gex_embedding=all_tcr_gex_embedding, 
        query_tcr_gex_embedding=all_tcr_gex_embedding, 
        tcr_dim=tcr_dim,
        gex_dim=gex_dim,
        label_key=label_key, 
        include_hla_keys=include_hla_keys,
        gpu=gpu,
        pure_label=pure_label,
        pure_criteria=pure_criteria,
        max_distance=max_distance,
        max_cluster_size=max_cluster_size,
        k=k,
        filter_intersection_fraction=filter_intersection_fraction,
    )
    return result

def cluster_tcr_from_reference(
    tcr_adata: sc.AnnData,
    tcr_reference_adata: sc.AnnData,
    label_key: str = None,  
    include_hla_keys: Optional[Iterable[str]] = None,
    gpu=0,
    layer_norm: bool = True,
    pure_label: bool = True,
    pure_criteria: Callable = default_pure_criteria,
    max_distance:float = 3.,
    max_cluster_size: int = 40,
    use_gex: bool = True,
    filter_intersection_fraction: float = 0.7,
    k: int = -1
) -> TDIResult:
    """
    Cluster TCRs from reference. 

    :param tcr_adata: AnnData object containing TCR data
    :param tcr_reference_adata: AnnData object containing reference TCR data
    :param label_key: Key of the label to use for clustering
    :param reference_path: Path to reference TCR data. Default: None
    :param gpu: GPU to use. Default: 0
    :param layer_norm: Whether to use LayerNorm. Default: True
    :param max_distance: Maximum distance of the most dissimilar TCR in a dTCR cluster . Default: 25
    :param max_cluster_size: Maximum cluster size for dTCR clusters. Default: 40
    :param use_gex: Whether to use GEX embedding. Default: True
    :param filter_intersection_fraction: Filter intersection fraction. Default: 0.7
    :param k: Number of nearest neighbors for background. Default: -1, which means background neighbors equal to cluster size

    :return: TDIResult object containing clustered TCRs.

    .. note::
        The `gpu` parameter indicates GPU to use for clustering. If `gpu` is 0, CPU is used.
    
    """
    all_tcr_gex_embedding_reference, tcr_dim, gex_dim = _prepare_tcr_embedding(
        tcr_reference_adata,
        layer_norm=layer_norm,
        use_gex=use_gex
    )
    all_tcr_gex_embedding_query, tcr_dim, gex_dim = _prepare_tcr_embedding(
        tcr_adata,
        layer_norm=layer_norm,
        use_gex=use_gex
    )

    df = pd.concat([
        tcr_reference_adata.obs,
        tcr_adata.obs
    ])
    all_tcr_gex_embedding_reference = np.vstack([
        all_tcr_gex_embedding_reference,
        all_tcr_gex_embedding_query
    ])
    result = _cluster_tcr_by_label_core(
        df, 
        all_tcr_gex_embedding=all_tcr_gex_embedding_reference, 
        query_tcr_gex_embedding=all_tcr_gex_embedding_query, 
        tcr_dim=tcr_dim,
        gex_dim=gex_dim,
        label_key=label_key, 
        include_hla_keys=include_hla_keys,
        gpu=gpu,
        pure_label=pure_label,
        pure_criteria=pure_criteria,
        max_distance=max_distance,
        max_cluster_size=max_cluster_size,
        k=k,
        filter_intersection_fraction=filter_intersection_fraction,
    )
    return result


@typed(
    {
        "df": pd.DataFrame,
        "all_tcr_gex_embedding": np.ndarray,
        "query_tcr_gex_embedding": np.ndarray,
        "tcr_dim": int,
        "gex_dim": int,
        "label_key": str,
        "gpu": int,
        "pure_label": bool,
        "pure_criteria": Callable,
        "calculate_tcr_gex_distance": bool,
        "max_distance": float,
        "max_cluster_size": int,
    }
)
def _cluster_tcr_by_label_core(
    df,
    *,
    all_tcr_gex_embedding,
    query_tcr_gex_embedding,
    tcr_dim,
    gex_dim,
    label_key=None,
    include_hla_keys: Optional[List[str]] = None,
    use_gpu: bool = False,
    gpu=0,
    pure_label: bool = True,
    pure_criteria: Callable = lambda x, z: Counter(x).most_common()[0][1] / len(x) > 0.7
    and Counter(x).most_common()[0][0] == z,
    max_distance=3,
    max_cluster_size=40,
    calculate_tcr_gex_distance: bool = False,
    n_jobs: int = 1,
    k: int = -1,
    filter_intersection_fraction: float = 0.7,
    low_memory: bool = False,
    faiss_index_backend: FAISS_INDEX_BACKEND = FAISS_INDEX_BACKEND.KMEANS
) -> TDIResult:
    mt("Building Faiss index")
    all_tcr_gex_embedding = all_tcr_gex_embedding.astype(np.float32)
    query_tcr_gex_embedding = query_tcr_gex_embedding.astype(np.float32)

    if faiss_index_backend == FAISS_INDEX_BACKEND.KMEANS:
        if use_gpu:
            kmeans = faiss.Kmeans(
                all_tcr_gex_embedding.shape[1],
                all_tcr_gex_embedding.shape[0],
                niter=20,
                verbose=True,
                gpu=gpu
            )
        else:
            kmeans = faiss.Kmeans(
                all_tcr_gex_embedding.shape[1],
                all_tcr_gex_embedding.shape[0],
                niter=20,
                verbose=True,
            )
            
        kmeans.train(all_tcr_gex_embedding)
        D, I = kmeans.index.search(query_tcr_gex_embedding, max_cluster_size)
        index = kmeans.index
    elif faiss_index_backend == FAISS_INDEX_BACKEND.FLAT:
        use_gpu = False
        quantizer = faiss.IndexFlatL2(all_tcr_gex_embedding.shape[1])
        index = faiss.IndexIVFFlat(quantizer, all_tcr_gex_embedding.shape[1], 100)
        if use_gpu:
            try:
                res = faiss.StandardGpuResources()
            except:
                use_gpu = False
                mw("Warning: faiss GPU resources not found. Using CPU.")
            
            mt("Moving Faiss index to GPU")
            index = faiss.index_cpu_to_gpu(res, gpu, index)

        mt("Training Faiss index")
        index.train(all_tcr_gex_embedding.astype(np.float32))
        mt("Adding TCRs to Faiss index. This may take a while.")
        index.add(all_tcr_gex_embedding.astype(np.float32))
        mt("Searching for nearest neighbors")
        offset = k if k > 0 else 0
        D, I = index.search(query_tcr_gex_embedding, max_cluster_size + offset)
    else:
        raise ValueError("Invalid faiss index backend")
    
    D = np.sqrt(D) # squared L2 distance to L1 distance


    CAT_STORAGE_SIZE = 10
    NUM_CAT_STORAGE = int(max_cluster_size / CAT_STORAGE_SIZE)

    _result = [[] for _ in range(NUM_CAT_STORAGE+1)]

    if label_key is not None:
        label_map = dict(zip(range(len(df)), df[label_key]))
    else:
        label_map = dict(zip(range(len(df)), ["undefined"] * len(df)))
        mt("Warning: No label key is provided.")
        pure_label = False

    hla_map = None
    if include_hla_keys is not None:
        hla_map = {key: dict(zip(range(len(df)), df.loc[:, val].to_numpy())) for key, val in include_hla_keys.items()}

    cell_number_map = dict(zip(range(len(df)), df["number_of_cell"]))

    mt("Iterative select TCRs as clustering anchors")

    def FLATTEN(x):
        return [i for s in x for i in s]

    def par_func(data, queue):
        ret = []
        for i in data:
            label = np.array([label_map[x] for x in I[i]])
            cell_number = np.array([cell_number_map[x] for x in I[i]])
            hla_allelles = {k: np.array([hla_map[k][x] for x in I[i]]) for k in hla_map.keys()}
            if max_distance > 0:
                mp = np.argwhere(D[i] > max_distance)
                if len(mp) > 0:
                    init_j = mp[0][0]
                else:
                    init_j = max_cluster_size
            else:
                init_j = max_cluster_size

            for j in list(range(2, max(2, init_j + 1)))[::-1]:
                pure_criteria_pass = pure_criteria(label[:j], label[0])
                if pure_criteria_pass or (not pure_label):
                    d = label[0]
                    if pure_label:
                        pomc = list(filter(lambda x: label[x] == d, range(0, j)))
                        comp = list(
                            filter(lambda x: label[x] != d, range(1, max_cluster_size))
                        )[:j]
                    else:
                        pomc = list(range(0, j))
                        comp = list(range(j, max_cluster_size))[:j]

                    cluster_size = len(pomc)
                    
                    same_hla = None
                    if hla_map is not None:
                        same_hla = []
                        for h in hla_allelles.keys():
                            for k in np.unique(FLATTEN(hla_allelles[h][:j])):
                                if all(list(map(lambda x: k in x or x[0] == '-', zip(*list(hla_allelles[h][:j]))))):
                                    same_hla.append(k)
                        same_hla = ','.join(same_hla)


                    comp_size = k if k > 0 else cluster_size
                    ret.append(
                        (
                            (cluster_size - 1) // CAT_STORAGE_SIZE,
                            [label[0]]
                            + list(I[i][pomc])
                            + [-1]
                            * (
                                (
                                    CAT_STORAGE_SIZE
                                    * (((cluster_size - 1) // CAT_STORAGE_SIZE) + 1)
                                )
                                - len(pomc)
                            )
                            + [cluster_size]
                            + [
                                D[i][pomc].mean(),  # for tcr/gex similarity measurement
                                D[i][comp].mean(),  # for background tcr/gex distance
                                # (D[i][comp][:k] * cell_number[comp][:k]).sum() / sum(cell_number[comp][:k]) - \
                                #    (D[i][pomc] * cell_number[pomc]).sum() / sum(cell_number[pomc]),  # for disease association measurement
                                D[i][comp][:comp_size].mean() - D[i][pomc].mean(),
                                i,
                                same_hla,
                            ],
                        )
                    )
                    break
            if queue is not None:
                queue.put(1)
        if queue is not None:
            queue.put(0)
        return ret
    
    mt(f"Clustering clonotypes using {n_jobs} threads")
    start = time.time()
    if n_jobs > 1:
        mt("The estimated time is based on the number of threads. The actual time should be much shorter.")
        p = Parallelizer(n_jobs=n_jobs,)
        ret = p.parallelize(
            map_func=par_func,
            map_data=list(range(I.shape[0])),
            reduce_func=FLATTEN,
            progress=True,
            backend='loky'
        )()

    else:
        ret = []
        pbar = get_tqdm()(
            total=I.shape[0], bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}"
        )
        for i in range(I.shape[0]):
            label = np.array([label_map[x] for x in I[i]])
            cell_number = np.array([cell_number_map[x] for x in I[i]])
            hla_allelles = {k: np.array([hla_map[k][x] for x in I[i]]) for k in hla_map.keys()}
            mp = np.argwhere(D[i] > max_distance)
            if len(mp) > 0:
                init_j = mp[0][0]
            else:
                init_j = max_cluster_size

            for j in list(range(2, max(2, init_j + 1)))[::-1]:
                pure_criteria_pass = pure_criteria(label[:j], label[0])
                if hla_map is not None:
                    same_hla = False
                    same_hla_keys = []
                    hla = [np.array([hla_map[k][x] for x in I[i]]) for k in range(len(hla_map))]
                    # check if any of the HLA is same
                    for k in range(len(hla_map)):
                        if len(set(hla[k])) == 1:
                            same_hla = True
                            same_hla_keys.append(hla[k][0])

                if pure_criteria_pass or (not pure_label):
                    d = label[0]
                    if pure_label:
                        pomc = list(filter(lambda x: label[x] == d, range(0, j)))
                        comp = list(
                            filter(lambda x: label[x] != d, range(1, max_cluster_size))
                        )[:j]
                    else:
                        pomc = list(range(0, j))
                        comp = list(range(j, max_cluster_size))[:j]
                    cluster_size = len(pomc)

                    same_hla = None
                    if hla_map is not None:
                        same_hla = []
                        for h in hla_allelles.keys():
                            for k in np.unique(FLATTEN(hla_allelles[h][:j])):
                                if all(list(map(lambda x: k in x or x[0] == '-', zip(*list(hla_allelles[h][:j]))))):
                                    same_hla.append(k)
                        same_hla = ','.join(same_hla)


                    comp_size = k if k > 0 else cluster_size
                    ret.append(
                        (
                            (cluster_size - 1) // CAT_STORAGE_SIZE,
                            [label[0]]
                            + list(I[i][pomc])
                            + [-1]
                            * (
                                (
                                    CAT_STORAGE_SIZE
                                    * (((cluster_size - 1) // CAT_STORAGE_SIZE) + 1)
                                )
                                - len(pomc)
                            )
                            + [cluster_size]
                            + [
                                D[i][pomc].mean(),  # for tcr/gex similarity measurement
                                D[i][comp].mean(),  # for background tcr/gex distance
                                # (D[i][comp][:k] * cell_number[comp][:k]).sum() / sum(cell_number[comp][:k]) - \
                                #    (D[i][pomc] * cell_number[pomc]).sum() / sum(cell_number[pomc]),  # for disease association measurement
                                D[i][comp][:comp_size].mean() - D[i][pomc].mean(),
                                i,
                                same_hla,
                            ],
                        )
                    )
                    break

            pbar.update(1)
        pbar.close()

    mt(f"Clustering finished in {time.time() - start:.2f} seconds")
    for i in ret:
        if i[0] >= NUM_CAT_STORAGE:
            _result[-1].append(i[1])
        else:
            _result[i[0]].append(i[1])

    all_result_tcr = []
    all_tcr_sequence = df["tcr"].to_numpy()

    for ii in range(len(_result)):

        result = pd.DataFrame(_result[ii])
        max_cluster_size_ = (ii + 1) * CAT_STORAGE_SIZE

        if result.shape[0] == 0:
            continue
        # remove singletons clusters

        result = result[result.iloc[:, max_cluster_size_ + 1] > 1]
        result = result.sort_values(max_cluster_size_ + 1, ascending=False)
        if result.shape[0] == 0:
            continue

        result.index = list(range(len(result)))
        selected_indices = set()
        appeared_clusters = set()
        appeared_clusters_2 = {}
        indices_mapping = {}

        for i in range(len(result)):
            # t is a tuple of TCR indices representing this cluster
            t = tuple(
                sorted(
                    list(
                        filter(
                            lambda x: x >= 0,
                            result.iloc[i, 1 : max_cluster_size_ + 1].to_numpy(),
                        )
                    )
                )
            )

            if t in appeared_clusters:
                # this cluster is a exactly of a already selected cluster, skip
                continue

            flag = False
            st = set(t)
            for j in t:
                # for every indice of clonotype
                if j in appeared_clusters_2.keys():
                    for k in appeared_clusters_2[j]:
                        sk = set(k)
                        if st.issubset(sk):
                            flag = True
                            break
                        if (
                            len(st.intersection(sk)) / len(st) > filter_intersection_fraction and \
                            len(st.intersection(sk)) / len(sk) > filter_intersection_fraction and \
                            k in appeared_clusters
                        ):
                            if result.iloc[
                                indices_mapping[k], max_cluster_size_ + 2 # similarity
                            ] > result.iloc[
                                i, max_cluster_size_ + 2
                            ]:
                                    # this cluster is a subset of a already selected cluster
                                    # and the selected cluster has a higher mean distance (i.e.
                                    # lower tcr/gex similarity), remove the selected cluster
                                    appeared_clusters.remove(k)
                                    if indices_mapping[k] in selected_indices:
                                        selected_indices.remove(indices_mapping[k])
                                    for j in k:
                                        if k in appeared_clusters_2[j]:
                                            appeared_clusters_2[j].remove(k)
                            else:
                                flag = True
                                break
                                
                        else:
                            flag = True
                            break

            if flag:
                continue

            appeared_clusters.add(t)
            for j in t:
                if j not in appeared_clusters_2:
                    appeared_clusters_2[j] = [t]
                else:
                    appeared_clusters_2[j].append(t)

            selected_indices.add(i)
            indices_mapping[t] = i

        result = result.iloc[list(selected_indices)]

        result_tcr = result.copy()
        result_cell_number = result.copy()
        
        if not low_memory:
            for i in list(range(1, max_cluster_size_ + 1))[::-1]:
                result_tcr.iloc[:, i] = result_tcr.iloc[:, i].apply(
                    lambda x: all_tcr_sequence[x] if x >= 0 else "-"
                )

        result_tcr["number_of_individuals"] = list(
            map(
                lambda x: len(
                    np.unique(
                        list(
                            map(
                                lambda z: z.split("=")[-1],
                                filter(lambda x: x != "-", x),
                            )
                        )
                    )
                ),
                result_tcr.iloc[:, 1 : max_cluster_size_ + 1].to_numpy(),
            )
        )

        result_tcr["number_of_unique_tcrs"] = list(
            map(
                lambda x: len(
                    np.unique(
                        list(
                            map(
                                lambda z: "=".join(z.split("=")[:-1]),
                                filter(lambda x: x != "-", x),
                            )
                        )
                    )
                ),
                result_tcr.iloc[:, 1 : max_cluster_size_ + 1].to_numpy(),
            )
        )
        offset=7
        if "number_of_cells" in df.columns:
            all_number_of_cell = df["number_of_cells"].to_numpy()
            for i in list(range(1, max_cluster_size_ + 1))[::-1]:
                result_cell_number.iloc[:, i] = result_cell_number.iloc[:, i].apply(
                    lambda x: all_number_of_cell[x] if x >= 0 else "-"
                )
            result_tcr["number_of_cells"] = result_cell_number.iloc[
                :, 1 : max_cluster_size_ + 1
            ].sum(axis=1)
            offset += 1


        result_tcr.columns = (
            [label_key]
            + [f"TCRab{x}" for x in range(1, result_tcr.shape[1] - offset)]
            + [
                "number_of_tcrs",
                "mean_distance",
                "mean_distance_other",
                "distance_difference",
                "cluster_index",
                "number_of_individuals",
                "number_of_unique_tcrs",
            ]
            + (["number_of_cells"] if "number_of_cells" in df.columns else [])
        )

        # result_tcr = result_tcr[result_tcr['number_of_unique_tcrs'] > 1]
        # result_tcr = result_tcr[result_tcr['mean_distance'] > 1e-3]

        result_tcr["disease_association_score"] = result_tcr["distance_difference"]
        result_tcr["tcr_similarity_score"] = result_tcr["mean_distance"]

        if calculate_tcr_gex_distance:
            a = time.time()
            D_gex = nearest_neighbor_eucliean_distances(
                query_tcr_gex_embedding[:, tcr_dim:], I, result_tcr["cluster_index"]
            )
            D_tcr = nearest_neighbor_eucliean_distances(
                query_tcr_gex_embedding[:, :tcr_dim], I, result_tcr["cluster_index"]
            )

            result_tcr["mean_distance_gex"] = D_gex.mean(axis=1)
            result_tcr["mean_distance_tcr"] = D_tcr.mean(axis=1)

        result_tcr = result_tcr.drop(
            columns=["mean_distance", "mean_distance_other", "distance_difference"]
        )

        if not pure_label:
            result_tcr = result_tcr.drop(
                columns=[label_key, "disease_association_score"]
            )

        if not low_memory:
            tcrs = list(
                map(
                    ",".join,
                    result_tcr.loc[
                        :, [f"TCRab{x}" for x in range(1, max_cluster_size_ + 1)]
                    ].to_numpy(),
                )
            )
        else:
            tcrs = list(filter(lambda x: x >= 0, result_tcr.loc[
                :, [f"TCRab{x}" for x in range(1, max_cluster_size_ + 1)]
            ].to_numpy()))

        result_tcr = result_tcr.loc[
            :, list(filter(lambda x: not x.startswith("TCRab"), result_tcr.columns))
        ]
        result_tcr.insert(0, "TCRab", tcrs)
        all_result_tcr.append(result_tcr)

    result_tcr = pd.concat(all_result_tcr)
    return TDIResult(
        _data=sc.AnnData(
            obs=result_tcr,
            uns={
                "I": I,
                "D": D,
                "max_cluster_size": max_cluster_size,
                "max_distance": max_distance,
            },
        ),
        _tcr_df=df,
        _cluster_label=label_key,
        low_memory=low_memory,
        faiss_index=index,
    )


@typed({
    "reference_adata": sc.AnnData,
    "tcr_cluster_adata": sc.AnnData,
    "label_key": str,
    "map_function": Callable
})
def inject_labels_for_tcr_cluster_adata(
    reference_data: Union[sc.AnnData, pd.DataFrame], 
    tcr_cluster_adata: sc.AnnData, 
    label_key: str, 
    map_function: Callable = default_aggrf
):
    """
    Inject labels for tcr_cluster_adata based on reference_adata

    :param reference_adata: sc.AnnData. Reference AnnData object containing labels
    :param tcr_cluster_adata: sc.AnnData
    :param label_key: str. Key of the label to use for clustering in reference_adata.obs.columns
    :param map_function: Callable. Default: function that returns the most frequent label otherwise "Ambigious"
    
    :note:
        This method modifies the `df` inplace and adds the following columns:
            - `label_key`: The most frequent label in the cluster

    """
    # Get a list of lists of TCRs
    if isinstance(reference_data, sc.AnnData):
        reference_data = reference_data.obs 
       
    if 'tcr' not in reference_data.columns:
        raise ValueError("tcr column not found in reference_data. Please run `tdi.pp.update_anndata` first.")
        
    tcr_list = list(map(lambda x: x.split(','), tcr_cluster_adata.obs['TCRab']))
    
    # For each list of TCRs, find the label of the most abundant TCR
    labels = []
    tcr2int = dict(zip(reference_data['tcr'], range(len(reference_data['tcr']))))

    for tcrs in tcr_list:
        tcrs = list(filter(lambda x: x != '-', tcrs))
        if len(tcrs) == 0:
            labels.append("NA")
        else:
            # tcrs contains a list of TCRs
            # reference_adata.obs.loc[] finds the rows of reference_adata.obs where the TCRs are listed
            # reference_adata.obs.loc[... , label_key] gets the label_key column for those rows
            # map_function gets the aggregate function
            # map_function(reference_adata.obs.loc[... , label_key]) calls the aggregate function on the label_key column
            # labels.append(...) stores the result in the labels list
            labels.append(map_function(reference_data.iloc[
                list(map(lambda x: tcr2int[x], tcrs)),
            ][label_key])) 

    tcr_cluster_adata.obs[label_key] = labels

@typed({
    "tcr_adata": sc.AnnData,
    "tcr_reference_adata": sc.AnnData,
    "use_gpu": bool,
    "gpu": int,
    "max_tcr_distance": float,
    "layer_norm": bool,
    "faiss_index_backend": str
})
def query_tcr_without_gex(
    tcr_query_adata: sc.AnnData,
    tcr_reference_adata: sc.AnnData,
    use_gpu: bool = False,
    gpu=0,
    max_tcr_distance: float = 20.,
    layer_norm: bool = True,
    faiss_index_backend: FAISS_INDEX_BACKEND = FAISS_INDEX_BACKEND.KMEANS
):
    if layer_norm:
        ln_1 = torch.nn.LayerNorm(tcr_reference_adata.obsm["X_tcr_pca"].shape[1])
        X_tcr_pca_reference = ln_1(torch.tensor(tcr_reference_adata.obsm["X_tcr_pca"], dtype=torch.float32)).detach().cpu().numpy()
        X_tcr_pca_query = ln_1(torch.tensor(tcr_query_adata.obsm["X_tcr_pca"], dtype=torch.float32)).detach().cpu().numpy()
    else:
        X_tcr_pca_reference = tcr_reference_adata.obsm["X_tcr_pca"]
        X_tcr_pca_query = tcr_query_adata.obsm["X_tcr_pca"]

    mt("computing pairwise disease of TCRs...")

    if faiss_index_backend == "kmeans":
        kmeans = faiss.Kmeans(
            X_tcr_pca_query.shape[1],
            X_tcr_pca_query.shape[0],
            niter=20,
            verbose=True,
            gpu=gpu
        )
        # ignore faiss warnings
        kmeans.cp.min_points_per_centroid = 1
        kmeans.cp.max_points_per_centroid = 1000000000
        kmeans.train(X_tcr_pca_reference.astype(np.float32))
        D, I = kmeans.index.search(X_tcr_pca_query.astype(np.float32), 40)
    else:
        use_gpu = False
        quantizer = faiss.IndexFlatL2(X_tcr_pca_query.shape[1])
        index = faiss.IndexIVFFlat(quantizer, X_tcr_pca_query.shape[1], 100)
        if use_gpu:
            try:
                res = faiss.StandardGpuResources()
            except:
                use_gpu = False
                mt("Warning: faiss GPU resources not found. Using CPU.")
            index = faiss.index_cpu_to_gpu(res, gpu, index)

        index.train(X_tcr_pca_reference.astype(np.float32))
        index.add(X_tcr_pca_reference.astype(np.float32))
        D, I = index.search(X_tcr_pca_query.astype(np.float32), 40)
    D = np.sqrt(D) # squared L2 distance to L1 distance
    I[D > max_tcr_distance] = -1
    tcr_col_index = tcr_reference_adata.obs.columns.get_loc("tcr")
    result = []
    
    for i in I:
        result.append(
            json.dumps(
                list(
                    map(
                        lambda z: tcr_reference_adata.obs.iloc[z, tcr_col_index],
                        filter(lambda x: x != -1, i),
                    )
                )
            )
        )
    tcr_query_adata.obs["tcr_query"] = result
