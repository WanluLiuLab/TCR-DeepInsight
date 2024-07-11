import scanpy as sc
import numpy as np
import pandas as pd
from typing import Callable, List, Optional, Union, Iterable

from scatlasvae.utils._decorators import typed

from ..utils._logger import mt
from ..utils._definitions import (
    TRA_DEFINITION_ORIG,
    TRAB_DEFINITION_ORIG,
    TRB_DEFINITION_ORIG,
    TRA_DEFINITION,
    TRAB_DEFINITION,
    TRB_DEFINITION,
)
from ..utils._utilities import default_aggrf

@typed({
    'gex_adata': sc.AnnData, 
    'gex_embedding_key': str,
    'tcr_embedding_key': str,
    'joint_embedding_key': str,
})
def update_anndata(
    gex_adata: sc.AnnData,
    gex_embedding_key: str = 'X_gex',
    tcr_embedding_key: str = 'X_tcr',
    joint_embedding_key: str = 'X_gex_tcr',
):
    """
    Update the adata with the embedding keys

    .. note::
        This method modifies the `gex_adata` inplace. 
        added columns in .obs: `tcr`, `CDR3a`, `CDR3b`, `TRAV`, `TRAJ`, `TRBV`, `TRBJ`

    .. note::
        `"tcr"` should be in `gex_adata.obs.columns`.


    :param gex_adata: AnnData object
    :param gex_embedding_key: embedding key for gex
    :param tcr_embedding_key: embedding key for tcr
    :param joint_embedding_key: embedding key for joint
    """
    gex_adata.uns['embedding_keys'] = {
        "gex": gex_embedding_key,
        "tcr": tcr_embedding_key,
        "joint": joint_embedding_key,
    }
    for i,j in zip([
        'IR_VJ_1_junction_aa',
        'IR_VDJ_1_junction_aa',
        'IR_VJ_1_v_call',
        'IR_VJ_1_j_call',
        'IR_VDJ_1_v_call',
        'IR_VDJ_1_j_call'
    ], ['CDR3a', 
        'CDR3b', 
        'TRAV', 
        'TRAJ', 
        'TRBV', 
        'TRBJ'
    ]):
        if not (i in gex_adata.obs.columns or j in gex_adata.obs.columns):
            raise ValueError(f"{i} or {j} is not in adata.obs.columns")
        if i in gex_adata.obs.columns:
            gex_adata.obs[j] = gex_adata.obs[i]

    if not 'individual' in gex_adata.obs.columns:
        raise ValueError("individual is not in adata.obs.columns.")
    mt("TCRDeepInsight: initializing dataset")
    mt("TCRDeepInsight: adding 'tcr' to adata.obs")
    gex_adata.obs['tcr'] = None
    gex_adata.obs.iloc[:, list(gex_adata.obs.columns).index("tcr")] = list(map(lambda x: '='.join(x), gex_adata.obs.loc[:,TRAB_DEFINITION + ['individual']].to_numpy()))

@typed({
    'gex_adata': sc.AnnData,
    'gex_embedding_key': str,
    'agg_index': pd.DataFrame,
})
def aggregated_gex_embedding_by_tcr(
    gex_adata: sc.AnnData,
    gex_embedding_key: str,
    agg_index: pd.DataFrame,
):
    """
    Aggregate GEX embedding by TCR

    :param gex_adata: AnnData object
    :param agg_index: Aggregated index

    :return: aggregated GEX embedding
    """
    all_gex_embedding = []
    for i in agg_index['index']:
        all_gex_embedding.append(
            gex_adata.obsm[gex_embedding_key][i].mean(0)
        )
    all_gex_embedding = np.vstack(all_gex_embedding)
    return all_gex_embedding

@typed({
    'gex_adata': sc.AnnData,
    'gex_embedding_key': str,
    'additional_label_keys': list,
    'map_function': Callable
})
def unique_tcr_by_individual(
    gex_adata: sc.AnnData,
    embedding_key: Union[str, Iterable[str]] = 'X_gex',
    label_key: Optional[str] = None,
    additional_label_keys: Iterable[str] = None,
    aggregate_func: Callable = default_aggrf
):
    """
    Unique TCRs by individual and aggregate GEX embedding by TCR. Unique TCR is defined by the combination of TRAV,TRAJ,TRBV,TRBJ,CDR3α,CDR3β and individual. 
    Also aggregate GEX embedding by TCR, and add the aggregated GEX embedding to the tcr_adata.obsm[gex_embedding_key].
    
    :param gex_adata: AnnData object of gene expression data
    :param embedding_key: Key(s) in adata.obsm where GEX embedding is stored. Default: 'X_gex'
    :param label_key: Key in adata.obs where TCR type abels are stored. Default: 'cell_type', where 'cell_type' should be included in adata.obs.columns
    :param additional_label_keys: Additional keys in adata.obs where TCR type labels are stored. Default: None
    :param map_function: Function to map TCR type labels. Default: default_aggrf

    :return: TCR adata
    """
    if "individual" not in gex_adata.obs.columns:
        raise ValueError("individual is not in adata.obs.columns")
    
    gex_adata.obs['tcr'] = None
    gex_adata.obs.iloc[:, list(gex_adata.obs.columns).index("tcr")]=list(map(lambda x: '='.join(x), gex_adata.obs.loc[:,TRAB_DEFINITION + ['individual']].to_numpy()))
    gex_adata.obs['index'] = list(range(len(gex_adata)))

    if label_key is not None:
        agg_index = gex_adata.obs.groupby("tcr").agg({
            "index": list,
            label_key: aggregate_func
        })
    else:
        agg_index = gex_adata.obs.groupby("tcr").agg({
            "index": list
        })

    if isinstance(embedding_key, str) and embedding_key in gex_adata.obsm.keys():
        tcr_adata = sc.AnnData(
            obs = pd.DataFrame(
                list(map(lambda x: x.split("="), agg_index.index)),
                columns = TRAB_DEFINITION + ['individual']
            ),
            obsm={
                embedding_key: aggregated_gex_embedding_by_tcr(gex_adata, embedding_key, agg_index)
            }
        )
    elif isinstance(embedding_key, Iterable):
        tcr_adata = sc.AnnData(
            obs = pd.DataFrame(
                list(map(lambda x: x.split("="), agg_index.index)),
                columns = TRAB_DEFINITION + ['individual']
            ),
            obsm={
                k: aggregated_gex_embedding_by_tcr(gex_adata, k, agg_index)
                for k in embedding_key if k in gex_adata.obsm.keys()
            }
        )

    if label_key:
        tcr_adata.obs[label_key] = list(agg_index[label_key])

    tcr_adata.obs['number_of_cell'] = list(map(len, agg_index['index']))
        
    if additional_label_keys is not None:
        for k in additional_label_keys:
            agg_index = gex_adata.obs.groupby("tcr").agg({
                k: aggregate_func
            })
            tcr_adata.obs[k] = list(agg_index[k])

    tcr_adata.obs['tcr'] = None
    tcr_adata.obs.iloc[:, list(tcr_adata.obs.columns).index("tcr")] = list(
        map(
            lambda x: '='.join(x), 
            tcr_adata.obs.loc[:,TRAB_DEFINITION + ['individual']].to_numpy()
        )
    )
    return tcr_adata
