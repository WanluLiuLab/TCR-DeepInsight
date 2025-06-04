import os
from typing import Dict, List, Optional, Union, Mapping
import pandas as pd
import scanpy as sc 
import numpy as np 

import faiss
from faiss import IndexFlatL2

from ._constants import (
    FAISS_INDEX_BACKEND,
    TDI_RESULT_FIELD
)

from ..utils._compat import Literal
from ..utils._logger import Colors, get_tqdm
from ..utils._utilities import FLATTEN
from ..utils._tcr_definitions import TRAB_DEFINITION
from ..utils._tcr import TCR


def check_is_parquet_serializable(obj: pd.DataFrame):
    for c in obj.columns:
        if obj[c].dtype == 'O':
            obj[c] = obj[c].astype(str)
        elif isinstance(obj[c].dtype, pd.CategoricalDtype):
            obj[c] = obj[c].astype(str)
    return obj

class TDICluster(list):
    @property
    def number_of_unique_tcrs(self):
        return len(list(map(lambda x: x.to_tcr_string(), self)))
    
    @property
    def number_of_individuals(self):
        return len(np.unique(list(map(lambda x: x.individual, self))))
                   
    def __repr__(self):
        return f'<TDICluster at {hex(id(self))} with {len(self)} TCRs>'


class TDIResult:
    """
    A class to store the TDI clustering result

    :param _data: the clustering result
    :param _tcr_df: the tcr dataframe
    :param _tcr_adata: the tcr anndata
    :param _gex_adata: the gex anndata
    :param _cluster_label: the cluster label
    :param faiss_index: the faiss index
    :param low_memory: whether to use low memory mode
    """

    def __init__(
        self,
        _data: sc.AnnData,
        _tcr_df: pd.DataFrame = None,
        _tcr_adata: sc.AnnData = None,
        _gex_adata: Optional[sc.AnnData] = None,
        _cluster_label: Optional[str] = None,
        faiss_index: Optional[IndexFlatL2] = None,
        low_memory: bool = False
    ):
        self._data = _data
        self._data.obs['TCRab'] = list(map(TDICluster, self._data.obs['TCRab']))
        self._tcr_df = _tcr_df

        self._tcr_adata = _tcr_adata
        self._gex_adata = _gex_adata

        self._cluster_label = _cluster_label
        self._data.uns["_cluster_label"] = _cluster_label
        self.faiss_index = faiss_index
        self.low_memory = low_memory

    @property
    def data(self):
        return self._data

    @property
    def cluster_label(self):
        return self._cluster_label

    @property
    def I(self):
        return self._data.uns['I']

    @property
    def D(self):
        return self._data.uns['D']

    @property
    def tcr_df(self):
        return self._tcr_df
    
    def __getattribute__(self, name):
        if name in object.__getattribute__(self, '_data').obs.columns:
            return object.__getattribute__(self, '_data').obs[name].to_numpy()
        else:
            return super().__getattribute__(name)
        
    def select(self, indices):
        return TDIResult(
            self._data[indices],
            self._tcr_df,
            self._tcr_adata,
            self._gex_adata,
            self._cluster_label,
            self.faiss_index,
            self.low_memory
        )

    @cluster_label.setter
    def cluster_label(self, cluster_label: str):
        assert(cluster_label in self._data.obs.columns)
        self._cluster_label = cluster_label

    @property
    def tcr_adata(self):
        return self._tcr_adata

    @tcr_adata.setter
    def tcr_adata(self, tcr_adata: sc.AnnData):
        self._tcr_adata = tcr_adata

    @property
    def gex_adata(self):
        return self._gex_adata

    @gex_adata.setter
    def gex_adata(self, gex_adata: sc.AnnData):
        self._gex_adata = gex_adata

    def calculate_cluster_additional_information(self):

        self.data.obs['unique_tra'] = list(map(lambda x: 
            len(np.unique(list(map(lambda x: 
                    '='.join([x.cdr3a, x.trav, x.traj]),
                    x
            ))))==1, self.data.obs['TCRab']
        ))
        self.data.obs['unique_trb'] = list(map(lambda x: 
            len(np.unique(list(map(lambda x: 
                    '='.join([x.cdr3b, x.trbv, x.trbj]),
                    x
            ))))==1, self.data.obs['TCRab']
        ))

        self.data.obs['MAIT_TCR'] = list(map(lambda x: 
            all(list(map(lambda x: 
                    x.trav == 'TRAV1-2' and x.traj in ['TRAJ33','TRAJ20', 'TRAJ12'], 
                    x
            ))), self.data.obs['TCRab']
        ))


    def __repr__(self) -> str:
        base_string = f'{Colors.GREEN}TDIResult{Colors.NC} object containing {Colors.CYAN}{self.data.shape[0]}{Colors.NC} clusters\n'

        if self._cluster_label is not None:
            base_string += f'     Cluster label: {self._cluster_label}\n'

        if self._gex_adata is not None:
            base_string += f'     GEX data: \n' +\
                 '     ' + self._gex_adata.__repr__().replace('\n', '\n     ')
        if self._tcr_adata is not None:
            base_string += f'     TCR data: \n' +\
                 '     ' + self._tcr_adata.__repr__().replace('\n', '\n     ') 
        return base_string


    def save_to_disk(self, 
        save_path, 
        save_cluster_result_as_csv=True,
        save_tcr_data=True, 
        save_gex_data=True
    ):
        """
        Save the cluster result to disk

        :param save_path: the path to save the cluster result
        save_cluster_result_as_csv: whether to save the cluster result as csv files
        :param save_tcr_data: whether to save the tcr data
        :param save_gex_data: whether to save the gex data
        """
        if not os.path.exists(save_path):
            os.mkdir(save_path)
        column0 = None
        if isinstance(self._data.obs['TCRab'].iloc[0], list):
            column0 = list(self._data.obs["TCRab"])
            self._data.obs["TCRab"] = list(
                map(
                    lambda x: ",".join(list(map(lambda x: x.to_string(), x))),
                    column0,
                )
            )
        self._data.write_h5ad(os.path.join(save_path, 'cluster_data.h5ad'))
        if column0 is not None:
            self._data.obs["TCRab"] = column0
        if self._tcr_df is not None:
            self._tcr_df = check_is_parquet_serializable(self._tcr_df)
            self._tcr_df.to_parquet(os.path.join(save_path, 'tcr_data.parquet'))
        if save_tcr_data and self._tcr_adata is not None:
            self._tcr_adata.write_h5ad(os.path.join(save_path, 'tcr_data.h5ad'))
        if save_gex_data and self._gex_adata is not None:
            self._gex_adata.write_h5ad(os.path.join(save_path, 'gex_data.h5ad'))
        if self.faiss_index is not None:
            try:
                faiss.write_index(self.faiss_index, os.path.join(save_path, 'faiss_index.faiss'))
            except:
                pass
                # print('Failed to save faiss index')
        if save_cluster_result_as_csv:
            cluster_tcr = self.to_pandas_dataframe_tcr()
            cluster_index = self.to_pandas_dataframe_cluster_index()
            cluster_tcr.to_csv(os.path.join(save_path, "cluster_tcr.csv"), index=False)
            cluster_index.to_csv(
                os.path.join(save_path, "cluster_index.csv"), index=False
            )

    @classmethod
    def load_from_disk(
        cls, 
        save_path: str, 
        tcr_data_path: Optional[str] = None,
        gex_adata_path: Optional[str] = None
    ):
        """
        Load the cluster result from disk
        
        :param save_path: the path to load the cluster result
        :param tcr_data_path: the path to load the tcr data
        :param gex_adata_path: the path to load the gex data
        """
        data = sc.read_h5ad(os.path.join(save_path, 'cluster_data.h5ad'))
        if type(data.obs['TCRab'].iloc[0]) == str:
            data.obs["TCRab"] = list(
                map(
                    lambda x: list(
                        map(lambda z: TCR.from_string(z), filter(lambda tcr: tcr != '-', x.split(",")))
                    ),
                    data.obs["TCRab"],
                )
            )
        tcr_df = None
        if os.path.exists(os.path.join(save_path, 'tcr_data.parquet')):
            tcr_df = pd.read_parquet(os.path.join(save_path, 'tcr_data.parquet'))
        if tcr_data_path is not None:
            print('Loading tcr data from {}'.format(tcr_data_path))
            tcr_adata = sc.read_h5ad(tcr_data_path)
        else: 
            if os.path.exists(os.path.join(save_path, 'tcr_data.h5ad')):
                tcr_adata = sc.read_h5ad(os.path.join(save_path, 'tcr_data.h5ad'))
            else:
                tcr_adata = None

        if gex_adata_path is not None:
            print('Loading gex data from {}'.format(gex_adata_path))
            gex_adata = sc.read_h5ad(gex_adata_path)
        else:
            if os.path.exists(os.path.join(save_path, 'gex_data.h5ad')):
                gex_adata = sc.read_h5ad(os.path.join(save_path, 'gex_data.h5ad'))
            else:
                gex_adata = None

        cluster_label = None
        if '_cluster_label' in data.uns:
            cluster_label = data.uns['_cluster_label']

        faiss_index = None
        
        if os.path.exists(os.path.join(save_path, 'faiss_index.faiss')):
            try:
                faiss_index = faiss.read_index(os.path.join(save_path, 'faiss_index.faiss'))
            except:
                pass

        return cls(data, tcr_df, tcr_adata, gex_adata, cluster_label, faiss_index=faiss_index)

    def get_tcrs_for_cluster(
        self,
        label: Optional[Union[Mapping[str, str], Mapping[str, List[str]]]] = None,
        rank: int = 0, 
        rank_by: Literal["convergence", "disease_association"] = 'convergence',
        min_unique_tcr_number: int = 4,
        min_individual_number: int = 2,
        min_cell_number: int = 10,
        min_tcr_convergence_score: Optional[float] = None,
        min_disease_association_score: Optional[float] = None,
        return_background_tcrs: bool = False,
        additional_label_key_values: Optional[Dict[str, List[str]]] = None
    ):
        """
        Get the tcrs for a specific cluster

        :param label: the cluster label
        :param rank: the rank of the tcrs to return
        :param rank_by: the metric to rank the tcrs
        :param min_unique_tcr_number: the minimum number of unique tcrs in the cluster
        :param min_individual_number: the minimum number of individuals in the cluster
        :param min_cell_number: the minimum number of cells in the cluster
        :param min_tcr_convergence_score: the minimum convergence score
        :param min_disease_association_score: the minimum disease specificity score
        :param return_background_tcrs: whether to return other tcrs in the cluster
        :param additional_label_key_values: additional label key values to filter the cluster

        :return: a dictionary containing the tcrs and their metadata

        """
        return self._get_tcrs_for_cluster(
            label, 
            rank_by,
            rank, 
            min_unique_tcr_number=min_unique_tcr_number, 
            min_individual_number=min_individual_number,
            min_cell_number=min_cell_number,
            min_tcr_convergence_score = min_tcr_convergence_score,
            min_disease_association_score = min_disease_association_score,
            return_background_tcrs=return_background_tcrs,                
            additional_label_key_values=additional_label_key_values
        )

    def to_pandas_dataframe_tcr(
        self, 
        rank_by: Literal["convergence", "disease_association", "individual", "unique_tcr"] = 'convergence',
        return_background_tcrs: bool = False
    ):
        """
        Convert the cluster result to a pandas dataframe.

        :param rank_by: the metric to rank the tcrs
        :param return_background_tcrs: whether to return background tcrs for each cluster

        :return: a pandas dataframe containing the cluster result
        """
        ret = []
        cluster_indices = []
        cluster_labels = []
        if rank_by == 'convergence':
            self.data.obs.sort_values(TDI_RESULT_FIELD.CONVERGENCE.value, ascending=False, inplace=True)
        elif rank_by == 'disease_association':
            self.data.obs.sort_values(TDI_RESULT_FIELD.DISEASE_ASSOCIATION.value, ascending=True, inplace=True)
        elif rank_by == 'individual':
            self.data.obs.sort_values(TDI_RESULT_FIELD.NUMBER_OF_INDIVIDUAL.value, ascending=False, inplace=True)
        elif rank_by == 'unique_tcr':
            self.data.obs.sort_values(TDI_RESULT_FIELD.NUMBER_OF_UNIQUE_TCR.value, ascending=False, inplace=True)
        if self.cluster_label in self.data.obs.columns:
            for i, l, j in zip(
                self.data.obs.iloc[:, 0],
                self.data.obs[self.cluster_label],
                self.data.obs["cluster_index"],
            ):
                if isinstance(i, str):
                    ret.append(
                        list(
                            map(
                                lambda x: x.split("=")[:7] + ['foreground'],
                                filter(lambda x: x != "-", i.split(",")),
                            )
                        )
                    )
                elif isinstance(i, list):
                    ret.append(
                        list(
                            map(
                                lambda x: (
                                    x.cdr3a,
                                    x.cdr3b,
                                    x.trav,
                                    x.traj,
                                    x.trbv,
                                    x.trbj,
                                    x.individual,
                                    'foreground'
                                ),
                                i,
                            )
                        )
                    )
                cluster_indices.append([j] * len(ret[-1]))
                cluster_labels.append([l] * len(ret[-1]))
                size = len(ret[-1])
                if return_background_tcrs:
                    cluster_labels.append([])
                    ret.append([])
                    tcr_clusters = self.get_tcrs_by_cluster_index(j)
                    c = 0
                    for i, background_l in zip(tcr_clusters['tcrs'], tcr_clusters['labels']):
                        if background_l != l:
                            ret[-1].append(i.split("=")[:7] + ['background'])
                            c += 1
                            cluster_labels[-1].append(background_l)
                            if c == size:
                                break
                    cluster_indices.append([j] * c)

            df = pd.DataFrame(FLATTEN(ret), columns=TRAB_DEFINITION + ["individual", "type"])
            df["cluster_index"] = FLATTEN(cluster_indices)
            df["cluster_label"] = FLATTEN(cluster_labels)
        else:
            for i, j in zip(self.data.obs.iloc[:, 0], self.data.obs["cluster_index"]):
                if isinstance(i, str):
                    ret.append(
                        list(
                            map(
                                lambda x: x.split("=")[:7],
                                filter(lambda x: x != "-", i.split(",")),
                            )
                        )
                    )
                elif isinstance(i, list):
                    ret.append(
                        list(
                            map(
                                lambda x: (
                                    x.cdr3a,
                                    x.cdr3b,
                                    x.trav,
                                    x.traj,
                                    x.trbv,
                                    x.trbj,
                                    x.individual,
                                ),
                                i,
                            )
                        )
                    )
                cluster_indices.append([j] * len(ret[-1]))
            df = pd.DataFrame(FLATTEN(ret), columns=TRAB_DEFINITION + ["individual"])
            df["cluster_index"] = FLATTEN(cluster_indices)
        return df

    def to_pandas_dataframe_cluster_index(self):
        return self.data.obs.loc[
            :,
            ["cluster_index"]
            + list(
                filter(
                    lambda x: x not in ["cluster_index", "TCRab"], self.data.obs.columns
                )
            ),
        ]

    def _get_tcrs_for_cluster(
        self, 
        label: Optional[str] = None, 
        rank_by: Literal["convergence", "disease_association", "individual", "unique_tcr"] = 'convergence',
        rank: int = 0,
        min_unique_tcr_number: int = 4,
        min_individual_number: int = 2,
        min_cell_number: int = 10,
        min_tcr_convergence_score: Optional[float] = None,
        min_disease_association_score: Optional[float] = None,
        return_background_tcrs: bool = False,
        additional_label_key_values: Optional[Dict[str, List[str]]] = None
    ):
        if rank_by not in ['convergence', 'disease_association', "individual", "unique_tcr"]:
            raise ValueError('rank_by must be one of ["convergence", "disease_association"], got {}'.format(rank_by))

        min_tcr_convergence_score = min_tcr_convergence_score if min_tcr_convergence_score is not None else -1
        min_disease_association_score = min_disease_association_score if min_disease_association_score is not None else -1
        if additional_label_key_values is not None:
            additional_indices = np.bitwise_and.reduce(
                [
                    np.array(self.data.obs[key] == value) for key, value in additional_label_key_values.items()
                ]
            )
        else:
            additional_indices = np.ones(len(self.data.obs), dtype=bool)

        if label is not None:
            label_criteria = np.bitwise_and.reduce([
                np.array(self.data.obs[k] == v)
                if type(v) == str
                else np.array(self.data.obs[k].isin(v)) for k, v in label.items()])

            result_tcr = self.data.obs[
                label_criteria & 
                np.array(self.data.obs[TDI_RESULT_FIELD.NUMBER_OF_UNIQUE_TCR] >= min_unique_tcr_number) &
                np.array(self.data.obs[TDI_RESULT_FIELD.NUMBER_OF_INDIVIDUAL] >= min_individual_number) &
                np.array(self.data.obs[TDI_RESULT_FIELD.CONVERGENCE] >= min_tcr_convergence_score) &
                np.array(self.data.obs[TDI_RESULT_FIELD.DISEASE_ASSOCIATION] >= min_disease_association_score) &
                additional_indices
            ]
        else:
            result_tcr = self.data.obs[
                np.array(self.data.obs[TDI_RESULT_FIELD.NUMBER_OF_UNIQUE_TCR] >= min_unique_tcr_number) &
                np.array(self.data.obs[TDI_RESULT_FIELD.NUMBER_OF_INDIVIDUAL] >= min_individual_number) &
                np.array(self.data.obs[TDI_RESULT_FIELD.CONVERGENCE] >= min_tcr_convergence_score) &
                np.array(self.data.obs[TDI_RESULT_FIELD.DISEASE_ASSOCIATION] >= min_disease_association_score) &
                additional_indices
            ]

        if TDI_RESULT_FIELD.NUMBER_OF_CELL in result_tcr.columns:
            result_tcr = result_tcr[
                np.array(result_tcr[TDI_RESULT_FIELD.NUMBER_OF_CELL] >= min_cell_number)
            ]

        if rank_by == 'convergence':
            result = result_tcr.sort_values(
                TDI_RESULT_FIELD.CONVERGENCE, 
                ascending=False
            )
        elif rank_by == 'disease_association':
            result = result_tcr.sort_values(
                'disease_association_score', 
                ascending=False
            )
        elif rank_by == 'individual':
            result = result_tcr.sort_values(
                'number_of_individuals', 
                ascending=False
            )
        elif rank_by == 'unique_tcr':
            result = result_tcr.sort_values(
                'number_of_unique_tcrs', 
                ascending=False
            )

        tcrs = result.iloc[rank,0]
        cdr3a, cdr3b, trav, traj, trbv, trbj, individual = list(map(
            list, 
            np.array(list(map(lambda x: (x.cdr3a, x.cdr3b, x.trav, x.traj, x.trbv, x.trbj, x.individual), tcrs))).T
        ))

        if return_background_tcrs:
            return {
                'cluster_index': result['cluster_index'][rank],
                'tcrs': tcrs,
                'cdr3a': cdr3a,
                'cdr3b': cdr3b,
                'trav': trav,
                'traj': traj,
                'trbv': trbv,
                'trbj': trbj,
                'individual': individual
            }, self.get_tcrs_by_cluster_index(result['cluster_index'][rank], len(tcrs))
        else:
            return {
                'cluster_index': result['cluster_index'][rank],
                'tcrs': tcrs,
                'cdr3a': cdr3a,
                'cdr3b': cdr3b,
                'trav': trav,
                'traj': traj,
                'trbv': trbv,
                'trbj': trbj,
                'individual': individual
            }

    def _get_tcrs_gex_embedding_coordinates(self, use_rep: str = 'X_umap') -> Dict[str, np.ndarray]:
        """
        Get the gex embedding coordinates for each tcr

        :param use_rep: the representation to use. Default: 'X_umap'

        :return: a dictionary containing the gex embedding coordinates for each tcr
        """
        result_tcr = self.data.obs
        coordinates = {}
        tcrs2int = dict(zip(self.gex_adata.obs['tcr'], range(len(self.gex_adata))))
        pbar = get_tqdm()(total=len(result_tcr), bar_format='{l_bar}{bar:10}{r_bar}{bar:-10b}')
        for i in range(len(result_tcr)):
            tcrs = list(result_tcr.iloc[i,0].split(","))
            coordinates[result_tcr.index[i]] = self.gex_adata.obsm[use_rep][
                list(map(tcrs2int.get, list(filter(lambda x: x != '-', tcrs))))
            ]
            pbar.update(1)
        pbar.close()

        return coordinates

    def get_tcrs_by_cluster_index(
        self, 
        cluster_index: int,
        _n_after: int = 0
    ) -> List[str]:
        """
        Get the tcrs for a specific cluster

        :param cluster_index: the cluster index
        :param _n_after: the number of tcrs to skip
        """
        _result = self.tcr_df.iloc[
            self.I[int(cluster_index)]
        ]
        cdr3a = list(_result['CDR3a'])[_n_after:]
        cdr3b = list(_result['CDR3b'])[_n_after:]
        trav = list(_result['TRAV'])[_n_after:]
        traj = list(_result['TRAJ'])[_n_after:]
        trbv = list(_result['TRBV'])[_n_after:]
        trbj = list(_result['TRBJ'])[_n_after:]
        tcrs = list(_result['tcr'])[_n_after:]
        labels = None
        if self.cluster_label:
            labels = list(_result[self.cluster_label])[_n_after:]

        individual = list(_result['individual'])[_n_after:]

        return {
            'cluster_index': cluster_index,
            'tcrs': tcrs,
            'cdr3a': cdr3a,
            'cdr3b': cdr3b,
            'trav': trav,
            'traj': traj,
            'trbv': trbv,
            'trbj': trbj,
            'individual': individual,
            'labels': labels
        }
