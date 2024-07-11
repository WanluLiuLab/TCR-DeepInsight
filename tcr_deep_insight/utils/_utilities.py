import os
import sys
import numpy as np
import pandas as pd
import torch
from collections import Counter
from umap.distances import euclidean
import warnings
from ._logger import Colors
from typing import List
import torch
from functools import partial


from scatlasvae.utils._parallelizer import Parallelizer
from scatlasvae.utils._decorators import deprecated


def FLATTEN(x): 
    return [i for s in x for i in s]


def highlight_iterable(arr, print_fn, highlight_fn, highlight_condition, sep=', '):
    return sep.join([highlight_fn(print_fn(x)) if highlight_condition(x) else print_fn(x) for x in arr])


def _nearest_neighbor_eucliean_distances(
    X,
    queue
):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)    
        out = []
        for d in X:
            x1, x2 = d
            if queue is not None:
                queue.put(None)
            ret = []
            for x in x2:
                ret.append(euclidean(x1, x))
            out.append(np.array(ret, dtype=np.float32))
        return np.vstack(out)

def nearest_neighbor_eucliean_distances_parallel(
    X: np.ndarray, 
    neigh_indices: np.ndarray,
    sel_indices: np.ndarray = None,
    n_jobs: int = os.cpu_count(),
    backend='threading'
):
    data = np.array(list(zip(X, X[neigh_indices])))
    if sel_indices is not None:
        data = data[sel_indices]
    p = Parallelizer(n_jobs=n_jobs)
    result = p.parallelize(
        map_func=_nearest_neighbor_eucliean_distances, 
        map_data=data, 
        reduce_func=FLATTEN,
        backend=backend
    )()
    return np.vstack(result)

def nearest_neighbor_eucliean_distances(
    X: np.ndarray, 
    neigh_indices: np.ndarray,
    sel_indices: np.ndarray = None,
    backend='threading',
    n_jobs: int = os.cpu_count(),
):
    np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)    
    return nearest_neighbor_eucliean_distances_parallel(
        X, 
        neigh_indices, 
        sel_indices,
        n_jobs=n_jobs,
        backend=backend
    )

def multi_values_dict(keys, values):
    ret = {}
    for k,v in zip(keys, values):
        if k not in ret.keys():
            ret[k] = [v]
        else:
            ret[k].append(v)
    return ret

def mafft_alignment(sequences):
    result = []
    import sys
    from Bio.Align.Applications import MafftCommandline
    import tempfile
    with tempfile.NamedTemporaryFile() as temp:
        temp.write('\n'.join(list(map(lambda x: '>seq{}\n'.format(x[0]) + x[1], enumerate(sequences)))).encode())
        temp.seek(0)
        mafft_cline = MafftCommandline(input=temp.name)
        stdout,stderr=mafft_cline()
    for i,j in enumerate(stdout.split("\n")):
        if i % 2 != 0:
            result.append(j.replace("-","."))
    return result

def seqs2mat(sequences, char_set = list('ACDEFGHIKLMNPQRSTVWY'), gap_character = '.'):
    mat = np.zeros((len(sequences[0]), len(char_set)))
    for i in range(len(sequences[0])):
        count = Counter(list(map(lambda x:x[i], sequences)))
        for k,v in count.items():
            if k != gap_character:
                mat[i][char_set.index(k)] = v
    mat = pd.DataFrame(mat, columns = char_set)
    return mat

def random_subset_by_key(adata, key, n):
    from collections import Counter
    counts = {k:v/len(adata) for k,v in Counter(adata.obs[key]).items()}
    ns = [(k,int(v*n)) for k,v in counts.items()]
    all_indices = []
    for k,v in ns:
        indices = np.argwhere(adata.obs[key] == k).flatten()
        if len(indices) > 0:
            indices = np.random.choice(indices, v, replace=False)
            all_indices.append(indices)
    all_indices = np.hstack(all_indices)
    return adata[all_indices]

def exists(x):
    return x != None

def sliceSimutaneuously(a, index):
    return pd.DataFrame(a).iloc[index,index].to_numpy()

def mask_split(tensor, indices):
    sorter = torch.argsort(indices)
    _, counts = torch.unique(indices, return_counts=True)
    return torch.split(tensor[sorter], counts.tolist())

def print_version():
    print(Colors.YELLOW)
    print('Python VERSION:{}\n'.format(Colors.NC), sys.version)
    print(Colors.YELLOW)
    print('PyTorch VERSION:{}\n'.format(Colors.NC), torch.__version__)
    print(Colors.YELLOW)
    print('CUDA VERSION{}\n'.format(Colors.NC))
    from subprocess import call
    try: call(["nvcc", "--version"])
    except: pass
    print(Colors.YELLOW)
    print('CUDNN VERSION:{}\n'.format(Colors.NC), torch.backends.cudnn.version())
    print(Colors.YELLOW)
    print('Number CUDA Devices:{}\n'.format(Colors.NC), torch.cuda.device_count())
    try:
        print('Devices             ')
        call(["nvidia-smi", "--format=csv", "--query-gpu=index,name,driver_version,memory.total,memory.used,memory.free"])
    except FileNotFoundError:
        # There is no nvidia-smi in this machine
        pass
    if torch.cuda.is_available():
        print('Active CUDA Device: GPU', torch.cuda.current_device())
        print ('Available devices     ', torch.cuda.device_count())
        print ('Current cuda device   ', torch.cuda.current_device())
    else:
        # cuda not available
        pass

def read_tsv(path, header:bool = True, skip_first_line: bool = False, return_pandas: bool = True):
    result = []
    if os.path.exists(path):
        f = open(path)
        if skip_first_line:
            line = f.readline()
        header_length = None
        if header:
            header = f.readline().strip().split('\t')
            header_length = len(header)

        while 1:
            line = f.readline()
            if not line:
                break
            line = line.strip().split('\t')
            if not header_length:
                header_length = len(line)
            result.append(line[:header_length])
        f.close()
        if return_pandas:
            if header:
                return pd.DataFrame(result, columns = header)
            else:
                return pd.DataFrame(result)
        else:
            return result
    else:
        it = iter(path.split('\n'))
        if skip_first_line:
            line = next(it)
        header_length = None
        if header:
            header = next(it).strip().split('\t')
            header_length = len(header)

        while 1:
            try:
                line = next(it)
                if not line:
                    break
                line = line.strip().split('\t')
                if not header_length:
                    header_length = len(line)
                result.append(line[:header_length])
            except:
                break 
        if return_pandas:
            if header:
                return pd.DataFrame(list(filter(lambda x: len(x) == 125, result)), columns = header)
            else:
                return pd.DataFrame(list(filter(lambda x: len(x) == 125, result)))
        else:
            return result

def default_aggrf(i, ambiguous="Ambiguous"):
    if len(np.unique(i)) == 1:
        return i[0]
    if len(i) == 2:
        return ambiguous
    else:
        c = Counter(i)
        return sorted(c.items(), key=lambda x: -x[1])[0][0]
            
def default_pure_criteria(x,y):
    return (Counter(x).most_common()[0][1] / len(x) > 0.8) and Counter(x).most_common()[0][0] == y

def exists(x):
    return x is not None

def absent(x):
    return x is None
