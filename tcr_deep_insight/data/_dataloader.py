import os
from pathlib import Path
import scanpy as sc
import warnings
MODULE_PATH = Path(__file__).parent
warnings.filterwarnings("ignore")

zenodo_accession = '12176396'


################
#     Human    #
################


def human_gex_reference_v2():
    """
    Load the human gex reference v2
    """
    default_path = os.path.join(MODULE_PATH, '../data/datasets/human_gex_reference_v2.h5ad')
    if os.path.exists(
        default_path
    ):
        return sc.read_h5ad(default_path)
    else:
        os.system(f"curl -o {default_path} https://zenodo.org/record/{zenodo_accession}/files/human_gex_reference_v2.h5ad?download=1")
    

def human_tcr_reference_v2():
    """
    Load the human tcr reference v2. Can be generated from the human gex reference v2 via `tdi.pp.unique_tcr_by_individual`
    """
    default_path = os.path.join(MODULE_PATH, '../data/datasets/human_tcr_reference_v2.h5ad')
    if os.path.exists(
        default_path
    ):
        return sc.read_h5ad(default_path)
    else:
        os.system(f"curl -o {default_path} https://zenodo.org/record/{zenodo_accession}/files/human_tcr_reference_v2.h5ad?download=1")


################
#     Mouse    #
################

def mouse_gex_reference_v1():
    """
    Load the mouse gex reference v1
    """
    default_path = os.path.join(MODULE_PATH, '../data/datasets/mouse_gex_reference_v1.h5ad')
    if os.path.exists(
        default_path
    ):
        return sc.read_h5ad(default_path)
    else:
        raise FileNotFoundError("Dataset has not been published online. Please contact the author for more information.")

def mouse_tcr_reference_v1():
    """
    Load the mouse tcr reference v1. Can be generated from the mouse gex reference v1 via `tdi.pp.unique_tcr_by_individual`
    """
    default_path = os.path.join(MODULE_PATH, '../data/datasets/mouse_tcr_reference_v1.h5ad')
    if os.path.exists(
        default_path
    ):
        return sc.read_h5ad(default_path)
    else:
        raise FileNotFoundError("Dataset has not been published online. Please contact the author for more information.")



################
#     Model    #
################


def download_model_weights():
    os.system(f"curl -o ../data/pretrained_weights/human_bert_tcr_768_v2.ckpt https://zenodo.org/record/{zenodo_accession}/files/human_bert_tcr_768_v2.ckpt?download=1")
    os.system(f"curl -o ../data/pretrained_weights/human_bert_tcr_pca_v2.pkl https://zenodo.org/record/{zenodo_accession}/files/human_bert_tcr_pca_v2.pkl?download=1")

