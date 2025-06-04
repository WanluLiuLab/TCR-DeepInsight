import os
from pathlib import Path
import scanpy as sc
import warnings
import subprocess
MODULE_PATH = Path(__file__).parent
warnings.filterwarnings("ignore")

zenodo_accession = '12741480'


################
#     Human    #
################


def human_gex_reference_v2():
    """
    Load the human gex reference v2. If the dataset is not found, it will be downloaded from Zenodo.
    """
    default_path = os.path.join(MODULE_PATH, '../data/datasets/human_gex_reference_v2.h5ad')
    if os.path.exists(
        default_path
    ):
        return sc.read_h5ad(default_path)
    else:
        import subprocess
        result = subprocess.run(
            ["curl", "-L", "-o", default_path, f"https://zenodo.org/records/{zenodo_accession}/files/human_gex_reference_v2.h5ad?download=1"],
            check=True
        )
        if result.returncode == 0:
            return sc.read_h5ad(default_path)
        else:
            raise RuntimeError("Failed to download the dataset.")


def human_tcr_reference_v2():
    """
    Load the human tcr reference v2. Can be generated from the human gex reference v2 via `tdi.pp.unique_tcr_by_individual`. 
    If the dataset is not found, it will be downloaded from Zenodo.
    """
    default_path = os.path.join(MODULE_PATH, '../data/datasets/human_tcr_reference_v2.h5ad')
    if os.path.exists(
        default_path
    ):
        return sc.read_h5ad(default_path)
    else:
        result = subprocess.run(
            ["curl", "-L", "-o", default_path, f"https://zenodo.org/records/{zenodo_accession}/files/human_tcr_reference_v2.h5ad?download=1"],
            check=True
        )
        if result.returncode == 0:
            return sc.read_h5ad(default_path)
        else:
            raise RuntimeError("Failed to download the dataset.")


################
#     Mouse    #
################

def mouse_gex_reference_v1():
    """
    Load the mouse gex reference v1. 
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
    """
    Download the pretrained weights for the TCR-DeepInsight model. 
    """
    urls = [
        f"https://zenodo.org/records/{zenodo_accession}/files/human_gex_reference_v2.scatlasvae.ckpt?download=1",
        f"https://zenodo.org/records/{zenodo_accession}/files/human_bert_pseudosequence.tcr_v2.ckpt?download=1",
        f"https://zenodo.org/records/{zenodo_accession}/files/human_bert_pseudosequence_pca.tcr_v2.pkl?download=1"
    ]

    paths = [
        "../data/pretrained_weights/human_gex_reference_v2.scatlasvae.ckpt",
        "../data/pretrained_weights/human_bert_pseudosequence.tcr_v2.ckpt",
        "../data/pretrained_weights/human_bert_pseudosequence_pca.tcr_v2.pkl"
    ]

    for url, path in zip(urls, paths):
        result = subprocess.run(["curl", "-L", "-o", path, url], check=True)
        if result.returncode != 0:
            raise RuntimeError(f"Failed to download {url}")
