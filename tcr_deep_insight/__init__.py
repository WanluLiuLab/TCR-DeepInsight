from ._metadata import __version__, within_flit
import importlib
import subprocess
import warnings

try:
    subprocess.run('mafft --version'.split(), capture_output=True)
    mafft_is_installed = True
except:
    mafft_is_installed = False

faiss_is_installed = importlib.util.find_spec("faiss") is not None
# throw warnings
if not faiss_is_installed:
    warnings.warn(
        "The package `faiss` is not installed. "
        "You'll therefore not be able to use the `tl.cluster_tcr` and `tl.cluster_tcr_from_reference` function."
    )
if not mafft_is_installed:
    warnings.warn(
        "`mafft` is not installed. "
        "You'll therefore not be able to use the `pl.plot_gex_tcr_selected_tcrs` function."
    )

if not within_flit():  # see function docstring on why this is there
    # the actual API
    # (start with settings as several tools are using it)
    from . import model as model
    from . import tool as tl
    from . import preprocessing as pp
    from . import plotting as pl
    from . import utils as ut
    from . import data as data

    from .tool._deep_insight_result import TDIResult

    load_tdi_result = TDIResult.load_from_disk

    from anndata import AnnData, concat
    from anndata import (
        read_h5ad,
        read_csv,
        read_excel,
        read_hdf,
        read_loom,
        read_mtx,
        read_text,
        read_umi_tools,
    )


    # has to be done at the end, after everything has been imported
    import sys
    sys.modules.update({f'{__name__}.{m}': globals()[m] for m in ['model', 'tl', 'pp', 'pl', 'ut', 'data']})
    del sys