from collections import Counter
from typing import Any, List, Mapping, Optional, Tuple
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib.patches import Rectangle
import numpy as np 
import logomaker
import os
from ..utils._utilities import seqs2mat, mafft_alignment
from ._palette import godsnot_102
from ..tool._deep_insight_result import TDIResult

import matplotlib.font_manager as fm
import warnings
from pathlib import Path

MODULE_PATH = Path(__file__).parent
warnings.filterwarnings("ignore")

# Path to your OTF font file
font_path =  os.path.join(MODULE_PATH, "fonts", "Arial.ttf")

# Create a FontProperties object
# Load the custom font
arial_font = fm.FontProperties(fname=font_path)
plt.rcParams['font.family'] = arial_font.get_name()


def set_plotting_params(
    dpi: int,
    fontsize: int = 12,
    fontfamily: str = "Arial",
    linewidth: float = 0.5,
):
    """
    Set default plotting parameters for matplotlib.

    :param dpi: dpi for saving figures
    :param fontsize: default fontsize
    :param fontfamily: default fontfamily
    :param linewidth: default linewidth
    
    """
    plt.rcParams['figure.dpi'] = dpi 
    plt.rcParams['savefig.dpi'] = dpi 
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['axes.linewidth'] = linewidth
    plt.rcParams['font.family'] = fontfamily
    mpl.rcParams['pdf.fonttype'] = 42 # saving text-editable pdfs
    mpl.rcParams['ps.fonttype'] = 42 # saving text-editable pdfs

amino_acids_color_scale = {
    'R': "#4363AE",
    'H': "#8282D1",
    'K': "#4164AE",
    'D': "#E61D26",
    'E': "#E61D26",
    'S': "#F8971D",
    'T': "#F8971D",
    'N': "#4FC4CC",
    'Q': "#4FC4CC",
    'C': "#E5E515",
    'G': "#ECEDEE",
    'P': "#DC9682",
    'A': "#C8C7C7",
    'V': "#148340",
    'I': "#148340",
    'L': "#148340",
    'M': "#E5E515",
    'F': "#3B5CAA",
    'Y': "#4164AE",
    'W': "#148340",
    '-': "#F7F7F7",
    '.': "#F7F7F7",
}
"""
Amino acids color scale for logo plots.
reference link: http://yulab-smu.top/ggmsa/articles/guides/Color_schemes_And_Font_Families.html
"""


def create_fig(figsize=(8, 4)) -> Tuple[matplotlib.figure.Figure, plt.Axes] :
    """
    Create a figure with a single axis.

    :param figsize: figure size. Default: (8, 4)
    :return matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    fig,ax=plt.subplots()           
    ax.spines['right'].set_color('none')     
    ax.spines['top'].set_color('none')
    #ax.spines['bottom'].set_color('none')     
    #ax.spines['left'].set_color('none')
    for line in ax.yaxis.get_ticklines():
        line.set_markersize(5)
        line.set_color("#585958")
        line.set_markeredgewidth(0.5)
    for line in ax.xaxis.get_ticklines():
        line.set_markersize(5)
        line.set_markeredgewidth(0.5)
        line.set_color("#585958")
    ax.set_xbound(0,10)
    ax.set_ybound(0,10)
    fig.set_size_inches(figsize)
    return fig,ax

def create_subplots(
    nrow, ncol, figsize=(8,8),gridspec_kw={}
) -> Tuple[matplotlib.figure.Figure, np.ndarray]:
    """
    Create a figure with multiple axes.

    :param nrow: number of rows
    :param ncol: number of columns
    :param figsize: figure size. Default: (8, 8)
    :param gridspec_kw: gridspec_kw. Default: {}
    :return matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    fig,axes=plt.subplots(nrow, ncol, gridspec_kw=gridspec_kw)
    for ax in axes.flatten():
        ax.spines['right'].set_color('none')     
        ax.spines['top'].set_color('none')
        for line in ax.yaxis.get_ticklines():
            line.set_markersize(5)
            line.set_color("#585958")
            line.set_markeredgewidth(0.5)
        for line in ax.xaxis.get_ticklines():
            line.set_markersize(5)
            line.set_markeredgewidth(0.5)
            line.set_color("#585958")
    fig.set_size_inches(figsize)
    return fig,axes


def piechart(
    ax: matplotlib.axes.Axes,
    annotation: Mapping[str, int], 
    cm_dict: Mapping[str, str],
    radius=1, 
    width=1, 
    setp=False, 
    show_annotation=False
):
    """
    Draw a pie chart.
    :param ax: matplotlib.axes.Axes
    :param annotation: annotationtation
    :param cm_dict: color map
    :param radius: radius
    :param width: width
    :param setp: whether to set width and edgecolor
    :param show_annotation: whether to show annotationtation
    """
    kw = dict(arrowprops=dict(arrowstyle="-"), zorder=0, va="center")
    pie, _ = ax.pie(annotation.values(),
                    radius=radius,
                    colors=[cm_dict[p] for p in annotation.keys()],
                    # wedgeprops=dict(width=width, edgecolor='w')
                    )
    if setp:
        plt.setp(pie, width=width, edgecolor='w')
    for i, p in enumerate(pie):
        theta1, theta2 = p.theta1, p.theta2
        center, r = p.center, p.r
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        x = r * np.cos(np.pi / 180 * (theta1+theta2)/2) + center[0]
        y = r * np.sin(np.pi / 180 * (theta1+theta2)/2) + center[1]
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "arc3, rad=0"
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        percentage = annotation[list(annotation.keys())[i]] / sum(list(annotation.values()))
        if show_annotation:
            if annotation[list(annotation.keys())[i]] / sum(list(annotation.values())) > 0.005:
                ax.annotate(list(annotation.keys())[i] + ", " + str(round(percentage * 100, 2)) + "%", xy=(x, y), xytext=(
                    x*1.2, y*1.2), horizontalalignment=horizontalalignment, size=6, fontweight=100)
    return pie



def plot_cdr3_sequence(
    sequences: List[str],
    alignment: bool = False, 
    labels: Optional[List[str]] = None,
    labels_palette: Optional[Mapping[str, Any]] = None,
    labels_postfix: Optional[str] = None,
    labels_postfix_palette: Optional[Mapping[str, Any]] = None,
    ax: Optional[plt.Axes] = None,
) -> Tuple[matplotlib.figure.Figure, plt.Axes]:
    """
    Plot CDR3 sequences.

    :param sequences: a list of CDR3 sequences
    :param alignment: whether to align the sequences. Default: False
    :param labels: a list of labels. Default: None
    :param labels_palette: a dictionary of labels and colors. Default: None
    :param ax: matplotlib.axes.Axes. Default: None

    :return: matplotlib.figure.Figure, matplotlib.axes.Axes
    """
    if alignment:
        sequences = mafft_alignment(sequences)
    if ax is None:
        fig,ax=create_fig((8, 0.36 * len(sequences)))
    else: 
        fig = ax.get_figure()
    renderer = fig.canvas.get_renderer()
    default_offset = 0

    if labels is not None:
        t = ax.text(x=0,y=0,s=sorted(labels, key=lambda x: len(x))[-1],fontfamily='arial', size=12)
        if hasattr(t.get_window_extent(renderer=renderer), 'inverse_transformed'):
            # matplotlib 3.3.0+
            bb = t.get_window_extent(renderer=renderer).inverse_transformed(ax.transData)
        else: 
            # matplotlib < 3.3.0
            bb = t.get_window_extent(renderer=renderer).transformed(ax.transData.inverted())

        default_offset = bb.width
        t.remove()
    max_x = 0
    max_y = 0.7 * len(sequences) + .9
    for i,s in enumerate(sequences):
        y = 0.7 * (len(sequences) - i) + .9
        for j,c in enumerate(s):
            ax.add_patch(
                Rectangle(
                    (1+j*0.24+default_offset, y), 
                    0.24, 
                    0.6, 
                    facecolor = amino_acids_color_scale[c]
                )
            )
            if 1+(1+j*0.24)+default_offset > max_x:
                max_x = 1+(1+j*0.24)+default_offset
            t = ax.text(x=0,y=0,s=c,fontfamily='Droid Sans Mono for Powerline', size=12)
            if hasattr(t.get_window_extent(renderer=renderer), 'inverse_transformed'):
                # matplotlib 3.3.0+
                bb = t.get_window_extent(renderer=renderer).inverse_transformed(ax.transData)
            else: 
                # matplotlib < 3.3.0
                bb = t.get_window_extent(renderer=renderer).transformed(ax.transData.inverted())

            offset = (0.24 - bb.width) / 2
            t.remove()
            t = ax.text(
                x = 1+j*0.24 + offset + default_offset,
                y= y + .1,
                s=c,
                fontfamily='Droid Sans Mono for Powerline', 
                size=12
            )
        if labels is not None:
            ax.text(
                x = 0,
                y = y + .1,
                s = labels[i],
                fontfamily='arial', 
                size=12,
                color = labels_palette[labels[i]] if labels_palette is not None else 'black'
            )

        if labels_postfix is not None:
            ax.text(
                x = max_x + 0.1,
                y = y + .1,
                s = labels_postfix[i],
                fontfamily='arial', 
                size=12,
                color = labels_postfix_palette[labels_postfix[i]] if labels_postfix_palette is not None else 'black'
            )

    ax.spines['right'].set_color('none')     
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')     
    ax.spines['left'].set_color('none')
    ax.set_xticks([])
    ax.set_yticks([])
    default_x_bound = ax.get_xbound()[1]
    max_x += 5
    if max_x > default_x_bound:
        ax.set_xbound(0, max_x)
        fig.set_size_inches(8 * (max_x / default_x_bound), 0.36 * len(sequences))
    ax.set_ybound(0, 0.7*i+2.5)
    return fig, ax 

def _plot_gex_selected_tcrs(
    adata,
    tcrs,
    color,
    palette,
    **kwargs
):
    """
    Plot the tcrs on the umap of the gex data

    :param gex_adata: sc.AnnData
    :param color: str
    :param tcrs: list
    :param palette: dict
    :return: fig, ax

    .. note::
        You should have `mafft` installed in your system to use this function
    """
    fig,ax=create_fig()
    fig.set_size_inches(3,3)

    ax.scatter(
        adata.obsm["X_umap"][:,0],
        adata.obsm["X_umap"][:,1],
        s=0.1,
        color=list(map(lambda x: palette[x], adata.obs[color])),
        linewidths=0,
        alpha=0.2,
    )

    obsm = adata[
        np.array(list(map(lambda x: x in tcrs,adata.obs['tcr'])))
    ].obsm["X_umap"]

    ax.scatter(
        obsm[:,0],
        obsm[:,1],
        s=10,
        marker='*',
        color='red',
    )

def _plot_gex_tcr_selected_tcrs(
    gex_adata: sc.AnnData,
    color: str,
    tcrs: list,
    tcrs_background: Optional[List[str]] = None,
    palette: Optional[dict] = None
):
    """
    Plot the tcrs on the umap of the gex data, with the TCRs as a pie chart and logo plot

    :param gex_adata: sc.AnnData
    :param color: str
    :param tcrs: list
    :param palette: dict (optional) 
    :return: fig, ax

    .. note::
        You should have `mafft` installed in your system to use this function
    """
    if palette is None:
        if len(set(gex_adata.obs[color])) <= 20:
            palette = sc.pl.palettes.default_20
        elif len(set(gex_adata.obs[color])) <= 28:
            palette = sc.pl.palettes.default_28
        else:
            palette = sc.pl.palettes.default_102
    
        palette = dict(zip(set(gex_adata.obs[color]), palette))

        
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['font.family'] = "arial"

    if tcrs_background is not None:
        gs_kw = dict(width_ratios=[2,1,1], height_ratios=[1,1,1,1,1])
        fig, axes = plt.subplot_mosaic(
            [[0,3,1],[0,4,2],[0,9,7],[0,10,8],[0,5,6]],
            gridspec_kw=gs_kw, 
            figsize=(7,3.2),
            #layout="constrained"
        )
    else: 
        gs_kw = dict(width_ratios=[2,1,1], height_ratios=[1,1,1])
        fig, axes = plt.subplot_mosaic(
            [[0,3,1],[0,4,2],[0,5,6]],
            gridspec_kw=gs_kw, 
            figsize=(7,3.2),
            #layout="constrained"
        )

    logomaker.Logo(seqs2mat(mafft_alignment(list(map(lambda x: x.split("=")[0], tcrs)))), ax=axes[1])
    logomaker.Logo(seqs2mat(mafft_alignment(list(map(lambda x: x.split("=")[1], tcrs)))), ax=axes[2])


    obsm = gex_adata.obsm["X_umap"]

    axes[0].scatter(
        gex_adata.obsm["X_umap"][:,0],
        gex_adata.obsm["X_umap"][:,1],
        s=0.5,
        color=list(map(lambda x: palette[x], gex_adata.obs[color])),
        linewidths=0,
        alpha=0.2,
    )

    
    if tcrs_background is not None:
        obsm = gex_adata[
            np.array(list(map(lambda x: x in tcrs_background, gex_adata.obs['tcr'])))
        ].obsm["X_umap"]

        axes[0].scatter(
            obsm[:,0],
            obsm[:,1],
            s=10,
            color='grey'
        )


    obsm = gex_adata[
            np.array(list(map(lambda x: x in tcrs, gex_adata.obs['tcr'])))
    ].obsm["X_umap"]

    axes[0].scatter(
        obsm[:,0],
        obsm[:,1],
        s=10,
        color='red'
    )
    
    obs = gex_adata[
        np.array(list(map(lambda x: x in tcrs, gex_adata.obs['tcr'])))
    ].obs

    piechart(
        axes[3],
        Counter(obs['TRAV']),
        show_annotation=True,
        cm_dict=dict(zip(Counter(obs['TRAV']).keys(), godsnot_102))
    )


    piechart(
        axes[4],
        Counter(obs['TRBV']),
        show_annotation=True,
        cm_dict=dict(zip(Counter(obs['TRBV']).keys(), godsnot_102))
    )

    axes[1].set_title("CDR3a dTCR", fontdict={'fontsize': 6})
    axes[2].set_title("CDR3b dTCR", fontdict={'fontsize': 6})
    axes[3].set_title("TRAV dTCR", fontdict={'fontsize': 6})
    axes[4].set_title("TRBV dTCR", fontdict={'fontsize': 6})

    if tcrs_background is not None:
        logomaker.Logo(seqs2mat(mafft_alignment(list(map(lambda x: x.split("=")[0], tcrs_background)))), ax=axes[7])
        logomaker.Logo(seqs2mat(mafft_alignment(list(map(lambda x: x.split("=")[1], tcrs_background)))), ax=axes[8])

        obs = gex_adata[
            np.array(list(map(lambda x: x in tcrs_background, gex_adata.obs['tcr'])))
        ].obs

        piechart(
            axes[9],
            Counter(obs['TRAV']),
            show_annotation=True,
            cm_dict=dict(zip(Counter(obs['TRAV']).keys(), godsnot_102))
        )


        piechart(
            axes[10],
            Counter(obs['TRBV']),
            show_annotation=True,
            cm_dict=dict(zip(Counter(obs['TRBV']).keys(), godsnot_102))
        )

        axes[7].set_title("CDR3a bTCR", fontdict={'fontsize': 6})
        axes[8].set_title("CDR3b bTCR", fontdict={'fontsize': 6})
        axes[9].set_title("TRAV bTCR", fontdict={'fontsize': 6})
        axes[10].set_title("TRBV bTCR", fontdict={'fontsize': 6})

    
    for ax in axes.values():
        ax.spines['right'].set_color('none') 
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.set_xticks([])
        ax.set_yticks([])

        for line in ax.yaxis.get_ticklines():
            line.set_markersize(5)
            line.set_color("#585958")
            line.set_markeredgewidth(0.5)
        for line in ax.xaxis.get_ticklines():
            line.set_markersize(5)
            line.set_markeredgewidth(0.5)
            line.set_color("#585958")


def plot_selected_tcrs(
    tcr_cluster_result: TDIResult,
    color: str,
    tcrs: List[str],
    tcrs_background: List[str],
    palette: Optional[dict] = None
):
    """
    Plot the tcrs on the umap of the gex data, with the TCRs as a pie chart and logo plot

    :param tcr_cluster_result: TDIResult
    :param color: str
    :param tcrs: list
    :param palette: dict (optional) 
    :return: fig, ax

    .. note::
        You should have `mafft` installed in your system to use this function
    """
    assert(tcr_cluster_result.gex_adata is not None)
    _plot_gex_tcr_selected_tcrs(
        tcr_cluster_result.gex_adata,
        color,
        tcrs,
        tcrs_background,
        palette
    )



