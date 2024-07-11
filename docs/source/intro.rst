Introduction
============

The emergence of single-cell immune profiling technology has led to the production of a large amount of data on single-cell gene expression (GEX) and T cell receptor (TCR), which has great potential for studying TCR biology and identifying effective TCRs. However, one of the major challenges is the lack of a reference atlas that provides easy access to these datasets. On the other hand, the use of TCR engineering in disease immunotherapy is rapidly advancing, and single-cell immune profiling data can be a valuable resource for identifying functional TCRs. Nevertheless, the lack of efficient computational tools to integrate and identify functional TCRs is a significant obstacle in this field.

To robustly identify potential disease associated TCRαβ clonotypes considering both TCR sequence similarity and transcriptome features from million-level paired TCRα/β repertoire, we developed a deep-learning based framework named TCR-DeepInsight. 

Rationale
---------

For large-scale and heterogeneous scRNA-seq gene expression (GEX) data, we use an Variational Autoencoder (VAE) to capture the biological signal and regress out technical or biological batch efffects.

.. image:: static/images/gex_vae.svg
   :width: 600
   :align: center

For full-length TCR repertoire data, we use a transformer-based model, BERT, to learn the TCR sequence features and TCR sequence similarity. 
We use the CDR1α, CDR2α, CDR3α, CDR1β, CDR2β, and CDR3β to represent the TCR sequence features. 

.. image:: static/images/tcr_bert.svg
   :width: 600
   :align: center

By concatenating the TCR sequence features and the GEX features, we use a faiss-based similarity search to identify the potential disease associated TCRαβ clonotypes 
with convergene gene expression profile and TCR sequence (cGxTr-TCRαβ), accelerated by GPU computing.

.. image:: static/images/trgx_embedding.svg
   :width: 600
   :align: center

We defined a TCR/GEX (TrGx) convergence score to measure the similarity of GEX and TCR within a TCR cluster. We also defined a disease-association score to measure the disease association of a TCR cluster.


Installation
------------

Hardware requirement for TCR-DeepInsight includes

1. CPU: >= 1 cores. Recommended >= 8 cores for large-scale dataset
2. RAM: >=1 Gb. Recommended >=64 Gb for large-scale dataset
3. VRAM >= 1Gb of CUDA-enabled GPU. Recommended >= 8 Gb for large-scale dataset
4. Disk space >= 1Gb. Recommended >= 100Gb for large-scale dataset


Operation System requirements for running TCR-DeepInsight include the installation of Python3 (Python3.9 used for development) and several PyPI packages. You can create a running environment using `conda`.

.. code-block:: shell
  :linenos:

    conda create -n t-deep-insight -f environment.yml
    conda activate t-deep-insight
    git clone git@github.com:WanluLiuLab/TCR-DeepInsight.git


Usage
-----

In IPython, simply import the package to get started:

.. code-block:: python
  :linenos:
    
    import tcr_deep_insight as tdi 
    
For more details, please refer to the tutorials.

Package Features 
----------------

Update Plan  
-----------
