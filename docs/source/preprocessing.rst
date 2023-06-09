Preprocessing Immune Profiling Data 
===================================

We fork the preprocessing codes from the huARdb paper for convenience in introducing 
the dataset.

We use CellRanger (https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome) to 10x process single-cell GEX/TCR library and generate the raw data.


.. code-block::shell
  :linenos 

  mkdir EXAMPLE_SAMPLE
  cellranger count \
        --id=ProcessedData \
        --transcriptome=/path/to/transcriptome/reference \
        --jobmode=local \
        --localmem=64 \
        --localcores=32 \
        --sample=GEX \
        --fastqs=/path/to/fastq/
  
  cellranger-6.1.2 vdj \
        --id=ProcessedData \
        --reference=/path/to/vdj/reference \
        --jobmode=local \
        --localmem=64 \
        --localcores=16 \
        --sample=TCR \
        --chain=TR \
        --fastqs=/path/to/fastq/


The raw output file from CellRanger for RNA and VDJ library should be placed in the separated folders.
For example, the RNA library should be placed in the folder named "RNA" and the VDJ library should be placed in the folder named "VDJ".
The directory strucutre should look like:
    
.. code-block:: 
  :linenos:

    ├── EXAMPLE_SAMPLE
    │   ├── GEX/ProcessedData/outs
    │   │   ├── filtered_feature_bc_matrix
    │   │   │   ├── barcodes.tsv
    │   │   │   ├── features.tsv
    │   │   │   └── matrix.mtx
    │   │   ├── ...
    │   ├── VDJ/ProcessedData/outs
    │   │   ├── filtered_contig_annotations.csv
    │   │   ├── all_contig_annotations.csv
    │   │   ├── ...
    │   ├── ...

You can simply use our preprocessing pipeline to preprocess the data.

.. note::
  The preprocessing pipeline is only tested on Linux and Mac OS. 
  It is not tested on Windows. 


.. code-block:: python
  :linenos:

    from tcr_deep_insight.preprocessing._preprocess import *
    # VDJPreprocessingV1Human is for human data
    # VDJPreprocessingV1Mouse is for mouse data
    pp = VDJPreprocessingV1Human(
        cellranger_gex_output_path = "./EXAMPLE_SAMPLE/RNA/",
        cellranger_vdj_output_path = "./EXAMPLE_SAMPLE/VDJ/",
        output_path = "./EXAMPLE_SAMPLE_OUTPUT_DIR/"
    )
    pp.process(
        r_path = "/opt/anaconda3/envs/r403/bin/Rscript", 
        ref_data_path = HSAP_REF_DATA["Tonly"]
    )
    

.. note::
  You should have R installed and the Rscript executable should be in the PATH. 
  You can also specify the path to the Rscript executable by the argument "r_path". Make 
  sure your R version is >= 4.0.3 and the following R packages are installed:

  * Seurat.
  * DoubletFinder.
  * dplyr.
  * SingleR .
  * zellkonverter.
  * SeuratDisk.

The output files should include the following files:

.. code-block:: 
  :linenos:

    ├── EXAMPLE_SAMPLE_OUTPUT_DIR
    │   ├── all_contig_annotations.json
    │   ├── filterDoublet_prediction_sce.h5ad
    │   ├── results_preprocessed.h5ad
    │   ├── results_raw.h5ad
    │   ├── results_single_chain_preprocessed.h5ad
    │   ├── results_single_chain_raw.h5ad
    │   ├── rscript.R

* The `results_raw.h5ad` is the raw data file. 
* The `results_preprocessed.h5ad` is the preprocessed data file.
* The `results_single_chain_raw.h5ad` is the raw data file for single chain TCRs (or BCRs).
* The `results_single_chain_preprocessed.h5ad` is the preprocessed data file for single chain TCRs.
* The `filterDoublet_prediction_sce.h5ad` is the doublet prediction file.
* The `all_contig_annotations.json` is the json file for the VDJ annotation.
* The `rscript.R` is the R script for the preprocessing.


