Transfer and Clustering TCRαβ and GEX
=====================================

This notebook demonstrates how to transfer and cluster TCRαβ from new (query) datasets to reference dataset. We use the human huARdb v2 reference dataset as the reference dataset, and a
new dataset as the query dataset.


Load the reference dataset
--------------------------

.. code-block:: python
  :linenos:
  
  import tcr_deep_insight as tdi
  import torch 

  tcr_reference_adata = tdi.data.human_tcr_reference_v2()
  
  # Assume that both the GEX and TCR data are stored in the same anndata object
  assert("X_gex" in gex_reference_adata.obsm.keys())
  assert("X_tcr" in tcr_reference_adata.obsm.keys())


Load the Query dataset
----------------------

.. code-block:: python
  :linenos:
  
  gex_query_adata = sc.read_h5ad('path/to/query/adata.h5ad')

  # Update the anndata object
  tdi.pp.update_anndata(gex_query_adata)
  
  # Get unique TCRs by individual
  tcr_query_adata = tdi.pp.unique_tcr_by_individual(gex_query_adata)

  # Remove TRBV20OR9-2 as it is not in the IMGT reference
  tcr_query_adata = tcr_query_adata[tcr_query_adata.obs['TRBV']  != 'TRBV20OR9-2']

  # Add pseudo sequence to the dataframe
  tdi.tl.add_tcr_pseudosequence_to_dataframe(
    tcr_query_adata.obs
  )

  # We have trained a BERT model on the TCR sequence, 
  # which we can use to get the TCR embeddings
  tdi.tl.get_pretrained_tcr_embedding(
    tcr_adata=tcr_query_adata,
    bert_config=tdi.model.modeling_bert.get_human_config(),
    checkpoint_path='./tcr_deep_insight/data/pretrained_weights/human_bert_pseudosequence.tcr_v2.ckpt',
    pca_path='./tcr_deep_insight/data/pretrained_weights/human_bert_pseudosequence_pca.tcr_v2.pkl',
    use_pca=True,
    encoding='cdr123',
  )

  assert("X_tcr" in tcr_query_adata.obsm.keys())


Cluster the TCR and GEX data together using both the reference and query datasets
---------------------------------------------------------------------------------


.. code-block:: python
  :linenos:
  
  import os
  # We can now cluster the TCR and GEX data together
  tdi_cluster_result = tdi.tl.cluster_tcr_from_reference(
    tcr_query_adata,
    tcr_reference_adata,
    'disease_type',
    max_distance=4.,
    max_cluster_size=100,
    n_jobs=os.cpu_count(),
  )

  tdi_cluster_result.save_to_disk('path/to/save/cluster_result')