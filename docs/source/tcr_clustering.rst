Clustering TCRαβ and GEX
========================================

This example shows how to cluster the joint information from 
TCRαβ and GEX data. We use the human huARdb v2 reference dataset as an example.

.. code-block:: python
  :linenos:
  
  import tcr_deep_insight as tdi
  import torch 
  
  # We have trained a VAE model on the GEX data
  # which we can use to get the GEX embeddings
  gex_reference_adata = tdi.data.human_gex_reference_v2()

  assert("X_gex" in gex_reference_adata.obsm.keys())

  # Update the anndata object
  tdi.pp.update_anndata(gex_reference_adata)
  
  # Get unique TCRs by individual
  tcr_reference_adata = tdi.pp.unique_tcr_by_individual(gex_reference_adata)

  # Remove TRBV20OR9-2 as it is not in the IMGT reference
  tcr_reference_adata = tcr_reference_adata[tcr_reference_adata.obs['TRBV']  != 'TRBV20OR9-2']

  # We have trained a BERT model on the TCR sequence, 
  # which we can use to get the TCR embeddings
  tdi.tl.get_pretrained_tcr_embedding(
    tcr_adata=tcr_reference_adata,
    bert_config=tdi.model.modeling_bert.get_human_config(),
    checkpoint_path='./tcr_deep_insight/data/pretrained_weights/human_bert_pseudosequence.tcr_v2.ckpt',
    pca_path='./tcr_deep_insight/data/pretrained_weights/human_bert_pseudosequence_pca.tcr_v2.pkl',
    use_pca=True,
    encoding='cdr123',
  )

  assert("X_tcr" in tcr_reference_adata.obsm.keys())

  # We can now cluster the TCR and GEX data together
  tdi_cluster_result = tdi.tl.cluster_tcr(
    tcr_reference_adata,
    'disease_type',
    max_distance=4.,
    max_cluster_size=100,
  )





