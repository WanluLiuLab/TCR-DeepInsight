Training Gene Expression (GEX) data by VAE
==========================================

This is a tutorial for training GEX data by VAE. We use scAtlasVAE package to integrate the multi-batch GEX data. We use the human huARdb v2 reference dataset as an example.

Load the reference dataset
--------------------------

.. code-block:: python
  :linenos:
  
  import tcr_deep_insight as tdi
  import torch 
  
  gex_reference_adata = tdi.data.human_gex_reference_v2()


Construct and train the GEX model
---------------------------------

.. code-block:: python
  :linenos:

  # GEXModelingVAE is an alias of scatlasvae.model.scAtlasVAE

  model = tdi.model.GEXModelingVAE(
    gex_reference_adata,
    batch_key=['study_name','sample_name'],
    n_latent=10,
    batch_hidden_dim=24
  )

  model.fit()


Extract the GEX embedding and Save the trained model
----------------------------------------------------

.. code-block:: python
  :linenos:

  gex_reference_adata.obsm['X_gex'] = model.get_latent_representation()

  torch.save(
    model.state_dict(), 
    "/PATH/TO//tcr_deep_insight/data/pretrained_weights/human_scatlasvae_gex_v2.ckpt"
  )

For more detailed information of `scAtlasVAE`, please refer to the `scAtlasVAE documentation <https://scatlasvae.readthedocs.io/en/latest/gex_integration.htmls>`_.


.. note::
  The trained model is available at `Zenodo <https://zenodo.org/records/12741480>`_.


Clustering and UMAP visualization
---------------------------------

The downstream analysis can be performed using `scanpy` package's standard workflow.

.. code-block:: python
  :linenos:

  import scanpy as sc

  sc.pp.neighbors(gex_reference_adata, use_rep='X_gex', n_neighbors=15)
  sc.tl.umap(gex_reference_adata)
  sc.tl.leiden(gex_reference_adata)

  sc.pl.umap(gex_reference_adata, color='leiden')