Training Gene Expression (GEX) data by VAE
==========================================

This is a tutorial for training GEX data by VAE. We use scAtlasVAE package to integrate the multi-batch GEX data. We use the human huARdb v2 reference dataset as an example.

.. code-block:: python
  :linenos:
  
  import tcr_deep_insight as tdi
  import torch 
  
  gex_reference_adata = tdi.data.human_gex_reference_v2()

  # Update the anndata object

  model = tdi.model.GEXModelingVAE(
    gex_reference_adata,
    batch_key=['study_name','sample_name'],
    n_latent=10,
    batch_hidden_dim=24
  )

  model.fit()

  gex_reference_adata.obsm['X_gex'] = model.get_latent_representation()

  torch.save(
    model.state_dict(), 
    "/PATH/TO//tcr_deep_insight/data/pretrained_weights/human_scatlasvae_gex_v2.ckpt"
  )

For more detailed information, please refer to the `scAtlasVAE documentation <https://scatlasvae.readthedocs.io/en/latest/gex_integration.htmls>`_.