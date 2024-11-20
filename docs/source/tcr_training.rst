Training T cell receptor alpha and beta chains (TCRαβ) seqeunce by BERT
===============================================

This is a tutorial for training TCR sequence by BERT. We use the huARdb v2 reference dataset as an example.


Load the reference dataset
--------------------------

.. code-block:: python
  :linenos:
  
  import tcr_deep_insight as tdi
  import torch 
  
  # Load the human gex reference dataset
  gex_reference_adata = tdi.data.human_gex_reference_v2()

  # Update the anndata object
  tdi.pp.update_anndata(gex_reference_adata)
  
  # Get unique TCRs by individual
  tcr_reference_adata = tdi.unique_tcr_by_individual(gex_reference_adata)

  # Remove TRBV20OR9-2 as it is not in the IMGT reference
  tcr_reference_adata = tcr_reference_adata[tcr_reference_adata.obs['TRBV']  != 'TRBV20OR9-2']

  # Add pseudo sequence to the dataframe
  tdi.tl.add_tcr_pseudosequence_to_dataframe(
    tcr_reference_adata.obs
  )

  # convert the adata object to dataset.Dataset object
  tcr_dataset = tdi.tl.tcr_adata_to_datasets(tcr_reference_adata)


Training the TCR sequence by BERT
---------------------------------

.. code-block:: python
  :linenos:

  # create a TCR BERT model for training
  tcr_model = tdi.model.TRabModelingBertForPseudoSequence(
    tdi.model.modeling_bert.get_human_config(),
  ).to("cuda")

  # create a collator for tcr sequence to mask the input sequence
  tcr_collator = tdi.model.default_collator(
    tcr_model,
    species='human'
  )

  tcr_model_optimizer = tdi.model.default_optimizer(
    tcr_model,
    lr=1e-4,
  )
  
  tcr_model_trainer = tdi.model.TRabModelingBertForPseudoSequenceTrainer(
    model = tcr_model, 
    collator=tcr_collator, 
    train_dataset=tcr_dataset['train'], 
    test_dataset=tcr_dataset['test'], 
    optimizers=(tcr_model_optimizer, tcr_model_scheduler), 
    device='cuda'
  )

  tcr_model_trainer.fit(
    max_epoch=5, 
    show_progress=True, 
    n_per_batch=128, 
    shuffle=True
  )

  torch.save(
    tcr_model.state_dict(), 
    "/PATH/TO/tcr_deep_insight/data/pretrained_weights/human_bert_tcr_v2.ckpt"
  )
  
Obtaining TCR sequence embedding
--------------------------------

We can obtain the TCR sequence embedding by the following code:

.. code-block:: python
  :linenos:
  
  # This will add 'X_tcr' and 'X_tcr_pca' to the adata object
  # If the code is run in first time, it will save the pca model to the pca_path
  # If the pca model is already saved, it will load the pca model from the pca_path
  tdi.tl.get_pretrained_tcr_embedding(
    tcr_adata=tcr_reference_adata,
    bert_config=tdi.model.config.get_human_config(),
    checkpoint_path='./tcr_deep_insight/data/pretrained_weights/human_bert_pseudosequence.tcr_v2.ckpt',
    pca_path='./tcr_deep_insight/data/pretrained_weights/human_bert_pseudosequence_pca.tcr_v2.pkl',
    use_pca=True
  )

