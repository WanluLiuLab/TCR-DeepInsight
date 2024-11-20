API
===


.. autoclass:: tcr_deep_insight.SPECIES
    :members:
    :undoc-members:

.. autoclass:: tcr_deep_insight.TCR
    :members:
    :undoc-members:

Preprocessing
-------------

.. autofunction:: tcr_deep_insight.pp.update_anndata
.. autofunction:: tcr_deep_insight.pp.unique_tcr_by_individual

Data
----

.. autofunction:: tcr_deep_insight.data.human_gex_reference_v2
.. autofunction:: tcr_deep_insight.data.human_tcr_reference_v2
.. autofunction:: tcr_deep_insight.data.mouse_gex_reference_v1
.. autofunction:: tcr_deep_insight.data.mouse_tcr_reference_v1

Tool
----

.. autofunction:: tcr_deep_insight.tl.add_tcr_pseudosequence_to_dataframe
.. autofunction:: tcr_deep_insight.tl.add_pmhc_pseudosequence_to_dataframe

.. autofunction:: tcr_deep_insight.tl.tcr_adata_to_datasets
.. autofunction:: tcr_deep_insight.tl.tcr_dataframe_to_datasets
.. autofunction:: tcr_deep_insight.tl.to_embedding_tcr_only
.. autofunction:: tcr_deep_insight.tl.to_embedding_tcr_only_from_pandas
.. autofunction:: tcr_deep_insight.tl.get_pretrained_tcr_embedding

Clustering
~~~~~~~~~~

.. autofunction:: tcr_deep_insight.tl.cluster_tcr
.. autofunction:: tcr_deep_insight.tl.cluster_tcr_from_reference
.. autofunction:: tcr_deep_insight.tl.inject_labels_for_tcr_cluster_adata

.. autoclass:: tcr_deep_insight.tl.TDIResult
    :members:
    :undoc-members:
    :show-inheritance:

Clustering options
~~~~~~~~~~~~~~~~~~

.. autoclass:: tcr_deep_insight.tl.FAISS_INDEX_BACKEND
    :members:
    :undoc-members:


Models
------

.. autofunction:: tcr_deep_insight.model.modeling_bert.get_human_config

.. autoclass:: tcr_deep_insight.model.TRabModelingBertForPseudoSequence
    :members:
    :undoc-members:
    :show-inheritance:
    
.. autoclass:: tcr_deep_insight.model.TRabModelingBertForVJCDR3
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: tcr_deep_insight.model.GEXModelingVAE
    :members:
    :undoc-members:
    :show-inheritance:

Model options 
~~~~~~~~~~~~~

.. autoclass:: tcr_deep_insight.model.TCR_BERT_ENCODING
    :members:
    :undoc-members:

.. autoclass:: tcr_deep_insight.model.TCR_BERT_POOLING
    :members:
    :undoc-members:

Training Utilities
------------------

.. autoclass:: tcr_deep_insight.model.TRabModelingBertForVJCDR3Trainer
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: tcr_deep_insight.model.TRabTokenizerForVJCDR3
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: tcr_deep_insight.model.TRabCollatorForVJCDR3
    :members:
    :undoc-members:
    :show-inheritance:


Plotting
--------

.. autofunction:: tcr_deep_insight.pl.set_plotting_params
.. autofunction:: tcr_deep_insight.pl.create_fig
.. autofunction:: tcr_deep_insight.pl.create_subplots
.. autofunction:: tcr_deep_insight.pl.plot_cdr3_sequence
.. autofunction:: tcr_deep_insight.pl.plot_selected_tcrs

Additional Utilities
--------------------

