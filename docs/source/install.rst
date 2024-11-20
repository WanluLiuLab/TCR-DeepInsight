Installation
------------

.. note::

   Hardware requirement for TCR-DeepInsight includes

   1. CPU: >= 1 cores. Recommended >= 8 cores for large-scale dataset
   2. RAM: >=1 Gb. Recommended >=64 Gb for large-scale dataset
   3. VRAM >= 1Gb of CUDA-enabled GPU. Recommended >= 8 Gb for large-scale dataset
   4. Disk space >= 1Gb. Recommended >= 100Gb for large-scale dataset


Operation System requirements for running TCR-DeepInsight include the installation of Python3 (Python3.9 used for development) and several PyPI packages. You can create a running environment using `conda`.


Install from PyPI
~~~~~~~~~~~~~~~~~

.. code-block:: shell
  :linenos:

    pip install tcr-deep-insight


Install from source
~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell
  :linenos:

    conda create -n tcr-deep-insight -f environment.yml
    conda activate tcr-deep-insight
    git clone git@github.com:WanluLiuLab/TCR-DeepInsight.git
    python3 setup.py install


Usage
-----

In IPython, simply import the package to get started:

.. code-block:: python
  :linenos:
    
    import tcr_deep_insight as tdi 
    
For more details, please refer to the :doc:`tutorial`.

Example dataset
---------------

.. warning::
   TCR-DeepInsight require AnnData objects as input.

Example processed datasets are available at `Zenodo <https://zenodo.org/records/12741480>`_.

or you can use the following code to download the example dataset:

.. code-block::
  :linenos:

  import tcr_deep_insight as tdi
  gex_adata = tdi.data.human_gex_reference_v2()
  tcr_adata = tdi.data.human_tcr_reference_v2()