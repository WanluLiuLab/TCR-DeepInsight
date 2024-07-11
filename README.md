# TCR-DeepInsight

## Aims

We have previously developed [huARdb](https://huarc.net/database) and an [updated version (v2)](https://huarc.net/v2/) currently in developing which collects single-cells immune profiling datasets including linked transcriptome and full-length TCR informations. However, one of the main obstacles in using single-cell immune profiling datasets for disease immunotherapy is the absence of a convenient reference atlas to access this information. Despite the growing potential of TCR engineering in this area, there is a significant challenge in identifying functional TCRs due to the lack of efficient computational tools to integrate the data. Nonetheless, these datasets offer a valuable resource to identify such TCRs and further advance TCR-T technology.

<img src="./imgs/img1.png" alt="TCRDeepInsight" style="zoom:150%;" />

### Detailed collection of datasets

We aims to build the most comprehensive atlas containing matched transcriptome and full-length V(D)J sequence of T and B cells in **Cancer**, **Autoimmune diseases**, and **Infections**. 


In 2023, we make the first release of an integrated datasets containing more than 1,000,000 hcT cells with full length TCR sequence, including the following studies. 



| **Study name**            | **Number of T cells** | **Disease**                                                  |
| ------------------------- | --------------------- | ------------------------------------------------------------ |
| [Abbas et al., 2021](https://doi.org/10.1038/s41467-021-26282-z)        | 20302                 | Acute Myeloid Leukimia                                       |
| [Azizi et al., 2018](https://doi.org/10.1016/j.cell.2018.05.060)       | 20851                 | Breast Cancer                                                |
| [Bacher et al., 2020](https://doi.org/10.1016/j.immuni.2020.11.016.)      | 39950                 | COVID-19                                                     |
| [Boland et al., 2020](https://doi.org/10.1126/sciimmunol.abb4432)       | 84076                 | Ulcerative colitis and Healthy controls                      |
| [Borcherding et al., 2021](https://doi.org/10.1038/s42003-020-01625-6)  | 8494                  | Clear cell renal cell carcinoma                              |
| [Cheon et al., 2021](https://doi.org/10.1126/sciimmunol.abk1741)      | 16271                 | COVID-19                                                     |
| [Corridoni et al., 2020](https://doi.org/10.1038/s41591-020-1003-4)    | 6915                  | Ulcerative colitis and Healthy controls                      |
| [Gao et al., 2020](https://doi.org/10.1038/s41467-022-29175-x)          | 210216                | Large granular lymphocyte leukemia and healthy controls      |
| [Gate et al., 2020](https://doi.org/10.1038/s41586-019-1895-7)         | 871                   | Healthy controls                                             |
| [He et al., 2020](https://doi.org/10.1186/s13059-020-02210-14)       | 34058                 | Healthy                                                      |
| [Kim et al., 2022](https://doi.org/10.1038/s41467-022-29539-3)          | 51343                 | Checkpoint inhibitor associated arthritis                    |
| [Krishna et al., 2021](https://doi.org/10.1016/j.ccell.2021.03.007)      | 3635                  | Clear cell renal cell carcinoma                              |
| [Liao et al., 2020](https://doi.org/10.1038/s41591-020-0901-9)         | 76477                 | COVID-19                                                     |
| [Liu et al., 2021](https://doi.org/10.1038/s41467-021-21043-4)         | 34767                 | Nasopharyngeal carcinoma                                     |
| [Lu et al., 2019](https://doi.org/10.1038/s41467-022-29539-3)           | 57220                 | Metastatic colorectal cancer                                 |
| [Luoma et al., 2020](https://doi.org/10.1016/j.cell.2020.06.001)        | 4201                  | Checkpoint inhibitor associated colitis and no colitis       |
| [Mahuron et al., 2020](https://doi.org/10.1084/jem.20192080)      | 4315                  | Melanoma                                                     |
| [Neal et al., 2018](https://doi.org/10.1016/j.cell.2018.11.021)         | 31568                 | Clear cell renal cell carcinoma                              |
| [Notarbartolo et al., 2021](https://doi.org/10.1126/sciimmunol.abg502) | 36558                 | COVID-19 and healthy controls                                |
| [Penkava et al., 2020](https://doi.org/10.1038/s41467-020-18513-6)      | 56281                 | Psoriatic arthritis                                          |
| [Ramaswamy et al., 2021](https://doi.org/10.1016/j.immuni.2021.04.003)    | 29207                 | SARS-CoV-2-associated  multisystem inflammatory syndrome and healthy controls |
| [Simone et al., 2021](https://doi.org/10.1038/s42003-021-02931-3)       | 4633                  | Ankylosing spondylitis                                       |
| [Suo et al., 2022](https://doi.org/10.1126/science.abo0516)          | 40115                 | Healthy                                                      |
| [Wang et al., 2021](https://doi.org/10.1038/s41467-021-25771-5)         | 54474                 | Kawasaki disease and healthy controls                        |
| [Wang et al., 2022](https://doi.org/10.3389/fimmu.2022.812514)         | 21664                 | COVID-19 and healthy controls                                |
| [Wen et al., 2020](https://doi.org/10.1038/s41421-020-0168-9)          | 44855                 | COVID-19 and healthy controls                                |
| [Yost et al., 2019](https://doi.org/10.1038/s41591-019-0522-3)         | 24560                 | Basal cell carinoma and squamous cell carcinoma              |
| [Zheng et al., 2020](https://doi.org/10.1038/s41467-020-20019-0)        | 20302                 | Esophagus squamous cell carcinoma                            |

### Study/Dataset being processed for the database and reference dataset

| **Study name**            | **Number of T cells** | **Disease**                                                  |
| ------------------------- | --------------------- | ------------------------------------------------------------ |
| [Minervina et al.](https://doi.org/10.1038/s41590-022-01184-4) | TBD | COVID-19 |
| [Poon et al., 2023](https://doi.org/10.1038/s41590-022-01395-9) | TBD | Healthy Barrier Sites |
| [Bieberich et al., 2021](https://www.frontiersin.org/articles/10.3389/fimmu.2021.701085/full) | TBD | COVID-19  |
| [Schalck et al.](https://doi.org/10.1158/2159-8290.CD-21-1248) | TBD | Pancreatic Cancer |
| [Eberhardt et al.]() | TBD | Head and neck squamous cell carcinomas  |
| [Leader et al.]() | TBD | Non-Small Cell Lung Cancer |
| [Ren et al.]() | TBD | Ovarian Cancer |
| [Shi et al.]() | TBD | Biliary Cancer | 
| [Tong et al.](https://doi.org/10.1038/s41467-022-34581-2) | TBD | Breast Cancer |
| [Pai et al., 2023](https://doi.org/10.1016/j.ccell.2023.03.009) | TBD | Multiple solid tumors |
| [Rahim et al., 2023](https://doi.org/10.1016/j.cell.2023.02.021) | TBD | HNSCC |
| [Ogino et al., 2022](https://www.jci.org/articles/view/151239) | TBD | Low-grade Glioma |
| [Suo et al., 2023](10.1038/s41587-023-01734-7) | TBD | Nonsense-mediated decay (NMD)  |
| [Moon et al., 2023](10.1038/s41467-022-35264-8) | TBD | Rheumatoid arthritis | TBD |
| [Zeng et al., 2023](10.1016/j.chom.2023.02.001) | TBD | Alcohol-associated liver disease | TBD |
| [Xiao et al., 2023](https://www.nature.com/articles/s43587-023-00379-0) | TBD | COVID-19 | 
| [Sonigra et al., 2023](https://insight.jci.org/articles/view/160964) | TBD | Rheumatoid arthritis | 
| [Jiang et al., 2021](https://insight.jci.org/articles/view/148035) | TBD  | skin lesion erythema migrans |
| [Buggert et al., 2022](10.1016/j.cell.2020.11.019) | TBD | HIV |
| [Xu et al., 2023](https://www.nature.com/articles/s41590-022-01367-z) | TBD | COVID-19 |
| [Ogino et al., 2022](https://www.jci.org/articles/view/151239) | TBD | Low-grade Glioma |
| [Friedrich et al., 2023](https://doi.org/10.1016/j.ccell.2023.02.008) | TBD | Multiple Myeloma |
| [Ali et al., 2023](https://rupress.org/jem/article/220/4/e20220729/213819/PD-1-blockade-and-CDK4-6-inhibition-augment) |  TBD | Breast and Ovarian Cancer |


## Introduction to TCR-DeepInsight

To robustly identify potential disease associated TCRα/β pairs considering both TCR sequence similarity and transcriptome features from million-level paired TCRα/β repertoire, we developed a deep-learning based framework named TCR-DeepInsight. 


<img src="./imgs/TCRDeepInsight.png" alt="TCRDeepInsight" style="zoom:150%;" />

## Installation

The installation of TCR-DeepInsight needs about 5 minutes.


**Hardware requirement for TCR-DeepInsight includes**
1. RAM: >16Gb for larger dataset
2. VRAM of CUDA-enabled GPU: >8Gb 


Operation System requirements for running TCR-DeepInsight include the installation of Python3 (Python3.8 used for development) and several PyPI packages. You can create a running environment using CONDA ([Anaconda](https://www.anaconda.com/download#Downloads) or [Miniconda](https://docs.conda.io/en/main/miniconda.html))

```shell
conda create -n tcr-deep-insight -f environment.yml
conda activate tcr-deep-insight
git clone git@github.com:WanluLiuLab/TCR-DeepInsight.git
cd TCR-DeepInsight
```

In IPython, simply import the package to get started:


```python
import tcr_deep_insight as tdi 
```

**RAPIDS**

See [RAPIDS] https://docs.rapids.ai/install if using cuml in TDI package.


## Usage


For detailed usage, Please see the [Documentation for TCR-DeepInsight](https://tcr-deepinsight.readthedocs.io/en/latest/index.html).


The entire process for CD8Mapper will take less than 30 minutes on a dataset of approximately 50,000 cells. Please see [Tutorial for CD8Mapper](https://tcr-deepinsight.readthedocs.io/en/latest/notebooks/huARdb_CD8_Mapper.html) for more details. 