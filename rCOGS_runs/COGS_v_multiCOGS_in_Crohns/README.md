## COGS and multiCOGS runs in Crohn's Disease, using 3D chromatin data from hILC3s and CD4+ T cells

In these runs, we compared different iterations of COGS, running with Crohn's Disease GWAS data and pCHI-C/ABCC data in hILC3s and CD4+ T cells. The results are presented in our paper on Biorxiv [Malysheva/Ray-Jones/Lakes et al., 2026](https://www.biorxiv.org/content/10.1101/2022.10.19.512842v4)

The contents of the scripts are as follows:

[01_preparing_PMs_for_COGS.ipynb](https://pages.github.com/) - the peakmatrices (containing pCHi-C interactions and ABCC pairings) required some modification in order to be compatible with rCOGS. This script details these modifications, as well as some downstream comparisons post-COGS. Please note, the ILC3 and CD4+ T cell peakmatrices provided in our [OSF repository](https://osf.io/aq9fb/wiki?wiki=4j6ea) have already been formatted for COGS, so this step is not necessary for them.

[02_wrapper_ILCs_updateABC_deLange_classicCOGS_combinedInteractions_Extended.ipynb](https://pages.github.com/) - wrapper script for classic (univariate fine-mapping) COGS, run without summary statistics imputation, on Crohn's Disease with hILC3 data.

[03_wrapper_ILCs_updateABC_deLange_classicCOGS_imputed_combinedInteractions_Extended.ipynb](https://pages.github.com/) - wrapper script for classic (univariate fine-mapping) COGS, run with summary statistics imputation, on on Crohn's Disease with hILC3 data.

[04_wrapper_ILCs_updateABC_deLange_SuSIE_fix_combinedInteractions_Extended.ipynb](https://pages.github.com/) - wrapper script for multiCOGS (multivariate fine-mapping with SuSIE), which includes summary statistics imputation, on on Crohn's Disease with hILC3 data.

[05_wrapper_CD4s_updateABC_deLange_classicCOGS_combinedInteractions_Extended.ipynb](https://pages.github.com/) - wrapper script for classic (univariate fine-mapping) COGS, run without summary statistics imputation, on Crohn's Disease with CD4+ T cell data.

[06_wrapper_CD4s_updateABC_deLange_classicCOGS_imputed_combinedInteractions_Extended.ipynb](https://pages.github.com/) - wrapper script for classic (univariate fine-mapping) COGS, run with summary statistics imputation, on on Crohn's Disease with CD4+ T cell data.

[07_wrapper_CD4s_updateABC_deLange_SuSIE_fix_combinedInteractions_Extended.ipynb](https://pages.github.com/) - wrapper script for multiCOGS (multivariate fine-mapping with SuSIE), which includes summary statistics imputation, on on Crohn's Disease with CD4+ T cell data.
