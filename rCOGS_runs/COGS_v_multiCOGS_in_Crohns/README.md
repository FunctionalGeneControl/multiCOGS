## COGS and multiCOGS runs in Crohn's Disease (hILC3s and CD4+ T cells)

In these runs, we compared different iterations of COGS, running with Crohn's Disease GWAS data and pCHI-C/ABCC data in hILC3s and CD4+ T cells. The results are presented in our paper on Biorxiv [Malysheva/Ray-Jones/Lakes et al., 2026](https://www.biorxiv.org/content/10.1101/2022.10.19.512842v4)

### Data sources

The input data can be found in our large supplementary data files in our [OSF repository](https://osf.io/aq9fb/overview):
- The peakmatrices for these runs are available as files Data_S5 (ILC3s) and Data_S9 (CD4+ T cells).
- The rmap and baitmap (CHi-C design), here used at DpnII-fragment level, are found in Data_S12.
- The SuSIE fine-mapped files are available in Data_S13.

### Description of the scripts

The contents of the scripts are as follows:

[01_preparing_PMs_for_COGS.ipynb](/rCOGS_runs/COGS_v_multiCOGS_in_Crohns/01_preparing_PMs_for_COGS.ipynb) - the peakmatrices (containing pCHi-C interactions and ABCC pairings) required some modification in order to be compatible with rCOGS. This script details these modifications. 
*Please note, the ILC3 and CD4+ T cell peakmatrices provided in our [OSF repository](https://osf.io/aq9fb/wiki?wiki=4j6ea) have already been formatted for COGS, so this step is not necessary for them.*

[02_wrapper_ILCs_updateABC_deLange_classicCOGS_combinedInteractions_Extended.sh](/rCOGS_runs/COGS_v_multiCOGS_in_Crohns/02_wrapper_ILCs_updateABC_deLange_classicCOGS_combinedInteractions_Extended.sh) - wrapper script for classic (univariate fine-mapping) COGS, run without summary statistics imputation, on Crohn's Disease with hILC3 data.

[03_wrapper_ILCs_updateABC_deLange_classicCOGS_imputed_combinedInteractions_Extended.sh](/rCOGS_runs/COGS_v_multiCOGS_in_Crohns/03_wrapper_ILCs_updateABC_deLange_classicCOGS_imputed_combinedInteractions_Extended.sh) - wrapper script for classic (univariate fine-mapping) COGS, run with summary statistics imputation, on on Crohn's Disease with hILC3 data.

[04_wrapper_ILCs_updateABC_deLange_SuSIE_fix_combinedInteractions_Extended.sh](/rCOGS_runs/COGS_v_multiCOGS_in_Crohns/04_wrapper_ILCs_updateABC_deLange_SuSIE_fix_combinedInteractions_Extended.sh) - wrapper script for multiCOGS (multivariate fine-mapping with SuSIE), which includes summary statistics imputation, on on Crohn's Disease with hILC3 data. As well as running multiCOGS with all annotation features, we also run it separately on each feature for comparison (CHi-C at fragment-level or 5kb interactions, ABCC, VProm and coding SNPs).

[05_wrapper_CD4s_updateABC_deLange_classicCOGS_combinedInteractions_Extended.sh](/rCOGS_runs/COGS_v_multiCOGS_in_Crohns/05_wrapper_CD4s_updateABC_deLange_classicCOGS_combinedInteractions_Extended.sh) - wrapper script for classic (univariate fine-mapping) COGS, run without summary statistics imputation, on Crohn's Disease with CD4+ T cell data.

[06_wrapper_CD4s_updateABC_deLange_classicCOGS_imputed_combinedInteractions_Extended.sh](/rCOGS_runs/COGS_v_multiCOGS_in_Crohns/06_wrapper_CD4s_updateABC_deLange_classicCOGS_imputed_combinedInteractions_Extended.sh) - wrapper script for classic (univariate fine-mapping) COGS, run with summary statistics imputation, on on Crohn's Disease with CD4+ T cell data.

[07_wrapper_CD4s_updateABC_deLange_SuSIE_fix_combinedInteractions_Extended.sh](/rCOGS_runs/COGS_v_multiCOGS_in_Crohns/07_wrapper_CD4s_updateABC_deLange_SuSIE_fix_combinedInteractions_Extended.sh) - wrapper script for multiCOGS (multivariate fine-mapping with SuSIE), which includes summary statistics imputation, on on Crohn's Disease with CD4+ T cell data. As well as running multiCOGS with all annotation features, we also run it separately on each feature for comparison (CHi-C at fragment-level or 5kb interactions, ABCC, VProm and coding SNPs).

[08_downstream_checks_on_COGS.ipynb](/rCOGS_runs/COGS_v_multiCOGS_in_Crohns/08_downstream_checks_on_COGS.ipynb) - code during the project for checking the input peakmatrices and comparing COGS runs downstream.
