## multiCOGS runs in autoimmune traits

In these runs, we ran multiCOGS for autoimmune traits, using pCHI-C/ABCC data in hILC3s and CD4+ T cells. The results are presented in our paper on Biorxiv [Malysheva/Ray-Jones/Lakes et al., 2026](https://www.biorxiv.org/content/10.1101/2022.10.19.512842v4)

### Data sources

- The peakmatrices for these runs are available as files Data_S5 (ILC3s) and Data_S9 (CD4+ T cells) in our [OSF repository](https://osf.io/aq9fb/overview)
- The SuSIE fine-mapped data are available at: [XXX]()
- The rmap and baitmap (CHi-C design) at DpnII-fragment level are available in Data_S12 in our [OSF repository](https://osf.io/aq9fb/files/osfstorage)


### The following traits were analysed in the wrapper scripts:

- Asthma, adult onset (ASTAO) (PMID [30929738](https://pubmed.ncbi.nlm.nih.gov/30929738/))
  - multiCOGS in ILC3s: [01_wrapper_ILCs_ASTAO_Ferreira_30929738.sh](/rCOGS_runs/multiCOGS_in_6_autoimmune_traits/01_wrapper_ILCs_ASTAO_Ferreira_30929738.sh)
  - multiCOGS in CD4+ T cells: [07_wrapper_CD4s_ASTAO_Ferreira_30929738.sh](/rCOGS_runs/multiCOGS_in_6_autoimmune_traits/07_wrapper_CD4s_ASTAO_Ferreira_30929738.sh)

- Ulcerative Colitis (UC) (PMID [28067908](https://pubmed.ncbi.nlm.nih.gov/28067908/))
  - multiCOGS in ILC3s: [02_wrapper_ILCs_UC_DeLange_28067908.sh](/rCOGS_runs/multiCOGS_in_autoimmune_traits/02_wrapper_ILCs_UC_DeLange_28067908.sh)
  - multiCOGS in CD4+ T cells: [07_wrapper_CD4s_UC_DeLange_28067908.sh](/rCOGS_runs/multiCOGS_in_autoimmune_traits/06_wrapper_CD4s_ASTAO_Ferreira_30929738.sh)

- Inflammatory Bowel Disease (IBD) (PMID [28067908](https://pubmed.ncbi.nlm.nih.gov/28067908/))
  - multiCOGS in ILC3s: [03_wrapper_ILCs_IBD_DeLange_28067908](/rCOGS_runs/multiCOGS_in_autoimmune_traits/03_wrapper_ILCs_IBD_DeLange_28067908.sh)
  - multiCOGS in CD4+ T cells: [08_wrapper_CD4s_IBD_DeLange_28067908.sh](/rCOGS_runs/multiCOGS_in_autoimmune_traits/08_wrapper_CD4s_IBD_DeLange_28067908.sh)

- Celiac Disease (CEL) (PMID [20190752](https://pubmed.ncbi.nlm.nih.gov/20190752/))
  - multiCOGS in ILC3s: [04_wrapper_ILCs_CEL_Dubois_20190752.sh](/rCOGS_runs/multiCOGS_in_autoimmune_traits/04_wrapper_ILCs_CEL_Dubois_20190752.sh)
  - multiCOGS in CD4+ T cells: [09_wrapper_CD4s_CEL_Dubois_20190752](/rCOGS_runs/multiCOGS_in_autoimmune_traits/09_wrapper_CD4s_CEL_Dubois_20190752.sh)

- Primary Sclerosing Cholangitis (PSC) (PMID [27992413](https://pubmed.ncbi.nlm.nih.gov/27992413/))
  - multiCOGS in ILC3s: [05_wrapper_ILCs_PSC_Ji_27992413.sh](/rCOGS_runs/multiCOGS_in_autoimmune_traits/05_wrapper_ILCs_PSC_Ji_27992413.sh)
  - multiCOGS in CD4+ T cells: [10_wrapper_CD4s_PSC_Ji_27992413.sh](/rCOGS_runs/multiCOGS_in_autoimmune_traits/10_wrapper_CD4s_PSC_Ji_27992413.sh)
    


