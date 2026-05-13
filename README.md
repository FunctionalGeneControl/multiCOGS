# multiCOGS: COGS for multivariate finemapping data

These R scripts implement an updated version of the rCOGS package as used in Malysheva/Ray-Jones et al., bioRxiv 2022. 

## The key differences from the original rCOGS package (https://github.com/ollyburren/rCOGS) are as follows:

* The algorithm has been modified to enable use with multivariate fine-mapped GWAS data. Note that currently read_gwas() can only perform Wakefield synthesis based on a single causal variant assumption and multivariate fine-mapped data need to be provided directly to the compute_cogs() function

* All operations on genomic intervals are performed using data.tables directly. GenomicRanges are no longer used, avoiding the need for interconversion

* The number of promoter-proximal fragments included in the analysis can be specified by the user. Also note that make_vprom() is now called from within make_pchic() that now gains the vPromLen parameter

* Several minor bug fixes

Note that multiROGS R package is still in progress, and currently all scripts should be sourced directly into the analysis notebooks. 

The folders rCOGS_runs and rCOGS_scripts contain analysis code used in the Malysheva/Ray-Jones et al., bioRxiv 2022 that utilises multiCOGS functionality and can be seen as use examples for the new package.

## Explanation on the scripts for running multiCOGS in rCOGS_scripts

The scripts should be run in the following order. For **full examples** (wrapper scripts) and **data sources**, see [rCOGS_runs/multiCOGS_in_autoimmune_traits](/rCOGS_runs/multiCOGS_in_autoimmune_traits)

*Please see also original [rCOGS github repo](https://ollyburren.github.io/rCOGS/articles/Quickstart.html) for an explanation of the main input files for COGS.*

### 1. Make the COGS input files using [Make_rCOGS_input_files.sh](/rCOGS_scripts/Make_rCOGS_input_files.sh)

This script generates input directories containing the following formatted files needed to run COGS:
- Approx. LD independent region file
- MAF file
- Restriction fragment digest file
- pCHi-C design/annotation file

Note, [annot_CHi-C_files.R](/rCOGS_scripts/annot_CHi-C_files.R) will be run from within this script. Please ensure that annot_CHi-C_files.R is located within the same directory as Make_rCOGS_input_files.sh, and that R is available with libraries *argparser* and *data.table* installed.

**The following inputs are needed to run Make_rCOGS_input_files.sh:**
- `--inputDir` - full path to where you want preliminary files to go. These are interim files, not the final formatted files for COGS.
- `--outDir` - full path to where you want the final, formatted files for COGS to go.
- `--pchic` - full path to the pchic peak matrix formatted for COGS. The following column headers are required at this stage: 
    - baitChr - Bait (A captured restriction fragment) chromosome
    - baitStart - Bait fragment start
    - baitEnd - Bait fragment end
    - baitID - Bait fragment ID should match ID in digest file
    - baitName - Free text of what is captured - deprecated.
    - oeChr - 'Other end' (Promoter interacting restriction fragment) chromosome
    - oeStart - 'Other end' fragment start
    - oeEnd - 'Other end' fragment end
    - oeID - 'Other end' fragment ID should match ID in digest file
    - oeName - Free text of what is interacting - deprecated.
    - dist - Distance between bait and other end - deprecated.
    - Additional columns contain CHiCAGO scores for one or more analyses, e.g. "ILC3_chicago_score_fres" and "ILC3_ABCC_score"
- `--baitmap` - full path to the (fragment-level) baitmap file containing captured fragments (no header; columns contain: chr, start, end, fragid, fragname)
- `--rmap` - full path to the (fragment-level) digest file, i.e. fragments across the whole genome (no header; columns contain: chr, start, end, fragid)

Additionally:
- <ins>LD blocks (currently hardcoded!)</ins>: these are files of approximately independent LD regions in the genome (column headers: chr, start, end). We use the files from Berisa and Pickrell:
    - [GRCh37](https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/fourier_ls-all.bed)
    - [GRCh38](https://github.com/jmacdon/LDblocks_GRCh38)
- <ins>MAF file (currently hardcoded!)</ins>: This file contains MAF for all SNPs in your population of interest. We used 1000 Genomes in European individuals (GRCh37 or GRCh38, as needed). Column headers: chr, pos, maf.
- <ins>TSS file (currently hardcoded!)</ins>: This file, containing TSS locations (in GRCh37 or GRCh38, as needed), is used to annotated the peakmatrix with genes and biotypes (e.g. protein coding). We used TSS information collated from Ensembl and Havana. Columns (in this order): ensg, genename, chr, TSS, strand, type

**Explanation of the optional arguments for Make_rCOGS_input_files.sh:**
- `--GRCh37` - Flag that the required assembly is GRCh37 (this currently affects which of the harcoded LD block, MAF and TSS files are used)
- `--GRCh38` - Flag that the required assembly is GRCh38 (this currently affects which of the harcoded LD block, MAF and TSS files are used)
- `--gwas` - The full path to the trait data (summary statistics file), obtained from the GWAS catalog, required if running univate (classic) COGS. Must have cols named chromosome, base_pair_location, p_value and variant_id. However, for multiCOGS we instead use a pre-made table of variants fine-mapped with SuSIE, which already contains SNPs and associated posterior probabilities (PPIs). If you already have such fine-mapped data, a GWAS file does not need to be supplied here. The SNPs with their associated PPIs are used in steps 2 and 3 below.
- `--expandBaitmap` - Flag to expand the baitmap from binned to fragment level. This is in case you want to analyse the data at fragment level resolution, but the other-ends in the peak matrix are at binned resolution. Note that, in our case, we already expanded our peakmatrices so that each row is at the level of a single fragment, such that this flag was not needed. If this flag is used, you must also supply the binned rmap with solitary baits using `--rmapSolBaits`.
- `--rmapSolBaits` - Full path to the rmap in binned setting with solitary baits, if `--expandBaitmap` used.
- `--tempdir` - Directory for temporary files, default = current wd
- `--help` - prints help message

### 2. Make a list of coding SNPs from within your GWAS or SuSIE fine-mapped data using [run_vep.sh](/rCOGS_scripts/run_vep.sh)

This script is used to identify SNPs within protein-coding regions in your input GWAS or SuSIE fine-mapped data. You need a working version of the [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/script/index.html) to run this script.
 
Note! If you are running COGS using a GWAS file that includes rsids, the previous script should have generated a list of these IDs in the preliminary input files folder, ready to compare against VEP information. However, our fine-mapped input files did not have rsids. Therefore, we used VEP to find coding SNPs based on chr, pos instead. The chr/pos information is first extracted from the SuSIE fine-mapped data as follows: 

Go to the preliminary files folder

`cd ${DIR}/prelim_files_${MYNAME}`

Extract a file of chr, pos from the SuSIE data. *Note, this line accounts for the fact that some SNPs have _ref_alt added after the position, but not all of them!*

`cat ${SUSIE} | tail -n +2 | cut -d',' -f1 | awk -F '_' '{print $1}' | awk -F':' '{print $1 "\t" $2}' > SNP_positions.txt`

Please see a [wrapper file](/rCOGS_runs/multiCOGS_in_autoimmune_traits/01_wrapper_ILCs_ASTAO_Ferreira_30929738.sh) for a full example of this.

**The following inputs are needed to run run_vep.sh:**
- `--vep` - The full path to the vep executable
- `--inputdir` - The full path to preliminary COGs input files (i.e. the same inputDir as above, in Make_rCOGS_input_files.sh). Note, if you ran Make_rCOGS_input_files.sh with a GWAS summary statistics file, there should already be a list of rsids in this same folder.
- `--assembly` - The required assembly: GRCh37 or GRCh38 (will use the hardcoded SNP sites file in hg19 or hg38; see below for explanation on the sites file)
- `--outdir` - full path to where you want the final, formatted files for COGS to go (i.e. the same outDir as above, in Make_rCOGS_input_files.sh).

Additionally:
- <ins>SNP sites file (currently hardcoded!)</ins>: This is a SNP sites vcf that we obtained from 1000 Genomes (in hg19 or hg38) with the following columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO

**The following options can be used with run_vep.sh:**
- `--method` - The method with which to run VEP. rsid: runs VEP using rsids. varVCF: gets vcf using variant IDs. posVCF: gets vcf using variant positions. posVCF requires --positions option. Default = rsids.
- `--positions` - A file with chr, position. Needs to be supplied for method = posVCF. 
- `--help` - prints help message

The VEP script also edits the final coding SNPS file ("coding.txt") to make it compatible for rCOGS. Please double check that this formatting is working correctly in your case. The final coding SNPs file is named "coding.format.txt".

### 3. Run rCOGS (or multiCOGS) using [run_rCOGS.R](/rCOGS_scripts/run_rCOGS.R)
If steps 1 and 2 were followed, the input files should now be ready to run with COGS/multiCOGS. The run_COGS.R script utilises the scripts in the present folder:
- cogs.R
- data.R
- gwas.R
- pchic.R
- sCVPP.R

Please see the wrapper scripts in [COGS_v_multiCOGS_in_Crohns](/rCOGS_runs/COGS_v_multiCOGS_in_Crohns) for examples of how to use run_rCOGS.R in different fine-mapping settings (univariate vs multivariate) and with different features (e.g. chicago, ABCC, VProm, coding).

**The following inputs are needed to run run_COGS.R:**
- `--cogsIn` - This is the directory with the prepared files from steps 1 and 2 above, i.e. the `--outDir` from Make_rCOGS_input_files.sh and run_vep.sh above. In this folder, there should be the following files, which will be detected by run_COGS.R (these can also be supplied separately to the run_COGS.R script):
    - LD regions file with columns named: chr, start, end", default="_ld.format.bed$", can also be supplied with `--ld`
    - MAF file with columns named: chr, pos, maf.", default="formatted.maf.txt$", can also be supplied with `--maf`
    - GWAS data with columns named: chr, pos, p.", default="_gwas.format.txt$", can also be supplied with `--gwas`. NOTE - if using a mutivariate fine-mapped SuSIE table, supply the path to this file in place of the default. You also need to use the `--susie` flag (optional argument, below). See a wrapper script such as [01_wrapper_ILCs_ASTAO_Ferreira_30929738.sh](/rCOGS_runs/multiCOGS_in_autoimmune_traits/01_wrapper_ILCs_ASTAO_Ferreira_30929738.sh) as an example.
    - Formatted peak matrix with biotypes", default = "_pm.format.txt$", can also be supplied with `--pmFormat`
    - RMAP with columns named: chr, start, end, fragid", default = ".rmap_wHeader.txt$", can also be supplied with `--rmap`
    - Annotated baits with columns named: fragid, ensg and biotype", default = "PCHiC_design_annotation_plusUnbaited_with_geneType.txt", can also be supplied with `--bannot`
    - Coding SNPs with columns named: chr, pos, ensg", default = "coding.format.txt", can also be supplied with `--coding`
- `--ncases` - The number of cases in the GWAS
- `--ncontrols` - The number of controls in the GWAS
- `--assembly` - Assembly as either GRCh37 or GRCh38, this is used to remove the MHC region. GRCh37 region defined as in rCOGS vignette: https://ollyburren.github.io/rCOGS/articles/Quickstart.html; GRCh38 region defined as in NCBI: https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13  
- `--cogsOut` - Directory for COGS output files

Additionally:
- <ins>TSS/gene information (currently hardcoded!)</ins> - The output ranked gene list of COGS has ENSG IDs but not gene names. The supplied gene information is used to further annotated the output COGS results. This table requires the column "ensg" and any additional annotation columns, such as gene name and TSS.

**Explanation of the optional arguments for run_COGS.R:**
- `--featureNames` - This is a comma separated list of which score columns COGS should consider. These columns should have been included in the input peakmatrix and can include, for example, chicago_score or ABCC_score. If using specific features, make sure to add the "VProm" and "coding_snp" features. Default is that all features are used, i.e. all score columns, plus VProm and coding_snp. For examples of runs using different features, see a wrapper script such as [04_wrapper_ILCs_updateABC_deLange_SuSIE_fix_combinedInteractions_Extended.sh](/rCOGS_runs/COGS_v_multiCOGS_in_Crohns/04_wrapper_ILCs_updateABC_deLange_SuSIE_fix_combinedInteractions_Extended.sh)
- `--vProm` - This is the number of fragments to use when creating virtual promoter regions, default = 1. Note, we used 5 fragments in our analyses.
- `--chicThresh` - The hard threshold for CHiC interactions. Scores will only be considered ABOVE this value. This is why we set all ABCC scores to 5.1 in our input peakmatrices.
- `--susie` - flag to run on the SuSIE (multiCOGS) setting. If so, provide the SuSIE .csv file, containing PPIs, in place of the GWAS file above. Do not need to supply LD blocks. Please note, for SNPs where SuSIE's data aren't available or were filtered out, this script will currently ensure that rCOGS uses the single.pp value from a single causal variant model (Wakefield synthesis), provided in a separate column of the SuSIE input file. Please see the input SuSIE files included in our paper: **LINK TO PROVIDE**







