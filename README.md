## multiCOGS: COGS for multivariate finemapping data

These R scripts implement an updated version of the rCOGS package as used in Malysheva/Ray-Jones et al., bioRxiv 2022. 

### The key differences from the original rCOGS package (https://github.com/ollyburren/rCOGS) are as follows:

* The algorithm has been modified to enable use with multivariate fine-mapped GWAS data. Note that currently read_gwas() can only perform Wakefield synthesis based on a single causal variant assumption and multivariate fine-mapped data need to be provided directly to the compute_cogs() function

* All operations on genomic intervals are performed using data.tables directly. GenomicRanges are no longer used, avoiding the need for interconversion

* The number of promoter-proximal fragments included in the analysis can be specified by the user. Also note that make_vprom() is now called from within make_pchic() that now gains the vPromLen parameter

* Several minor bug fixes

Note that multiROGS R package is still in progress, and currently all scripts should be sourced directly into the analysis notebooks. 

The folders rCOGS_runs and rCOGS_scripts contain analysis code used in the Malysheva/Ray-Jones et al., bioRxiv 2022 that utilises multiCOGS functionality and can be seen as use examples for the new package.

### Explanation on the scripts for running multiCOGS in rCOGS_scripts

The scripts should be run in the following order. For full examples, see the wrapper scripts in [rCOGS_runs/multiCOGS_in_autoimmune_traits](/rCOGS_runs/multiCOGS_in_autoimmune_traits)

*Please see also original [rCOGS github repo](https://ollyburren.github.io/rCOGS/articles/Quickstart.html) for an explanation of the main input files for COGS.*

1. Make the COGS input files using [Make_rCOGS_input_files.sh](/rCOGS_scripts/Make_rCOGS_input_files.sh)
This script generates input directories containing the following formatted files needed to run COGS:
- Approx. LD independent region file, 
- MAF file,
- Restriction fragment digest file,
- pCHi-C design/annotation file

Note, [annot_CHi-C_files.R](/rCOGS_scripts/annot_CHi-C_files.R) will be run from within this script. Please ensure that annot_CHi-C_files.R is located within the same directory as Make_rCOGS_input_files.sh, and that R is available with libraries *argparser* and *data.table* installed.

The following inputs are needed to run Make_rCOGS_input_files.sh:
- --inputDir - full path to where you want preliminary files to go. These are interim files, not the final formatted files for COGS.
- --outDir - full path to where you want the final, formatted files for COGS to go.
- --pchic - full path to the pchic peak matrix formatted for COGS. The following column headers are required at this stage: 
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
- --baitmap - full path to the baitmap file containing captured fragments (column headers: )
- --rmap - full path to the digest file, i.e. fragments across the whole genome (column headers: chr, start, end, fragid)
- LD blocks (currently hardcoded!): these are files of approximately independent LD regions in the genome (column headers: chr, start, end). We use the files from Berisa and Pickrell:
    - [GRCh37](https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/fourier_ls-all.bed)
    - [GRCh38] (https://github.com/jmacdon/LDblocks_GRCh38)
- MAF file (currently hardcoded!): This file contains MAF for all SNPs in your population of interest. We used 1000 Genomes in European individuals (GRCh37 or GRCh38, as needed). Column headers: chr, pos, maf.
- TSS file (currently hardcoded!): This file, containing TSS locations (in GRCh37 or GRCh38, as needed), is used to annotated the peakmatrix. We used TSS information collated from Ensembl and Havana. Columns (in this order): ensg, genename, chr, TSS, strand, type

Explanation of the optional arguments for Make_rCOGS_input_files.sh:
- --GRCh37 - Flag that the required assembly is GRCh37 (this currently affects which of the harcoded LD block, MAF and TSS files are used)
- --GRCh38 - Flag that the required assembly is GRCh38 (this currently affects which of the harcoded LD block, MAF and TSS files are used)
- --gwas - The full path to the trait data (summary statistics file), obtained from the GWAS catalog, required if running univate (classic) COGS. Must have cols named chromosome, base_pair_location, p_value and variant_id. However, for multiCOGS we instead use a pre-made SuSIE table, which already contains SNPs and associated PPIs. In this case, a GWAS file does not need to be supplied here. The SuSIE table will be required instead in steps 2 and 3 below.
- --expandBaitmap - Flag to expand the baitmap from binned to fragment level. This is in case you want to analyse the data at fragment level resolution, but the other-ends in the peak matrix are at binned resolution. Note that, in our case, we already expanded our peakmatrices such that each row is at the level of a single fragment, so this flag was not needed in most of our analyses. If this flag is used, you must also supply binned rmap with solitary baits using --rmapSolBaits.
- --rmapSolBaits - Full path to the rmap in binned setting with solitary baits, if used.
- --tempdir - Directory for temporary files, default = current wd
- --help - prints help message

2. Make a list of coding SNPs from within your GWAS or SuSIE fine-mapped data using [run_vep.sh](/rCOGS_scripts/run_vep.sh)

This script is used to identify SNPs within protein-coding regions in your inpur GWAS or SuSIE fine-mapped data. You need a working version of the [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/script/index.html) to run this script.
 
Note! If you are running COGS using a GWAS file that includes rsids, the previous script should have generated a list of these IDs in the preliminary input files folder, ready to compare against VEP information. However, when running with SuSIE input files, we did not have rsids. Therefore, VEP is used to find coding SNPs based on chr, pos instead. The chr/pos information is first extracted from the SuSIE fine-mapped data as follows: 

Go to the preliminary files folder
`cd ${DIR}/prelim_files_${MYNAME}`

Extract a file of chr, pos from the SuSIE data. *Note, this line accounts for the fact that some SNPs have _ref_alt added after the position, but not all of them!*
`cat ${SUSIE} | tail -n +2 | cut -d',' -f1 | awk -F '_' '{print $1}' | awk -F':' '{print $1 "\t" $2}' > SNP_positions.txt`

Please see a [wrapper file](/rCOGS_runs/multiCOGS_in_autoimmune_traits/01_wrapper_ILCs_ASTAO_Ferreira_30929738.sh) for a full example of this.

The following inputs are needed to run run_vep.sh:
--vep               The full path to the vep executable"
    echo -e "  -i,  --inputdir          The full path to preliminary COGs input files. If you ran Make_rCOGS_input_files.sh there will already be a list of rsids in this same folder."
    echo -e "  -a,  --assembly          The required assembly: GRCh37 or GRCh38"
    echo -e "  -o,  --outdir            The full path to directory where desired COGS input file will go\n"
    echo -e "OPTIONAL:"
    echo -e "  -m,  --method            The method with which to run VEP. rsid: runs VEP using rsids. varVCF: gets vcf using variant IDs. posVCF: gets vcf using variant positions. posVCF requires --positions option. Default = rsids."
    echo -e "  -p,  --positions         Use with posVCF. file with chr, position."
    echo -e "  -h,  --help              Prints this help\n"


3. Run rCOGS (or multiCOGS) using [run_rCOGS.R](/rCOGS_scripts/run_rCOGS.R)







