#PBS -S /bin/bash
#PBS -N deLange_SuSIE_Jav_hg38
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=32:mem=96gb
#PBS -o ~/rCOGS_in/OU
#PBS -e ~/rCOGS_in/ER

###### Specify all input data and genome assembly.
###### Run for all required GWAS.

##########################################################################################################################
### Requirements for rCOGS input files are described here: https://ollyburren.github.io/rCOGS/articles/Quickstart.html ###
##########################################################################################################################

DIR=~/rCOGS_in
SCRIPTS=~/scripts/helen_scripts_for_rCOGS_in

################ hILCS, deLange, GRCh38
INDIR=${DIR}/prelim_files_deLange_Javierre_SuSIE_fix_hg38
OUTDIR=${DIR}/COGS_input_deLange_Javierre_SuSIE_fix_hg38

###### Now running having lifted over Javierre peak matrix to hg38.
###### THEN intersected the lifted over PM with the design files in Design/Human_hg38
###### I required a min overlap of 20 bp using genomic ranges. If a frag was split into two, we keep it.
###### See the script in ~/scripts/liftover_javierre_peakm_and_Design_files.R

###### Make input files and copy coding snps over from hILCs.
source activate DT_DPLYR
### Run this script from the scripts directory, because it needs to find the R script:
cd $SCRIPTS
./Make_rCOGS_input_files.sh \
        --inputDir ${INDIR} \
        --outDir ${OUTDIR} \
        --gwas ~/external_data/gwas/deLange/28067908-GCST004132-EFO_0000384.h.tsv \
        --pchic ~/external_data/Javierre/PCHiC_peak_matrix_cutoff5_hg38Liftover_overlapped.txt \
        --baitmap ~/Design/Human_hg38/Human_HindIII_hg38.baitmap \
        --rmap ~/Design/Human_hg38/Human_HindIII_hg38.rmap \
        --GRCh38

conda deactivate


### Copy over coding SNPs in hg38 from hILCs with SuSIE CD.
cd ${OUTDIR}
cp ~/rCOGS_in/COGS_input_deLange_ILCs_hg38_SuSIE_fix/coding.format.txt ./


######### RUNNING RCOGS

### Looking at VProm/coding will be the same across cell types, so we can just look in "all cells"

source activate DT_DPLYR
cd ${OUTDIR}
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
      	--susie \
      	--gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_allCells_VProm_only \
	--featureNames VProm
conda deactivate

source activate DT_DPLYR
cd ${OUTDIR}
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
      	--susie \
      	--gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_allCells_VProm_coding_only \
	--featureNames VProm,coding_snp

#### Run "all cells" and we will eventually store the combined dataset in here.
#cd ${OUTDIR}
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_allCells


### Then run for all cells one by one.
# Run for Mon
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_Mon \
        --featureNames Mon,VProm,coding_snp


# Run for Mac0
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_Mac0 \
        --featureNames Mac0,VProm,coding_snp

# Run for Mac1
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_Mac1 \
        --featureNames Mac1,VProm,coding_snp

# Run for Mac2
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_Mac2 \
        --featureNames Mac2,VProm,coding_snp

# Run for Neu
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_Neu \
        --featureNames Neu,VProm,coding_snp

# Run for EP
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_EP \
        --featureNames EP,VProm,coding_snp

# Run for Ery
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_Ery \
        --featureNames Ery,VProm,coding_snp

# Run for MK
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_MK \
        --featureNames MK,VProm,coding_snp

# Run for FoeT
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_FoeT \
        --featureNames FoeT,VProm,coding_snp

# Run for nCD4
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_nCD4 \
        --featureNames nCD4,VProm,coding_snp

# Run for tCD4
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_tCD4 \
        --featureNames tCD4,VProm,coding_snp

# Run for aCD4
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_aCD4 \
        --featureNames aCD4,VProm,coding_snp

# Run for naCD4
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_naCD4 \
        --featureNames naCD4,VProm,coding_snp

# Run for nCD8
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_nCD8 \
        --featureNames nCD8,VProm,coding_snp

# Run for tCD8
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_tCD8 \
        --featureNames tCD8,VProm,coding_snp

# Run for nB
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_nB \
        --featureNames nB,VProm,coding_snp

# Run for tB
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --susie \
        --gwas ~/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_Javierre_hg38_SuSIE_tB \
        --featureNames tB,VProm,coding_snp

conda deactivate




