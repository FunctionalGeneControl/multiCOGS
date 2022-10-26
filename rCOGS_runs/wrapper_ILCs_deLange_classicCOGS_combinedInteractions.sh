#PBS -S /bin/bash
#PBS -N run_rCOGS_deLange_SuSIE
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=32:mem=62gb
#PBS -o ~/rCOGS_in/OU
#PBS -e ~/rCOGS_in/ER


###### Specify all input data and genome assembly.
###### Run for all required GWAS.

##########################################################################################################################
### Requirements for rCOGS input files are described here: https://ollyburren.github.io/rCOGS/articles/Quickstart.html ###
##########################################################################################################################

DIR=~/rCOGS_in
SCRIPTS=~/scripts/helen_scripts_for_rCOGS_in

################ hILCS, deLange, GRCh38, SuSIE, frag res and 5Kb and ABC interactions.
INDIR=${DIR}/prelim_files_deLange_ILCs_hg38_classicCOGS_combinedInteractions
OUTDIR=${DIR}/COGS_input_deLange_ILCs_hg38_classicCOGS_combinedInteractions

### Our new interaction PM uses frag res, 5kb and ABC interactions.

source activate DT_DPLYR
cd $SCRIPTS
./Make_rCOGS_input_files.sh \
        --inputDir ${INDIR} \
        --outDir ${OUTDIR} \
        --gwas ~/external_data/gwas/deLange/28067908-GCST004132-EFO_0000384.h.tsv \
	--pchic ~/chicago/ILC3_chicago_fres_5kb_abc_fres_consensus_peakm.txt \
        --baitmap ~/Design/Human_hg38_DpnII_75_1200/hg38_dpnII.baitmap \
        --rmap ~/Design/Human_hg38_DpnII_75_1200/hg38_dpnII.rmap \
        --GRCh38 

conda deactivate

### Run VEP
cd ${OUTDIR}
source activate VEP
./run_vep.sh \
        -v ~/anaconda3/envs/VEP/bin/ensembl-vep/vep \
        -i ${INDIR} \
        -a GRCh38 \
        -m posVCF \
        -p ~/external_data/gwas/deLange/28067908-GCST004132-EFO_0000384.sites.txt \
        -o ${OUTDIR}
conda deactivate

# Now run rCOGS on: frag res, 5Kb, ABC, All, Vprom/coding.

source activate DT_DPLYR
cd ${OUTDIR}
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_ILCs_hg38_ClassicCOGS_combinedInteractions_ALL \
        --vProm 5 

#### Run for 5Kb interactions only
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_ILCs_hg38_ClassicCOGS_combinedInteractions_5KbOnly \
        --vProm 5 \
	--featureNames chicago_5kb

#### Run for fres only
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_ILCs_hg38_ClassicCOGS_combinedInteractions_fresOnly \
        --vProm 5 \
        --featureNames chicago_fres

#### Run for ABC only
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_ILCs_hg38_ClassicCOGS_combinedInteractions_ABCOnly \
        --vProm 5 \
        --featureNames ABC.Score

#### Run for VProm and coding only
Rscript ${SCRIPTS}/run_Mikhails_rCOGS.R \
        --ncases 12194 \
        --ncontrols 28072 \
        --cogsIn ${OUTDIR} \
        --assembly GRCh38 \
        --cogsOut ~/COGS_results/COGS_out/CD_deLange_ILCs_hg38_ClassicCOGS_combinedInteractions_VProm_coding_only \
        --vProm 5 \
        --featureNames VProm,coding_snp

conda deactivate



