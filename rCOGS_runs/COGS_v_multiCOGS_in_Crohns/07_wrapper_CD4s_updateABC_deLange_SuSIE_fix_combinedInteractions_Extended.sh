#PBS -S /bin/bash
#PBS -N SuSIE_CD4
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=32:mem=62gb
#PBS -o /rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/hILCs/rCOGS_in/OU
#PBS -e /rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/hILCs/rCOGS_in/ER

###### Specify all input data and genome assembly.
###### Run for all required GWAS.
###### We've discussed ABC thresholds now and have decided to use 0.023 for both ILC and CD4, and use ABC in 5kb bins only.
###### This is because:
###### (a) at this cutoff, both ILC and CD4 are close to max cor(ABCnumerator, GE)
###### (b) both datasets have similar numbers of peaks
###### (c) for 5kb, the absolute values of cor(ABCnumerator, GE) are about twice as high as for single frag.


##### This is the folder where all required scripts can be found (Make_rCOGS_input_files.sh, run_rCOGS.R, and dependent scripts)
SCRIPTS=~/HRJ_monocytes/hILCs/scripts/helen_scripts_for_rCOGS_in
#####

##########################################################################################################################
### Requirements for rCOGS input files are described here: https://ollyburren.github.io/rCOGS/articles/Quickstart.html ###
##########################################################################################################################

### Our new interaction PM uses frag res, 5kb and ABC interactions. However, the first three cols correspond to fragment level ints. 
### Only using the extended PM.
### Now using the peakmatrix with increased ABC cutoff. But there are more ABC interactions than before (other improvements were made to the ABC pipeline.)
# the extended (hopefully used for the final paper version): ~/spivakov/miniPCHiC/hILCs/ILC3/PCHiC/data/ILC3_chicago_fres_bin_5kb_abc_023_fres_extended_peakm_13012025.txt
# The PM was modified in the first script in this folder to make it suitable to run with COGS scripts.

######## MODIFY THE FOLLOWING PATHS.

############ PARENT FOLDER WHERE COGS INPUT FILES SHOULD GO ############
DIR=~/HRJ_monocytes/hILCs/rCOGS_in/Version3_revision2
########################################################################

############ PARENT FOLDER WHERE COGS OUTPUT FILES SHOULD GO ###########
COGS_OUT=~/HRJ_monocytes/hILCs/COGS_results/COGS_out/Version3_revision2
########################################################################

################### NAME FOR THE OUTPUT FILES - THE EXPERIMENT SETTINGS ##################
MYNAME=revision_deLange_CD4s_hg38_SuSIE_fix_combinedInteractions_extended_ABC023
##########################################################################################

############ LOCATION OF THE (MODIFIED) PM, SUITABLE FOR COGS INPUT ######################
PM_MOD=${DIR}/peakmatrices/CD4_chicago_fres_5kb_abc_023_fres_extended_peakm_13012025_modified.txt
##########################################################################################

############################## THE LOCATION OF THE GWAS FILE #############################
GWAS=~/HRJ_monocytes/external_data/gwas/deLange/28067908-GCST004132-EFO_0000384.h.tsv
##########################################################################################

#################### THE LOCATION OF THE BAITMAP AND RMAP (AT FRAG RES) ###################
BAITS=~/spivakov/Design/Human_hg38_DpnII_75_1200/hg38_dpnII.baitmap
RMAP=~/spivakov/Design/Human_hg38_DpnII_75_1200/hg38_dpnII.rmap
###########################################################################################

#################### N CASES AND CONTROLS FOR THE GWAS ###################
NCASES=12194
NCONTROLS=28072
##########################################################################

### CHECK THE PATHS IN THE COMMANDS.

################ 1. Make COGS input files.
source activate DT_DPLYR
cd $SCRIPTS
./Make_rCOGS_input_files.sh \
        --inputDir ${DIR}/prelim_files_${MYNAME} \
        --outDir ${DIR}/COGS_input_${MYNAME} \
        --gwas ${GWAS} \
	--pchic ${PM_MOD} \
        --baitmap ${BAITS} \
        --rmap ${RMAP} \
        --GRCh38 

conda deactivate

### 2. Copy SNPs already run using VEP on GRCh38, using positions.
cd ${DIR}/COGS_input_${MYNAME}
cp ~/HRJ_monocytes/hILCs/rCOGS_in/Archive/COGS_input_deLange_ILCs_hg38_newVEP_pos/coding.format.txt ./

### 3. Now run rCOGS on: frag res, 5Kb, ABC, All, Vprom/coding.
source activate DT_DPLYR
cd ${DIR}/COGS_input_${MYNAME}

################# Provide the SuSIE input file, if using ##################
SUSIE=~/HRJ_monocytes/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv
###########################################################################

#### Here the SuSIE input file is used, in place of GWAS. 
#### Run on everything
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME} \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
        --featureNames chicago_score_fres,chicago_score_5kb,ABC.Score,VProm,coding_snp

##### Run for all CHiC interactions and ABC (no VProm or coding)
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME}_CHiCandABCOnly \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
        --featureNames chicago_score_fres,chicago_score_5kb,ABC.Score
# I want to double check that the result of this run is the same as putting all chicago interactions and ABC into one column.

##### Run for 5Kb interactions only
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME}_5kbOnly \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
        --featureNames chicago_score_5kb


##### Run for fres only
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME}_fresOnly \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
        --featureNames chicago_score_fres
#
##### Run for ABC only
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME}_ABConly \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
        --featureNames ABC.Score
#
##### Run for VProm and coding only
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME}_VPromCodingOnly \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
        --featureNames VProm,coding_snp
#
##### Run for everything except coding
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME}_NoCoding \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
        --featureNames VProm,chicago_score_fres,chicago_score_5kb,ABC.Score
#
##### Run for everything except ABC
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME}_NoABC \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
	--featureNames VProm,coding_snp,chicago_score_fres,chicago_score_5kb
#
#### Run for VProm only
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME}_VPromOnly \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
        --featureNames VProm
#
##### Run for coding only
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME}_CodingOnly \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
        --featureNames coding_snp

##### Run for all CHiC interactions
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --cogsOut ${COGS_OUT}/${MYNAME}_CHiCOnly \
        --susie \
        --gwas ${SUSIE} \
        --vProm 5 \
        --featureNames chicago_score_fres,chicago_score_5kb

conda deactivate


