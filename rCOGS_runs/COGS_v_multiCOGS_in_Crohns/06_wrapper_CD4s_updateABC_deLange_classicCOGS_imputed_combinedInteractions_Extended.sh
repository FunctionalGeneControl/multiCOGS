#PBS -S /bin/bash
#PBS -N CD4_classicImp
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=32:mem=62gb
#PBS -o /rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/hILCs/rCOGS_in/OU
#PBS -e /rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/hILCs/rCOGS_in/ER


###### Specify all input data and genome assembly.
###### Run for all required GWAS.
##### This is the folder where all required scripts can be found (Make_rCOGS_input_files.sh, run_rCOGS.R, and dependent scripts)
SCRIPTS=~/HRJ_monocytes/hILCs/scripts/helen_scripts_for_rCOGS_in
#####

##########################################################################################################################
### Requirements for rCOGS input files are described here: https://ollyburren.github.io/rCOGS/articles/Quickstart.html ###
##########################################################################################################################

### Our new interaction PM uses frag res, 5kb and ABC interactions. However, the first three cols correspond to fragment level ints. 
# The PM was modified in the first script in this folder to make it suitable to run with COGS scripts.

######## MODIFY THE FOLLOWING PATHS.

############ PARENT FOLDER WHERE COGS INPUT FILES SHOULD GO ############
DIR=~/HRJ_monocytes/hILCs/rCOGS_in/Version3_revision2
########################################################################

############ PARENT FOLDER WHERE COGS OUTPUT FILES SHOULD GO ###########
COGS_OUT=~/HRJ_monocytes/hILCs/COGS_results/COGS_out/Version3_revision2
########################################################################

################### NAME FOR THE OUTPUT FILES - THE EXPERIMENT SETTINGS ##################
MYNAME=revision_deLange_CD4s_hg38_classicCOGS_imputed_combinedInteractions_extended_ABC023
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

### 3. Now run rCOGS - on everything.

################# Provide the SuSIE input file, if using ##################
SUSIE=~/HRJ_monocytes/external_data/gwas/SuSIE/cd_for_mikhail_SuSIE_fix.csv
###########################################################################

### Note, the PM contains additional columns N_fres, N_5kb and N_abc. These should be ignored, so required feature names are given.
source activate DT_DPLYR
cd ${DIR}/COGS_input_${MYNAME}
Rscript ${SCRIPTS}/run_rCOGS.R \
        --ncases ${NCASES} \
        --ncontrols ${NCONTROLS} \
        --cogsIn ${DIR}/COGS_input_${MYNAME} \
        --assembly GRCh38 \
        --susie.single \
        --gwas ${SUSIE} \
        --cogsOut ${COGS_OUT}/${MYNAME} \
        --vProm 5 \
	--featureNames chicago_score_fres,chicago_score_5kb,ABC.Score,VProm,coding_snp
conda deactivate


