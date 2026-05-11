#PBS -S /bin/bash
#PBS -N UC_DeLange_CD4
#PBS -l walltime=6:00:00
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

######## MODIFY THE FOLLOWING PATHS.

############ PARENT FOLDER WHERE COGS INPUT FILES SHOULD GO ############
DIR=~/HRJ_monocytes/hILCs/rCOGS_in/UC_DeLange_28067908
########################################################################

############ PARENT FOLDER WHERE COGS OUTPUT FILES SHOULD GO ###########
COGS_OUT=~/HRJ_monocytes/hILCs/COGS_results/COGS_out/UC_DeLange_28067908
########################################################################

################### NAME FOR THE OUTPUT FILES - THE EXPERIMENT SETTINGS ##################
MYNAME=UC_DeLange_28067908_CD4_SuSIE
##########################################################################################

############ LOCATION OF THE (MODIFIED) PM, SUITABLE FOR COGS INPUT ######################
PM_MOD=~/HRJ_monocytes/hILCs/rCOGS_in/Version3_revision2/peakmatrices/CD4_chicago_fres_5kb_abc_023_fres_extended_peakm_13012025_modified.txt
##########################################################################################

#################### THE LOCATION OF THE BAITMAP AND RMAP (AT FRAG RES) ###################
BAITS=~/spivakov/Design/Human_hg38_DpnII_75_1200/hg38_dpnII.baitmap
RMAP=~/spivakov/Design/Human_hg38_DpnII_75_1200/hg38_dpnII.rmap
###########################################################################################

################# Provide the SuSIE input file, if using ##################
SUSIE=~/HRJ_monocytes/external_data/gwas/SuSIE/UC_DeLange_28067908_1-hg38.csv
###########################################################################

#################### N CASES AND CONTROLS FOR THE GWAS ###################
# https://www.ebi.ac.uk/gwas/studies/GCST004133
NCASES=12366
NCONTROLS=33609
##########################################################################

### CHECK THE PATHS IN THE COMMANDS.

################ 1. Make COGS input files.
source activate DT_DPLYR
cd $SCRIPTS
./Make_rCOGS_input_files.sh \
        --inputDir ${DIR}/prelim_files_${MYNAME} \
        --outDir ${DIR}/COGS_input_${MYNAME} \
    	--pchic ${PM_MOD} \
        --baitmap ${BAITS} \
        --rmap ${RMAP} \
        --GRCh38 

conda deactivate

### NOT DOING - coding SNPs extraction
### We can copy paste the coding SNPs from the file previously run using the ILC3 pipeline.

### 2. Get coding SNPs included in the GWAS (here SuSIE) using VEP on GRCh38, using positions.
### Note that, if the GWAS had been included and contained rsids, there would be a list of these in the prelim folder.
### As it is, we need to extract the cols chr, pos from the SuSIE file and use as input to run VEP using singularity.
#cd ${DIR}/prelim_files_${MYNAME}

# extract a file of chr, pos. This line accounts for the fact that some SNPs have _ref_alt added after the position, but not all of them!
#cat ${SUSIE} | tail -n +2 | cut -d',' -f1 | awk -F '_' '{print $1}' | awk -F':' '{print $1 "\t" $2}' > SNP_positions.txt

# Run VEP using this file as input.
#source activate VEP
#cd ${DIR}/COGS_input_${MYNAME}
#~/HRJ_monocytes/hILCs/scripts/helen_scripts_for_rCOGS_in/run_vep.sh \
#	-v /rds/general/user/hrayjone/home/anaconda3/envs/VEP/bin/ensembl-vep/vep \
#	-i ${DIR}/prelim_files_${MYNAME} \
#	-a GRCh38 \
#	-o ${DIR}/COGS_input_${MYNAME} \
#	-m posVCF \
#	-p ${DIR}/prelim_files_${MYNAME}/SNP_positions.txt
#conda deactivate

cd ${DIR}/COGS_input_${MYNAME}
cp ${DIR}/COGS_input_UC_DeLange_28067908_ILC3_SuSIE/coding.format.txt ./

### 3. Now run rCOGS on: frag res, 5Kb, ABC, All, Vprom/coding.
source activate DT_DPLYR
cd ${DIR}/COGS_input_${MYNAME}

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

conda deactivate



