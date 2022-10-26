
#!/bin/bash

script="run_vep.sh"


#### prior to running this script, you must activate a VEP conda env and have path to VEP ready.
#### This script has been written to run VEP and then modify the output so that the coding SNPs are ready to submit to rCOGS.

#Declare the number of mandatory args
margs=4

# Common functions - BEGIN
function example {
    echo -e "example: $script -v path/to/vep -i path/to/inDir -a assembly -o path/to/outDir"
}

function usage {
    echo -e "usage: $script MANDATORY args, OPTIONAL args \n"
}

function help {
  usage
    echo -e "MANDATORY:"
    echo -e "  -v,  --vep               The full path to the vep executable"
    echo -e "  -i,  --inputdir          The full path to preliminary COGs input files. If you ran Make_rCOGS_input_files.sh there will already be a list of rsids in this same folder."
    echo -e "  -a,  --assembly          The required assembly: GRCh37 or GRCh38"
    echo -e "  -o,  --outdir            The full path to directory where desired COGS input file will go\n"
    echo -e "OPTIONAL:"
    echo -e "  -m,  --method            The method with which to run VEP. rsid: runs VEP using rsids. varVCF: gets vcf using variant IDs. posVCF: gets vcf using variant positions. posVCF requires --positions option. Default = rsids."
    echo -e "  -p,  --positions         Use with posVCF. file with chr, position."
    echo -e "  -h,  --help              Prints this help\n"
  example
}

# Ensures that the number of passed args are at least equals
# to the declared number of mandatory args.
# It also handles the special case of the -h or --help arg.
function margs_precheck {
        if [ $2 ] && [ $1 -lt $margs ]; then
                if [ $2 == "--help" ] || [ $2 == "-h" ]; then
                        help
                        exit
                else
                usage
                        example
                exit 1 # error
                fi
        fi
}

# Ensures that all the mandatory args are not empty
function margs_check {
        if [ $# -lt $margs ]; then
            usage
                example
            exit 1 # error
        fi
}

# Exit if subcommands error
function check_error {
retval=$1
  if [ $retval -ne 0 ]; then
    exit $retval;
  fi
}
# Common functions - END

# Main
margs_precheck $# $1

vep=
inputDir=
assembly=
outDir=
method=rsid

# Args while-loop
while [ "$1" != "" ];
do
   case $1 in
   -v   | --vep )     shift
                            vep=$1
                                  ;;
   -i   | --inputdir )     shift
                            inputDir=$1
                                  ;;
   -a   | --assembly )        shift
                            assembly=$1
                                  ;;
   -m   | --method )          shift
	                    method=$1
			          ;;
   -p   | --positions )       shift
	                    positions=$1
			          ;;
   -o   | --outdir )   shift
                            outDir=$1
                                  ;;
   -t   | --tempdir )   shift
                            tempdir=$1
                                  ;;
   -h   | --help )        help
                            exit
                                  ;;
   *)
                            echo "$script: illegal option $1"
                            usage
                                                  example
                                                  exit 1 # error
                          ;;
    esac
    shift
done

# Pass here your mandatory args for check
margs_check $vep $inputDir $assembly $outDir

command -v ${vep} >/dev/null 2>&1 || { echo >&2 "Error: Cannot execute VEP. Check that the path is correct. Aborting."; exit 1; }

# Direct output to log file
exec > ${inputDir}/vep_output.log 2>&1

################ Run processes
echo "running VEP on ${inputDir}. Will output all coding SNPs for genome assembly ${assembly}"
echo "using VEP cache... make sure you have it available!"

### NOTE you cannot use the ID format in offline mode
### I have preselected the options for VEP in order to get coding SNPs for rCOGS. Please check for your purposes.
cd ${tempdir}

if [ $assembly == "GRCh37" ]; then
	echo "Using VEP options for GRCh37"
	if [ $method == "rsids" ]; then
               ${vep} -i ${inputDir}/rsids.txt --fork 4 --cache --merged --coding_only --no_stats --port 3337 --assembly ${assembly} -o ${inputDir}/coding.txt -v --tab --force_overwrite
		check_error $?

	       elif [ $method == "varVCF" ]; then
               echo "extracting SNPs from 1KG VCF using rsids"
               module load vcftools
               command -v vcftools >/dev/null 2>&1 || { echo >&2 "Error: Cannot execute vcftools. Check that it's installed and added to PATH. Aborting."; exit 1; }
               vcftools --gzvcf ~/1000Genomes/GRCh37/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --snps ${inputDir}/rsids.txt --out ${inputDir}/snps_1KG --recode
               echo "Running VEP from VCF"
               ${vep} -i ${inputDir}/snps_1KG.recode.vcf --use_given_ref --offline --fork 4 --cache --merged --coding_only --no_stats --port 3337 --assembly ${assembly} -o ${inputDir}/coding.txt -v --tab --force_overwrite
               check_error $?

	       elif [ $method == "posVCF" ]; then
               echo "extracting SNPs from 1KG VCF using positions"
               module load vcftools
               command -v vcftools >/dev/null 2>&1 || { echo >&2 "Error: Cannot execute vcftools. Check that it's installed and added to PATH. Aborting."; exit 1; }
               vcftools --gzvcf ~/1000Genomes/GRCh37/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz --positions $positions --recode --out ${inputDir}/snps_1KG

               echo "Running VEP from VCF"
               ${vep} -i ${inputDir}/snps_1KG.recode.vcf --use_given_ref --offline --fork 4 --cache --merged --coding_only --no_stats --port 3337 --assembly ${assembly} -o ${inputDir}/coding.txt -v --tab --force_overwrite
               check_error $?
       fi

########

########

elif [ $assembly == "GRCh38" ]; then
	echo "Using VEP options for GRCh38"
       if [ $method == "rsids" ]; then	
	       ${vep} -i ${inputDir}/rsids.txt --fork 4 --cache --merged --coding_only --no_stats --assembly ${assembly} -o ${inputDir}/coding.txt -v --tab --force_overwrite
	       check_error $?

       elif [ $method == "varVCF" ]; then 
	       echo "extracting SNPs from 1KG VCF using rsids" 
	       module load vcftools
	       command -v vcftools >/dev/null 2>&1 || { echo >&2 "Error: Cannot execute vcftools. Check that it's installed and added to PATH. Aborting."; exit 1; }
	       vcftools --vcf ~/1000Genomes/GRCh38/ALL_GRCh38_sites.vcf --snps ${inputDir}/rsids.txt --out ${inputDir}/snps_1KG --recode
	       echo "Running VEP from VCF"
	       ${vep} -i ${inputDir}/snps_1KG.recode.vcf --use_given_ref --offline --fork 4 --cache --merged --coding_only --no_stats --assembly ${assembly} -o ${inputDir}/coding.txt -v --tab --force_overwrite
	       check_error $?


       elif [ $method == "posVCF" ]; then
	       echo "extracting SNPs from 1KG VCF using positions"
	       module load vcftools
	       command -v vcftools >/dev/null 2>&1 || { echo >&2 "Error: Cannot execute vcftools. Check that it's installed and added to PATH. Aborting."; exit 1; }
	       vcftools --vcf ~/1000Genomes/GRCh38/ALL_GRCh38_sites.vcf --positions $positions --recode --out ${inputDir}/snps_1KG

	       echo "Running VEP from VCF"
               ${vep} -i ${inputDir}/snps_1KG.recode.vcf --use_given_ref --offline --fork 4 --cache --merged --coding_only --no_stats --assembly ${assembly} -o ${inputDir}/coding.txt -v --tab --force_overwrite
	       check_error $?
       fi
fi


echo "finished VEP run! The full results can be found in ${inputDir}/coding.txt."
echo "Now modifying the results file for rCOGs input"


cd ${inputDir}
# remove all lines in the file beginning with a double hash
sed '/^##/d' coding.txt | cut -f 2,4 | uniq > temp
# Split chr:pos and remove non standard chromosomes
sed 's/:/\t/g' temp | sed 1d | egrep "^[0-9XY]" > temp2
# May still need cleaning up. Sometimes snps have positions with "-number" and some genes have a weird id
awk '{$2+=0}1' temp2 > temp3
grep ENSG* temp3 | uniq > temp4
# Add header
echo -e "chr\tpos\tensg" > header
cat header temp4 > temp5
sed 's/ /\t/g' temp5 > ${outDir}/coding.format.txt
check_error $?
rm temp*
rm header

echo "finished formatting the VEP file! The COGS input can be found in ${outDir}/coding.format.txt."

