#!/bin/bash

script="Make_rCOGS_input_files.sh"

#### This script prepares input files for rCOGS.
#### Made this for easy switching between GWAS, CHiC data and assemblies.
#### Please run this script from the scripts directory.

#Declare the number of mandatory args
margs=6

# Common functions - BEGIN
function example {
    echo -e "example: $script -i path/to/input/files -o path/to/output/files -g gwas.txt -p pchic.peakmatrix -b baitmap -r rmap -hg38 -expa --rmapSolBaits path/to/rmapSolBaits"
}

function usage {
    echo -e "usage: $script MANDATORY args, OPTIONAL args \n"
}

function help {
  usage
    echo -e "MANDATORY:"
    echo -e "  -i,  --inputDir               The full path to where you want preliminary input files to go, like non formatted snps etc"
    echo -e "  -o,  --outDir                 The full path to where you want COGS input files to go"
    echo -e "  -g,  --gwas                   The full path to the trait data from GWAS catalog. Must have cols named chromosome, base_pair_location, p_value and variant_id"
    echo -e "  -p,  --pchic                  The full path to the pchic peak matrix"
    echo -e "  -b,  --baitmap                The full path to the baitmap"  
    echo -e "  -r,  --rmap                   The full path to the rmap. If bin+solBaits was used, this should be the fragment level rmap.\n"
    echo -e "OPTIONAL:"
    echo -e "  -hg19, --GRCh37                Flag that the required assembly is GRCh37"
    echo -e "  -hg38, --GRCh38                Flag that the required assembly is GRCh38"
    echo -e "  -expa, --expandBaitmap         Flag to expand the baitmap from binned to fragment level. Must also supply binned rmap with solitary baits using rmapSolBaits." 
    echo -e "  -sol,  --rmapSolBaits          Full path to the rmap in binned setting with solitary baits, if used"
    echo -e "  -t,  --tempdir                Directory for temporary files, default = current wd"
    echo -e "  -h,  --help                   Prints this help\n"
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

inputDir=
outDir=
tempDir=.
gwas=
pchic=
baitmap=
rmap=
GRCh37="FALSE"
GRCh38="FALSE"
rmapSolBaits=
expandBaitmap="FALSE"

# Args while-loop
while [ "$1" != "" ];
do
   case $1 in
   -i    | --inputDir )     shift
                            inputDir=$1
                                  ;;
   -o    | --outDir )       shift
                            outDir=$1
                                  ;;
   -t    | --tempDir )      shift
                            tempDir=$1
                                  ;;
   -g    | --gwas )         shift
                            gwas=$1
                                  ;;
   -p    | --pchic )        shift
	                    pchic=$1
			          ;;
   -b    | --baitmap )      shift
	                    baitmap=$1
			          ;;
   -r    | --rmap )         shift
	                    rmap=$1
			          ;;
   -hg19 | --GRCh37 )       GRCh37="TRUE"
                                  ;;
   -hg38 | --GRCh38 )       GRCh38="TRUE"
			          ;;
   -sol  | --rmapSolBaits ) shift
	                    rmapSolBaits=$1
			          ;;
   -expa | --expandBaitmap ) expandBaitmap="TRUE"
			          ;;
   -h    | --help )         help
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
margs_check $inputDir $outDir $gwas $rmap $baitmap $pchic


################ Run processes

mkdir -p $inputDir
mkdir -p $outDir

# Send output to log file
exec > ${outDir}/Make_rCOGS_input_files.log 2>&1

echo "Making rCOGS input files using the data in ${inputDir}".

### Requirements for rCOGS input files are described here: https://ollyburren.github.io/rCOGS/articles/Quickstart.html

########## 1. Make the approx. LD independent region file

## should be non-zero bed in the format chr, start, end
#### Using Berisa and Pickrell file in EUR from Bitbucket:
# https://bitbucket.org/nygcresearch/ldetect-data/src/master/EUR/fourier_ls-all.bed

## For GRCh38, using this new paper which uses the same method as Berisa and Pickrell:
## https://www.biorxiv.org/content/10.1101/2022.03.04.483057v1.full.pdf
## And downloaded from https://github.com/jmacdon/LDblocks_GRCh38

# input LD blocks are currently hard coded. Just need to specify GRCh37 or GRCh38.
echo "Finding the input LD files and formatting for rCOGS."

if [ $GRCh37 == "TRUE" ]; then 
	LD=~/external_data/fourier/fourier_ls-all_GRCh37.bed 
	echo "LD file used is from ~/external_data/fourier/fourier_ls-all_GRCh37.bed"
elif [ $GRCh38 == "TRUE" ]; then
	LD=~/external_data/fourier/LDblocks_GRCh38
	echo "LD file used is from ~/external_data/fourier/LDblocks_GRCh38"
else
	echo "Please specify GRCh37 or GRCh38"
fi

# Edit the header 
cd ${outDir}
echo -e "chr\tstart\tend" > header
sed 1d ${LD} > temp
cat header temp > temp2
mv temp2 LDblocks.bed
rm header temp
## The end of first block = beginning of second block, which is counted as overlap. 
## Take away one from the end.
awk -vOFS="\t" 'NR>1{print $1,$2,$3-1}' LDblocks.bed > temp
sed 's/chr//g' temp > temp2
LD_NAME="${LD##*/}" 
echo -e "chr\tstart\tend" | cat - temp2 > ${LD_NAME}_ld.format.bed
rm temp* LDblocks.bed

########## 2. Get the MAF for our SNPs
### using 1KG.
### eventually will need to be: chr, pos, maf
echo "Linking the MAF data to the COGS directory"
if [ $GRCh37 == "TRUE" ]; then
        ln -s ~/1000Genomes/GRCh37/EUR_HRJ/ALL.EUR.GRCh37.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.MAF01.formatted.maf.txt ${outDir}/ALL.EUR.GRCh37.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.MAF01.formatted.maf.txt
	MAF_FILE=ALL.EUR.GRCh37.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.MAF01.formatted.maf.txt
	echo "MAF file is from ~/1000Genomes/GRCh37/EUR_HRJ/ALL.EUR.GRCh37.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.MAF01.formatted.maf.txt"

elif [ $GRCh38 == "TRUE" ]; then
	ln -s ~/1000Genomes/GRCh38/EUR/COGS_michiel/MINORallelefrequencies_EUR_GRCh38.tsv ${outDir}/allelefrequencies_EUR_GRCh38.tsv.formatted.maf.txt
	MAF_FILE=allelefrequencies_EUR_GRCh38.tsv.formatted.maf.txt
	echo "MAF file is from ~/1000Genomes/GRCh38/EUR/COGS_michiel/MINORallelefrequencies_EUR_GRCh38.tsv"
else
        echo "Please specify GRCh37 or GRCh38"
fi


########## 3. Get the restriction fragment digest file - needs to be chr, start, end, frag ID
echo "Adding header to the restriction fragment digest file"

cd ${outDir}
echo -e "chr\tstart\tend\tfragid" > header
RMAP="${rmap##*/}" # this gets the name of the file.
cat header ${rmap} > ${RMAP}.rmap_wHeader.txt
check_error $?
rm header

########## 4. Get the pcHi-C design/annotation file. Columns are: fragID, ENSG
echo "Running Rscript to annotate the pCHiC data with gene names"

# hard coded TSS files from ensembl.
# This was after discussion and comparison between various sources for TSS. Currently using these sources:
if [ $GRCh37 == "TRUE" ]; then
        ln -s ~/external_data/ensembl/for_rCOGS_testing/EnsemblHavanaMerge_plusHavana_Jul2022_V107_TSS_GRCh37.txt ${inputDir}/
	TSS=${inputDir}/EnsemblHavanaMerge_plusHavana_Jul2022_V107_TSS_GRCh37.txt
	echo "TSS file used is ~/external_data/ensembl/for_rCOGS_testing/EnsemblHavanaMerge_plusHavana_Jul2022_V107_TSS_GRCh37.txt"

elif [ $GRCh38 == "TRUE" ]; then
	ln -s ~/external_data/ensembl/for_rCOGS_testing/EnsemblHavanaMerge_plusHavana_Mar2017_V88_TSS_GRCh38.txt ${inputDir}/
	TSS=${inputDir}/EnsemblHavanaMerge_plusHavana_Mar2017_V88_TSS_GRCh38.txt
	echo "TSS file used is ~/external_data/ensembl/for_rCOGS_testing/EnsemblHavanaMerge_plusHavana_Mar2017_V88_TSS_GRCh38.txt"
else
        echo "Please specify GRCh37 or GRCh38"
fi

############## Annotate PCHi-C and design files using annot_CHi-C_files.R from the scripts folder (working directory)
cd $PBS_O_WORKDIR
### This script will annotate pchic and design files. It will add unbaited promoters to the baitmap.

echo "Annotating the CHiC files using annot_CHi-C_files.R"

if [ $expandBaitmap == "TRUE" ]; then
	Rscript annot_CHi-C_files.R \
		${pchic} \
		--rmap ${rmap} \
		--baitmap ${baitmap} \
		--TSS ${TSS} \
		--outDir ${outDir} \
		--inputDir ${inputDir} \
		--gwas ${gwas} \
		--expandBaitmap \
		--rmapSolBaits ${rmapSolBaits}
	check_error $?
else
	Rscript annot_CHi-C_files.R \
		${pchic} \
                --rmap ${rmap} \
                --baitmap ${baitmap} \
                --TSS ${TSS} \
                --outDir ${outDir} \
                --inputDir ${inputDir} \
                --gwas ${gwas}
	check_error $?
fi

echo "...done!"












