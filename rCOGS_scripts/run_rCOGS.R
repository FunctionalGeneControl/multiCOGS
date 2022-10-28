#!/usr/bin/env Rscript
suppressMessages(library(argparser))
args = commandArgs(trailingOnly=T)
suppressMessages(library(data.table))

########################################################
#### Run rCOGS using the data.table version. 
#### Run this from within folder with rCOGS input files
########################################################

p <- arg_parser("Running rCOGS the data.table way", name="Rscript run_Mikhails_rCOGS.R")
p <- add_argument(p, arg="--ld",
                  help="LD regions file with columns named: chr, start, end", default="_ld.format.bed$")
p <- add_argument(p, arg = "--maf", 
                  help = "MAF file with columns named: chr, pos, maf.", default="formatted.maf.txt$")
p <- add_argument(p, arg = "--gwas", 
                  help = "GWAS data with columns named: chr, pos, p.", default="_gwas.format.txt$")
p <- add_argument(p, arg="--ncases", 
                  help = "Integer for the number of cases in the gwas", type= "numeric")
p <- add_argument(p, arg="--ncontrols",
                  help="Integer for the number of controls in the gwas", type= "numeric")
p <- add_argument(p, arg="--pmFormat",
                  help="Formatted peak matrix with biotypes", default = "_pm.format.txt$")
p <- add_argument(p, arg = "--rmap",
                  help = "RMAP with columns named: chr, start, end, fragid", default = ".rmap_wHeader.txt$")
p <- add_argument(p, arg="--bannot", 
                  help = "Annotated baits with columns named: fragid, ensg and biotype", default = "PCHiC_design_annotation_plusUnbaited_with_geneType.txt")
p <- add_argument(p, arg = "--coding", 
                  help = "Coding SNPs with columns named: chr, pos, ensg", default = "coding.format.txt")
p <- add_argument(p, arg = "--assembly", 
                  help = "Assembly as either GRCh37 or GRCh38")
p <- add_argument(p, arg="--cogsIn", 
                  help="Directory containing the COGS formatted input files mentioned above", default = ".")
p <- add_argument(p, arg = "--cogsOut", 
                  help = "Directory for COGS output files")
p <- add_argument(p, arg="--verbose",
                  help = "Flag specifying whether to print process steps.", flag=TRUE)
p <- add_argument(p, arg = "--featureNames", 
                  help = "Comma separated list of cell types to use. If using, make sure to add VProm and coding_snp",
                  default = NA)
p <- add_argument(p, arg = "--vProm", 
                  help = "Number of fragments to use when creating virtual promoter regions", default = 1, type = "numeric")
p <- add_argument(p, arg = "--chicThresh", type = "numeric", default = 5, 
                  help = "The hard threshold for CHiC interactions. Scores will only be considered ABOVE this value.")
p <- add_argument(p, arg = "--susie", flag = TRUE,
                  help = "Flag to run with SuSIE. If so, provide the SuSIE .csv file in place of the GWAS file above. Do not need to supply LD blocks.")

opts = parse_args(p, args)

ld = opts[["ld"]]
maf= opts[["maf"]]
gwas=opts[["gwas"]]
ncases=opts$ncases
ncontrols=opts$ncontrols
pmFormat=opts[["pmFormat"]]
rmap = opts[["rmap"]]
bannot = opts[["bannot"]]
coding = opts[["coding"]]
assembly = opts[["assembly"]]
cogsIn = opts[["cogsIn"]]
cogsOut = opts[["cogsOut"]]
verb = opts[["verbose"]]
featureNames = opts[["featureNames"]]
vPromLen = opts[["vProm"]]
chic_thresh = opts[["chicThresh"]]
susie = opts[["susie"]]


############################# Load sources #################################
source("~/rCOGS/rCOGS_SUSIE/rCOGS/R/cogs.R")
source("~/rCOGS/rCOGS_SUSIE/rCOGS/R/data.R")
source("~/rCOGS/rCOGS_SUSIE/rCOGS/R/gwas.R")
source("~/rCOGS/rCOGS_SUSIE/rCOGS/R/pchic.R")
source("~/rCOGS/rCOGS_SUSIE/rCOGS/R/sCVPP.R")
############################################################################

setDTthreads(threads = 4)
setwd(cogsIn)

dir.create(cogsOut, showWarnings = FALSE, recursive = TRUE)
sink(file = paste0(cogsOut, "/rCOGS.log"))

############################# Load input files  ###############################

cwd <- as.character(getwd())
print(paste0("Loading input files for COGS run, using the files in ", cwd))

if(susie == FALSE) {
  my_ld <- load_ld_regions(dir(pattern=ld)[1])
  if(exists("my_ld") == FALSE) {
    stop("LD file not found") 
  }
}

my_maf <- load_ref_maf(dir(pattern=maf)[1], min.maf=0.05)
if(exists("my_maf") == FALSE) {
	   stop("MAF file not found")
}

if(susie == TRUE) {
  print(paste0("Running with SuSIE results from ", gwas))
  sus <- fread(gwas)
  
  sus[, chr:=as.numeric(gsub("(\\S+)\\:\\d+", "\\1", snp))]
  sus[, pos:=as.numeric(gsub("\\S+\\:(\\d+)", "\\1", snp))]
  
  # include all available pip_sets
  pip_sets = paste0("pip_set_", 1:100)
  pip_sets = c(pip_sets[pip_sets %in% names(sus)])
  
  # where SuSIE's data aren't available or were filtered out, use single.pp from a single causal variant model
  # Require that all pip_sets are NA, and put single.pp into a new column.
  sus[, sum_pip := rowSums(.SD, na.rm = TRUE), .SDcols = pip_sets]
  sus[sum_pip == 0, only.single.pp := single.pp]
  
  pip_sets_wSingle <- c("only.single.pp", pip_sets)
  
  sus1 = sus[, c("chr", "pos", "block", pip_sets_wSingle), with=F]
  sus.melt = melt(sus1, id.vars = c("chr", "pos", "block"), variable.name = "cred_set", value.name = "ppi")
  my_gwas = sus.melt[!is.na(ppi)]
  setnames(my_gwas, "block", "ld")
} else {
  print("Loading GWAS file")
  my_gwas <- load_gwas(dir(pattern=gwas)[1], my_maf, my_ld, n.cases=ncases, n.controls=ncontrols, 
                       chrcol = "chr", poscol = 'pos', pcol='p')
}

if(exists("my_gwas") == FALSE) {
  stop("GWAS file not found")
}

print(paste0("Making feature sets with a CHi-C threshold of ", chic_thresh))
feature.sets <- make_pchic(dir(pattern=pmFormat)[1], biotype.filter='protein_coding', 
                           vPromLen = vPromLen, f.digest = dir(pattern=rmap)[1],
                           f.design=dir(pattern=bannot)[1], chic_thresh=chic_thresh)
if(exists("feature.sets") == FALSE) {
	stop("Something is wrong with the peak matrix, rmap or annotated baitmap! Please check the sources.")
}

csnps <- make_csnps(dir(pattern = coding)[1])
no_csnps <- nrow(csnps)
print(paste0("No. coding snps is ", no_csnps))

digest <- load_digest(dir(pattern = rmap)[1])

###############################################################################

############################# Remove the MHC  #################################
# GRCh37 coords from RCOGS vignette: https://ollyburren.github.io/rCOGS/articles/Quickstart.html 
# GRCh38 coords from NCBI: https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13
print(paste0("Removing MHC region using ", assembly, " coordinates"))
if(assembly == "GRCh37") {
  mhcStart <- 25000000
  mhcEnd <- 35000000 
} else {
  if(assembly == "GRCh38") {
    mhcStart <- 28510120
    mhcEnd <- 33480577
  } else {
    stop("Please specify assembly as GRCh37 or GRCh38")
  }
}

mhc.idx <- which(my_gwas$chr==6 & between(my_gwas$pos,mhcStart,mhcEnd)) 
if(length(mhc.idx)>0) {
  my_gwas <- my_gwas[-mhc.idx,]
}

##############################################################################

############################# Run COGS  ######################################
print("Running COGS...")

# parse feature names, if they have been specified

if(!is.na(featureNames)) {
  print(paste0("You have chosen the following features: ", featureNames))
  feature.names <- as.character(strsplit(featureNames, ",")[[1]])
  overall.scores <- compute_cogs(my_gwas,csnps,digest,feature.sets,feature.names)
} else {
  print(paste0("You have chosen all features in the peak matrix"))
  overall.scores <- compute_cogs(my_gwas,csnps,digest,feature.sets)
}


overall.scores <- overall.scores[order(cogs,decreasing = TRUE),]

#### Create an annotated version
# Already made a TSS file for the ensemblHavana_Havana TSS in GRCh37 and GRCh38. See make_TSS_annotation_table.R
TSS <- fread("~/external_data/ensembl/for_rCOGS_testing/ANNOTATIONS_EnsemblHavanaMerge_plusHavana_V88hg38_V107hg19.txt")
print("Annotating the COGS results using TSS from ~/external_data/ensembl/for_rCOGS_testing/ANNOTATIONS_EnsemblHavanaMerge_plusHavana_V88hg38_V107hg19.txt")

###
annot.scores <- TSS[overall.scores, on = "ensg"]


fwrite(overall.scores, file = paste0(cogsOut, "/COGS_scores_data.table.txt"), 
       sep="\t", quote = F, row.names = F, col.names = T)

fwrite(annot.scores, file = paste0(cogsOut, "/Annotated_COGS_scores_data.table.txt"), 
       sep="\t", quote = F, row.names = F, col.names = T)

### Save the ppi results.
fwrite(my_gwas, file = paste0(cogsOut, "/COGS_PPIs_data.table.txt"), 
       sep="\t", quote = F, row.names = F, col.names = T)


print(paste0("...finished running COGS! The results are saved in: ", cogsOut))
###########################################################################




