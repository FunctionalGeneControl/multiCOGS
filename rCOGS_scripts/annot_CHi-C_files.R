#!/usr/bin/env Rscript
suppressMessages(library(argparser))
args = commandArgs(trailingOnly=T)
######## 

p <- arg_parser("Adding ENSG IDs to the pCHiC data", name="Rscript annot_CHi-C_files.R")
p <- add_argument(p, arg="<PeakMatrix>",
                  help="Full path to the peak matrix")
p <- add_argument(p, arg = "--rmap", 
                  help = "Full path to the rmap at fragment level. Will be used to make a new promoter baitmap that includes unbaited promoters.")
p <- add_argument(p, arg = "--baitmap", 
                  help = "Full path to the original baitmap at fragment level.")
p <- add_argument(p, arg="--TSS", 
                  help = "Full path to TSS file with at least these first 6 cols: ensg, genename, Chr, TSS, strand, type" )
p = add_argument(p, arg="--outDir",
                 help="Required output directory", default = ".")
p <- add_argument(p, arg="--inputDir", 
		  help="Directory for preliminary, non formatted input files")
p <- add_argument(p, arg="--gwas", 
		  help = "Full path to GWAS file with required columns named chromosome, base_pair_location, p_value and variant_id")
p <- add_argument(p, arg="--expandBaitmap", 
                  help = "Flag for expanding the baitmap to individual restriction fragments. Use if wanting to work with rmap which has unbinned baits and binned other ends. If true, must supply rmap_solBaits", 
                  flag = TRUE)
p <- add_argument(p, arg ="--rmapSolBaits", 
                  help = "Full path to the rmap with solitary baits and binned other ends. Required if interactions were called in this setting.")
p <- add_argument(p, arg="--verbose",
                 help = "Flag specifying whether to print process steps.", flag=TRUE)

opts = parse_args(p, args)

my_pm = opts[["<PeakMatrix>"]]
my_pm_name <- basename(my_pm)
rmap = opts[["rmap"]]
baits = opts[["baitmap"]]
tss = opts[["TSS"]]
outdir = opts[["outDir"]]
inputdir = opts[["inputDir"]]
gwas = opts[["gwas"]]
expa = opts[["expandBaitmap"]]
rmapB = opts[["rmapSolBaits"]]
verb = opts[["verbose"]]

###############

suppressMessages(library(data.table))
suppressMessages(library(tidyr))

sink(file = paste0(outdir, "/annot_CHiC_files.log"), type = c("output", "message")) 

## Format the GWAS file
print(paste0("Formatting GWAS file: ", gwas))
g <- fread(gwas)
g1 <- g[, .(chromosome, base_pair_location, p_value)]
names(g1) = c("chr", "pos", "p")
gName <- basename(gwas)
fwrite(g1, file = paste0(outdir, "/", gName, "_gwas.format.txt"), sep = "\t", 
       quote = F, row.names = F, col.names = T)
g2 <- as.data.table(unique(g[, variant_id]))
fwrite(g2, file = paste0(inputdir, "/rsids.txt"), sep = "\t", quote = F, 
       row.names = F, col.names = F)

## Annotate pCHi-C files
setwd(outdir)
cwd <- getwd()
print(paste0("Formatting the CHi-C files found in ", cwd))

# Make the pCHiC design/annotation file: fragID, ensg. This needs to include unbaited promoters.

baitmap_dt <- fread(baits)
rmap_dt <- fread(rmap)
names(baitmap_dt) = c("Chr", "start", "end", "fragID", "bait_gene")
names(rmap_dt) = c("Chr", "start", "end", "fragID")

get_new_baitmap <- function(tss_dt) {
  mytss <- tss_dt[, 1:6]
  names(mytss) = c("ensg", "genename", "Chr", "TSS", "strand", "type")
  mytss[, TSS2 := TSS]
  setkey(mytss,  Chr, TSS, TSS2)
  v <- foverlaps(rmap_dt, mytss, by.x = c("Chr", "start", "end"), nomatch = NULL)
  
  return(v)
}
eh_h_tss <- fread(tss)
print("Formatting baitmap file")
my_new_baitmap <- get_new_baitmap(tss_dt = eh_h_tss)
check <- nrow(my_new_baitmap)
if(check == 0) {
	stop("Something went wrong with baitmap annotation, please check!")
}

############

my_new_baitmap2 <- unique(my_new_baitmap[, .(fragID, ensg, type)])
names(my_new_baitmap2) = c("fragid", "ensg", "biotype")
my_new_baitmap3 <- my_new_baitmap2[!is.na(ensg)] # FragIDs in this baitmap correspond to unbinned fragments.
check2 <- nrow(my_new_baitmap3)
if(check2 == 0) {
        stop("Something went wrong with baitmap annotation, please check!")
}
fwrite(my_new_baitmap3, file = "PCHiC_design_annotation_plusUnbaited_with_geneType.txt", 
       row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
print("Annotated baitmap saved as PCHiC_design_annotation_plusUnbaited_with_geneType.txt")

# Make the pCHi-C data file: add ensg, name, biotype and strand to the pCHiC peak matrix
# using the baitmap + ensg + biotype.
pm <- fread(my_pm)
check3 <- nrow(pm)
if(check3 == 0) {
	stop("Peakmatrix file is empty, please check!")
}

if("clusterID" %in% names(pm)) {
  pm[, "clusterID" := NULL]
}

if("clusterPostProb" %in% names(pm)) {
  pm[, "clusterPostProb" := NULL]
}

# We have to get the ENSG IDs for the original baits, based on fragID.
baitmap_dt_nl <- unique(baitmap_dt[, .(bait_gene, fragID)])
my_new_baitmap_small <- my_new_baitmap[, c("ensg", "genename", "strand", "type", "fragID")]
baitmap_dt_nl_genes <- my_new_baitmap_small[baitmap_dt_nl, on = c("fragID"), nomatch = NULL]

h <- unique(baitmap_dt_nl_genes[, .(ensg, genename, type, strand, fragID)])
check4 <- nrow(h)
if(check4 == 0) {
	stop("Something went wrong with annotating the baitmap, please check!")
}



#### In case the data had been produced using sol baits/binning:
if(expa == TRUE) {
  print("Expanding the binned other ends in peak matrix to individual fragments")
  binnedOE <- fread(rmapB)
  names(binnedOE) = c("Chr", "start", "end", "fragID")
  # Do foverlaps between this and rmap_dt to get the corresponding frags
  setkey(rmap_dt, Chr, start, end)
  
  ### Firstly, for the baits.
  both_rmaps <- foverlaps(binnedOE, rmap_dt, by.x=c("Chr", "start", "end"), nomatch = NULL) # they should overlap perfectly
  names(both_rmaps) = c("Chr", "DpnStart", "DpnEnd", "DpnFragID", "binStart", "binEnd", "binFragID")
  setkey(pm, baitID)
  annot1 <- pm[both_rmaps, on = c(baitID = "binFragID"), allow.cartesian=TRUE, nomatch = NULL]
  # PM may or may not have bait coordinates, but the frag level baitmap must have bait coordinates - which we now use.
  annot1[, c("baitID", "binStart", "binEnd") := NULL]
  if("baitChr" %in% names(annot1)) {
    annot1[, "baitChr" := NULL]
  }
  if("baitStart" %in% names(annot1)) {
    annot1[, "baitStart" := NULL]
  }
  if("baitEnd" %in% names(annot1)) {
    annot1[, "baitEnd" := NULL]
  }
  # Name the new DpnII coords as the baits.
  setnames(annot1, c("Chr", "DpnStart", "DpnEnd", "DpnFragID"), c("baitChr", "baitStart", "baitEnd", "baitID"))
  
  ### Secondly, for the other ends.
  ### Modify the name of oeID if needed.
  if("otherEndID" %in% names(annot1)) {
    setnames(annot1, "otherEndID", "oeID")
  }
  setkey(annot1, oeID)
  pm <- annot1[both_rmaps, on = c(oeID = "binFragID"), allow.cartesian=TRUE, nomatch = NULL]
  # PM may or may not have OE coordinates, but the frag level baitmap must have OE coordinates - which we now use.
  pm[, c("oeID", "binStart", "binEnd") := NULL]
  if("oeChr" %in% names(pm)) {
    pm[, "baitChr" := NULL]
  }
  if("oeStart" %in% names(pm)) {
    pm[, "oeStart" := NULL]
  }
  if("oeEnd" %in% names(pm)) {
    pm[, "oeEnd" := NULL]
  }
  setnames(pm, c("Chr", "DpnStart", "DpnEnd", "DpnFragID"), c("oeChr", "oeStart", "oeEnd", "oeID"))
  # The coordinates in the new pm correspond to DpnII fragments.
  
  ### In case there is no "dist" column in the peak matrix. This is required for rCOGS input, but not used.
  ### "The dist column lists the linear distance between each bait and other end, with the positive sign indicating that the bait is 
  ###  located upstream of the other end (based on chromosomal coordinates), and the negative sign indicating otherwise. NA in the dist
  ###  field indicates a trans-chromosomal interaction". From https://osf.io/cn4k8
  ### add names too. Deprecated in COGS function; we use TSS annotations of the baitmap instead.
  if(!"dist" %in% names(pm)) {
    pm[, "baitMid" := (baitStart + baitEnd)/2]
    pm[, "oeMid" := (oeStart + oeEnd)/2]
    pm[, "dist" := oeMid - baitMid]
    pm[, c("baitMid", "oeMid") := NULL]
  }
  if(!"baitName" %in% names(pm)) {
    pm[, "baitName" := NA]
  }
  if(!"oeName" %in% names(pm)) {
    pm[, "oeName" := NA]
  }
  
}

print("Annotating the peak matrix with gene information")
setkey(pm, baitID)
annot <- pm[h, on = c(baitID = "fragID"), allow.cartesian=TRUE]

moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

setcolorder(annot, moveme(names(annot), "ensg first; genename after ensg; type after genename; strand after type; baitName after baitID; oeName after oeID"))
names(annot)[2] = "name"
names(annot)[3] = "biotype"

annot2 <- annot[!is.na(baitChr)]
annot3 <- annot2[!is.na(oeChr)]
annot4 <- annot3[oeChr != "MT"]

check5 <- nrow(annot4)
if(check5 == 0) {
	stop("The final formatted peak matrix was empty, something went wrong!")
}

fwrite(annot4, file = paste0(my_pm_name, "_pm.format.txt"), 
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

print("Finished annotating PCHiC files")














