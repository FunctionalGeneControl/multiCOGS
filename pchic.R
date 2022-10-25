## routine to load promoter pcHi-C DATA

library(data.table)
library(magrittr)
# library(GenomicRanges)

PCHIC_COLS <- c("ensg","name","biotype","strand","baitChr","baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName","dist")
DIGEST_COLS <- c('chr','start','end','fragid')
DESIGN_COLS <- c('fragid','ensg','biotype')
CSNP_COLS <- c('chr','pos','ensg')

#' This loads a restriction digest file with columns chr,start,end and fragid
#' @param f a scalar representing a file name of digest file
#' @return a data.table representing restriction digest
#' @export

### GRanges to data.table conversion tested
load_digest <- function(f){
  DT <- fread(f)
  if(!identical(names(DT),DIGEST_COLS)){
    stop(sprintf("COGSR:load_digest file format error, expected columns ''%s' got '%s'",
      paste(DIGEST_COLS,sep=',',collapse=',',),
      paste(names(DT),sep=',',collapse=',',))
    )
  }
  DT[,chr:=fix_chr(chr,f)]
  DT <- DT[order(chr,start),][chr!=25,]
  setkey(DT, chr, start, end)
  if(nrow(foverlaps(DT,DT))!=nrow(DT)){
    stop(sprintf("OGSR:load_digest_regions file format error, digest regions appear to overlap, please check."))
  }
  DT
}

#' This loads a capture design file with columns fragid and ensg. fragid should match digest file.
#' @param f a scalar representing a file name of design file
#' @return data.table object representing capture design

load_capture_design <- function(f){
  DT <- fread(f)
  if(!identical(names(DT),DESIGN_COLS)){
    stop(sprintf("COGSR:load_capture_design file format error, expected columns ''%s' got '%s'",
      paste(DESIGN_COLS,sep=',',collapse=',',),
      paste(names(DT),sep=',',collapse=',',))
    )
  }
  DT
}

#' This function computes 'Virtual' promoter regions. These are baited fragments and the fragments either side.
#' @param f.digest a scalar representing a file name of digest file
#' @param f.design a scalar representing a file name of design file
#' @param vPromLen number of restriction fragments on either side to include into vprom 
#' @param keep.ensg a vector of ensembl id's to include. If ommited then all ensembl ids regardless of further filtering are included
#' @return data.table with two columns mapping ensg to virtual promoter fragment ids
#' @export

### New N parameter and major algorithm rewrite tested

make_vprom <- function(f.digest,f.design,vPromLen=1, keep.ensg){
  frags <- load_digest(f.digest)
  if(!missing(keep.ensg)){
    des.DT <- load_capture_design(f.design)[ensg %in% keep.ensg]
  }else{
    des.DT <- load_capture_design(f.design)
  }
  
  baitIDs = unique(des.DT$fragid)
  vprom.DT = data.table(baitID = baitIDs, vpromID = baitIDs)
  if(vPromLen>0){
	  for(i in 1:vPromLen){
    		vprom.DT = rbindlist(list(vprom.DT,
                              data.table(baitID=baitIDs, vpromID=baitIDs-i),
                              data.table(baitID=baitIDs, vpromID=baitIDs+i)))
  	  }
  }
  vprom.DT<-merge(des.DT,vprom.DT,by.x='fragid',by.y='baitID', allow.cartesian=TRUE) 
  # we can have the same Vprom fragment mapping to multiple genes and vice versa
  out<-vprom.DT[,.(ensg,fragid=vpromID)]
  setkey(out, ensg,fragid)
  # but from allow.cartesian it seems we are also potentially getting duplicates?
  return(unique(out, by=c("ensg", "fragid")))
}

PCHIC_COLS <- c("ensg","name","biotype","strand","baitChr","baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName","dist")

#' This function computes loads in a set of promoter capture Hi-C files (pcHi-C). This is a tab delimited file that contains the following mandatory columns
#' \enumerate{
#' \item ensg - Ensembl gene id for gene
#' \item name - Gene name
#' \item biotype - Gene biotype from ensembl e.g. protein_coding
#' \item strand - Gene strand
#' \item baitChr - Bait (A captured restriction fragment) chromosome
#' \item baitStart - Bait fragment start
#' \item baitEnd - Bait fragment end
#' \item baitID - Bait fragment ID should match ID in digest file
#' \item baitName - Free text of what is captured - deprecated.
#' \item oeChr - 'Other end' (Promoter interacting restriction fragment) chromosome
#' \item oeStart - 'Other end' fragment start
#' \item oeEnd - 'Other end' fragment end
#' \item oeID - 'Other end' fragment ID should match ID in digest file
#' \item oeName - Free text of what is interacting - deprecated.
#' }
#' Additional columns contain CHiCAGO scores for one or more columns.
#' @param f a scalar representing file name of pcHi-C file
#' @param chic_thresh scalar representing CHiCAGO threshold for calling a significant interaction. Default is set to 5.
#' @param biotype.filter scalar/vector of characters representing a biotype filter. Default is set to 'protein_coding'
#' @param remove.trans.contacts a boolean on whether to remove interaction between chromosomes. Default is set to TRUE.
#' @param vPromLen - the number of fragments around the baits to be scored as virtual promoters
#' @param f.digest a scalar representing a file name of digest file
#' @param f.design a scalar representing a file name of design file
#' @param vprom.keep.ensg a vector of ensembl id's to include into virtual promoter generation. If ommited then all ensembl ids regardless of further filtering are included
#' @return list of data.tables where each data.table has three columns ensg,fragid and chic.score
#' @export

#### Note vpromoter fragments are now added from within make_pchic, so no need to execute make_vprom
#### If not desirable, set vPromLen to zero.
make_pchic <- function(f,chic_thresh=5,biotype.filter='protein_coding',remove.trans.contacts=TRUE, 
                       vPromLen = 1, f.digest, f.design, vprom.keep.ensg, vprom.filter.biotype=TRUE){
  DT <- fread(f)
  cnames <- names(DT)
  idx.ct.cols <- which(cnames %in% PCHIC_COLS)
  if(length(idx.ct.cols) != length(PCHIC_COLS))
    stop(sprintf("COGSR:load_pchic Missing required columns in pcHi-C file %s",f))
  if(!missing(biotype.filter)){
      message(sprintf("Filtering by biotypes %s",paste(biotype.filter,sep=",",collapse=',')))
      DT <- DT[biotype %in% biotype.filter,]
  }
  trans.idx <- which(DT$baitChr != DT$oeChr)
  if(remove.trans.contacts & length(trans.idx)>0){
    message(sprintf("Removing %d trans contacts",length(trans.idx)))
    DT <- DT[-trans.idx]
  }
  cell.types <- cnames[!cnames %in% PCHIC_COLS]
  DT <- DT[,c('baitChr','oeChr'):=list(fix_chr(baitChr),fix_chr(oeChr))]
  lenBefore <- nrow(DT)
  DT <- DT[abs(baitID-oeID)>vPromLen] # so we don't count vPromoters twice
  lenAfter <- nrow(DT)
  if (lenBefore-lenAfter>0){
    message(sprintf("Removed %d contacts whose otherEnds map to +/- %d fragment(s) around the bait that will be scored separately as virtual promoters", 
                  lenBefore-lenAfter, vPromLen))
  }
  
  message(sprintf("Processing %d cell types",length(cell.types)))
  cts <- lapply(cell.types,function(ct){
    message(sprintf("Processing %s",ct))
    tDT <- DT[get(`ct`)>chic_thresh,.(ensg,fragid=oeID,chic.score=get(`ct`))]
  })
  names(cts) <- cell.types
  
  if(!is.null(f.digest) & !is.null(f.design)){ 
   # note this is where the actual baited fragment is added to the analysis space, 
   # even if vPromLen=0 and irrespectively of whether or not there are any contacts detected for it 
      message("Generating virtual promoters...")
      if(vprom.filter.biotype & !missing(biotype.filter)){
              message("Filtering capture design for biotype...")
              design = load_capture_design(f.design)
              # a tiny fraction of genes will have biotype listed in DT but not in baitmap, so include it as well
              incl.ensg = unique(c(design[biotype %in% biotype.filter]$ensg, DT$ensg))
              message("Including the vpromoters of ", length(incl.ensg), " genes (make sure all genes you want included are listed in capture design file even if unbaited)")
              if(!missing(vprom.keep.ensg)){
                  vprom.keep.ensg = c(vprom.keep.ensg, incl.ensg)
              }else{
                vprom.keep.ensg = incl.ensg
              }
      }
      cts[['VProm']] = make_vprom(f.digest = f.digest, f.design=f.design,
                                  vPromLen = vPromLen, keep.ensg = vprom.keep.ensg)
  }else{
      message("f.digest and/or f.design not defined, make_vprom() cannot be executed internally.")
  }
  
  return(cts)
}

### Example:
### feature.sets <- make_pchic('peakmatrix_GRCh38.tsv',biotype.filter='protein_coding', 
###  vPromLen = 1, f.digest="rmap_GRCh38.tsv",f.design="bannot_GRCh38.tsv")
### Note no need to run feature.sets[['Vprom']] <- make_vprom(<...>) separately


#' This function loads in a coding SNPs matched to a given ensembl id, columns are chr,pos and ensg
#' @param f a scalar representing file name of coding SNP file
#' @param keep.ensg a vector of ensembl id's to include. If ommited then all ensembl ids regardless of further filtering are included
#' @return data.table object of coding SNPs
#' @export

make_csnps <- function(f,keep.ensg){
  if(!missing(keep.ensg)){
    DT <- fread(f)[ensg %in% keep.ensg]
  }else{
    DT <- fread(f)
  }
  if(!identical(names(DT),CSNP_COLS)){
    stop(sprintf("COGSR:load_capture_design file format error, expected columns ''%s' got '%s'",
      paste(CSNP_COLS,sep=',',collapse=',',),
      paste(names(DT),sep=',',collapse=',',))
    )
  }
  DT[,chr:=fix_chr(chr,f)]
  DT <- DT[order(chr,pos,ensg),]
  
  DT
}
