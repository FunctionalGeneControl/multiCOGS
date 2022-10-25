#' This function computes COGS scores for a given trait across a set of pchic regions
#' @param gwas - data.table obtained from load_gwas
#' @param csnps = data.table object obtained from make_csnps
#' @param digest data.table object obtained from load_digest
#' @param regions list object obtained from make_pchic
#' @param feature.names - a character vector of regions to include these must be present in pcHi-C file
#' @return data.table of results
#' @export

compute_cogs <- function(gwas,csnps,digest,regions,feature.names){
  ## annotate coding SNPs - we remove these and add back in for gene of interest

  message("Extracting coding SNPs")
  
  # cred_set is a parameter defined by SUSIE reflecting the potentially multiple credible sets of causal SNPs per LD block, 
  # so if it's listed then we are working with single-variant fine-mapping and should just account for LD blocks
  if(!"cred_set" %in% names(gwas)){
    gwas[, cred_set:=1]
  }
  
  gwas[, pos1:=pos]
  setnames(csnps, c(1,2), c("chr", "csnp.pos"))
  csnps[, csnp.pos1:=csnp.pos]
  
  setkey(csnps, chr, csnp.pos, csnp.pos1)
  
  # As a result of this, gwas will contain, csnp.pos, ensg and csnp.pos1, which will be non-NA only for coding SNPs
  gwas <- foverlaps(gwas, csnps, by.x = c("chr", "pos", "pos1"), by.y=c("chr", "csnp.pos", "csnp.pos1"), nomatch=NA)
  # Now we are computing sum Ppi from coding SNPs only. Note that only coding SNPs will have ensg!=NA at this stage  
  csnps.DT <- gwas[,list(csnp.ppi=sum(ppi)),by=c('ensg','ld', 'cred_set')][!is.na(ensg),][,.(ensg=ensg,ld=ld,cred_set=cred_set,frag.ppi=csnp.ppi)]

  gwas <- gwas[!is.na(pos)]
  
  message("Mapping SNPs to restriction fragments")
  setkey(digest, chr, start, end)
  ol.DT <- foverlaps(gwas, digest, by.x=c("chr", "pos", "pos1"), by.y=c("chr", "start", "end"), nomatch=0)
  gwas[, pos1:=NULL]
  gwas[, csnp.pos1:=NULL]
  ol.DT[, pos1:=NULL]
  ol.DT[, csnp.pos1:=NULL]
  
  #next we remove cSNPs and compute sum for each fragment, LD block and credible set
  ## note that we can have ld block boundaries going across fragments, in which case the fragments will essentially be split into two
  ## because the primary key is always "fragid, ld, (cred_set)" 
  ## note also that cred_sets can overlap, and so in this case we count the same SNPs/fragids multiple times for different cred_sets
  ol.DT <- ol.DT[is.na(csnp.pos),][,list(frag.ppi=sum(ppi)),by=c('cred_set', 'ld','fragid')]

  ## to compute scores we want for each lu[['ensg']][['tissue']] is a vector of fragids
  message("Gathering fragments of interest")
  all <- lapply(names(regions),function(ct){
    tDT <- regions[[ct]]
    tDT[,.(ensg,fragid,tissue=ct)]
  }) %>% rbindlist
  mDT <- merge(all,ol.DT,by.x='fragid',by.y='fragid')
  mDT[,tissue:=factor(tissue)]
  ## we can compute any score by selecting one or more tissues
  if(missing(feature.names)){
      feature.names <- c(names(regions),'coding_snp')
      message(sprintf("feature.names argument missing - will compute the score using %s",paste(feature.names,collapse=",",sep=",")))
  }
  fs.DT <- mDT[tissue %in% feature.names,]
  ## remove duplicate fragments by gene, LD block and credible set so don't count twice
  fs.DT <- fs.DT[!duplicated(fs.DT[,.(cred_set,ld, fragid,ensg)]),][,.(ensg,ld,cred_set,frag.ppi,fragid)]
  ## check to see if coding snps need to be included if so add
  if(any(feature.names=='coding_snp')){
    message("Adding coding SNPs")
    # note that the coding SNP's sum(ppi) have already been counted separately per gene/cred_set/LD block, 
    # irrespectively of what fragids they map to
    # csnps.DT doesn't have a fragid column, so add it, setting it to zero (value doesn't matter - just so rbindlist works)
    fs.DT <- rbindlist(list(fs.DT,csnps.DT[, fragid:=0]))
  }
  message("Computing COGS scores")
  fs.ld.score <- fs.DT[,list(ld.score=sum(frag.ppi)),by=c('cred_set','ld','ensg')]
  fs.ld.score[,list(cogs=1-prod(1-ld.score)),by='ensg']
}
