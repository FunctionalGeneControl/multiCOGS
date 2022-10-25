These R scripts implement an updated version of the rCOGS package with the following key differences:
* The algorithm has been modified to enable use with multivariate fine-mapped GWAS data (note that currently read_gwas() can only perform Wakefield synthesis based on a single causal variant assumption and multivariate fine-mapped data need to be provided directly to the compute_cogs() function
* All operations on genomic intervals are performed using data.tables directly. GenomicRanges are no longer used, avoiding the need for interconversion
* The number of promoter-proximal fragments included in the analysis can be specified by the user. Also note that make_vprom() is now called from within make_pchic() that now gains the vPromLen parameter
* Several minor bug fixes

Note that rROGS2 R package is still in progress, and currently the scripts should be are sourced directly into the analysis notebooks as in the example provided.
