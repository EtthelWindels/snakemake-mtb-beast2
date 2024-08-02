## ------------------------------------------------------------------------
## Subsample alignment
## 2023-07-26 Etthel Windels
## ------------------------------------------------------------------------


# Load libraries ----------------------------------------------------------

library(dplyr)
library(seqinr)



# Read files --------------------------------------------------------------

alignment_in <- seqinr::read.fasta(file = snakemake@input[["alignment_in"]], seqtype="DNA", forceDNAtolower = F)
full_metadata <- read.table(file = snakemake@input[["meta_all"]], header=T,sep= '\t')


  
# Create subset -----------------------------------------------------------

subsample <- function(fasta, n, seed){
  set.seed(seed)
  id <- sample(names(fasta)[names(fasta) != "Mycobacterium_canettii"], n, replace=F)  
  fasta_new <- fasta[names(fasta) %in% id]
  return(fasta_new)
  }

if (length(alignment_in)>400){
  alignment_subsampled <- subsample(fasta=alignment_in, n=400, seed=1234)
} else {
  alignment_subsampled <- alignment_in
}
seqinr::write.fasta(alignment_subsampled,names=names(alignment_subsampled),file.out=snakemake@output[["alignment_out"]])


