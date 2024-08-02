## ------------------------------------------------------------------------
## Adjust metadata and alignment sequence names.
## Replace names of sequences with complete names for 
## analysis (GNUMBER/date/type) 
## 2022-02-02 Etthel Windels
## ------------------------------------------------------------------------


# Load libraries ----------------------------------------------------------

library(seqinr)



# Read files --------------------------------------------------------------

alignment <- seqinr::read.fasta(file = snakemake@input[["alignment_in"]], seqtype="DNA", forceDNAtolower = F)
metadata <- read.table(file = snakemake@input[["meta_all"]], header=T, sep= '\t')



# Adjust alignment --------------------------------------------------------

# 1. Filter alignment to remove outgroup G00157 (which is not present in metadata)
alignment <- 
  alignment[names(alignment) %in% metadata$G_NUMBER]

# 2. Adjust alignment names
gnumber <- names(alignment)
full_names <- metadata$new_id[match(gnumber, metadata$G_NUMBER)]
print(full_names)
if (any(is.na(full_names))) {
  stop("Not all strains have full names in metadata.")
}
names(alignment) <- full_names



# Save output files -------------------------------------------------------

seqinr::write.fasta(alignment, names = names(alignment), file.out = snakemake@output[["alignment_out"]])