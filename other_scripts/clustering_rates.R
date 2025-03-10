## ------------------------------------------------------------------------
## Calculate clustering rates
## NOTE: This script is not part of the Snakemake workflow

## 2024-12-20 Etthel Windels
## ------------------------------------------------------------------------



# Load libraries ----------------------------------------------------------

library(ape)
library(cluster)
library(seqinr)
library(adegenet)
library(ggplot2)
library(cowplot)


# Read files --------------------------------------------------------------

Ma_L1 <- "alignments/Malawi/Malawi_L1_outgr_var_aln.fasta"
Ma_L2 <- "alignments/Malawi/Malawi_L2_outgr_var_aln.fasta"
Ma_L3 <- "alignments/Malawi/Malawi_L3_outgr_var_aln.fasta"
Ma_L4 <- "alignments/Malawi/Malawi_L4_outgr_var_aln.fasta"
Ta_L1 <- "alignments/Tanzania/Tanzania_L1_outgr_var_aln.fasta"
Ta_L2 <- "alignments/Tanzania/Tanzania_L2_outgr_var_aln.fasta"
Ta_L3 <- "alignments/Tanzania/Tanzania_L3_outgr_var_aln.fasta"
Ta_L4 <- "alignments/Tanzania/Tanzania_L4_outgr_var_aln.fasta"
TG_L2 <- "alignments/TheGambia/TheGambia_L2_outgr_var_aln.fasta"
TG_L4 <- "alignments/TheGambia/TheGambia_L4_outgr_var_aln.fasta"
TG_L6 <- "alignments/TheGambia/TheGambia_L6_outgr_var_aln.fasta"
VN_L1 <- "alignments/Vietnam/Vietnam_L1_unique.fasta"
VN_L2 <- "alignments/Vietnam/Vietnam_L2_outgr_var_aln.fasta"
VN_L4 <- "alignments/Vietnam/Vietnam_L4_outgr_var_aln.fasta"



# Calculate clustering rate per alignment ---------------------------------

# Calculate SNP distance matrix

get_distance_matrix <- function(path_to_alignment){
  fasta <- seqinr::read.fasta(path_to_alignment, seqtype='DNA',forceDNAtolower=F)
  length_align <- getLength(fasta)[1]
  DNAbin <- fasta2DNAbin(file=path_to_alignment,chunkSize=100)
  d_TN93 <- dist.dna(DNAbin, model="TN93",as.matrix=T)*length_align
  return(d_TN93)
}

# Define transmission clusters

get_clusters <- function(distance_matrix, snp_threshold){
  dataset_clustering <- agnes(distance_matrix, diss=T, method='average') # UPGMA method
  clusters <- as.data.frame(cutree(as.hclust(dataset_clustering), h = snp_threshold))
  colnames(clusters) <- "cluster_id"
  return(clusters)
}

# Calculate clustering rate

get_clustering_rate <- function(clusters){
  n <- dim(clusters)[1] - 1
  # Get cluster_ids that are "duplicated" (i.e. they have more than one patient)
  clusterIDs <- unique(subset(clusters, duplicated(clusters$cluster_id))$cluster_id)
  # Get samples within these "duplicated" clusters
  clusters_new <- subset(clusters, cluster_id %in% clusterIDs)
  clustering_rate <- dim(clusters_new)[1]/n
  return(clustering_rate)
}

# Clustering rate per alignment

alignment_to_clustering_rate <- function(path_to_alignment, snp_threshold){
  dist <- get_distance_matrix(path_to_alignment)
  clusters <- get_clusters(dist,snp_threshold)
  clust_rate <- get_clustering_rate(clusters)
  return(clust_rate)
}



# Calculate all clustering rates ------------------------------------------

# Malawi
alignment_to_clustering_rate(Ma_L1, 5)
alignment_to_clustering_rate(Ma_L2, 5)
alignment_to_clustering_rate(Ma_L3, 5)
alignment_to_clustering_rate(Ma_L4, 5)

# Tanzania
alignment_to_clustering_rate(Ta_L1, 5)
alignment_to_clustering_rate(Ta_L2, 5)
alignment_to_clustering_rate(Ta_L3, 5)
alignment_to_clustering_rate(Ta_L4, 5)

# The Gambia
alignment_to_clustering_rate(TG_L2, 5)
alignment_to_clustering_rate(TG_L4, 5)
alignment_to_clustering_rate(TG_L6, 5)

# Vietnam
alignment_to_clustering_rate(VN_L1, 5)
alignment_to_clustering_rate(VN_L2, 5)
alignment_to_clustering_rate(VN_L4, 5)

clustering_data <- data.frame(location = c(rep("Malawi",4),rep("Tanzania",4),rep("The Gambia",3),rep("Vietnam",3)),
                              lineage = c("L1", "L2", "L3", "L4", "L1", "L2", "L3", "L4", "L2", "L4", "L6", "L1", "L2", "L4"),
                              clust_rate = c(alignment_to_clustering_rate(Ma_L1,5),alignment_to_clustering_rate(Ma_L2,5),alignment_to_clustering_rate(Ma_L3,5),alignment_to_clustering_rate(Ma_L4,5),
                                             alignment_to_clustering_rate(Ta_L1,5),alignment_to_clustering_rate(Ta_L2,5),alignment_to_clustering_rate(Ta_L3,5),alignment_to_clustering_rate(Ta_L4,5),
                                             alignment_to_clustering_rate(TG_L2,5),alignment_to_clustering_rate(TG_L4,5),alignment_to_clustering_rate(TG_L6,5),
                                             alignment_to_clustering_rate(VN_L1,5),alignment_to_clustering_rate(VN_L2,5),alignment_to_clustering_rate(VN_L4,5)))

ggplot(clustering_data, x=lineage, y=clust_rate) +
  geom_point(aes(x=lineage, y=clust_rate, col=location))
     
L1_col <- 'darksalmon'
L2_col <- '#8a96a3ff'
L3_col <- '#bfcbdbff'
L4_col <- 'slategray1'
L6_col <- 'indianred4'


Ma_plot <- ggplot(clustering_data[clustering_data$location=="Malawi",], aes(x=lineage, y=clust_rate, fill=lineage)) +
  geom_bar(stat='identity') +
  labs(title="Malawi", x="Lineage", y="Proportion of clustered cases") +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  theme_classic() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text.x = element_text(hjust=0.5, size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.2,"cm"),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none',
        plot.title = element_text(size=25, hjust=0.5, face='bold')) +
  scale_fill_manual(values=c(L1_col,L2_col,L3_col,L4_col))
  

Ta_plot <- ggplot(clustering_data[clustering_data$location=="Tanzania",], aes(x=lineage, y=clust_rate, fill=lineage)) +
  geom_bar(stat='identity') +
  labs(title="Tanzania", x="Lineage", y="Proportion of clustered cases") +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  theme_classic() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text.x = element_text(hjust=0.5, size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.2,"cm"),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none',
        plot.title = element_text(size=25, hjust=0.5, face='bold')) +
  scale_fill_manual(values=c(L1_col,L2_col,L3_col,L4_col))


TG_plot <- ggplot(clustering_data[clustering_data$location=="The Gambia",], aes(x=lineage, y=clust_rate, fill=lineage)) +
  geom_bar(stat='identity') +
  labs(title="The Gambia", x="Lineage", y="Proportion of clustered cases") +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  theme_classic() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text.x = element_text(hjust=0.5, size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.2,"cm"),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none',
        plot.title = element_text(size=25, hjust=0.5, face='bold')) +
  scale_fill_manual(values=c(L2_col,L4_col,L6_col))


VN_plot <- ggplot(clustering_data[clustering_data$location=="Vietnam",], aes(x=lineage, y=clust_rate, fill=lineage)) +
  geom_bar(stat='identity') +
  labs(title="Vietnam", x="Lineage", y="Proportion of clustered cases") +
  scale_y_continuous(limits=c(0,1), expand = c(0, 0)) +
  theme_classic() +
  theme(axis.line = element_line(color="black"),
        axis.ticks = element_line(color="black"),
        axis.text.x = element_text(hjust=0.5, size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title.x = element_blank(),
        axis.ticks.length = unit(0.2,"cm"),
        axis.title.y = element_text(size=18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position = 'none',
        plot.title = element_text(size=25, hjust=0.5, face='bold')) + 
  scale_fill_manual(values=c(L1_col,L2_col,L4_col))


all_plot <- plot_grid(Ma_plot, Ta_plot, TG_plot, VN_plot, ncol=4) 
all_plot

ggsave("clustering.png", plot=all_plot, path="figures/", width=500, height=144, units="mm")
ggsave("clustering.svg", plot=all_plot, path="figures/", width=500, height=144, units="mm")
