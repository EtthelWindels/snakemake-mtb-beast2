## ------------------------------------------------------------------------
## Calculate terminal branch lengths
## NOTE: This script is not part of the Snakemake workflow

## 2024-12-20 Etthel Windels
## ------------------------------------------------------------------------



# Load libraries ----------------------------------------------------------

library(ape)
library(seqinr)
library(ggplot2)
library(cowplot)
library(ggtree)


# Read files --------------------------------------------------------------

Ma_L1 <- read.tree("MLtrees/Malawi_L1.treefile")
Ma_L2 <- read.tree("MLtrees/Malawi_L2.treefile")
Ma_L3 <- read.tree("MLtrees/Malawi_L3.treefile")
Ma_L4 <- read.tree("MLtrees/Malawi_L4.treefile")
Ta_L1 <- read.tree("MLtrees/Tanzania_L1.treefile")
Ta_L2 <- read.tree("MLtrees/Tanzania_L2.treefile")
Ta_L3 <- read.tree("MLtrees/Tanzania_L3.treefile")
Ta_L4 <- read.tree("MLtrees/Tanzania_L4.treefile")
TG_L2 <- read.tree("MLtrees/TheGambia_L2.treefile")
TG_L4 <- read.tree("MLtrees/TheGambia_L4.treefile")
TG_L6 <- read.tree("MLtrees/TheGambia_L6.treefile")
VN_L1 <- read.tree("MLtrees/Vietnam_L1.treefile")
VN_L2 <- read.tree("MLtrees/Vietnam_L2.treefile")
VN_L4 <- read.tree("MLtrees/Vietnam_L4.treefile")



# Calculate terminal branch lenghts per alignment ---------------------------------

genome_length = 4411532
Ma_L1 <- drop.tip(Ma_L1, "Mycobacterium_canettii")
Ma_L2 <- drop.tip(Ma_L2, "Mycobacterium_canettii")
Ma_L3 <- drop.tip(Ma_L3, "Mycobacterium_canettii")
Ma_L4 <- drop.tip(Ma_L4, "Mycobacterium_canettii")
Ta_L1 <- drop.tip(Ta_L1, "Mycobacterium_canettii")
Ta_L2 <- drop.tip(Ta_L2, "Mycobacterium_canettii")
Ta_L3 <- drop.tip(Ta_L3, "Mycobacterium_canettii")
Ta_L4 <- drop.tip(Ta_L4, "Mycobacterium_canettii")
TG_L2 <- drop.tip(TG_L2, "Mycobacterium_canettii")
TG_L4 <- drop.tip(TG_L4, "Mycobacterium_canettii")
TG_L6 <- drop.tip(TG_L6, "Mycobacterium_canettii")
VN_L1 <- drop.tip(VN_L1, "Mycobacterium_canettii")
VN_L2 <- drop.tip(VN_L2, "Mycobacterium_canettii")
VN_L4 <- drop.tip(VN_L4, "Mycobacterium_canettii")

get_tbl <- function(tree){
  num_tips <- Ntip(tree)
  terminal_edges <- tree$edge[, 2] <= num_tips
  terminal_branch_lengths <- tree$edge.length[terminal_edges]*genome_length # branch lengths from IQTREE expressed as subst/site (entire genome)
  return(TBL=terminal_branch_lengths)
}

tbl_data <- data.frame(location = c(rep("Malawi",Ntip(Ma_L1)+Ntip(Ma_L2)+Ntip(Ma_L3)+Ntip(Ma_L4)), rep("Tanzania",Ntip(Ta_L1)+Ntip(Ta_L2)+Ntip(Ta_L3)+Ntip(Ta_L4)), rep("TheGambia",Ntip(TG_L2)+Ntip(TG_L4)+Ntip(TG_L6)), rep("Vietnam",Ntip(VN_L1)+Ntip(VN_L2)+Ntip(VN_L4))),
                              lineage = c(rep("L1",Ntip(Ma_L1)), rep("L2",Ntip(Ma_L2)), rep("L3",Ntip(Ma_L3)), rep("L4",Ntip(Ma_L4)), rep("L1",Ntip(Ta_L1)), rep("L2",Ntip(Ta_L2)), rep("L3",Ntip(Ta_L3)), rep("L4",Ntip(Ta_L4)), rep("L2",Ntip(TG_L2)), rep("L4",Ntip(TG_L4)), rep("L6",Ntip(TG_L6)), rep("L1",Ntip(VN_L1)), rep("L2",Ntip(VN_L2)), rep("L4",Ntip(VN_L4))),
                              tbl = c(get_tbl(Ma_L1), get_tbl(Ma_L2), get_tbl(Ma_L3), get_tbl(Ma_L4), get_tbl(Ta_L1), get_tbl(Ta_L2), get_tbl(Ta_L3), get_tbl(Ta_L4), get_tbl(TG_L2), get_tbl(TG_L4), get_tbl(TG_L6), get_tbl(VN_L1), get_tbl(VN_L2), get_tbl(VN_L4)))

L1_col <- 'darksalmon'
L2_col <- '#8a96a3ff'
L3_col <- '#bfcbdbff'
L4_col <- 'slategray1'
L6_col <- 'indianred4'


Ma_plot <- ggplot(tbl_data[tbl_data$location=="Malawi",], aes(x=lineage, y=tbl, fill=lineage)) +
  geom_violin(aes(col=lineage)) +
  labs(title="Malawi", x="Lineage", y="Terminal branch length (substitutions)") +
  scale_y_continuous(limits=c(0,300), expand = c(0.01, 0.01)) +
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
  scale_fill_manual(values=c(L1_col,L2_col,L3_col,L4_col)) +
  scale_colour_manual(values=c(L1_col,L2_col,L3_col,L4_col))

Ta_plot <- ggplot(tbl_data[tbl_data$location=="Tanzania",], aes(x=lineage, y=tbl, fill=lineage)) +
  geom_violin(aes(col=lineage)) +
  labs(title="Tanzania", x="Lineage", y="Terminal branch length (substitutions)") +
  scale_y_continuous(limits=c(0,300), expand = c(0.01, 0.01)) +
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
  scale_fill_manual(values=c(L1_col,L2_col,L3_col,L4_col)) +
  scale_colour_manual(values=c(L1_col,L2_col,L3_col,L4_col))

TG_plot <- ggplot(tbl_data[tbl_data$location=="TheGambia",], aes(x=lineage, y=tbl, fill=lineage)) +
  geom_violin(aes(col=lineage)) +
  labs(title="The Gambia", x="Lineage", y="Terminal branch length (substitutions)") +
  scale_y_continuous(limits=c(0,300), expand = c(0.01, 0.01)) +
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
  scale_fill_manual(values=c(L2_col,L4_col,L6_col)) +
  scale_colour_manual(values=c(L2_col,L4_col,L6_col))

VN_plot <- ggplot(tbl_data[tbl_data$location=="Vietnam",], aes(x=lineage, y=tbl, fill=lineage)) +
  geom_violin(aes(col=lineage)) +
  labs(title="Vietnam", x="Lineage", y="Terminal branch length (substitutions)") +
  scale_y_continuous(limits=c(0,300), expand = c(0.01, 0.01)) +
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
  scale_fill_manual(values=c(L1_col,L2_col,L4_col)) +
  scale_colour_manual(values=c(L1_col,L2_col,L4_col))


all_plot <- plot_grid(Ma_plot, Ta_plot, TG_plot, VN_plot, ncol=4) 
all_plot

median(tbl_data$tbl[tbl_data$location=="Malawi" & tbl_data$lineage=="L1"])
median(tbl_data$tbl[tbl_data$location=="Malawi" & tbl_data$lineage=="L2"])
median(tbl_data$tbl[tbl_data$location=="Malawi" & tbl_data$lineage=="L3"])
median(tbl_data$tbl[tbl_data$location=="Malawi" & tbl_data$lineage=="L4"])

median(tbl_data$tbl[tbl_data$location=="Tanzania" & tbl_data$lineage=="L1"])
median(tbl_data$tbl[tbl_data$location=="Tanzania" & tbl_data$lineage=="L2"])
median(tbl_data$tbl[tbl_data$location=="Tanzania" & tbl_data$lineage=="L3"])
median(tbl_data$tbl[tbl_data$location=="Tanzania" & tbl_data$lineage=="L4"])

median(tbl_data$tbl[tbl_data$location=="TheGambia" & tbl_data$lineage=="L2"])
median(tbl_data$tbl[tbl_data$location=="TheGambia" & tbl_data$lineage=="L4"])
median(tbl_data$tbl[tbl_data$location=="TheGambia" & tbl_data$lineage=="L6"])

median(tbl_data$tbl[tbl_data$location=="Vietnam" & tbl_data$lineage=="L1"])
median(tbl_data$tbl[tbl_data$location=="Vietnam" & tbl_data$lineage=="L2"])
median(tbl_data$tbl[tbl_data$location=="Vietnam" & tbl_data$lineage=="L4"])



ggsave("tbl.png", plot=all_plot, path="figures/", width=500, height=150, units="mm")
ggsave("tbl.svg", plot=all_plot, path="figures/", width=500, height=150, units="mm")


