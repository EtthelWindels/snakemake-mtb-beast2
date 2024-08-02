## ------------------------------------------------------------------------
## Plot BDMM results
## 2024-01-16 Etthel Windels
## ------------------------------------------------------------------------


# Load libraries ----------------------------------------------------------

library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(hrbrthemes)
library(ggridges)



# Read input --------------------------------------------------------------

# 1 : L1 10-8 - L2 6.4x10-7
LOGFILE_Ta_1_L1 <- snakemake@input[[1]]
LOGFILE_Ta_1_L2 <- snakemake@input[[2]]

# 2 : L1 2x10-8 - L2 3.2x10-7
LOGFILE_Ta_2_L1 <- snakemake@input[[3]]
LOGFILE_Ta_2_L2 <- snakemake@input[[4]]

# 3 : L1 4x10-8 - L2 1.6x10-7
LOGFILE_Ta_3_L1 <- snakemake@input[[5]]
LOGFILE_Ta_3_L2 <- snakemake@input[[6]]

# 4 : L1 8x10-8 - L2 8x10-8
LOGFILE_Ta_4_L1 <- snakemake@input[[7]]
LOGFILE_Ta_4_L2 <- snakemake@input[[8]]

# 5 : L1 1.6x10-7 - L2 4x10-8
LOGFILE_Ta_5_L1 <- snakemake@input[[9]]
LOGFILE_Ta_5_L2 <- snakemake@input[[10]]

# 6 : L1 3.2x10-8 - L2 2x10-7
LOGFILE_Ta_6_L1 <- snakemake@input[[11]]
LOGFILE_Ta_6_L2 <- snakemake@input[[12]]

# 7 : L1 6.4x10-7 - L2 10-8
LOGFILE_Ta_7_L1 <- snakemake@input[[13]]
LOGFILE_Ta_7_L2 <- snakemake@input[[14]]



# Data loading functions --------------------------------------------------

load_trajectories <- function(filename, burninFrac=0, subsample=NA) {
  df_in <- as.matrix(read.table(filename, header=T))
  if (burninFrac>0) {
    n <- dim(df_in)[1]
    df_in <- df_in[-(1:ceiling(burninFrac*n)),]
  }
  if (!is.na(subsample)) {
    indices <- unique(round(seq(1, dim(df_in)[1], length.out=subsample)))
    df_in <- df_in[indices,]
  }
  return(df_in)
}

get_extended_log <- function(path_to_logfile) {
  log <- load_trajectories(path_to_logfile, burninFrac=0) %>%
    as.data.frame()
  log$bUR <- log$becomeUninfectiousRate.2
  log$inf_period <- 1/log$bUR
  log$Re <- log$R0AmongDemes.2
  log$lambda <- log$Re*log$bUR 
  log$migr <- log$rateMatrix.1
  log$latent_period <- 1/log$migr
  log$inf_period_tot <- log$latent_period + log$inf_period
  log$time_until_transm <- log$latent_period + 1/log$lambda
  return(log)
}

get_migr_mean <- function(extended_logfile){
  post_mean <- mean(extended_logfile$migr)
  return(post_mean)
}



# Get posterior densities -------------------------------------------------

# 1
log_Ta_1_L1 <- get_extended_log(LOGFILE_Ta_1_L1)
log_Ta_1_L2 <- get_extended_log(LOGFILE_Ta_1_L2)
# 2
log_Ta_2_L1 <- get_extended_log(LOGFILE_Ta_2_L1)
log_Ta_2_L2 <- get_extended_log(LOGFILE_Ta_2_L2)
# 3
log_Ta_3_L1 <- get_extended_log(LOGFILE_Ta_3_L1)
log_Ta_3_L2 <- get_extended_log(LOGFILE_Ta_3_L2)
# 4
log_Ta_4_L1 <- get_extended_log(LOGFILE_Ta_4_L1)
log_Ta_4_L2 <- get_extended_log(LOGFILE_Ta_4_L2)
# 5
log_Ta_5_L1 <- get_extended_log(LOGFILE_Ta_5_L1)
log_Ta_5_L2 <- get_extended_log(LOGFILE_Ta_5_L2)
# 6
log_Ta_6_L1 <- get_extended_log(LOGFILE_Ta_6_L1)
log_Ta_6_L2 <- get_extended_log(LOGFILE_Ta_6_L2)
# 7
log_Ta_7_L1 <- get_extended_log(LOGFILE_Ta_7_L1)
log_Ta_7_L2 <- get_extended_log(LOGFILE_Ta_7_L2)



# Generate posterior density datasets -------------------------------------

# Posterior estimates

get_log_length <- function(logfile){
  length <- dim(logfile)[1]
  return(length)
}

n_1 <- get_log_length(log_Ta_1_L1) + get_log_length(log_Ta_1_L2)
n_2 <- get_log_length(log_Ta_2_L1) + get_log_length(log_Ta_2_L2)
n_3 <- get_log_length(log_Ta_3_L1) + get_log_length(log_Ta_3_L2)
n_4 <- get_log_length(log_Ta_4_L1) + get_log_length(log_Ta_4_L2)
n_5 <- get_log_length(log_Ta_5_L1) + get_log_length(log_Ta_5_L2)
n_6 <- get_log_length(log_Ta_6_L1) + get_log_length(log_Ta_6_L2)
n_7 <- get_log_length(log_Ta_7_L1) + get_log_length(log_Ta_7_L2)

all_logs <- as.data.frame(rbind(log_Ta_1_L1, log_Ta_1_L2, log_Ta_2_L1, log_Ta_2_L2, log_Ta_3_L1, log_Ta_3_L2, log_Ta_4_L1, log_Ta_4_L2,
                                log_Ta_5_L1, log_Ta_5_L2, log_Ta_6_L1, log_Ta_6_L2, log_Ta_7_L1, log_Ta_7_L2))
all_logs$number <- c(rep("1",n_1),rep("2",n_2),rep("3",n_3),rep("4",n_4),rep("5",n_5),rep("6",n_6),rep("7",n_7))
all_logs$lineage <- as.factor(c(rep("L1",get_log_length(log_Ta_1_L1)), rep("L2",get_log_length(log_Ta_1_L2)),
                                rep("L1",get_log_length(log_Ta_2_L1)), rep("L2",get_log_length(log_Ta_2_L2)),
                                rep("L1",get_log_length(log_Ta_3_L1)), rep("L2",get_log_length(log_Ta_3_L2)),
                                rep("L1",get_log_length(log_Ta_4_L1)), rep("L2",get_log_length(log_Ta_4_L2)),
                                rep("L1",get_log_length(log_Ta_5_L1)), rep("L2",get_log_length(log_Ta_5_L2)),
                                rep("L1",get_log_length(log_Ta_6_L1)), rep("L2",get_log_length(log_Ta_6_L2)),
                                rep("L1",get_log_length(log_Ta_7_L1)), rep("L2",get_log_length(log_Ta_7_L2))))
all_logs$lineage <- factor(all_logs$lineage, levels=c("L2","L1"))



# Plot posterior densities ------------------------------------------------

plot_latentperiod <- function(number, x_min=0, x_max=15, lineage_labels, lineage_colours){
  plot <- ggplot(all_logs[all_logs$number==number,], aes(lineage, x=latent_period, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend=F, scale=3) +
    xlab("Infected non-infectious period (years)") +
    ylab("") +
    theme_ipsum(axis_title_just = 'm') +
    theme(axis.text.y = element_text(size=10, colour='black'),
          axis.text.x = element_text(size=10, colour='black'),
          axis.title.y = element_text(size=12, vjust=0.5),
          axis.title.x = element_blank(),
          axis.ticks.length.x = unit(0.1, "cm"),
          panel.margin.x = NULL,
          panel.margin.y = NULL) +
    xlim(x_min,x_max) +
    scale_fill_manual(labels=lineage_labels, values=lineage_colours) +
    scale_colour_manual(labels=lineage_labels, values=lineage_colours)
  return(plot)  
}

# Colour code

L1_col <- 'darksalmon'
L2_col <- '#8a96a3ff'

# Infected uninfectious period plot

uninfperiod_1 <- plot_latentperiod("1", lineage_labels=c("L1","L2"), lineage_colours=c(L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.2, l=0.5, r=0.5, unit="cm"))
uninfperiod_2 <- plot_latentperiod("2", lineage_labels=c("L1","L2"), lineage_colours=c(L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.2, l=0.5, r=0.5, unit="cm"))
uninfperiod_3 <- plot_latentperiod("3", lineage_labels=c("L1","L2"), lineage_colours=c(L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.2, l=0.5, r=0.5, unit="cm"))
uninfperiod_4 <- plot_latentperiod("4", lineage_labels=c("L1","L2"), lineage_colours=c(L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.2, l=0.5, r=0.5, unit="cm"))
uninfperiod_5 <- plot_latentperiod("5", lineage_labels=c("L1","L2"), lineage_colours=c(L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.2, l=0.5, r=0.5, unit="cm"))
uninfperiod_6 <- plot_latentperiod("6", lineage_labels=c("L1","L2"), lineage_colours=c(L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.2, l=0.5, r=0.5, unit="cm"))
uninfperiod_7 <- plot_latentperiod("7", lineage_labels=c("L1","L2"), lineage_colours=c(L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, l=0.5, r=0.5, unit="cm"))

# Combined plot

plot_all <- plot_grid(uninfperiod_1, uninfperiod_2, uninfperiod_3, uninfperiod_4, uninfperiod_5, uninfperiod_6, uninfperiod_7, ncol=1, align='v', axis='l', rel_heights=c(1,1,1,1,1,1,1))



# Save output plots -------------------------------------------------------

ggsave("bdmm_clockgrid.png", plot=plot_all, path=dirname(snakemake@output[["png"]]), width=75, height=200, units="mm", bg='white')
ggsave("bdmm_clockgrid.svg", plot=plot_all, path=dirname(snakemake@output[["svg"]]), width=75, height=200, units="mm", bg='white')
