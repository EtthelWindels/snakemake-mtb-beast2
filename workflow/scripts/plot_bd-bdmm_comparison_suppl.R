## ------------------------------------------------------------------------
## Plot BD-BDMM comparisons
## 2024-02-23 Etthel Windels
## ------------------------------------------------------------------------



# Load libraries ----------------------------------------------------------

library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(hrbrthemes)
library(ggridges)



# Read input --------------------------------------------------------------

# BD

# Malawi
LOGFILE_BD_Ma_1 <- snakemake@input[["input_logs_bd"]][[1]]
LOGFILE_BD_Ma_2 <- snakemake@input[["input_logs_bd"]][[2]]
LOGFILE_BD_Ma_3 <- snakemake@input[["input_logs_bd"]][[3]]
LOGFILE_BD_Ma_4 <- snakemake@input[["input_logs_bd"]][[4]]

# Tanzania
LOGFILE_BD_Ta_1 <- snakemake@input[["input_logs_bd"]][[5]]
LOGFILE_BD_Ta_2 <- snakemake@input[["input_logs_bd"]][[6]]
LOGFILE_BD_Ta_3 <- snakemake@input[["input_logs_bd"]][[7]]
LOGFILE_BD_Ta_4 <- snakemake@input[["input_logs_bd"]][[8]]

# The Gambia
LOGFILE_BD_TG_2 <- snakemake@input[["input_logs_bd"]][[9]]
LOGFILE_BD_TG_4 <- snakemake@input[["input_logs_bd"]][[10]]
LOGFILE_BD_TG_6 <- snakemake@input[["input_logs_bd"]][[11]]

# Vietnam
LOGFILE_BD_VN_1 <- snakemake@input[["input_logs_bd"]][[12]]
LOGFILE_BD_VN_2 <- snakemake@input[["input_logs_bd"]][[13]]
LOGFILE_BD_VN_4 <- snakemake@input[["input_logs_bd"]][[14]]

# BDMM

# Malawi
LOGFILE_BDMM_Ma_1 <- snakemake@input[["input_logs_bdmm"]][[1]]
LOGFILE_BDMM_Ma_2 <- snakemake@input[["input_logs_bdmm"]][[2]]
LOGFILE_BDMM_Ma_3 <- snakemake@input[["input_logs_bdmm"]][[3]]
LOGFILE_BDMM_Ma_4 <- snakemake@input[["input_logs_bdmm"]][[4]]

# Tanzania
LOGFILE_BDMM_Ta_1 <- snakemake@input[["input_logs_bdmm"]][[5]]
LOGFILE_BDMM_Ta_2 <- snakemake@input[["input_logs_bdmm"]][[6]]
LOGFILE_BDMM_Ta_3 <- snakemake@input[["input_logs_bdmm"]][[7]]
LOGFILE_BDMM_Ta_4 <- snakemake@input[["input_logs_bdmm"]][[8]]

# The Gambia
LOGFILE_BDMM_TG_2 <- snakemake@input[["input_logs_bdmm"]][[9]]
LOGFILE_BDMM_TG_4 <- snakemake@input[["input_logs_bdmm"]][[10]]
LOGFILE_BDMM_TG_6 <- snakemake@input[["input_logs_bdmm"]][[11]]

# Vietnam
LOGFILE_BDMM_VN_1 <- snakemake@input[["input_logs_bdmm"]][[12]]
LOGFILE_BDMM_VN_2 <- snakemake@input[["input_logs_bdmm"]][[13]]
LOGFILE_BDMM_VN_4 <- snakemake@input[["input_logs_bdmm"]][[14]]



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

get_extended_log_BD <- function(path_to_logfile) {
  log <- load_trajectories(path_to_logfile, burninFrac=0) %>%
    as.data.frame()
  log$bUR <- log$becomeUninfectiousRate.2
  log$inf_period_tot <- 1/log$bUR
  log$inf_period <- 1/log$bUR
  log$Re <- log$R0WithinDemes.2
  log$lambda <- log$Re*log$bUR
  log$time_until_transm <- 1/log$lambda
  return(log)
}

get_extended_log_BDMM <- function(path_to_logfile) {
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



# Get posterior densities -------------------------------------------------

# BD - Malawi
log_BD_Ma_1 <- get_extended_log_BD(LOGFILE_BD_Ma_1)
log_BD_Ma_2 <- get_extended_log_BD(LOGFILE_BD_Ma_2)
log_BD_Ma_3 <- get_extended_log_BD(LOGFILE_BD_Ma_3)
log_BD_Ma_4 <- get_extended_log_BD(LOGFILE_BD_Ma_4)

# BD - Tanzania
log_BD_Ta_1 <- get_extended_log_BD(LOGFILE_BD_Ta_1)
log_BD_Ta_2 <- get_extended_log_BD(LOGFILE_BD_Ta_2)
log_BD_Ta_3 <- get_extended_log_BD(LOGFILE_BD_Ta_3)
log_BD_Ta_4 <- get_extended_log_BD(LOGFILE_BD_Ta_4)

# BD - The Gambia
log_BD_TG_2 <- get_extended_log_BD(LOGFILE_BD_TG_2)
log_BD_TG_4 <- get_extended_log_BD(LOGFILE_BD_TG_4)
log_BD_TG_6 <- get_extended_log_BD(LOGFILE_BD_TG_6)

# BD - Vietnam
log_BD_VN_1 <- get_extended_log_BD(LOGFILE_BD_VN_1)
log_BD_VN_2 <- get_extended_log_BD(LOGFILE_BD_VN_2)
log_BD_VN_4 <- get_extended_log_BD(LOGFILE_BD_VN_4)

# BDMM - Malawi
log_BDMM_Ma_1 <- get_extended_log_BDMM(LOGFILE_BDMM_Ma_1)
log_BDMM_Ma_2 <- get_extended_log_BDMM(LOGFILE_BDMM_Ma_2)
log_BDMM_Ma_3 <- get_extended_log_BDMM(LOGFILE_BDMM_Ma_3)
log_BDMM_Ma_4 <- get_extended_log_BDMM(LOGFILE_BDMM_Ma_4)

# BDMM - Tanzania
log_BDMM_Ta_1 <- get_extended_log_BDMM(LOGFILE_BDMM_Ta_1)
log_BDMM_Ta_2 <- get_extended_log_BDMM(LOGFILE_BDMM_Ta_2)
log_BDMM_Ta_3 <- get_extended_log_BDMM(LOGFILE_BDMM_Ta_3)
log_BDMM_Ta_4 <- get_extended_log_BDMM(LOGFILE_BDMM_Ta_4)

# BDMM - The Gambia
log_BDMM_TG_2 <- get_extended_log_BDMM(LOGFILE_BDMM_TG_2)
log_BDMM_TG_4 <- get_extended_log_BDMM(LOGFILE_BDMM_TG_4)
log_BDMM_TG_6 <- get_extended_log_BDMM(LOGFILE_BDMM_TG_6)

# BDMM - Vietnam
log_BDMM_VN_1 <- get_extended_log_BDMM(LOGFILE_BDMM_VN_1)
log_BDMM_VN_2 <- get_extended_log_BDMM(LOGFILE_BDMM_VN_2)
log_BDMM_VN_4 <- get_extended_log_BDMM(LOGFILE_BDMM_VN_4)



# Generate posterior density datasets -------------------------------------

# BD posterior estimates

get_log_length <- function(logfile){
  length <- dim(logfile)[1]
  return(length)
}

n_BD_Ma <- get_log_length(log_BD_Ma_1) + get_log_length(log_BD_Ma_2) + get_log_length(log_BD_Ma_3) + get_log_length(log_BD_Ma_4)
n_BD_Ta <- get_log_length(log_BD_Ta_1) + get_log_length(log_BD_Ta_2) + get_log_length(log_BD_Ta_3) + get_log_length(log_BD_Ta_4)
n_BD_TG <- get_log_length(log_BD_TG_2) + get_log_length(log_BD_TG_4) + get_log_length(log_BD_TG_6)
n_BD_VN <- get_log_length(log_BD_VN_1) + get_log_length(log_BD_VN_2) + get_log_length(log_BD_VN_4)

all_logs_BD <- as.data.frame(rbind(log_BD_Ma_1, log_BD_Ma_2, log_BD_Ma_3, log_BD_Ma_4, 
                                log_BD_Ta_1, log_BD_Ta_2, log_BD_Ta_3, log_BD_Ta_4, 
                                log_BD_TG_2, log_BD_TG_4, log_BD_TG_6,
                                log_BD_VN_1, log_BD_VN_2, log_BD_VN_4))
all_logs_BD$country <- c(rep("Malawi",n_BD_Ma),rep("Tanzania",n_BD_Ta),rep("The Gambia",n_BD_TG),rep("Vietnam",n_BD_VN))
all_logs_BD$lineage <- as.factor(c(rep("L1",get_log_length(log_BD_Ma_1)), rep("L2",get_log_length(log_BD_Ma_2)), rep("L3",get_log_length(log_BD_Ma_3)), rep("L4",get_log_length(log_BD_Ma_4)),
                                rep("L1",get_log_length(log_BD_Ta_1)), rep("L2",get_log_length(log_BD_Ta_2)), rep("L3",get_log_length(log_BD_Ta_3)), rep("L4",get_log_length(log_BD_Ta_4)),
                                rep("L2",get_log_length(log_BD_TG_2)), rep("L4",get_log_length(log_BD_TG_4)), rep("L6",get_log_length(log_BD_TG_6)),
                                rep("L1",get_log_length(log_BD_VN_1)), rep("L2",get_log_length(log_BD_VN_2)), rep("L4",get_log_length(log_BD_VN_4))))
all_logs_BD$lineage <- factor(all_logs_BD$lineage, levels=c("L6","L4","L3","L2","L1"))


# BDMM posterior estimates

n_BDMM_Ma <- get_log_length(log_BDMM_Ma_1) + get_log_length(log_BDMM_Ma_2) + get_log_length(log_BDMM_Ma_3) + get_log_length(log_BDMM_Ma_4)
n_BDMM_Ta <- get_log_length(log_BDMM_Ta_1) + get_log_length(log_BDMM_Ta_2) + get_log_length(log_BDMM_Ta_3) + get_log_length(log_BDMM_Ta_4)
n_BDMM_TG <- get_log_length(log_BDMM_TG_2) + get_log_length(log_BDMM_TG_4) + get_log_length(log_BDMM_TG_6)
n_BDMM_VN <- get_log_length(log_BDMM_VN_1) + get_log_length(log_BDMM_VN_2) + get_log_length(log_BDMM_VN_4)

all_logs_BDMM <- as.data.frame(rbind(log_BDMM_Ma_1, log_BDMM_Ma_2, log_BDMM_Ma_3, log_BDMM_Ma_4, 
                                   log_BDMM_Ta_1, log_BDMM_Ta_2, log_BDMM_Ta_3, log_BDMM_Ta_4, 
                                   log_BDMM_TG_2, log_BDMM_TG_4, log_BDMM_TG_6,
                                   log_BDMM_VN_1, log_BDMM_VN_2, log_BDMM_VN_4))
all_logs_BDMM$country <- c(rep("Malawi",n_BDMM_Ma),rep("Tanzania",n_BDMM_Ta),rep("The Gambia",n_BDMM_TG),rep("Vietnam",n_BDMM_VN))
all_logs_BDMM$lineage <- as.factor(c(rep("L1",get_log_length(log_BDMM_Ma_1)), rep("L2",get_log_length(log_BDMM_Ma_2)), rep("L3",get_log_length(log_BDMM_Ma_3)), rep("L4",get_log_length(log_BDMM_Ma_4)),
                                   rep("L1",get_log_length(log_BDMM_Ta_1)), rep("L2",get_log_length(log_BDMM_Ta_2)), rep("L3",get_log_length(log_BDMM_Ta_3)), rep("L4",get_log_length(log_BDMM_Ta_4)),
                                   rep("L2",get_log_length(log_BDMM_TG_2)), rep("L4",get_log_length(log_BDMM_TG_4)), rep("L6",get_log_length(log_BDMM_TG_6)),
                                   rep("L1",get_log_length(log_BDMM_VN_1)), rep("L2",get_log_length(log_BDMM_VN_2)), rep("L4",get_log_length(log_BDMM_VN_4))))
all_logs_BDMM$lineage <- factor(all_logs_BDMM$lineage, levels=c("L6","L4","L3","L2","L1"))



# Plot posterior densities ------------------------------------------------

plot_Re <- function(dataset, country, x_min=0.5, x_max=1.5, lineage_labels, lineage_colours){
  plot <- ggplot(dataset[dataset$country==country,], aes(lineage, x=Re, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend=F) +
    coord_cartesian(clip = "off") +
    labs(title=country) + 
    xlab(expression(R[e])) +
    ylab("") +
    theme_ipsum(axis_title_just = 'm') +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.y = element_text(size=12, vjust=0.5),
          axis.title.x = element_text(size=12, vjust=0.5)) +
    xlim(x_min,x_max) +
    scale_fill_manual(labels=lineage_labels, values=lineage_colours) +
    scale_colour_manual(labels=lineage_labels, values=lineage_colours)
  return(plot)  
}

plot_time_until_transm <- function(dataset, country, x_min=0, x_max=5, lineage_labels, lineage_colours){
  plot <- ggplot(dataset[dataset$country==country,], aes(lineage, x=time_until_transm, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend = F) +
    coord_cartesian(clip = "off") +
    labs(title=country) + 
    xlab("Time until transmission (years)") +
    ylab("") +
    theme_ipsum(axis_title_just = 'm') +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.y = element_text(size=12, vjust=0.5),
          axis.title.x = element_text(size=12, vjust=0.5)) +
    xlim(x_min,x_max) +
    scale_fill_manual(labels=lineage_labels, values=lineage_colours) +
    scale_colour_manual(labels=lineage_labels, values=lineage_colours)
  return(plot)  
}

plot_infperiod_tot <- function(dataset, country, x_min=0, x_max=5, lineage_labels, lineage_colours){
  plot <- ggplot(dataset[dataset$country==country,], aes(lineage, x=inf_period_tot, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend = F) +
    coord_cartesian(clip = "off") +
    labs(title=country) + 
    xlab("Total infected period (years)") +
    ylab("")+
    theme_ipsum(axis_title_just = 'm') +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.y = element_text(size=12, vjust=0.5),
          axis.title.x = element_text(size=12, vjust=0.5),
          plot.title = element_text(size=14, hjust=0.5)) +
    xlim(x_min,x_max) +
    scale_fill_manual(labels=lineage_labels, values=lineage_colours) +
    scale_colour_manual(labels=lineage_labels, values=lineage_colours)
  return(plot)  
}

plot_transrate <- function(dataset, country, x_min=0, x_max=20, lineage_labels, lineage_colours){
  plot <- ggplot(dataset[dataset$country==country,], aes(lineage, x=lambda, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend=F) +
    coord_cartesian(clip = "off") +
    labs(title=country) + 
    xlab(expression(Transmission~rate~(year^"-1"))) +
    ylab("") +
    theme_ipsum(axis_title_just = 'm') +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.y = element_text(size=12, vjust=0.5),
          axis.title.x = element_text(size=12, vjust=0.5)) +
    xlim(x_min,x_max) +
    scale_fill_manual(labels=lineage_labels, values=lineage_colours) +
    scale_colour_manual(labels=lineage_labels, values=lineage_colours)
  return(plot)  
}

plot_infperiod <- function(dataset, country, x_min=0, x_max=5, lineage_labels, lineage_colours){
  plot <- ggplot(dataset[dataset$country==country,], aes(lineage, x=inf_period, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend=F) +
    coord_cartesian(clip = "off") +
    labs(title=country) + 
    xlab("Infectious period (years)") +
    ylab("") +
    theme_ipsum(axis_title_just = 'm') +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.y = element_text(size=12, vjust=0.5),
          axis.title.x = element_text(size=12, vjust=0.5)) +
    xlim(x_min,x_max) +
    scale_fill_manual(labels=lineage_labels, values=lineage_colours) +
    scale_colour_manual(labels=lineage_labels, values=lineage_colours)
  return(plot)  
}

# Colour code

L1_col <- 'darksalmon'
L2_col <- '#8a96a3ff'
L3_col <- '#bfcbdbff'
L4_col <- 'slategray1'
L6_col <- 'indianred4'

# Reproductive number plot

Re_Ma_BD <- plot_Re(all_logs_BD, "Malawi", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, unit="cm"),
  plot.title = element_text(face="bold", size=14, hjust=0.5))
Re_Ta_BD <- plot_Re(all_logs_BD, "Tanzania", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
Re_TG_BD <- plot_Re(all_logs_BD, "The Gambia", lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
Re_VN_BD <- plot_Re(all_logs_BD, "Vietnam", x_max=1.5, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, r=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))

Re_Ma_BDMM <- plot_Re(all_logs_BDMM, "Malawi", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, l=0.5, unit="cm"),
        plot.title = element_blank())
Re_Ta_BDMM <- plot_Re(all_logs_BDMM, "Tanzania", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"),
        plot.title = element_blank())
Re_TG_BDMM <- plot_Re(all_logs_BDMM, "The Gambia", lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"),
        plot.title = element_blank())
Re_VN_BDMM <- plot_Re(all_logs_BDMM, "Vietnam", x_max=4, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, r=0.5, unit="cm"),
        plot.title = element_blank())

plots_col1 <- align_plots(Re_Ma_BD, Re_Ma_BDMM, align='v', axis='l')
plots_col2 <- align_plots(Re_Ta_BD, Re_Ta_BDMM, align='v', axis='l')
plots_col3 <- align_plots(Re_TG_BD, Re_TG_BDMM, align='v', axis='l')
plots_col4 <- align_plots(Re_VN_BD, Re_VN_BDMM, align='v', axis='l')
plot_a <- plot_grid(plots_col1[[1]], plots_col2[[1]], plots_col3[[1]], plots_col4[[1]], ncol=4, align='h', axis='b', labels='a')
plot_b <- plot_grid(plots_col1[[2]], plots_col2[[2]], plots_col3[[2]], plots_col4[[2]], ncol=4, align='h', axis='b', labels='b')
plot_Re <- plot_grid(plot_a, plot_b, ncol=1, align='v', axis='l', rel_heights=c(1.3,1))

# Time until transmission plot

time_until_transm_Ma_BD <- plot_time_until_transm(all_logs_BD, "Malawi", x_max=3, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
time_until_transm_Ta_BD <- plot_time_until_transm(all_logs_BD, "Tanzania", x_max=4, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
time_until_transm_TG_BD <- plot_time_until_transm(all_logs_BD, "The Gambia", x_max=6, lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
time_until_transm_VN_BD <- plot_time_until_transm(all_logs_BD, "Vietnam", x_max=10, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, r=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))

time_until_transm_Ma_BDMM <- plot_time_until_transm(all_logs_BDMM, "Malawi", x_max=4, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, l=0.5, unit="cm"),
        plot.title = element_blank())
time_until_transm_Ta_BDMM <- plot_time_until_transm(all_logs_BDMM, "Tanzania", x_max=4, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"),
        plot.title = element_blank())
time_until_transm_TG_BDMM <- plot_time_until_transm(all_logs_BDMM, "The Gambia", x_max=6, lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"),
        plot.title = element_blank())
time_until_transm_VN_BDMM <- plot_time_until_transm(all_logs_BDMM, "Vietnam", x_max=18, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, r=0.5, unit="cm"),
        plot.title = element_blank())

plots_col1 <- align_plots(time_until_transm_Ma_BD, time_until_transm_Ma_BDMM, align='v', axis='l')
plots_col2 <- align_plots(time_until_transm_Ta_BD, time_until_transm_Ta_BDMM, align='v', axis='l')
plots_col3 <- align_plots(time_until_transm_TG_BD, time_until_transm_TG_BDMM, align='v', axis='l')
plots_col4 <- align_plots(time_until_transm_VN_BD, time_until_transm_VN_BDMM, align='v', axis='l')
plot_a <- plot_grid(plots_col1[[1]], plots_col2[[1]], plots_col3[[1]], plots_col4[[1]], ncol=4, align='h', axis='b', labels='a')
plot_b <- plot_grid(plots_col1[[2]], plots_col2[[2]], plots_col3[[2]], plots_col4[[2]], ncol=4, align='h', axis='b', labels='b')
plot_time_until_transm <- plot_grid(plot_a, plot_b, ncol=1, align='v', axis='l', rel_heights=c(1.3,1))

# Total infected period plot

infperiod_tot_Ma_BD <- plot_infperiod_tot(all_logs_BD, "Malawi", x_max=4, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
infperiod_tot_Ta_BD <- plot_infperiod_tot(all_logs_BD, "Tanzania", x_max=4, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
infperiod_tot_TG_BD <- plot_infperiod_tot(all_logs_BD, "The Gambia", x_max=6, lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
infperiod_tot_VN_BD <- plot_infperiod_tot(all_logs_BD, "Vietnam", x_max=10, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, r=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))

infperiod_tot_Ma_BDMM <- plot_infperiod_tot(all_logs_BDMM, "Malawi", x_max=4, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, l=0.5, unit="cm"),
        plot.title = element_blank())
infperiod_tot_Ta_BDMM <- plot_infperiod_tot(all_logs_BDMM, "Tanzania", x_max=4, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"),
        plot.title = element_blank())
infperiod_tot_TG_BDMM <- plot_infperiod_tot(all_logs_BDMM, "The Gambia", x_max=6, lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"),
        plot.title = element_blank())
infperiod_tot_VN_BDMM <- plot_infperiod_tot(all_logs_BDMM, "Vietnam", x_max=18, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, r=0.5, unit="cm"),
        plot.title = element_blank())

plots_col1 <- align_plots(infperiod_tot_Ma_BD, infperiod_tot_Ma_BDMM, align='v', axis='l')
plots_col2 <- align_plots(infperiod_tot_Ta_BD, infperiod_tot_Ta_BDMM, align='v', axis='l')
plots_col3 <- align_plots(infperiod_tot_TG_BD, infperiod_tot_TG_BDMM, align='v', axis='l')
plots_col4 <- align_plots(infperiod_tot_VN_BD, infperiod_tot_VN_BDMM, align='v', axis='l')
plot_a <- plot_grid(plots_col1[[1]], plots_col2[[1]], plots_col3[[1]], plots_col4[[1]], ncol=4, align='h', axis='b', labels='a')
plot_b <- plot_grid(plots_col1[[2]], plots_col2[[2]], plots_col3[[2]], plots_col4[[2]], ncol=4, align='h', axis='b', labels='b')
plot_infperiod_tot <- plot_grid(plot_a, plot_b, ncol=1, align='v', axis='l', rel_heights=c(1.3,1))

# Transmission rate and infectious period plot

transrate_Ma_BD <- plot_transrate(all_logs_BD, "Malawi", x_max=2, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
transrate_Ta_BD <- plot_transrate(all_logs_BD, "Tanzania", x_max=2, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
transrate_TG_BD <- plot_transrate(all_logs_BD, "The Gambia", x_max=2, lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))
transrate_VN_BD <- plot_transrate(all_logs_BD, "Vietnam", x_max=2, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, r=0.5, unit="cm"),
        plot.title = element_text(face="bold", size=14, hjust=0.5))

transrate_Ma_BDMM <- plot_transrate(all_logs_BDMM, "Malawi", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, l=0.5, unit="cm"),
        plot.title = element_blank())
transrate_Ta_BDMM <- plot_transrate(all_logs_BDMM, "Tanzania", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"),
        plot.title = element_blank())
transrate_TG_BDMM <- plot_transrate(all_logs_BDMM, "The Gambia", x_max=35, lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"),
        plot.title = element_blank())
transrate_VN_BDMM <- plot_transrate(all_logs_BDMM, "Vietnam", x_max=35, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, r=0.5, unit="cm"),
        plot.title = element_blank())

infperiod_Ma_BD <- plot_infperiod(all_logs_BD, "Malawi", x_max=3, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, unit="cm"),
        plot.title = element_blank())
infperiod_Ta_BD <- plot_infperiod(all_logs_BD, "Tanzania",  x_max=3, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"),
        plot.title = element_blank())
infperiod_TG_BD <- plot_infperiod(all_logs_BD, "The Gambia", x_max=6, lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"),
        plot.title = element_blank())
infperiod_VN_BD <- plot_infperiod(all_logs_BD, "Vietnam", x_max=10, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, r=0.5, unit="cm"),
        plot.title = element_blank())

infperiod_Ma_BDMM <- plot_infperiod(all_logs_BDMM, "Malawi", x_max=2, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, l=0.5, unit="cm"),
        plot.title = element_blank())
infperiod_Ta_BDMM <- plot_infperiod(all_logs_BDMM, "Tanzania", x_max=2, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"),
        plot.title = element_blank())
infperiod_TG_BDMM <- plot_infperiod(all_logs_BDMM, "The Gambia", x_max=2, lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"),
        plot.title = element_blank())
infperiod_VN_BDMM <- plot_infperiod(all_logs_BDMM, "Vietnam", x_max=2, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, r=0.5, unit="cm"),
        plot.title = element_blank())

plots_col1 <- align_plots(transrate_Ma_BD, transrate_Ma_BDMM, infperiod_Ma_BD, infperiod_Ma_BDMM, align='v', axis='l')
plots_col2 <- align_plots(transrate_Ta_BD, transrate_Ta_BDMM, infperiod_Ta_BD, infperiod_Ta_BDMM, align='v', axis='l')
plots_col3 <- align_plots(transrate_TG_BD, transrate_TG_BDMM, infperiod_TG_BD, infperiod_TG_BDMM, align='v', axis='l')
plots_col4 <- align_plots(transrate_VN_BD, transrate_VN_BDMM, infperiod_VN_BD, infperiod_VN_BDMM, align='v', axis='l')
plot_a <- plot_grid(plots_col1[[1]], plots_col2[[1]], plots_col3[[1]], plots_col4[[1]], ncol=4, align='h', axis='b', labels='a')
plot_b <- plot_grid(plots_col1[[2]], plots_col2[[2]], plots_col3[[2]], plots_col4[[2]], ncol=4, align='h', axis='b', labels='b')
plot_c <- plot_grid(plots_col1[[3]], plots_col2[[3]], plots_col3[[3]], plots_col4[[3]], ncol=4, align='h', axis='b', labels='c')
plot_d <- plot_grid(plots_col1[[4]], plots_col2[[4]], plots_col3[[4]], plots_col4[[4]], ncol=4, align='h', axis='b', labels='d')

plot_transrate_infperiod <- plot_grid(plot_a, plot_b, plot_c, plot_d, ncol=1, align='v', axis='l', rel_heights=c(1.3,1,1,1))



# Save output plots -------------------------------------------------------

ggsave("bd-bdmm_comparison_Re.png", plot=plot_Re, path=dirname(snakemake@output[["Re_png"]]), width=250, height=130, units="mm", bg='white')
ggsave("bd-bdmm_comparison_Re.svg", plot=plot_Re, path=dirname(snakemake@output[["Re_svg"]]), width=250, height=130, units="mm", bg='white')
ggsave("bd-bdmm_comparison_time_until_transm.png", plot=plot_Re, path=dirname(snakemake@output[["time_until_transm_png"]]), width=250, height=130, units="mm", bg='white')
ggsave("bd-bdmm_comparison_time_until_transm.svg", plot=plot_Re, path=dirname(snakemake@output[["time_until_transm_svg"]]), width=250, height=130, units="mm", bg='white')
ggsave("bd-bdmm_comparison_infperiod_tot.png", plot=plot_Re, path=dirname(snakemake@output[["infperiod_tot_png"]]), width=250, height=130, units="mm", bg='white')
ggsave("bd-bdmm_comparison_infperiod_tot.svg", plot=plot_Re, path=dirname(snakemake@output[["infperiod_tot_svg"]]), width=250, height=130, units="mm", bg='white')
ggsave("bd-bdmm_comparison_transrate_infperiod.png", plot=plot_Re, path=dirname(snakemake@output[["transrate_infperiod_png"]]), width=250, height=250, units="mm", bg='white')
ggsave("bd-bdmm_comparison_transrate_infperiod.svg", plot=plot_Re, path=dirname(snakemake@output[["transrate_infperiod_svg"]]), width=250, height=250, units="mm", bg='white')
