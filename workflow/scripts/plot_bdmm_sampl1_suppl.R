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

# Malawi
LOGFILE_Ma_1 <- snakemake@input[[1]]
LOGFILE_Ma_2 <- snakemake@input[[2]]
LOGFILE_Ma_3 <- snakemake@input[[3]]
LOGFILE_Ma_4 <- snakemake@input[[4]]

# Tanzania
LOGFILE_Ta_1 <- snakemake@input[[5]]
LOGFILE_Ta_2 <- snakemake@input[[6]]
LOGFILE_Ta_3 <- snakemake@input[[7]]
LOGFILE_Ta_4 <- snakemake@input[[8]]

# The Gambia
LOGFILE_TG_2 <- snakemake@input[[9]]
LOGFILE_TG_4 <- snakemake@input[[10]]
LOGFILE_TG_6 <- snakemake@input[[11]]

# Vietnam
LOGFILE_VN_1 <- snakemake@input[[12]]
LOGFILE_VN_2 <- snakemake@input[[13]]
LOGFILE_VN_4 <- snakemake@input[[14]]



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

# Malawi
log_Ma_1 <- get_extended_log(LOGFILE_Ma_1)
log_Ma_2 <- get_extended_log(LOGFILE_Ma_2)
log_Ma_3 <- get_extended_log(LOGFILE_Ma_3)
log_Ma_4 <- get_extended_log(LOGFILE_Ma_4)

# Tanzania
log_Ta_1 <- get_extended_log(LOGFILE_Ta_1)
log_Ta_2 <- get_extended_log(LOGFILE_Ta_2)
log_Ta_3 <- get_extended_log(LOGFILE_Ta_3)
log_Ta_4 <- get_extended_log(LOGFILE_Ta_4)

# The Gambia
log_TG_2 <- get_extended_log(LOGFILE_TG_2)
log_TG_4 <- get_extended_log(LOGFILE_TG_4)
log_TG_6 <- get_extended_log(LOGFILE_TG_6)

# Vietnam
log_VN_1 <- get_extended_log(LOGFILE_VN_1)
log_VN_2 <- get_extended_log(LOGFILE_VN_2)
log_VN_4 <- get_extended_log(LOGFILE_VN_4)



# Generate posterior density datasets -------------------------------------

# Posterior estimates

get_log_length <- function(logfile){
  length <- dim(logfile)[1]
  return(length)
}

n_Ma <- get_log_length(log_Ma_1) + get_log_length(log_Ma_2) + get_log_length(log_Ma_3) + get_log_length(log_Ma_4)
n_Ta <- get_log_length(log_Ta_1) + get_log_length(log_Ta_2) + get_log_length(log_Ta_3) + get_log_length(log_Ta_4)
n_TG <- get_log_length(log_TG_2) + get_log_length(log_TG_4) + get_log_length(log_TG_6)
n_VN <- get_log_length(log_VN_1) + get_log_length(log_VN_2) + get_log_length(log_VN_4)

all_logs <- as.data.frame(rbind(log_Ma_1, log_Ma_2, log_Ma_3, log_Ma_4, 
                                log_Ta_1, log_Ta_2, log_Ta_3, log_Ta_4, 
                                log_TG_2, log_TG_4, log_TG_6,
                                log_VN_1, log_VN_2, log_VN_4))
all_logs$country <- c(rep("Malawi",n_Ma),rep("Tanzania",n_Ta),rep("The Gambia",n_TG),rep("Vietnam",n_VN))
all_logs$lineage <- as.factor(c(rep("L1",get_log_length(log_Ma_1)), rep("L2",get_log_length(log_Ma_2)), rep("L3",get_log_length(log_Ma_3)), rep("L4",get_log_length(log_Ma_4)),
                                rep("L1",get_log_length(log_Ta_1)), rep("L2",get_log_length(log_Ta_2)), rep("L3",get_log_length(log_Ta_3)), rep("L4",get_log_length(log_Ta_4)),
                                rep("L2",get_log_length(log_TG_2)), rep("L4",get_log_length(log_TG_4)), rep("L6",get_log_length(log_TG_6)),
                                rep("L1",get_log_length(log_VN_1)), rep("L2",get_log_length(log_VN_2)), rep("L4",get_log_length(log_VN_4))))
all_logs$lineage <- factor(all_logs$lineage, levels=c("L6","L4","L3","L2","L1"))

# Distribution of times until infectiousness

x <- seq(0,10,by=0.01)
all_lin_uninf_distr <- data.frame(country=c(rep("Malawi",4*length(x)),rep("Tanzania",4*length(x)),rep("The Gambia",3*length(x)),rep("Vietnam",3*length(x))),
                                  lineage=c(rep("L1",length(x)),rep("L2",length(x)),rep("L3",length(x)),rep("L4",length(x)),rep("L1",length(x)),rep("L2",length(x)),rep("L3",length(x)),rep("L4",length(x)),rep("L2",length(x)),rep("L4",length(x)),rep("L6",length(x)),rep("L1",length(x)),rep("L2",length(x)),rep("L4",length(x))),
                                  exp_x=rep(x,times=14),
                                  exp_y=c(dexp(x,rate=get_migr_mean(log_Ma_1)),dexp(x,rate=get_migr_mean(log_Ma_2)),dexp(x,rate=get_migr_mean(log_Ma_3)),dexp(x,rate=get_migr_mean(log_Ma_4)),dexp(x,rate=get_migr_mean(log_Ta_1)),dexp(x,rate=get_migr_mean(log_Ta_2)),dexp(x,rate=get_migr_mean(log_Ta_3)),dexp(x,rate=get_migr_mean(log_Ta_4)),dexp(x,rate=get_migr_mean(log_TG_2)),dexp(x,rate=get_migr_mean(log_TG_4)),dexp(x,rate=get_migr_mean(log_TG_6)),dexp(x,rate=get_migr_mean(log_VN_1)),dexp(x,rate=get_migr_mean(log_VN_2)),dexp(x,rate=get_migr_mean(log_VN_4))))



# Plot posterior densities ------------------------------------------------

plot_Re <- function(country, x_min=0.5, x_max=1.5, lineage_labels, lineage_colours){
  plot <- ggplot(all_logs[all_logs$country==country,], aes(lineage, x=Re, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend=F) +
    coord_cartesian(clip = "off") +
    labs(title=country) + 
    xlab(expression(R[e])) +
    ylab("") +
    theme_ipsum(axis_title_just = 'm') +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.y = element_text(size=12, vjust=0.5),
          axis.title.x = element_text(size=12, vjust=0.5),
          plot.title = element_text(face="bold", size=14, hjust=0.5)) +
    xlim(x_min,x_max) +
    scale_fill_manual(labels=lineage_labels, values=lineage_colours) +
    scale_colour_manual(labels=lineage_labels, values=lineage_colours)
  return(plot)  
}

plot_latentperiod <- function(country, x_min=0, x_max=4, lineage_labels, lineage_colours){
  plot <- ggplot(all_logs[all_logs$country==country,], aes(lineage, x=latent_period, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend=F) +
    coord_cartesian(clip = "off") +
    xlab("Infected non-infectious period (years)") +
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

plot_distr <- function(country, x_min=0, x_max=6, lineage_labels, lineage_colours){
  plot <- ggplot(all_lin_uninf_distr[all_lin_uninf_distr$country==country,], aes(lineage, x=exp_x, y=exp_y, colour=lineage)) +
    geom_line(show.legend=T, linewidth=1) +
    xlab("Time until infectiousness (years)") +
    ylab("Density")+
    theme_ipsum(axis_title_just = 'm') +
    theme(axis.text.y = element_text(size=10),
          axis.text.x = element_text(size=10),
          axis.title.y = element_text(size=12, vjust=0.5),
          axis.title.x = element_text(size=12, vjust=0.5),
          legend.position = c(0.8,0.75),
          legend.background = element_rect(fill="white", colour="white"),
          legend.title = element_blank(),
          legend.text = element_text(size=10, vjust=0.5)) +
    xlim(x_min,x_max) +
    ylim(0,1) +
    scale_fill_manual(labels=lineage_labels, values=lineage_colours) +
    scale_colour_manual(labels=lineage_labels, values=lineage_colours)
  return(plot)  
}

plot_transrate <- function(country, x_min=0, x_max=20, lineage_labels, lineage_colours){
  plot <- ggplot(all_logs[all_logs$country==country,], aes(lineage, x=lambda, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend = F) +
    coord_cartesian(clip = "off") +
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

plot_infperiod <- function(country, x_min=0, x_max=1.5, lineage_labels, lineage_colours){
  plot <- ggplot(all_logs[all_logs$country==country,], aes(lineage, x=inf_period, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend = F) +
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

plot_time_until_transm <- function(country, x_min=0, x_max=5, lineage_labels, lineage_colours){
  plot <- ggplot(all_logs[all_logs$country==country,], aes(lineage, x=time_until_transm, fill=lineage, colour=lineage)) +
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

plot_infperiod_tot <- function(country, x_min=0, x_max=5, lineage_labels, lineage_colours){
  plot <- ggplot(all_logs[all_logs$country==country,], aes(lineage, x=inf_period_tot, fill=lineage, colour=lineage)) +
    geom_density_ridges(alpha=0.9, show.legend = F) +
    coord_cartesian(clip = "off") +
    labs(title=expression(atop(bold("Malawi")))) + 
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

# Colour code

L1_col <- 'darksalmon'
L2_col <- '#8a96a3ff'
L3_col <- '#bfcbdbff'
L4_col <- 'slategray1'
L6_col <- 'indianred4'

# Reproductive number plot

Re_Ma <- plot_Re("Malawi", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, l=0.5, unit="cm"))
Re_Ta <- plot_Re("Tanzania", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"))
Re_TG <- plot_Re("The Gambia", lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, unit="cm"))
Re_VN <- plot_Re("Vietnam", x_max=4, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(t=0.5, b=0.5, r=0.5, unit="cm"))

# Infected uninfectious period plot

uninfperiod_Ma <- plot_latentperiod("Malawi", x_max=7, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, l=0.5, unit="cm"))
uninfperiod_Ta <- plot_latentperiod("Tanzania", x_max=6, lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"))
uninfperiod_TG <- plot_latentperiod("The Gambia", x_max=10, lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm")) 
uninfperiod_VN <- plot_latentperiod("Vietnam", x_max=30, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, r=0.5, unit="cm"))

# Distribution of times until infectiousness plot

distr_Ma <- plot_distr("Malawi", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L1_col,L2_col,L3_col,L4_col)) +
  theme(plot.margin = margin(b=0.5, l=0.5, unit="cm"))
distr_Ta <- plot_distr("Tanzania", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L1_col,L2_col,L3_col,L4_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"))
distr_TG <- plot_distr("The Gambia", lineage_labels=c("L2","L4","L6"), lineage_colours=c(L2_col,L4_col,L6_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm")) 
distr_VN <- plot_distr("Vietnam", lineage_labels=c("L1","L2","L4"), lineage_colours=c(L1_col,L2_col,L4_col)) +
  theme(plot.margin = margin(b=0.5, r=0.5, unit="cm"))

# Infectious period plot

infperiod_Ma <- plot_infperiod("Malawi", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, l=0.5, unit="cm"))
infperiod_Ta <- plot_infperiod("Tanzania", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"))
infperiod_TG <- plot_infperiod("The Gambia", lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm")) 
infperiod_VN <- plot_infperiod("Vietnam", lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, r=0.5, unit="cm"))

# Transmission rate plot

transrate_Ma <- plot_transrate("Malawi", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, l=0.5, unit="cm"))
transrate_Ta <- plot_transrate("Tanzania", lineage_labels=c("L1","L2","L3","L4"), lineage_colours=c(L4_col,L3_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm"))
transrate_TG <- plot_transrate("The Gambia", x_max=50, lineage_labels=c("L2","L4","L6"), lineage_colours=c(L6_col,L4_col,L2_col)) +
  theme(plot.margin = margin(b=0.5, unit="cm")) 
transrate_VN <- plot_transrate("Vietnam", x_max=50, lineage_labels=c("L1","L2","L4"), lineage_colours=c(L4_col,L2_col,L1_col)) +
  theme(plot.margin = margin(b=0.5, r=0.5, unit="cm"))

# Combined plot

plots_col1 <- align_plots(Re_Ma, uninfperiod_Ma, distr_Ma, infperiod_Ma, transrate_Ma, align='v', axis='l')
plots_col2 <- align_plots(Re_Ta, uninfperiod_Ta, distr_Ta, infperiod_Ta, transrate_Ta, align='v', axis='l')
plots_col3 <- align_plots(Re_TG, uninfperiod_TG, distr_TG, infperiod_TG, transrate_TG, align='v', axis='l')
plots_col4 <- align_plots(Re_VN, uninfperiod_VN, distr_VN, infperiod_VN, transrate_VN, align='v', axis='l')
plot_a <- plot_grid(plots_col1[[1]], plots_col2[[1]], plots_col3[[1]], plots_col4[[1]], ncol=4, align='h', axis='b', labels='a')
plot_b <- plot_grid(plots_col1[[2]], plots_col2[[2]], plots_col3[[2]], plots_col4[[2]], ncol=4, align='h', axis='b', labels='b')
plot_c <- plot_grid(plots_col1[[3]], plots_col2[[3]], plots_col3[[3]], plots_col4[[3]], ncol=4, align='h', axis='b', labels='c')
plot_d <- plot_grid(plots_col1[[4]], plots_col2[[4]], plots_col3[[4]], plots_col4[[4]], ncol=4, align='h', axis='b', labels='d')
plot_e <- plot_grid(plots_col1[[5]], plots_col2[[5]], plots_col3[[5]], plots_col4[[5]], ncol=4, align='h', axis='b', labels='e')
plot_all <- plot_grid(plot_a, plot_b, plot_c, plot_d, plot_e, ncol=1, align='v', axis='l', rel_heights=c(1.3,1,1,1,1))



# Save output plots -------------------------------------------------------

ggsave("bdmm_sampl1.png", plot=plot_all, path=dirname(snakemake@output[["png"]]), width=250, height=300, units="mm", bg='white')
ggsave("bdmm_sampl1.svg", plot=plot_all, path=dirname(snakemake@output[["svg"]]), width=250, height=300, units="mm", bg='white')