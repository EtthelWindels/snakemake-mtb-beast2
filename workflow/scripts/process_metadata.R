## ------------------------------------------------------------------------
## Process metadata
## 2023-07-21 Etthel Windels
## ------------------------------------------------------------------------


# Load libraries ----------------------------------------------------------

library(dplyr)
library(stringr)



# Read files --------------------------------------------------------------

Malawi_raw <- read.table(file = snakemake@input[["meta_Ma"]], header=T,sep= '\t')
Tanzania_raw <- read.table(file = snakemake@input[["meta_Ta"]], header=T,sep= '\t')
TheGambia_raw <- read.table(file = snakemake@input[["meta_TG"]], header=T,sep= '\t')
Vietnam_raw <- read.table(file = snakemake@input[["meta_VN"]], header=T,sep= '\t')



# Process datasets --------------------------------------------------------

# Malawi

Malawi <- Malawi_raw %>%
  # 1. Add country column
  mutate(country=rep("Malawi",dim(Malawi_raw)[1])) %>%
  # 2. Rename (sub)lineage columns
  mutate(lineage=LINEAGE) %>%
  mutate(sublineage=SUBLINEAGE_COLL) %>%
  # 3. Adjust date format
  mutate(isolation_date = ifelse(str_split_fixed(Malawi_raw$date,pattern="\\.",3)[,3] > 90, 
                       paste(str_split_fixed(date,pattern="\\.",3)[,1],str_split_fixed(date,pattern="\\.",3)[,2],paste("19",str_split_fixed(date,pattern="\\.",3)[,3],sep="") , sep='-'),
                       paste(str_split_fixed(date,pattern="\\.",3)[,1],str_split_fixed(date,pattern="\\.",3)[,2],paste("20",str_split_fixed(date,pattern="\\.",3)[,3],sep="") , sep='-'))) %>%
    # 4. Add new_id column
  mutate(new_id = paste(G_NUMBER, isolation_date, rep("1",dim(Malawi_raw)[1]), sep = "/"))

# Tanzania

Tanzania <- Tanzania_raw %>%
  # 1. Add country column
  mutate(country=rep("Tanzania",dim(Tanzania_raw)[1])) %>%
  # 2. Rename (sub)lineage columns
  mutate(lineage=Lineage) %>%
  mutate(sublineage=Sublineage) %>%
  # 3. Adjust date format
  mutate(isolation_date = str_split_fixed(ids_date_interview,' ',2)[,1]) %>% 
  mutate(isolation_date=paste(str_split_fixed(isolation_date,pattern="\\.",3)[,1],str_split_fixed(isolation_date,pattern="\\.",3)[,2],paste("20",str_split_fixed(isolation_date,pattern="\\.",3)[,3],sep=""), sep='-')) %>%
    # 4. Add new_id column
  mutate(new_id = paste(G_NUMBER, isolation_date, rep("1",dim(Tanzania_raw)[1]), sep = "/"))

# The Gambia

TheGambia <- TheGambia_raw %>%
  # 1. Add country column
  mutate(country=rep("TheGambia",dim(TheGambia_raw)[1])) %>%
  # 2. Rename (sub)lineage columns
  mutate(lineage=LINEAGE_PROC) %>%
  mutate(sublineage=SUBLINEAGE) %>%
  # 3. Add new_id column
  mutate(new_id = paste(G_NUMBER, isolation_date, rep("1",dim(TheGambia_raw)[1]), sep = "/"))

# Vietnam

Vietnam <- Vietnam_raw %>%
  # 1. Add country column
  mutate(country=rep("Vietnam",dim(Vietnam_raw)[1])) %>%
  # 2. Rename (sub)lineage columns
  mutate(lineage=LINEAGE) %>%
  mutate(sublineage=SUBLINEAGE_COLL) %>%
  # 3. Adjust date format
  # Only month and year available -> set all dates to 15th of the month
  mutate(isolation_date = ifelse(Vietnam_raw$Year < 10,
                       paste(paste("15",ifelse(Vietnam_raw$Month<10,paste0("0",Month),Month),"200",sep="-"), Year, sep=""),
                       paste(paste("15",ifelse(Vietnam_raw$Month<10,paste0("0",Month),Month)  ,"20",sep="-"), Year, sep=""))) %>%
  # 4. Add new_id column
  mutate(new_id = paste(G_NUMBER, isolation_date, rep("1",dim(Vietnam_raw)[1]), sep = "/"))



# Subset datasets ---------------------------------------------------------

Malawi_subset <- subset(Malawi, select=c("G_NUMBER","country","lineage","sublineage","isolation_date","new_id"))
Tanzania_subset <- subset(Tanzania, select=c("G_NUMBER","country","lineage","sublineage","isolation_date","new_id"))
TheGambia_subset <- subset(TheGambia, select=c("G_NUMBER","country","lineage","sublineage","isolation_date","new_id"))
Vietnam_subset <- subset(Vietnam, select=c("G_NUMBER","country","lineage","sublineage","isolation_date","new_id"))
Vietnam_subset_unique <- distinct(Vietnam_subset)



# Combine datasets --------------------------------------------------------

full_dataset <- full_join(full_join(full_join(Malawi_subset,Tanzania_subset),TheGambia_subset),Vietnam_subset_unique)



# Save full dataset -------------------------------------------------------

write.table(full_dataset, file=snakemake@output[["meta_all"]], row.names=F, col.names=T, quote=F, sep= '\t')