#Authors: Timothy D. Berry, Thea L. Whitman, edits by Dana Johnson, 10/28,2022
#Associated works: An open-source, automated, gas sampling
#peripheral for laboratory incubation experiments

#Load necessary packages
library (tidyverse)
library (zoo)
library (stats)
library (RColorBrewer)
library (broom) 
library(dplyr)
library(plyr)

# November 9 - still need to run Oct and Nov through new processing (max instead 
# of end of curve peak picking) script.


#The working directory will need to be set by the user to point to where their
#data is stored, for example:

#setwd("C:/Users/User1/Desktop/Multiplexer/Demonstration_run")
setwd("../Box/WhitmanLab/Projects/WoodBuffalo/ModellingFires_DOE_TES/data/incubations/multiplexer")

#Changes overflow options to suppress automatic conversion to scientific notation
options(scipen=99)

#Changes from v9 - Moved various adjustable parameters up front, added additional calculations/graphs

################################################################################
#Parameter setup
################################################################################

#When smoothing data, how many points wide should the sliding window be? 
#Larger values = smoother
smooth_window = 50

#When checking for outliers, how many points should be analyzed (nslice) and what
#threshold constitutes an outlier (outliercutoff) smaller values will find more 
#outliers, but increase chances of false positives. It is recommended these values
#not be changed unless you have very choppy data even after smoothing
nslice = 10
outliercutoff = 99999

#Designates the valves used for flushing manifolds (flush_positions)
#and for moving gas through system without sampling jars (idle_positions).
#Anything not listed here will be treated like a sample!
flush_positions <-  c(1,17,33,49)
idle_positions <-  c(16,32,48,64)

#Which valves are connected to samples that are your soil treatment 
#(or your first C source) and which have had PyOM (or your second C source) added?
treatments_labeled<- c(3,5,6,7,8,9,
                       10,11,13,14,15,18,
                       19,20,21,22,23,24,
                       26,27,28,30,31,34,
                       35,36,37,38,39,40,
                       42,43,44,46,47,50)

trtmts_control <- c(3,7,10,15,19,21,30,34,39,40,47,50)
trtmts_low_burn <- c(5,6,13,18,20,24,26,31,35,38,42,44)
trtmts_high_burn <- c(8,9,11,14,22,23,27,28,36,37,43,46)
standards <- c(52,53,54,55,56)

#treatments_unlabeled<- c(11, 19, 25, 27) 
blanks <- c(2,4,12,25,29,41,45)

#Fraction of system volume that comes from the sample jar. In a system with no-deadspace 
#this would be 1. In reality, this depends on jar volume and
#must be determined for each specific sample container that is used
fj <- 0.7541  # 0.3131 = jam jars; 0.7541 = pint

#Parameters specific to incubation experimental design:
#Jar volume (L), dry mass of soil (grams), pyom added (mg), and
#PyOM d13C (d13C per mille vPDB)
jarvol_L <- 0.473 #0.1183 = jam; 0.473 = pint
#mass_soil_g <- 10
#label_added_mg <- 0.4
#label_d13C <- 100000 # this is a placeholder number (80 atom %)

#Sets the times during which standards are run. Measurements within these windows
#will not be analyzed as if they were samples
#standards_begin <- 1565634800
#standards_end <- 1565638000

################################################################################
#Data preprocessing
################################################################################

#Assembles a list of log files and Picarro data files from the chosen directory and combines them
#Verify that Demonstration_log_file.csv is in the working directory chosen above before running

#Aug
relay_log = rbind(read.csv("Aug/log_16.csv"),read.csv("Aug/log_17.csv"),
                 read.csv("Aug/log_18.csv"),read.csv("Aug/log_22.csv"),
                 read.csv("Aug/log_24.csv"),read.csv("Aug/log_25.csv"),
                  read.csv("Aug/log_30.csv"))

# Sept
relay_log = rbind(read.csv("Sept/log_30_Aug.csv"),
                  read.csv("Sept/log_09.csv"),
                  read.csv("Sept/log_16.csv"),
                  read.csv("Sept/log_23.csv"))

# Oct
relay_log <- rbind(read.csv("Oct/log_02.csv"),
                read.csv("Oct/log_10.csv"),
              read.csv("Oct/log_19.csv"),
            read.csv("Oct/log_24.csv"),
          read.csv("Oct/log_26.csv"),
        read.csv("Oct/log_27.csv"))

# Nov
#relay_log <- rbind(read.csv("Oct/log_27.csv"))


picarro_files <- list.files(path = "./Aug/", pattern='.*[.]dat$', recursive = T, full.names = TRUE)
picarro_files

all.pic.data <- as.data.frame(c())
total.pic.files <- length(picarro_files)

for(i in 1:total.pic.files){
  temp <- read_table(picarro_files[i])
  all.pic.data <- rbind.data.frame(temp, all.pic.data)
}

rm(i, picarro_files)

#Keep only data columns that are needed for this analysis
columns_to_keep <- c("TIME", "EPOCH_TIME", "12CO2_dry", "Delta_Raw_iCO2","Delta_5min_iCO2")

cleaned.pic.data <- all.pic.data[,columns_to_keep]
rm(all.pic.data, temp)


# # Archived files for September
# archived_files = rbind(read.csv("../Archived files/Archived Data/Sept/10/combined_data_file_20220910.csv"),
#                        read.csv("../Archived files/Archived Data/Sept/09/combined_filtered_data_20220909.csv"),
#                        read.csv("../Archived files/Archived Data/Sept/11/combined_filtered_data_20220911.csv"),
#                        read.csv("../Archived files/Archived Data/Sept/12/combined_filtered_data_20220912.csv"),
#                        read.csv("../Archived files/Archived Data/Sept/13/combined_filtered_data_20220913.csv"),
#                        read.csv("../Archived files/Archived Data/Sept/14/combined_filtered_data_20220914.csv"),
#                        read.csv("../Archived files/Archived Data/Sept/15/combined_filtered_data_20220915.csv"))
# 
# names(archived_files)[4] <- "dry_12CO2"
# 
# columns_to_keep <- c("EPOCH_TIME", "dry_12CO2", "Delta_Raw_iCO2","Delta_5min_iCO2")
# 
# cleaned.archived.data <- subset(archived_files, dry_12CO2 != '12CO2_dry')
# cleaned.archived.data <- archived_files[,columns_to_keep]
# 
# # Combine the archived files with the main dataset:
# names(cleaned.pic.data)[3] <- "dry_12CO2"
# 
# cleaned.pic.data <- cleaned.pic.data[,columns_to_keep]
# 
# head(cleaned.pic.data)
# head(cleaned.archived.data)
# head(data.frame(cleaned.pic.data))
# 
# 
# cleaned.archived.data$EPOCH_TIME <- as.numeric(cleaned.archived.data$EPOCH_TIME)
# cleaned.archived.data$dry_12CO2 <- as.numeric(cleaned.archived.data$dry_12CO2)
# cleaned.archived.data$Delta_Raw_iCO2 <- as.numeric(cleaned.archived.data$Delta_Raw_iCO2)
# cleaned.archived.data$Delta_5min_iCO2 <- as.numeric(cleaned.archived.data$Delta_5min_iCO2)
# 
# 
# cleaned.pic.data <- bind_rows(data.frame(cleaned.pic.data), cleaned.archived.data)
# 
# rm(archived_files)
# rm(cleaned.archived.data)
 
 



#Removes unneeded columns and times from data that are from 
#before or after the relay was active (since it's not experiment data)
#Also renames the CO2 concentration column to not start with
#a number so that GG plot doesn't encounter issues later

pic.trimmed <-  cleaned.pic.data[cleaned.pic.data$EPOCH_TIME > relay_log$Epoch_time[1],]
pic.trimmed <- pic.trimmed[pic.trimmed$EPOCH_TIME < relay_log$Epoch_time[length(relay_log$Epoch_time)],]
names(pic.trimmed)[3] <- "dry_12CO2" # packaged used later doesn't like strings starting with numbers


#Coverts epoch times into integers and removes duplicate timestamps caused by
#the Picarro and relay board not measuring at exact 1 second intervals. 
#This results in a loss of some data points, but there should still be plenty 
#left to redraw the curves and do our calculations with
pic.trimmed <- dplyr::mutate(pic.trimmed, EPOCH_TIME = as.integer(EPOCH_TIME))
pic.trimmed.deduped <- distinct(pic.trimmed, EPOCH_TIME, .keep_all = TRUE)


#Renames the Epoch_time in the relay file to EPOCH_TIME,
#removes duplicate time stamps, and creates a data frame for the
#relay logs wherein the epoch time matches with a picarro time
names(relay_log)[1] <- "EPOCH_TIME"
relay.deduped <- distinct(relay_log, EPOCH_TIME, .keep_all = TRUE)
relay.matched <- filter(relay.deduped, EPOCH_TIME %in% pic.trimmed.deduped$EPOCH_TIME)

merged <-  merge.data.frame(relay.matched, pic.trimmed.deduped, by.x = "EPOCH_TIME")


### Preliminary look at data
test <- merged %>%
  # filter(Active_relay1 =='3')
  filter(EPOCH_TIME < 1660800000 & EPOCH_TIME > 1660750000)
  filter(!(EPOCH_TIME %% 2))
dim(test)

test$Active_relay1 <- as.character(test$Active_relay1)
test$ID = 'TBD'
for (i in 1:nrow(test)) {
  if (test$Active_relay1[i] %in% c(blanks)) {
    test$ID[i] = 'blank' 
  } else if (test$Active_relay1[i] %in% c(flush_positions)) {
      test$ID[i] = 'flush'
  } else if (test$Active_relay1[i] %in% c(idle_positions)) {
    test$ID[i] = 'idle'
  } else if (test$Active_relay1[i] %in% c(trtmts_control)) {
    test$ID[i] = 'control'
  } else if (test$Active_relay1[i] %in% c(trtmts_low_burn)) {
    test$ID[i] = 'low'
  } else if (test$Active_relay1[i] %in% c(trtmts_high_burn)) {
    test$ID[i] = 'high'
  } else if (test$Active_relay1[i] %in% c(standards)) {
    test$ID[i] = 'standard'
  }
}

ggplot(subset(test, !(Active_relay1 %in% c(flush_positions, idle_positions)) & 
                EPOCH_TIME < 1660800000 & ID %in% c('blank','control')), 
       aes(x=EPOCH_TIME, y=dry_12CO2, color = Active_relay1)) + geom_point() + 
  facet_wrap(~Date, scales = 'fixed')


ggplot(subset(test, EPOCH_TIME < 1661000000), 
       aes(x=EPOCH_TIME, y=dry_12CO2, color = ID)) + geom_point() 




#Cleanup to salvage some memory
rm(cleaned.pic.data, pic.trimmed, relay.deduped, relay_log, logfile_list_read, logfile_list,
   pic.trimmed.deduped, relay.matched)


#Smoothing using the ZOO package. Takes a moving average (determined by the smooth_window parameter, set above). 
# For this dataset, 300 points roughly corresponds to a 5 minute moving average;
# 50 points corresponds to a 50 second moving average.
CO2_smooth <-rollapply(data = merged$dry_12CO2, FUN = "mean", width = smooth_window, partial = TRUE, align = "right")
merged <- cbind(merged,CO2_smooth)

iCO2_smooth <-rollapply(data = merged$Delta_Raw_iCO2, FUN = "mean", width = smooth_window, partial = TRUE, align = "right")
merged <- cbind(merged,iCO2_smooth)

rm(CO2_smooth, iCO2_smooth)

#Removes standard for calibration from merged data frame and creates a second 
#"standards dataframe". 

#standards <- filter(merged, EPOCH_TIME > standards_begin & EPOCH_TIME < standards_end)

#merged <- filter(merged, EPOCH_TIME > standards_end)


################################################################################
#Data Annotation (adds status and cycle info)
################################################################################

#Remove the redundant Time column in the merged dataframe and creates a status, 
#sample number,and cycle number, and cycle_bound columns which are then assigned 
#dummy values for now and are populated in the next steps.

merged$Time <- NULL
merged$status <- 0
merged$sample_num <- 0
merged$cycle <- NA
merged$cycle_bound <- FALSE

#Determines the status of the analysis (stored in the merged dataframe) at each
#time point and assigns sample numbers to jar-specific steps. The status column is 
#dplyr::mutated to replace the dummy value with a status value assigned based on
#the combination of valves open. For the sake of this process, it is assumed that the default 
#manifold configuration is used and the first position on each manifold is used to 
#flush while the last position on each is an "idle" position in which no jar is sampled.

# Tim's version - use this (11/09/2022)
merged <- merged %>% 
  dplyr::mutate(status = replace(status, is.element(Active_relay1, flush_positions) & 
                            is.element(Active_relay2, idle_positions),"flushing lines")) %>%
  dplyr::mutate(status = replace(status, is.element(Active_relay1, flush_positions) & 
                            is.na(Active_relay2),"flushing lines")) %>%
  dplyr::mutate(status = replace(status, is.element(Active_relay1,flush_positions) & 
                            !is.element(Active_relay2,idle_positions) & !is.na(Active_relay2),"flushing jar")) %>% 
  dplyr::mutate(status = replace(status, !is.element(Active_relay1,flush_positions) & 
                            !is.element(Active_relay1,idle_positions),"measuring jar")) %>%
  dplyr::mutate(status = replace(status, is.element(Active_relay1,idle_positions),"idle")) %>%
  dplyr::mutate(status = replace(status, status != lag(status,1),"boundary")) %>%
  dplyr::mutate(status = replace(status, is.na(lag(Date,1)),"boundary")) %>%
  dplyr::mutate(sample_num = case_when(.$status == "measuring jar" ~ .$Active_relay1,
                                .$status == "flushing jar" ~ .$Active_relay2,
                                .$status == "boundary" ~ as.integer(0),
                                .$status == "flushing lines" ~ as.integer(0),
                                .$status == "idle" ~ as.integer(0))) %>%
  dplyr::mutate(cycle_bound = 
           replace(cycle_bound,status == "boundary" &
                     lag(status) == "flushing jar" &
                     lead(status) == "idle", TRUE))



# my version?
# merged2 <- merged %>% 
#   dplyr::mutate(status = replace(status, is.element(Active_relay1, flush_positions) & 
#                             is.element(Active_relay2, idle_positions),"flushing lines")) %>%
#   dplyr::mutate(status = replace(status, is.element(Active_relay1, flush_positions) & 
#                             is.na(Active_relay2),"flushing lines")) %>%
#   dplyr::mutate(status = replace(status, is.element(Active_relay1,flush_positions) & 
#                             !is.element(Active_relay2,idle_positions) & !is.na(Active_relay2),"flushing jar")) %>% 
#   dplyr::mutate(status = replace(status, !is.element(Active_relay1,flush_positions) & 
#                             !is.element(Active_relay1,idle_positions),"measuring jar")) %>%
#   dplyr::mutate(status = replace(status, is.element(Active_relay1,idle_positions),"idle")) %>%
#   dplyr::mutate(status = replace(status, status != lag(status,1),"boundary")) %>%
#   dplyr::mutate(status = replace(status, is.na(lag(Date,1)),"boundary")) %>%
#   dplyr::mutate(sample_num = case_when(.$status == "measuring jar" ~ .$Active_relay1,
#                                 .$status == "flushing jar" ~ .$Active_relay2,
#                                 .$status == "boundary" ~ as.integer(0),
#                                 .$status == "flushing lines" ~ as.integer(0),
#                                 .$status == "idle" ~ as.integer(0))) %>%
#   dplyr::mutate(cycle_bound = 
#            replace(cycle_bound,status == "boundary" &
#                      lead(status) == "flushing lines" &
#                      lag(sample_num) == "9", TRUE))


head(merged,2)
# First look at data

#!!!NOTE: This function assigns cycle numbers based on boundaries between each 
#measurement cycle. When using this function it is assumed that cycles were run 
#to completion and following the suggested analysis timing template. If this
#Is not the case, this function must be modified to accurately label cycle numbers

cycle_counter <- function(x){
  boundary_vec <- which(x$cycle_bound == TRUE)
  boundary_vec <- c(1, boundary_vec)
  for(i in 2:length(boundary_vec)){
    x[c(boundary_vec[i-1]:boundary_vec[i]),]$cycle <- i -1 
  }
  return(x)
}

merged_complete<-cycle_counter(merged)

#rm(merged)
p1 = ggplot(subset(merged_complete, !(EPOCH_TIME %% 60)), aes(x=EPOCH_TIME, y=dry_12CO2)) + 
  geom_point() +
  labs(y='CO2 conc. (ppm)') +
  labs(title = 'respiration profiles')
p1

################################################################################
#Summarize Data and correct for dilution
################################################################################

#Remove boundary and idle statuses
jars <- merged_complete %>% group_by(sample_num, cycle) %>% 
  filter(cycle != 0) %>% filter(status != "boundary") %>% filter(status != "idle")

rm(merged_complete)

#To obtain the "dilute" value the minimum CO2 concentration during the "measuring jar" 
#phase is found and points within nslice/2 are taken as the representative area. Thus
#the dilute value CENTERS on the the most dilute value. Please note that this method
#works only when you have flat or positive peaks. "Negative" peaks where samples consume
#CO2 are not currently supported

most_dilute = jars %>% filter(status == "measuring jar") %>% 
  group_by(sample_num, cycle) %>% summarize(cycle_minimum = min(CO2_smooth))

min_time <- jars %>% filter(status == "measuring jar") %>%  
  group_by(sample_num,cycle) %>%  merge(.,most_dilute) %>%
  filter(CO2_smooth == cycle_minimum) %>% 
  select(sample_num, cycle, minimum_time = EPOCH_TIME)

iCO2_dilute = jars%>% 
  # Grab the datapoints right when the jar is being measured
  filter(status == "measuring jar") %>%
  # Create a series of columns identifying the difference between neighboring values,
  # This is not really necessary for the smoothed data but is not computationally
  #intensive so there is little reason to remove it yet.
  select(sample_num,cycle,iCO2_smooth,CO2_smooth, EPOCH_TIME,status) %>%
  dplyr::mutate(NeighborDist1 = c(NA,diff(iCO2_smooth,lag=1)),
         NeighborDist2 = c(NA,NA,diff(iCO2_smooth,lag=2)),
         NeighborDist3 = c(NA,NA,NA,diff(iCO2_smooth,lag=3)),
         NeighborDist4 = c(diff(iCO2_smooth,lead=1),NA),
         NeighborDist5 = c(diff(iCO2_smooth,lead=2),NA),
         NeighborDist6 = c(diff(iCO2_smooth,lead=3),NA)) %>%
  # Ungroup everything
  group_by(EPOCH_TIME)%>%
  # Create an outlier flag column, that pings when the distance between nearby points is above the outlier cutoff
  dplyr::mutate(Flag = ifelse(max(abs(NeighborDist1),abs(NeighborDist2),abs(NeighborDist3),abs(NeighborDist4),abs(NeighborDist5),abs(NeighborDist6),na.rm=TRUE)>outliercutoff,1,0))%>%
  # Group them up again
  group_by(sample_num,cycle) %>%  merge(.,min_time) %>%  
  dplyr::mutate(difference = minimum_time - EPOCH_TIME) %>% 
  filter(abs(difference) <= nslice/2) %>%
  filter(Flag != 1) %>% group_by(sample_num,cycle) %>%
  # Take the mean value of that raw iCO2 input
  # End result is a table with sample number, cycle number, and mean iCO2
  summarize(dilute_iCO2=mean(iCO2_smooth),dilute_CO2 = mean(CO2_smooth))

#Do the same thing at the end of the flush step to get the purge values (concentration)
#of CO2 in system after flushing. This is also where the summarized sample's time 
#entry comes from - the last point in the purge step is when accumulation begins.
iCO2_purged = jars%>%
  filter(status == "flushing jar") %>%
  select(sample_num,cycle,iCO2_smooth,CO2_smooth, EPOCH_TIME,status) %>%
  dplyr::mutate(NeighborDist1 = c(NA,diff(iCO2_smooth,lag=1)),
         NeighborDist2 = c(NA,NA,diff(iCO2_smooth,lag=2)),
         NeighborDist3 = c(NA,NA,NA,diff(iCO2_smooth,lag=3)),
         NeighborDist4 = c(diff(iCO2_smooth,lead=1),NA),
         NeighborDist5 = c(diff(iCO2_smooth,lead=2),NA),
         NeighborDist6 = c(diff(iCO2_smooth,lead=3),NA)) %>%
  group_by(EPOCH_TIME)%>%
  dplyr::mutate(Flag = ifelse(max(abs(NeighborDist1),abs(NeighborDist2),abs(NeighborDist3),
                           abs(NeighborDist4),abs(NeighborDist5),abs(NeighborDist6)
                           ,na.rm=TRUE)>outliercutoff,1,0)) %>%
  group_by(sample_num,cycle)%>%
  slice((n()-nslice+1):n())   %>%
  filter(Flag != 1) %>%
  summarize(purged_iCO2=mean(iCO2_smooth),
            purged_CO2 = mean(CO2_smooth), time = last(EPOCH_TIME))


#Original: "Peak selection takes the last point of the "measuring jar" status and points preceding
#this point by 10 seconds to be representative of the peak value: this tends 
#to represent the most stable part of the "peak" and makes the measurement
#much less likely to be erroneous when the previous sample was very concentrated"

# My notes: "But I'd like to use the max values as the measurement point because
## we were venting jars to atm during sampling. 

meas_end <- jars %>% group_by(sample_num, cycle) %>% 
  filter(status == "measuring jar") %>% summarize(last_point = last(EPOCH_TIME))

# Modified:
iCO2_peak <- jars %>%
  filter(status == "measuring jar") %>%
  select(sample_num,cycle,iCO2_smooth,CO2_smooth, EPOCH_TIME,status) %>%
  # make list of differences between given measurements and 3 preceeding and following measurements
  dplyr::mutate(NeighborDist1 = c(NA,diff(iCO2_smooth,lag=1)),
         NeighborDist2 = c(NA,NA,diff(iCO2_smooth,lag=2)),
         NeighborDist3 = c(NA,NA,NA,diff(iCO2_smooth,lag=3)),
         NeighborDist4 = c(diff(iCO2_smooth,lead=1),NA),
         NeighborDist5 = c(diff(iCO2_smooth,lead=2),NA),
         NeighborDist6 = c(diff(iCO2_smooth,lead=3),NA)) %>%
  # grouping by "EPOCH_TIME" undoes the previous grouping (since all EPOCH_TIMES 
  # here are unique)
  group_by(EPOCH_TIME) %>%
  # If any of the neighboring values are above our outlier, flag them (to be
  # removed later)
  dplyr::mutate(Flag = ifelse(max(abs(NeighborDist1),abs(NeighborDist2),abs(NeighborDist3),
                           abs(NeighborDist4),abs(NeighborDist5),abs(NeighborDist6),
                           na.rm=TRUE)>outliercutoff,1,0)) %>% 
  # Combine the EPOCH_TIME for the end of each measurement listed as "last_point"
  group_by(sample_num, cycle) %>% merge(.,meas_end) %>%
  # this creates a list "difference" of the time (s) before the end of the 
  # measuring period. E.g. at difference = 5, there are 5 s remaining in the 
  # measuring period. 
  dplyr::mutate(difference = last_point - EPOCH_TIME) %>%
  # This is set up to select for the last 10 (nslice = 10) points before the 
  # end of the measuring period
  #       filter(difference <= nslice) %>%
  # But we can switch it to filter for the beginning of the sampling or I could
  # just take the max value of the measurment period. 
  filter(Flag != 1) %>% group_by(sample_num,cycle) %>%
  # max value from sampling period 
  summarize(peak_iCO2 = max(iCO2_smooth), 
            peak_CO2 = max(CO2_smooth))





#Merge these 3 sets of values together and remove unnecessary intermediate files  
summary <- merge(iCO2_dilute, iCO2_peak) %>%
  merge(.,iCO2_purged) %>%
  dply::arrange(sample_num,cycle)

summary_long1 <- summary %>% rename("dilute" = dilute_CO2, "peak" = peak_CO2, "purged" = purged_CO2) %>%
  gather(key = "landmark", value = "CO2", dilute, peak, purged) %>% 
  select(-dilute_iCO2, -peak_iCO2, -purged_iCO2)

summary_long2 <- summary %>% rename("dilute" = dilute_iCO2, "peak" = peak_iCO2, "purged" = purged_iCO2) %>%
  gather(key = "landmark", value = "iCO2", dilute, peak, purged) %>% group_by(sample_num, landmark) %>%
  select(-dilute_CO2, -peak_CO2, -purged_CO2)

summary_long <- merge(summary_long1, summary_long2)






# quick look at data
test2 <- summary_long %>%
  filter(time < 1662000000)
dim(test2)

ggplot(subset(test2, sample_num %in% c(blanks) & time < 1662000000),
       aes(x=time, y=CO2, color=sample_num)) + 
  geom_point()

test3 <- final_summary %>%
  filter(time < 1662000000)

ggplot(subset(test3, sample_num==10 & time < 1662000000), 
       aes(x=time, y=sample_CO2, color = sample_num)) + 
  geom_point() 




ggplot(subset(test, EPOCH_TIME < 1662000000), 
       aes(x=EPOCH_TIME, y=dry_12CO2, color = ID)) + geom_point() 




#Plot figures to see the CO2 concentration and isotope profiles across cycles
#This is an easy way to visualize trends in respiration of samples and to assess
#whether your relay command file is adequate. Ideally the dilute and purge values
#should remain consistent throughout the course of the experiment. Changes in
#these values over the course of the experiment likely mean insufficient time is
#being given for the "flushing jars" step. The script can account for this deviations
#but to maximize data quality the flush time should be increased.
respiration_profiles <- ggplot(data = subset(summary_long, sample_num %in% c(50))) + 
  geom_col(aes(x = cycle, y = CO2, fill = landmark), position = "dodge") + 
  facet_wrap(~sample_num, scale = "free") + 
  labs(y='CO2 conc. (ppm)') +
  labs(title = 'respiration profiles')
respiration_profiles


isotope_profiles <- ggplot(data = subset(summary_long, sample_num %in% c(26:34) &
                                           time < 1666000000)) + 
  geom_point(aes(x=(time-1660000000)/100000, y=iCO2, color=landmark)) +
  # geom_col(aes(x = cycle, y = iCO2, fill = landmark), position = "dodge") + 
  facet_wrap(~sample_num, scale = "fixed") + 
  labs(y="delta 13C", title = 'isotope profiles')

isotope_profiles


#Remove unnecessary intermediate files to free up memory
rm(summary_long1, summary_long2)
rm(meas_end, min_time, most_dilute)
rm(iCO2_dilute, iCO2_peak, iCO2_purged, jars)

##########################################
#Calculations
##########################################

#Dilution parameters are determined by injecting known volumes of pure CO2 into
#purged vessels to determine the fraction of the contribution of the sample jar's gas (fj) to observed gas


#The total CO2 measured at the peak represents a combination of what's in the system
#before you start measuring and what's in the jar. The total volume of the system is jar
#+ gas lines, instrument, manifolds, etc. Thus, the CO2 in an equilibrated system will 
#be whatever  CO2 came from the system (dilute_CO2) multiplied by its fraction (1-fj),
#plus the CO2 that was in the jar multiplied by its fraction (fj).
#That results in the following equation:
corrected <- summary %>% 
  dplyr::mutate(peak_CO2_cor = (peak_CO2 - dilute_CO2 + (dilute_CO2 * fj))/fj)

#Next, we need to correct the iCO2 at the peak measurement to account for the fact
#that  our measurement includes mostly  CO2 from our sample, but also some
#from whatever gas remained in the system before initiating the measurement.
#Fortunately, we have those values. First, we calculate what fraction of all the CO2 
#molecules come from our actual sample (pj), vs. what was initially in the system (1-pj). 
#We can get to this by taking the amount of each source of CO2, multiplied by its 
#relative volume. Once we have that value, pj, we can use the same approach as we did
#for the total CO2, above.
i_corrected <- corrected %>% 
  dplyr::mutate(pj = (peak_CO2_cor * fj)/((peak_CO2_cor * fj + (1-fj) * dilute_CO2))) %>%
  dplyr::mutate(final_peak_iCO2 = (peak_iCO2 - dilute_iCO2 + (dilute_iCO2 * pj))/pj)

#**# Then, we just remove the columns we do not want to keep
final <- i_corrected %>% 
  select(-pj, -dilute_iCO2, -peak_iCO2, -peak_CO2, -peak_iCO2, -peak_iCO2)


#Now, we have the estimated values for the volume of CO2 in each jar at measurement (peak_CO2_cor)
#and the isotopic signature of that CO2 (final_peak_iCO2).
#Plus, we have all our baseline values of what the whole system was at when
#it was last purged: purged_CO2 and purged_iCO2_cor.Note, purged_CO2 does not need to
#be corrected, because the whole system was open and mixing at the time of its 
#measurement, so it should be a good estimate for what is in the jars when their
#valves close. The only correction applied to the purged_iCO2_cor was adjusting for
#any machine drift in isotopic measurement.

#To determine how much CO2 accumulated in the jar between the time it was last measured
#(the previous cycle's purge numbers), and each current measurement
#For total CO2, the calculations is simple subtraction.
#For the iCO2, we need to take into account the relative fractions in the jar when it last shut and at measurement.
final <- final %>% group_by(sample_num) %>% 
  dplyr::mutate(sample_CO2 = peak_CO2_cor - lag(purged_CO2)) %>%
  dplyr::mutate(sample_iCO2 = (final_peak_iCO2 - lag(purged_iCO2)*(lag(purged_CO2)/(lag(purged_CO2) + sample_CO2)))
         /(sample_CO2/(sample_CO2 + lag(purged_CO2))))

#Now we have all the values for the total CO2 emitted over each timepoint
#and its isotopic signature, for each jar and each cycle.
final_summary <- final %>% select(sample_num, cycle, sample_CO2, sample_iCO2, time)

head(final_summary)

p = ggplot(final_summary, aes(x=sample_num, y = time, color = sample_num)) + 
  geom_point()
p

#write.csv(final_summary, "../Compiled-Aug-CO2-and-iCO2-data.csv", row.names = FALSE)
#write.csv(final_summary, "../Compiled-Sept-CO2-and-iCO2-data.csv", row.names = FALSE)
#write.csv(final_summary, "../Compiled-Oct-CO2-and-iCO2-data.csv", row.names = FALSE)
#write.csv(final_summary, "../Compiled-Nov-CO2-and-iCO2-data.csv", row.names = FALSE)














##########################################
# Import processed respiration data:
Aug <- read.csv("Compiled-Aug-CO2-and-iCO2-data.csv")
Sept <- read.csv("Compiled-Sept-CO2-and-iCO2-data.csv")
Oct <- read.csv("Compiled-Oct-CO2-and-iCO2-data.csv")
Nov <- read.csv("Compiled-Nov-CO2-and-iCO2-data.csv")

Sept <- Sept %>% dplyr::mutate(cycle = cycle + max(Aug$cycle))
Oct <- Oct %>% dplyr::mutate(cycle = cycle + max(Sept$cycle))
Nov <- Nov %>% dplyr::mutate(cycle = cycle + max(Oct$cycle))

## Working to account for the blank readings within each set of measurements:
final_summary <- rbind(Aug, Sept, Oct, Nov)
head(final_summary)


#Add a treatment column to the finalized data summary
final_summary$treatment <- NA
final_summary <- final_summary %>% dplyr::mutate(treatment = case_when(
  is.element(sample_num, blanks) ~ "blank",
  is.element(sample_num, trtmts_control) ~ "control",
  is.element(sample_num, trtmts_low_burn) ~ "30 s burn",
  is.element(sample_num, trtmts_high_burn) ~ "120 s burn",
  is.element(sample_num, standards)~ "standards"))


# Import pre-burn mass data
mass_data <- read.csv("../../burn-simulations/wet-soil-mass-data.csv")
colnames(mass_data)[1] <- 'sample_num'
head(mass_data,2)

# Calculate post-burn wet and dry soil mass (g)
mass_data <- mass_data %>%
  dplyr::mutate(post.burn.wet.soil.mass.g = post.burn.wet.soil.and.foil.mass.g-foil.mass.g,
         estimated.post.burn.dry.soil.mass.g = post.burn.wet.soil.mass.g/(1+estimated.post.burn.moisture.fraction))


#Covert CO2 values from ppm to a mass basis. We make some simple assumptions that the 
#sample is at standard temperature and pressure

# Merge final_summary with core mass data
final_summary <- merge(final_summary, subset(mass_data, 
                                             select = c(sample_num,
                                                        jar.set,burn.trtmt,
                                                        O.hor.thickness.cm,
                                                        core.thickness.cm,
                                                        post.burn.wet.soil.mass.g,
                                                        estimated.post.burn.dry.soil.mass.g,
                                                        site.id,
                                                        epoch.wet.up.time)))

# I want to correct for the blank CO2 for a given set of jars and for each separate cycle. 
final_summary <- final_summary %>%
  group_by(jar.set,cycle) %>% 
  dplyr::mutate(blank_CO2 = mean(sample_CO2[treatment =="blank"])) %>%
  dplyr::mutate(sample_CO2_blank_corrected = sample_CO2 - blank_CO2)

# Calculate time of respiration (i.e. the time between measurement for each jar)
# and remove first measurement of each jar.
final_summary <- final_summary %>%
  group_by(sample_num) %>%
  dplyr::arrange(time, .by_group = TRUE) %>%
  dplyr::mutate(sampling_duration_s = time - lag(time)) %>%
  dplyr::mutate(first_cycle = min(cycle)) %>%
  subset(cycle != first_cycle)

final_summary <- final_summary %>%
  subset(!(treatment %in% c('blank','standards')))

#### WHAT IS THIS CSV FILE???
concat_final <- read.csv("sample-num-and-cycle-v2.csv")

# Pull out cycles to keep
concat_final <- concat_final %>%
  subset(select = c(sample_num, cycle, status)) %>%
  merge(final_summary) %>%
  subset(status == 'keep')

# Sanity check.
p = ggplot(subset(concat_final,sample_num %in% c(3,18,19,20,21,22,23,24,25,26,27,28) & time > 1662000000 ), 
           aes(x=time, y = sample_CO2_blank_corrected, color = sample_num)) + 
  geom_point() + 
  facet_wrap(~sample_num)
p

### ISSUES ####
# 1. What's going on with samples c(3,5,6,7,8,9,10,11,13,15...)


# Calculate g CO2-C respired, g CO2-C per g dry soil, and g CO2-C per g dry soil per unit time
#   To convert ppm --> divide by 1000?
#   CO2 = 44.01 g/mol
mineralization <- concat_final %>%
  dplyr::mutate(CO2C_g_respired = sample_CO2_blank_corrected/1000000/22.4*jarvol_L*12.01) %>% # volume of an ideal gas is 22.41 L/mol at STP
  filter(!is.na(sample_CO2)) %>%
  dplyr::mutate(CO2C_g_per_g_dry_soil = CO2C_g_respired / estimated.post.burn.dry.soil.mass.g,
         CO2C_g_per_s = CO2C_g_respired / sampling_duration_s) %>%
  dplyr::mutate(CO2C_g_per_g_dry_soil_per_s = CO2C_g_per_g_dry_soil/sampling_duration_s,
         CO2C_g_per_g_dry_soil_per_hr = CO2C_g_per_g_dry_soil/(sampling_duration_s/(60*60)),
         CO2C_g_per_g_dry_soil_per_day = CO2C_g_per_g_dry_soil_per_hr*24)

# I want to fill in gaps in data to make the cumulative C respired smoother

# Calculate cumulative gCO2-C respired per gram dry soil 
mineralization <- mineralization[order(mineralization$cycle),] %>%
  group_by(sample_num) %>%
  dplyr::arrange(time, .by_group = TRUE) %>%
  dplyr::mutate(cycle.min = min(cycle)) %>%
  subset(cycle != cycle.min) %>%
  dplyr::mutate(cumulative_g_CO2C_respired = cumsum(CO2C_g_respired))


# Sanity check.
ggplot(subset(mineralization))+
# ggplot(subset(mineralization, site.id==1 & time > 1664000000 & time < 1666500000 & sampling_duration_s<40500), 
#        aes(x=sampling_duration_s,y = CO2C_g)) + 
  geom_point(aes(x=time,y = cumulative_g_CO2C_respired, color=site.id), size=3) +
  # stat_smooth(method = "lm",formula = y ~ x,geom = "smooth")+
  facet_wrap(~burn.trtmt)



## 1. Create .csv with sample ID, dominate vegetation, initial wet mass, initial dry mass

## 2. Import percent total C data from FireMod2019 called "approximate-soil-C.csv"

## 3. Match everything up and calculate g CO2-C per estimated g initial C. 


# Import sample list to get vegetation and wet soil mass
sample_list <- read.csv("../../../methods/sample-list-for-incubations.csv")
head(sample_list)
colnames(sample_list)[1] <- 'core.id'
sample_list <- subset(sample_list, select = c(core.id, incubation.period.days, vegetation,Burn.set, wet.mass.g))

# Import data for sample id and corresponding relay numbers on multiplexer
relay_sample_id <- read.csv("../10wk-incubation-sample-order-for-picarro.csv")
head(relay_sample_id)
colnames(relay_sample_id)[1] <- 'sample_num'

# Combine mineralization dataframe with sample list and relay number
mineralization <- mineralization %>%
  merge(relay_sample_id, by = 'sample_num') %>%
  merge(sample_list, by = 'core.id') %>%
  dplyr::mutate(days.since.wet.up = (time-epoch.wet.up.time)/3600/24,
         whole.days.since.wet.up = plyr::round_any(days.since.wet.up,0.5))



# Create vegetation labels for later figures
VegLabels = c('Picea spp.','Pinus banksiana')
names(VegLabels) = c('Picea_spp.','Pinus_banksiana')

# Convert sample number and jar set to characters
mineralization$sample_num <- as.character(mineralization$sample_num)
mineralization$jar.set <- as.character(mineralization$jar.set)



### Remove outlying data generated by shutting off and re-starting multiplexer
#     during incubations. 
outlier.time = c(1666220000:1666260000, 1664750000:1664815000,
                 1666825000:1666875000, 1662715000:1662810000, 
                 1664750000:1664900000, 1663350000:1663405000,
                 1666640000:1666700000)

mineralization.cleaned <- mineralization %>%
  filter(!(time %in% outlier.time)) %>%
  filter(CO2C_g_per_g_dry_soil_per_day > 0)

cleaned.jar1 <- subset(mineralization, jar.set ==1) %>%
  filter(!(time %in% c(1665500000:1666400000,
                       1664750000:1664850000,
                       1664400000:1664430000,
                       1663350000:1663400000,
                       1662700000:1662740000,
                       1662000000:1662075000,
                       1661900000:1661940000)))
  
cleaned.jar2 <- subset(mineralization, jar.set ==2)%>%
  filter(!(time %in% c(1666000000:1666250000,
                       1664750000:1664850000)))

cleaned.jar3 <- subset(mineralization, jar.set ==3) %>%
  filter(!(time %in% c(1666600000:1668000000,
                       1666000000:1666200000,
                       1664500000:1665000000,
                       1661355000:1661400000,
                       1663350000:1663400000,
                       1662700000:1662755000,
                       1661900000:1661920000,
                       1661175000:1661250000,
                       1666000000:1666260000))) 


cleaned.jar4 <- subset(mineralization, jar.set ==4) %>%
  filter(!(time %in% c(1666850000:1666870000,
                       1666000000:1666250000,
                       1664700000:1664820000,
                       1663350000:1663440000,
                       1662750000:1662830000,
                       1666850000:1666950000)))


cleaned.jar5 <- subset(mineralization, jar.set ==5) %>%
  filter(!(time %in% c(1666850000:1666900000,
                       1666600000:1666660000,
                       1666000000:1666260000,
                       1664750000:1664860000,
                       1663350000:1663450000, 
                       1662750000:1662810000,
                       1661850000:1662000000)))

cleaned.jar6 <- subset(mineralization, jar.set ==6) %>%
  filter(!(time %in% c(1666830000:1666870000,
                       1666600000:1666700000,
                       1666000000:1666260000,
                       1664700000:1664820000)))

cleaned.combo <- rbind(cleaned.jar1, cleaned.jar2, 
                       cleaned.jar3, cleaned.jar4,
                       cleaned.jar5, cleaned.jar6)

# Sanity check.
ggplot(subset(cleaned.combo, burn.trtmt == 'control'))+
  geom_point(aes(x=whole.days.since.wet.up,y = CO2C_g_per_g_dry_soil_per_day, color=site.id), size=1) +
  # stat_smooth(method = "lm",formula = y ~ x,geom = "smooth")+
  facet_wrap(~burn.trtmt)






#### Normalizing C respired by a total initial C values:
# Import C and N data
df.CN <- read.csv('../../soil-properties/total-C-and-N.csv')
colnames(df.CN)[1] <- 'sample.id'

# Import metadata
df.meta <- read.csv('../../metadata.csv') 

colnames(df.meta)[1] = 'sample.id'

# Change site name and pull out only relevant columns
df.CN <- df.CN %>%
  dplyr::mutate(site.id = site) %>%
  merge(subset(df.meta, select = c(sample.id, horizon, vegetation,burn.trtmt.duration.seconds)))

# Create an O horizon-only dataframe and relabel C % column
df.O <- subset(df.CN, horizon == 'O') %>%
  dplyr::mutate(O.hor.C.percent = C.percent) %>%
  subset(select = c(site.id, burn.trtmt.duration.seconds, vegetation, O.hor.C.percent))

# Create an M horizon-only dataframe and relabel C% column
df.M <- subset(df.CN, horizon == 'M') %>%
  dplyr::mutate(M.hor.C.percent = C.percent) %>%
  subset(select = c(site.id, burn.trtmt.duration.seconds, M.hor.C.percent))

# Recombine O and M dataframes
df.horizons <- merge(df.O, df.M, by = c('site.id','burn.trtmt.duration.seconds'), all = TRUE)

# Create a burn.trtmt column that will match column in mineralization dataframe
for (i in 1:nrow(df.horizons)) {
  if (df.horizons$burn.trtmt.duration.seconds[i] == 0) {
    df.horizons$burn.trtmt[i] = 'control'
  } else if (df.horizons$burn.trtmt.duration.seconds[i] == 30) {
    df.horizons$burn.trtmt[i] = 'low severity'
  } else if (df.horizons$burn.trtmt.duration.seconds[i] == 120) {
    df.horizons$burn.trtmt[i] = 'high severity'
  }
}

# Merge the horizons dataframes and the cleaned up mineralization dataframe ("Cleaned.combo")
df.merge <- merge(df.horizons, cleaned.combo, by = c('site.id', 'burn.trtmt', 'vegetation'), all=TRUE) %>%
         # Calculate O and M horizon mass
  dplyr::mutate(O.hor.mass.g = O.hor.thickness.cm/core.thickness.cm*estimated.post.burn.dry.soil.mass.g,
         M.hor.mass.g = estimated.post.burn.dry.soil.mass.g-O.hor.mass.g,
         # Calculate O and M horizon mass of C
         O.hor.C.mass.g = O.hor.mass.g*(O.hor.C.percent/100),
         M.hor.C.mass.g = M.hor.mass.g*(M.hor.C.percent/100)) %>%
         # calculate initial total C
  group_by(core.id, cycle) %>%
  dplyr::mutate(initial.total.C.g = sum(O.hor.C.mass.g, M.hor.C.mass.g, na.rm=TRUE),
         # Calculate cumulative C respired per initial grams of C
         cumulative.g.CO2C.per.initial.g.C = cumulative_g_CO2C_respired/initial.total.C.g, 
         # And finally, calculate grams C respired per initial gram C per day
         g.CO2C.per.initial.g.C.per.day = CO2C_g_respired/initial.total.C.g*(60*60*24/sampling_duration_s)) %>%
  dplyr::mutate(g.CO2C.per.day = CO2C_g_respired/(60*60*24/sampling_duration_s))


df.merge$site.id <- as.character(df.merge$site.id)


# Sanity check
ggplot(subset(df.merge)) +
  geom_point(aes(x=whole.days.since.wet.up, y = g.CO2C.per.initial.g.C.per.day, color = site.id), size=1) + 
  geom_line(aes(x=whole.days.since.wet.up, y = g.CO2C.per.initial.g.C.per.day, color = site.id))+
  facet_wrap(~burn.trtmt, scales='fixed')



### Calculate fractional C remaining ###
remaining <- df.merge[order(df.merge$cycle),]  %>%
  group_by(core.id) %>%
  dplyr::mutate(cum_CO2C_g = cumsum(CO2C_g_respired),
         fractional.C.remaining = (initial.total.C.g - cum_CO2C_g)/initial.total.C.g)

# Sanity check...
ggplot(subset(remaining)) +
  geom_point(aes(x=time, y = fractional.C.remaining, color = site.id, shape=burn.trtmt), size=2) 


# Testing out filling in missing data to calculate fractional C remaining correctly.
df.test <-  remaining[order(remaining$cycle),] %>%
  group_by(core.id) %>%
  dplyr::mutate(days.between.measurements = whole.days.since.wet.up - lag(whole.days.since.wet.up)) %>% 
  subset(select = c(core.id, burn.trtmt, site.id, vegetation,cycle, sampling_duration_s,
                    CO2C_g_respired, whole.days.since.wet.up,time, days.between.measurements,
                    initial.total.C.g, g.CO2C.per.initial.g.C.per.day, g.CO2C.per.day))

# Create a column that is days 1-70 for every core.
core.ids <- subset(df.test, select = c(core.id, burn.trtmt, site.id, vegetation)) %>% unique()

time = data.frame(whole.days.since.wet.up=0, core.id = 'x', burn.trtmt = 'x', site.id = 0, vegetation='x')

for (i in 1:nrow(core.ids)) {
    x = data.frame(core.id = core.ids$core.id[i], 
                   whole.days.since.wet.up = seq(1,70,0.5), 
                   burn.trtmt = core.ids$burn.trtmt[i],
                   site.id = core.ids$site.id[i],
                   vegetation = core.ids$vegetation[i])
    time = rbind(time, x)
}

# Merge this with respired C dataframe:
df.test <- merge(df.test, time, by = c('core.id', 'whole.days.since.wet.up', 'burn.trtmt','site.id', 'vegetation'), all=TRUE)

# Create a new column for C respired
df.test$new.CO2C_g_respired = df.test$CO2C_g_respired

# Average CO2 values taken on the same day
df.test <- df.test[order(df.test$whole.days.since.wet.up),]  %>%
  group_by(whole.days.since.wet.up, core.id) %>%
  dplyr::mutate(mean.CO2C_g_respired = mean(CO2C_g_respired))

df.test <- df.test[order(df.test$whole.days.since.wet.up),] 
df.test <- df.test[order(df.test$core.id),]
df.test$data = 'missing'

# For days missing measurements, fill in C respired using the leading value
#    for C respired. 
for (i in 2:nrow(df.test)) {
  initial.total.C = df.test$initial.total.C.g[i-1]
  cycle = df.test$cycle[i]
  j=i
  k=i
  t=i
  leading = df.test$mean.CO2C_g_respired[i-1]
  previous.time = df.test$time[i-1]
  lagging = df.test$mean.CO2C_g_respired[i+1]
  
  if (is.na(df.test$cycle[i])==FALSE) {
    df.test$new.CO2C_g_respired[i] = df.test$CO2C_g_respired[i]
    df.test$data[i] = 'present'
  } else if (is.na(df.test$cycle[i])==TRUE & df.test$whole.days.since.wet.up[i]>=5) {
    # I'm only going to fill in data after the first 5 days of measuring because 
    #    the first 2-3 days see larger changes in respiration rates, so using avg.
    #    at this point to fill in data may be obscuring respiration trends
    df.test$data[i] = 'missing'
    while (j < 5795 & is.na(lagging)==TRUE) {
        j = j + 1
        lagging = df.test$mean.CO2C_g_respired[j]
        print(j)
    } 
    
    while (k>0 & is.na(leading)==TRUE) {
       k = k-1 
       leading = df.test$mean.CO2C_g_respired[k]
    }
    
    while (t>0 & is.na(previous.time)==TRUE) {
        t=t-1
        previous.time = df.test$time[t]
    }
    
    df.test$new.CO2C_g_respired[i] = (leading+lagging)/2
    df.test$initial.total.C.g[i] = initial.total.C
    df.test$cycle[i] = cycle + 1
    leading = df.test$mean.CO2C_g_respired[i]
    df.test$time[i] = previous.time+43200
  } else if (is.na(df.test$cycle[i])==TRUE) { 
    df.test$new.CO2C_g_respired[i] = df.test$new.CO2C_g_respired[i]
    df.test$initial.total.C.g[i] = df.test$initial.total.C.g[i]
    df.test$cycle[i] = df.test$cycle[i]
    df.test$time[i] = df.test$time[i]
    }
}


# Now arrange the dataset by cycle and calculate cumulative C respired and 
#    fractional C remaining. Also remove empty rows created earlier.
df.test <- df.test[order(df.test$time),]  %>%
  group_by(core.id) %>%
  dplyr::mutate(new.cum_CO2C_g = cumsum(new.CO2C_g_respired),
         new.fractional.C.remaining = (initial.total.C.g - new.cum_CO2C_g)/initial.total.C.g) %>%
  group_by(core.id, whole.days.since.wet.up) %>%
  dplyr::mutate(frac.C = mean(new.fractional.C.remaining)) %>%
  subset(site.id != 0)

# # Remove NAs from the first two days of the measuring period.
# df.test <- subset(df.test, is.na(new.fractional.C.remaining)==FALSE)


# Sanity check:
ggplot(subset(df.test, data != 'missing')) +
  geom_point(aes(x=whole.days.since.wet.up,
                 y=new.fractional.C.remaining,
                 color=burn.trtmt))+
  facet_wrap(~site.id)


ggplot(df.test) +
  geom_point(aes(x=whole.days.since.wet.up,
                 y=new.fractional.C.remaining,
                 color=data))+
  facet_wrap(~site.id)


### Save final dataframe:
# write.csv(df.test, 'processed-respiration-data.csv', row.names = FALSE)




# Extra... -----

#geom_smooth... can be enabled to plot the mean value for each treatment
#p = ggplot(data = mineralization) +
p = ggplot(data = subset(cleaned.combo
                         # jar.set %in% c(1,2,3,4,5,6)
                         # & treatment == 'control'
                         # & time > 1664300000
                         # & time < 1662000000
                         #time<1661253504  & time>1661150000
                         ), 
           aes(x = whole.days.since.wet.up, y= CO2C_g_per_g_dry_soil_per_day*1000, 
               color = treatment, shape = treatment), size=6) +
  geom_point()+
  #geom_line(aes(group = sample_num)) + 
  #geom_smooth(aes(x = wet_time, y = cumul_soilC_mg, color = treatment)) +
  #scale_color_manual(values=c("grey","black","firebrick1"), name = "Treatment")+
  #scale_shape_manual(values=c(19, 17, 15), name = "Treatment") + 
  facet_wrap(~vegetation, scales = 'free_y',labeller = labeller(vegetation = VegLabels)) +
  labs(x = "Incubation length (days)",
       y = "Soil Carbon Mineralized \n(mg CO2-C g initial dry soil)", 
       color = 'Burn treatment',
       shape = 'Burn treatment') +
  #labs(color = 'Burn duration',
  #     shape = 'Burn duration')+
  theme_bw() + 
  theme(legend.title = element_blank(),legend.text=element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title = element_text(size=14)) + 
  scale_color_manual(values = c('maroon','grey80', 'black'),
                     limits = c('120 s burn','30 s burn', 'control'),
                     labels = c('120 s burn','30 s burn', 'Unburned control')) +
  scale_shape_manual(values = c('circle','triangle', 'square'),
                     limits = c('120 s burn','30 s burn', 'control'),
                     labels = c('120 s burn','30 s burn', 'Unburned control')) +
  theme(axis.text.y = element_text(size = 14),
        legend.title = element_text(size=14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14, margin = margin(t=10, b=10)),
        axis.title.y = element_text(size = 14, margin = margin(r=10, l=10)),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.position = 'right',
        # legend.background = element_rect(fill = "darkgray"),
        legend.box = "vertical") #+ ylim(-0.1,0.1)

p




# Calculate average of lagging and leading readings.
clean.up <- mineralization %>%
  group_by(sample_num) %>%
  dplyr::mutate(lagging.reading = lag(sample_CO2_blank_corrected),
         leading.reading = lead(sample_CO2_blank_corrected),
         first.reading = min(cycle),
         last.reading = max(cycle)) %>%
  filter(cycle != first.reading, 
         cycle != last.reading, 
         is.na(lagging.reading) == FALSE, 
         is.na(leading.reading) == FALSE) %>%
  dplyr::mutate(avg.reading = (lagging.reading+leading.reading)/2) %>%
  filter(sample_CO2_blank_corrected < abs(avg.reading+(avg.reading*0.1)),
         sample_CO2_blank_corrected > avg.reading-avg.reading*0.1)

# Set range around average that read needs to fit into
# Remove readings outside of this range






p = ggplot(data = subset(normalization, ((time-epoch.wet.up.time)/3600/24) < 10), 
           aes(x = (time-epoch.wet.up.time)/3600/24, 
               y= CO2C_g_per_g_initial_C_per_day, 
               color = treatment, 
               shape = treatment)) +
  geom_point()+
  #geom_smooth(aes(x = wet_time, y = cumul_soilC_mg, color = treatment)) +
  #scale_color_manual(values=c("grey","black","firebrick1"), name = "Treatment")+
  #scale_shape_manual(values=c(19, 17, 15), name = "Treatment") + 
  facet_wrap(jar.set~vegetation, scales = "free_y",
             labeller = labeller(vegetation = VegLabels)) +
  xlab("Incubation length (days)") + 
  ylab("Soil Carbon Mineralized \n(g CO2-C g initial C)") +
  #xlim(0,15) + 
  theme_bw() + 
  labs(y = expression(Soil~respiration~rate~(mg~C-CO[2]~g^{"-1"}~day^{"-1"})),
       #y = expression(Soil~respiration~rate~(mg~C-CO[2]~g^{"-1"}~day^{"-1"})),
       color = 'Burn duration',
       shape = 'Burn duration')+
  theme(legend.title = element_blank(),legend.text=element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title = element_text(size=14)) + 
  scale_color_manual(values = c('maroon','grey80', 'black'),
                     limits = c('120 s burn','30 s burn', 'control'),
                     labels = c('120 s burn','30 s burn', 'control')) +
  scale_shape_manual(values = c('circle','triangle', 'square'),
                     limits = c('120 s burn','30 s burn', 'control'),
                     labels = c('120 s burn','30 s burn', 'control')) +
  theme(axis.text.y = element_text(size = 14),
        legend.title = element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t=10, b=10)),
        axis.title.y = element_text(size = 14, margin = margin(r=10, l=10)),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.position = 'right',
        # legend.background = element_rect(fill = "darkgray"),
        legend.box = "vertical")

p




p = ggplot(data = subset(cumulative, days.since.wet.up < 69), 
           aes(x = whole.days.since.wet.up, 
               y= cum_CO2C_mg_per_g_initial_C, 
               color = treatment, 
               shape = treatment),
           size=6) +
  geom_point()+
  #geom_line(aes(group = sample_num)) + 
  #geom_smooth(aes(x = wet_time, y = cumul_soilC_mg, color = treatment)) +
  #scale_color_manual(values=c("grey","black","firebrick1"), name = "Treatment")+
  #scale_shape_manual(values=c(19, 17, 15), name = "Treatment") + 
  facet_wrap(jar.set~vegetation, scales = "free",
             labeller = labeller(vegetation = VegLabels)) +
  xlab("Incubation length (days)") + 
  #xlim(0,15) + 
  labs(y = expression(Cumulative~C~respired~(mg~C-CO[2]~g^{"-1"})),
       shape = 'Burn duration',
       color = 'Burn duration')+
  theme_bw() + 
  theme(legend.title = element_blank(),legend.text=element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title = element_text(size=14)) + 
  scale_color_manual(values = c('maroon','grey80', 'black'),
                     limits = c('120 s burn','30 s burn', 'control'),
                     labels = c('120 s','30 s', '0 s')) +
  scale_shape_manual(values = c('circle','triangle', 'square'),
                     limits = c('120 s burn','30 s burn', 'control'),
                     labels = c('120 s','30 s', '0 s')) +
  theme(axis.text.y = element_text(size = 14),
        legend.title = element_text(size=14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14, margin = margin(t=10, b=10)),
        axis.title.y = element_text(size = 14, margin = margin(r=10, l=10)),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.position = 'right',
        # legend.background = element_rect(fill = "darkgray"),
        legend.box = "vertical")

p







p = ggplot(data = subset(remaining, days.since.wet.up < 69), 
           aes(x = whole.days.since.wet.up, 
               y= fractional.C.remaining, 
               color = treatment, 
               shape = treatment),
           size=6) +
  geom_point()+
  #geom_line(aes(group = sample_num)) + 
  #geom_smooth(aes(x = wet_time, y = cumul_soilC_mg, color = treatment)) +
  #scale_color_manual(values=c("grey","black","firebrick1"), name = "Treatment")+
  #scale_shape_manual(values=c(19, 17, 15), name = "Treatment") + 
  facet_wrap(jar.set~vegetation, scales = "free",
             labeller = labeller(vegetation = VegLabels)) +
  xlab("Incubation length (days)") + 
  #xlim(0,15) + 
  labs(y = 'Fractional C remaining',
       shape = 'Burn duration',
       color = 'Burn duration')+
  theme_bw() + 
  theme(legend.title = element_blank(),legend.text=element_text(size=14)) + 
  theme(axis.text=element_text(size=12),axis.title = element_text(size=14)) + 
  scale_color_manual(values = c('maroon','grey80', 'black'),
                     limits = c('120 s burn','30 s burn', 'control'),
                     labels = c('120 s','30 s', '0 s')) +
  scale_shape_manual(values = c('circle','triangle', 'square'),
                     limits = c('120 s burn','30 s burn', 'control'),
                     labels = c('120 s','30 s', '0 s')) +
  theme(axis.text.y = element_text(size = 14),
        legend.title = element_text(size=14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14, margin = margin(t=10, b=10)),
        axis.title.y = element_text(size = 14, margin = margin(r=10, l=10)),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.position = 'right',
        # legend.background = element_rect(fill = "darkgray"),
        legend.box = "vertical")

p
