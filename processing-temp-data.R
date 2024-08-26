##title: "Processing temp data"
##author: "Dana Johnson"
##date: "11/08/2022"


library(ggplot2)
library(dplyr)
library(tidyr)

setwd('../Box/WhitmanLab/Projects/WoodBuffalo/ModellingFires_DOE_TES')

ThermoPositions = c('O-Mnrl interface','Core base','Midpoint')
names(ThermoPositions) = c('U','L','M')

# Import raw data
df.Temp1 <- read.csv("data/burn-simulations/burn-data-Aug15-thermocouples_reformatted.csv")
df.Temp2 <- read.csv("data/burn-simulations/burn-data-Aug16-thermocouples_reformatted.csv")
df.Temp3 <- read.csv("data/burn-simulations/burn-data-Aug17-thermocouples_reformatted.csv")
df.Temp4 <- read.csv("data/burn-simulations/burn-data-Aug-23-thermocouples_corrected_reformatted.csv")
df.Temp5 <- read.csv("data/burn-simulations/burn-data-Aug24-thermocouples_reformatted.csv")
df.Temp6 <- read.csv("data/burn-simulations/burn-data-Aug-29-thermocouples_corrected_reformatted.csv")

head(df.Temp6)

# Correct heading in df.Temp1:
colnames(df.Temp1)[42] <- 'X12.12.M'
colnames(df.Temp1)[43] <- 'X12.13.M'

# Clean up all the data frames:
df.Temp1 <-  df.Temp1[,c(3,10:43)] %>%
  mutate(burn.set = '1') %>%
  pivot_longer(cols = c(2:35),names_to = c("site", 'core','thermocouple_position'), 
               names_pattern = "X(..).(..).(.)",  
               values_to = "temperature.C") %>%
  unite('core.id', c(site, core), sep = '-', remove = FALSE)%>%
  subset(temperature.C<1200) %>%
  group_by(core.id) %>%
  mutate(start.time.sec = min(Time..sec.),
         run.time.sec = Time..sec.-start.time.sec) 


df.Temp2 <-  df.Temp2[,c(3,10:43)] %>%
  mutate(burn.set = '2') %>%
  pivot_longer(cols = c(2:35),names_to = c("site", 'core','thermocouple_position'), 
               names_pattern = "X(..).(..).(.)",  
               values_to = "temperature.C") %>%
  unite('core.id', c(site, core), sep = '-', remove = FALSE) %>%
  subset(temperature.C<1200) %>%
  group_by(core.id) %>%
  mutate(start.time.sec = min(Time..sec.),
         run.time.sec = Time..sec.-start.time.sec) 



df.Temp3 <-  df.Temp3[,c(3,10:43)] %>%
  mutate(burn.set = '3') %>%
  pivot_longer(cols = c(2:35),names_to = c("site", 'core','thermocouple_position'), 
               names_pattern = "X(..).(..).(.)",  
               values_to = "temperature.C") %>%
  unite('core.id', c(site, core), sep = '-', remove = FALSE) %>%
  subset(temperature.C<1200) %>%
  group_by(core.id) %>%
  mutate(start.time.sec = min(Time..sec.),
         run.time.sec = Time..sec.-start.time.sec)




df.Temp4 <-  df.Temp4[,c(3,10:43)] %>%
  mutate(burn.set = '4') %>%
  pivot_longer(cols = c(2:35),names_to = c("site", 'core','thermocouple_position'), 
               names_pattern = "X(..).(..).(.)",  
               values_to = "temperature.C") %>%
  unite('core.id', c(site, core), sep = '-', remove = FALSE) %>%
  subset(temperature.C<1200) %>%
  group_by(core.id) %>%
  mutate(start.time.sec = min(Time..sec.),
         run.time.sec = Time..sec.-start.time.sec) 




df.Temp5 <-  df.Temp5[,c(3,10:43)] %>%
  mutate(burn.set = '5') %>%
  pivot_longer(cols = c(2:35),names_to = c("site", 'core','thermocouple_position'), 
               names_pattern = "X(..).(..).(.)",  
               values_to = "temperature.C") %>%
  unite('core.id', c(site, core), sep = '-', remove = FALSE) %>%
  subset(temperature.C<1200) %>%
  group_by(core.id) %>%
  mutate(start.time.sec = min(Time..sec.),
         run.time.sec = Time..sec.-start.time.sec) 




df.Temp6 <-  df.Temp6[,c(3,10:43)] %>%
  mutate(burn.set = '6') %>%
  pivot_longer(cols = c(2:35),names_to = c("site", 'core','thermocouple_position'), 
               names_pattern = "X(..).(..).(.)",  
               values_to = "temperature.C") %>%
  unite('core.id', c(site, core), sep = '-', remove = FALSE) %>%
  subset(temperature.C<1200) %>%
  group_by(core.id) %>%
  mutate(start.time.sec = min(Time..sec.),
         run.time.sec = Time..sec.-start.time.sec) 




# Clean up data frames
df.Temp <- rbind(df.Temp1, df.Temp2, df.Temp3, df.Temp4, df.Temp5, df.Temp6)

rm(df.Temp1, df.Temp2, df.Temp3, df.Temp4, df.Temp5, df.Temp6)

df.Temp.decimated <-df.Temp %>% 
  filter((run.time.sec/10) %% 1 == 0)

# Save output file:
# write.csv(df.Temp.decimated,'data/burn-simulations/Output-files/Compiled-decimated-temperature-all-cores.csv')


### PLOTS: 
# Take one data point per minute (to make it faster to graph)
# df.Temp.decimated <- read.csv('../Box/WhitmanLab/Projects/WoodBuffalo/ModellingFires_DOE_TES/data/burn-simulations/Output-files/Compiled-decimated-temperature-all-cores.csv')

df.Temp.filtered <- df.Temp.decimated %>%
  filter((run.time.sec/60) %% 1 == 0) %>%
  subset(run.time.sec/3600 < 7)

df.duration <- read.csv('data/burn-simulations/core-prep.csv')

df.duration <- df.duration %>%
  subset(select = c(core.id, vegetation, duration.s)) %>%
  separate(core.id, c('year','project','site','core'), sep = '-', remove = TRUE) %>%
  unite('core.id', c(site,core), sep='-', remove = TRUE)

df.Temp.filtered <- merge(df.Temp.filtered, df.duration, by = 'core.id')

ggplot(df.Temp.filtered, aes(x=run.time.sec/3600, y = temperature.C)) + 
  geom_point(aes(color = as.character(duration.s))) +
  facet_grid(thermocouple_position~vegetation, 
             scales = 'free',
             labeller = labeller(thermocouple_position = ThermoPositions)) + 
  theme_bw() + 
  scale_color_manual(values = c('orange','red4'),
                     limits = c('30','120'), 
                     labels = c('short duration burn','long duration burn'))+
  labs(x = 'Time (h)',
       y = 'Temperature (C)', 
       color = 'Treatment') + 
  theme(legend.position = 'bottom') 
