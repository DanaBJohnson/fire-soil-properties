##title: "Processing pH data"
##author: "Dana Johnson"
##date: "11/08/2022"


library(ggplot2)
library(dplyr)
library(tidyr)

# Import pH data
df.pH <- read.csv("../data/soil-properties/pH.csv")
head(df.pH)
colnames(df.pH)[1] <- 'core.id'

# Add vegetation and burn treatment to dataframe
sample_list <- read.csv("../methods/sample-list-for-incubations.csv")
head(sample_list)
colnames(sample_list)[1] <- 'core.id'

merged = merge(df.pH, sample_list)
merged$horizon <- factor(merged$horizon, levels = c('O','M'))


# Plot pH by burn treatment, horizon, incubation length, and vegetation
p = ggplot(subset(merged), aes(x=burn.trtmt, y = pH)) + 
  geom_boxplot(aes(fill=burn.trtmt), alpha=0.8)+
  #geom_point(aes(color = burn.trtmt), alpha = 0.5) + 
  facet_wrap(incubation.period.days+vegetation~horizon, scale = 'free') + 
  theme_bw() + 
  scale_fill_manual(values = c('black','grey80', 'darkorange'),
                    limits = c('control','low severity','high severity'),
                    labels = c('control','30 s','120 s')) +
  scale_x_discrete(limits = c('control','low severity','high severity'),
                   labels = c('0 s', '30 s','120 s')) +
  labs(x = 'Burn duration',
       y = 'Soil pH',
       fill = 'Burn duration') + 
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14)) + 
  ylim(4,9)
p

hist(subset(merged, incubation.period.days==2 & burn.trtmt == 'high severity')$pH)

test <- aov(pH~burn.trtmt, subset(merged, 
                                  incubation.period.days==2 & 
                                  vegetation != 'Pinus_banksiana' &
                                  horizon == 'O'))
test
summary(test)
TukeyHSD(test)

# Now calculate change in pH with burning:
df.control <- subset(merged, burn.trtmt == 'control') %>%
  subset(select = c(site.id, horizon, incubation.period.days, pH))
colnames(df.control)[4] <- 'pH.control'

df.burn <- subset(merged, burn.trtmt != 'control') %>%
  merge(df.control, by = c('site.id', 'horizon', 'incubation.period.days')) %>%
  mutate(change.pH = pH - pH.control)

df.burn$incubation.period.days <- as.character(df.burn$incubation.period.days)

# Plot pH by burn treatment, horizon, incubation length, and vegetation
p = ggplot(subset(df.burn, horizon == 'O'), aes(x=incubation.period.days, y = change.pH)) + 
  geom_boxplot(aes(fill=burn.trtmt), alpha=0.8)+
  #geom_point(aes(color = burn.trtmt), alpha = 0.5) + 
  facet_wrap(vegetation~horizon, scale = 'free') + 
  theme_bw() + 
  scale_fill_manual(values = c('gold1', 'darkorange2'),
                    limits = c('low severity','high severity'),
                    labels = c('30 s','120 s')) +
  #scale_x_discrete(limits = c('low severity','high severity'),
  #                 labels = c('30 s','120 s')) +
  labs(x = 'Days after burning',
       y = 'Delta pH with burning',
       fill = 'Burn duration') + 
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14))  
  
p


## Calculate change in pH with time since burn:
df.2days <- subset(merged,incubation.period.days == 2) %>%
  subset(select = c(site.id, horizon, pH))
colnames(df.2days)[3] <- 'pH.2days'

df.70days <- subset(merged, incubation.period.days == 70) %>%
  merge(df.2days, by = c('site.id','horizon')) %>%
  mutate(Delta.pH.with.time = pH-pH.2days)

# Plot pH by burn treatment, horizon, incubation length, and vegetation
p = ggplot(subset(df.70days, horizon == 'O' ), aes(x=burn.trtmt, y = Delta.pH.with.time)) + 
  geom_boxplot(aes(fill=burn.trtmt), alpha=0.8)+
  #geom_point(aes(color = burn.trtmt), alpha = 0.5) + 
  facet_wrap(vegetation~horizon, scale = 'free') + 
  theme_bw() + 
  scale_fill_manual(values = c('black','grey80', 'darkorange'),
                    limits = c('control','low severity','high severity'),
                    labels = c('control','30 s','120 s')) +
  scale_x_discrete(limits = c('control','low severity','high severity'),
                   labels = c('0 s', '30 s','120 s')) +
  labs(x = 'Burn duration',
       y = 'Change pH from 2 days to 10 wks post-burn',
       fill = 'Burn duration') + 
  theme(strip.text = element_text(size = 12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14)) 
p
