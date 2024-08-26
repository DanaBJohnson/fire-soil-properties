# --------------------
# Title: Fitting respiration to 2-pool decay model
# Author: TW modified by DBJ
# Date: 2022-Nov-14
#
# Purpose: Script for fitting two-pool models to exponential decay data
#           + loops

# --------------------
# Load package
library(minpack.lm)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(knitr)
# --------------------

### Running model on our own data ##

### Two-pool model ###
# Mt = M1*e^(-k1*t)+M2*e^(-k2*t)  The basic two-pool model
# Mt = 1 = M1+M2    Setting Mt as total initial C to 1, M1+M2 make up total
# Mt = M1*e^(-k1*t)+(1-M1)*e^(-k2*t)    Rearranging equation
# Mt = M1*e^(-k1*t)-M1*e^(-k2*t)+e^(-k2*t)  Rearranging equation

# Describe function form
TwoPoolDecay = function(params,t) params$M1*(exp(-params$k1*t))-params$M1*exp(-params$k2*t)+exp(-params$k2*t)

## Define residual function (our actual value minus the exponential function)
residualsFun.twopool = function(p, Mt, t) Mt - TwoPoolDecay(p, t)




# OBJECTIVE: set up input and output files for loop

df = read.csv('../data/incubations/multiplexer/processed-respiration-data.csv')  # This is a dataframe from the "processing-multiplexer-data.R"

df <- subset(df, select = c(core.id, site.id, burn.trtmt, whole.days.since.wet.up, new.fractional.C.remaining, data))
colnames(df)[5] = 'fractional.C.remaining'
head(df)

df <- df[order(df$core.id, df$whole.days.since.wet.up),]

df <- subset(df, is.na(fractional.C.remaining)==FALSE) # There can't be any NA's in the dataframe.

temp.dir <- df 
  
# Create a vector of Core IDs
# Optionally filter to pull out cores of interest
names.df <- df %>%
  subset(select = c(core.id)) %>%
  unique()

# Convert df to vector
names.vec = as.vector(names.df$core.id)
# Check if this worked. 
class(names.vec)
length(names.vec)

# --------------------
# Create an output df to populate with coefficients later
coefs.df <- data.frame("core.id" = "sample", "k1" = 0, "k2" = 0, "M1" = 0, "M2" = 0)

output.fit.df <- data.frame("core.id" = "sample", 
                            'r.squared' = 0,
                            "site.id" = 0,
                            "burn.trtmt" = "trtmt",
                            'data' = 'data',
                            "t" = 0, 
                            "Mt" = 0, "Mt.fit" = 0)



# List start parameters - can set these however
# StartPar.TwoPool = list(k1=0.005,k2=0.00001,M1=0.02)
# Dry: k1=0.005,k2=0.00001,M1=0.002
# Wet: k1=0.005,k2=0.00001,M1=0.02
# Control: k1=0.005,k2=0.00001,M1=0.02 --> SI: k1=0.00005,k2=0.0001,M1=0.0002


# Run the decay.models.R script in a for loop. 
for (i in seq_along(names.vec)) {
  
  df.x <- df %>%
    dplyr::rename(t = whole.days.since.wet.up, Mt = fractional.C.remaining) %>%
    subset(core.id == names.vec[[i]])
  
  if (df.x$burn.trtmt[i] == "high severity") {
    StartPar.TwoPool = list(k1=0.03,k2=0.000005,M1=0.0004)
  } else if (df.x$burn.trtmt[i] == "low severity") {
    StartPar.TwoPool = list(k1=0.02,k2=0.000002,M1=0.01)
  } else {
    StartPar.TwoPool = list(k1=0.01,k2=0.0001,M1=0.04)
  }
  
  TwoPool.nls.fit <- nls.lm(par=StartPar.TwoPool, fn = residualsFun.twopool, 
                            Mt = df.x$Mt, 
                            t = df.x$t, 
                            lower = c(0,0,0), upper = c(Inf, Inf, 1),
                            control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))
  
  df.x$Mt.fit <- TwoPoolDecay(as.list(coef(TwoPool.nls.fit)),df.x$t)
  
  df.lm <- lm(Mt ~ Mt.fit, data = df.x)
  df.x$r.squared <- summary(df.lm)$r.squared

  ## summary information on parameter estimates
  coefs <- summary(TwoPool.nls.fit) 
  k1 <- coefs$coefficients[[1]]
  k2 <- coefs$coefficients[[2]]
  M1 <- coefs$coefficients[[3]]
  M2 <- df.x$Mt[1] - M1
  
  coefs.df <- coefs.df %>%
    dplyr::add_row("core.id" = names.vec[[i]], "k1" = k1, "k2" = k2, "M1" = M1, "M2" = M2)
  coefs.df
  
  output.fit.df <- rbind(output.fit.df, df.x)
}



# --------------------
### Second iteration:

# Set cutoff for R squared
rsquared = 0.9

# If R-squared is above cut-off, it's good. Otherwise, redo parameters. 
output.fit.df <- merge(output.fit.df, coefs.df)
df.rsquared <- subset(output.fit.df, r.squared >rsquared)
df = temp.dir

df.redo <- subset(df, !(core.id %in% c(df.rsquared$core.id)))
#df.redo <- subset(output.fit.df, r.squared < rsquared)

names.vec <- df.redo %>%
  subset(select = c(core.id)) %>%
  unique()

#df <- read.csv("fractional-c-remaining.csv", row.names = 1)
df <- merge(names.vec, df) 


names.vec = as.vector(names.vec$core.id)
names.vec

# Create an output df to populate with coefficients later
coefs.df <- data.frame("core.id" = "sample", "k1" = 0, "k2" = 0, "M1" = 0, "M2" = 0)

output.fit.df <- data.frame("core.id" = "sample", 
                            'r.squared' = 0,
                            "site.id" = 0,
                            "burn.trtmt" = "trtmt",
                            'data' = 'data',
                            "t" = 0, 
                            "Mt" = 0, "Mt.fit" = 0)


# List start parameters - can set these however

StartPar.TwoPool = list(k1=0.001,k2=0.00001,M1=0.0001)

# StartPar.TwoPool = list(k1=0.008,k2=0.0000001,M1=0.0001)
head(df)


# Run the decay.models.R script in a for loop. 
for (i in seq_along(names.vec)) {
  
  df.x <- df %>%
    dplyr::rename(t = whole.days.since.wet.up, Mt = fractional.C.remaining) %>%
    subset(core.id == names.vec[[i]])
  
  TwoPool.nls.fit <- nls.lm(par=StartPar.TwoPool, fn = residualsFun.twopool, 
                            Mt = df.x$Mt, 
                            t = df.x$t, 
                            lower = c(0,0,0), upper = c(Inf, Inf, 1),
                            control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))
  
  df.x$Mt.fit <- TwoPoolDecay(as.list(coef(TwoPool.nls.fit)),df.x$t)
  
  df.lm <- lm(Mt ~ Mt.fit, data = df.x)
  df.x$r.squared <- summary(df.lm)$r.squared
  
  ## summary information on parameter estimates
  coefs <- summary(TwoPool.nls.fit) 
  k1 <- coefs$coefficients[[1]]
  k2 <- coefs$coefficients[[2]]
  M1 <- coefs$coefficients[[3]]
  M2 <- df.x$Mt[1] - M1
  
  coefs.df <- coefs.df %>%
    add_row("core.id" = names.vec[[i]], "k1" = k1, "k2" = k2, "M1" = M1, "M2" = M2)
  coefs.df
  
  output.fit.df <- rbind(output.fit.df, df.x)
}


# If R-squared is above cut-off, it's good. Otherwise, redo parameters.
output.fit.df <- merge(output.fit.df, coefs.df)
df.rsquared <- df.rsquared %>%
  rbind(output.fit.df) %>%
  subset(r.squared >rsquared)






# --------------------
### Third iteration


#df.redo <- subset(output.fit.df, r.squared < rsquared)
df <- temp.dir
df.redo <- subset(df, !(core.id %in% c(df.rsquared$core.id)))

names.vec <- df.redo %>%
  subset(select = c(core.id)) %>%
  unique()

#df <- read.csv("fractional-c-remaining.csv", row.names = 1)
df <- merge(names.vec, df)

names.vec = as.vector(names.vec$core.id)
names.vec

# Create an output df to populate with coefficients later
coefs.df <- data.frame("core.id" = "sample", "k1" = 0, "k2" = 0, "M1" = 0, "M2" = 0)

output.fit.df <- data.frame("core.id" = "sample", 
                            'r.squared' = 0,
                            "site.id" = 0,
                            "burn.trtmt" = "trtmt",
                            'data' = 'data',
                            "t" = 0, 
                            "Mt" = 0, "Mt.fit" = 0)

# List start parameters - can set these however
# StartPar.TwoPool = list(k1=0.005,k2=0.00001,M1=0.0001)

StartPar.TwoPool = list(k1=0.001,k2=0.00001,M1=0.0001) # this one works

# StartPar.TwoPool = list(k1=0.005,k2=0.00001,M1=0.0001)

# Run the decay.models.R script in a for loop. 
for (i in seq_along(names.vec)) {
  
  df.x <- df %>%
    dplyr::rename(t = whole.days.since.wet.up, Mt = fractional.C.remaining) %>%
    subset(core.id == names.vec[[i]])
  
  TwoPool.nls.fit <- nls.lm(par=StartPar.TwoPool, fn = residualsFun.twopool, 
                            Mt = df.x$Mt, 
                            t = df.x$t, 
                            lower = c(0,0,0), upper = c(Inf, Inf, 1),
                            control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))
  
  df.x$Mt.fit <- TwoPoolDecay(as.list(coef(TwoPool.nls.fit)),df.x$t)
  
  df.lm <- lm(Mt ~ Mt.fit, data = df.x)
  df.x$r.squared <- summary(df.lm)$r.squared
  
  ## summary information on parameter estimates
  coefs <- summary(TwoPool.nls.fit) 
  k1 <- coefs$coefficients[[1]]
  k2 <- coefs$coefficients[[2]]
  M1 <- coefs$coefficients[[3]]
  M2 <- df.x$Mt[1] - M1
  
  coefs.df <- coefs.df %>%
    add_row("core.id" = names.vec[[i]], "k1" = k1, "k2" = k2, "M1" = M1, "M2" = M2)
  coefs.df
  
  output.fit.df <- rbind(output.fit.df, df.x)
}


# If R-squared is above cut-off, it's good. Otherwise, redo parameters. 
output.fit.df <- merge(output.fit.df, coefs.df)

df.rsquared <- df.rsquared %>%
  rbind(output.fit.df) %>%
  subset(r.squared >rsquared)
#df.redo <- subset(output.fit.df, r.squared < rsquared)

df <- temp.dir
df.redo <- subset(df, !(core.id %in% c(df.rsquared$core.id)))

names.vec <- df.redo %>%
  subset(select = c(core.id)) %>%
  unique()
names.vec

# --------------------
### Fourth iteration

#df <- read.csv("fractional-c-remaining.csv", row.names = 1)
df <- merge(names.vec, df)


names.vec = as.vector(names.vec$core.id)
names.vec
# Create an output df to populate with coefficients later
coefs.df <- data.frame("core.id" = "sample", "k1" = 0, "k2" = 0, "M1" = 0, "M2" = 0)

output.fit.df <- data.frame("core.id" = "sample", 
                            'r.squared' = 0,
                            "site.id" = 0,
                            "burn.trtmt" = "trtmt",
                            'data' = 'data',
                            "t" = 0, 
                            "Mt" = 0, "Mt.fit" = 0)


# List start parameters - can set these however
StartPar.TwoPool = list(k1=0.003,k2=0.000005,M1=0.0004)

# StartPar.TwoPool = list(k1=0.001,k2=0.0000001,M1=0.0001) # this one works

# Run the decay.models.R script in a for loop. 
for (i in seq_along(names.vec)) {
  df.x <- df %>%
    dplyr::rename(t = whole.days.since.wet.up, Mt = fractional.C.remaining) %>%
    subset(core.id == names.vec[[i]])
  
  TwoPool.nls.fit <- nls.lm(par=StartPar.TwoPool, fn = residualsFun.twopool, 
                            Mt = df.x$Mt, 
                            t = df.x$t, 
                            lower = c(0,0,0), upper = c(Inf, Inf, 1),
                            control = nls.lm.control(nprint=1,maxiter=500,ftol=sqrt(.Machine$double.eps)/10))
  df.x$Mt.fit <- TwoPoolDecay(as.list(coef(TwoPool.nls.fit)),df.x$t)
  df.lm <- lm(Mt ~ Mt.fit, data = df.x)
  df.x$r.squared <- summary(df.lm)$r.squared
  ## summary information on parameter estimates
  coefs <- summary(TwoPool.nls.fit) 
  k1 <- coefs$coefficients[[1]]
  k2 <- coefs$coefficients[[2]]
  M1 <- coefs$coefficients[[3]]
  M2 <- df.x$Mt[1] - M1
  coefs.df <- coefs.df %>%
    add_row("core.id" = names.vec[[i]], "k1" = k1, "k2" = k2, "M1" = M1, "M2" = M2)
  coefs.df
  output.fit.df <- rbind(output.fit.df, df.x)
}


# If R-squared is above cut-off, it's good. Otherwise, redo parameters. 
output.fit.df <- merge(output.fit.df, coefs.df)






df.rsquared <- df.rsquared %>%
  rbind(output.fit.df)

#df.redo <- subset(output.fit.df, r.squared < rsquared)
df.redo <- subset(df, !(core.id %in% c(df.rsquared$core.id)))

names.vec <- df.redo %>%
  subset(select = c(core.id)) %>%
  unique()
names.vec
# Complete at R^2 = 0.98





# --------------------
# Clean up output files

head(df.rsquared)
unique(df.rsquared$site.id)
summary(df.rsquared$k1)
summary(df.rsquared$k2)
df.rsquared = subset(df.rsquared, r.squared != 0)
summary(df.rsquared$r.squared)

# Sanity check:
ggplot(subset(df.rsquared, site.id %in% c(5,7,4))) + 
  geom_point(aes(x=t, y = Mt, color=burn.trtmt))+
  geom_line(aes(x=t, y = Mt.fit, color=burn.trtmt))+
  facet_wrap(~site.id)
 
# write.csv(df.rsquared,'decay-model-coefs.csv', row.names = FALSE)

# df.rsquared.prelim <- write.csv(df.rsquared,'../../../incubations/2-pool-decay-model-output.csv', row.names = FALSE)

df.coefs <- df.rsquared %>%
  subset(select = c(core.id,r.squared,burn.trtmt,
                    k1,k2,M1,M2)) %>%
  unique()

df.fit <- df.rsquared %>%
  subset(select = c(core.id,site.id,r.squared,
                    burn.trtmt,Mt,Mt.fit))
  
# Look at the coefficient dataframe
hist(df.rsquared$k1,breaks=20)
hist(df.rsquared$k2,breaks=20)
hist(df.rsquared$M1,breaks=20)
hist(df.rsquared$M2,breaks=20)


#PLOT: Fractional C remaining at one site (13).
colnames(df.rsquared)
head(df.rsquared)



p = ggplot(df.rsquared, aes(x = t, y = Mt, 
                     color = burn.trtmt,
                     shape = burn.trtmt)) + 
  geom_point() +
  geom_line(aes(x=t, y = Mt.fit),color='black') +
  facet_wrap(~site.id) +
  theme_bw() +
  # ylim(0.997,1)+
  labs(x = "Day", 
       y = "C remaining as fraction of initial total C",
       title = "Two pool: Fractional C remaining during short-term incubation treatments")
p


values = subset(df.rsquared, select = c(core.id, r.squared, burn.trtmt, k1, k2, M1, M2)) %>%
  unique() %>%
  subset(burn.trtmt != 'trtmt')

p = ggplot(values, aes(x = burn.trtmt, 
                            y = k1)) + 
  geom_boxplot() +
  theme_bw() + 
  scale_x_discrete(limits = c('control','low severity','high severity'),
                   labels = c('0 s', '30 s','120 s')) +
  labs(x = "Burn duration", 
       y = "Decay rate coefficient \nfor fast C pool (k1)")  + 
  #ylim(0,0.025) + 
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





### Figures for poster...
# Respiration coefficients ----
df.coefs <- read.csv('../data/incubations/df.rsquared.with.estimated.initial.C.csv')

df.coefs <- mutate(df.coefs, Day = t) %>%
  subset(select = c(core.id, k1, k2, M1, M2)) %>%
  unique()

df.coefs <- merge(df.coefs, subset(df_C, select = c(core.id, vegetation, burn.trtmt)))

df.stack <- gather(df.coefs, 
                   key = 'coefs', 
                   value = value,M1,M2,k1,k2)

df.stack$coefs <- factor(df.stack$coefs, levels = c('M1', 'M2', 'k1', 'k2'))

pResp2 = ggplot(subset(df.stack, coefs %in% c('M1','M2')), aes(x=burn.trtmt)) + 
  geom_boxplot(aes(y=value, fill=burn.trtmt), alpha=0.7,varwidth = FALSE)+
  theme_bw()+
  labs(x= 'Burn treatment',y = 'Fractional pool size',
       fill = 'Burn treatment', color = 'Burn treatment') + 
  scale_fill_manual(values = c('black','grey80', 'darkorange'),
                    limits = c('control','low severity','high severity'),
                    labels = c('control','30 s burn','120 s burn')) +
  scale_x_discrete(limits = c('control','low severity','high severity'),
                   labels = c('Unburned\ncontrol','30 s burn','120 s burn')) +
  #scale_x_discrete(limits = c('O','A'), 
  #                 labels = c('Organic','Mineral'))
  facet_wrap(~coefs, ncol=1,scales = 'free_y')+
             #, labeller = labeller(coefs = CoefsLabel
  theme(axis.text.y = element_text(size = 14),
        legend.title = element_text(size=14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14, margin = margin(t=10, b=10)),
        axis.title.y = element_text(size = 14, margin = margin(r=25, l=10)),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        # legend.background = element_rect(fill = "darkgray"),
        legend.box = "vertical", 
        legend.position = '',
        strip.background = element_rect(colour="black", fill="grey92"))


pResp1 = ggplot(subset(df.stack, coefs %in% c('k1','k2')), aes(x=burn.trtmt)) + 
  geom_boxplot(aes(y=value, fill=burn.trtmt), alpha=0.7,varwidth = FALSE)+
  theme_bw()+
  labs(x= 'Burn treatment',y = 'Decay rate',
       fill = 'Burn treatment', color = 'Burn treatment') + 
  scale_fill_manual(values = c('black','grey80', 'darkorange'),
                    limits = c('control','low severity','high severity'),
                    labels = c('control','30 s burn','120 s burn')) +
  scale_x_discrete(limits = c('control','low severity','high severity'),
                   labels = c('Unburned\ncontrol','30 s burn','120 s burn')) +
  #scale_x_discrete(limits = c('O','A'), 
  #                 labels = c('Organic','Mineral'))
  facet_wrap(~coefs, ncol=1,scales = 'free_y')+
  #, labeller = labeller(coefs = CoefsLabel
  theme(axis.text.y = element_text(size = 14),
        legend.title = element_text(size=14),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14, margin = margin(t=10, b=10)),
        axis.title.y = element_text(size = 14, margin = margin(r=25, l=10)),
        strip.text = element_text(size=14),
        legend.text = element_text(size=14),
        # legend.background = element_rect(fill = "darkgray"),
        legend.box = "vertical", 
        legend.position = '',
        strip.background = element_rect(colour="black", fill="grey92"))

cowplot::plot_grid(pResp2, pResp1, ncol=2)

test = aov(k2 ~ burn.trtmt, df.coefs)
summary(test)
TukeyHSD(test)


for (i in 1:nrow(df.coefs)) {
  if (df.coefs$burn.trtmt[i] == 'control') {
    df.coefs$burned[i] = 'no'
  } else {
    df.coefs$burned[i] = 'yes'
  }
}

test = wilcox.test(k2~burned, data = df.coefs)
test


test = wilcox.test(k2~burn.trtmt, data = subset(df.coefs, burn.trtmt != 'control'))
test

test = aov(k2~burn.trtmt, data = df.coefs)
test
summary(test)
TukeyHSD(test)