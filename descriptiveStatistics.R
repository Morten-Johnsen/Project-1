rm(list = ls())
library(lubridate)
library(tidyverse)
library(reshape)
#setwd("Replace with path to directory containing project files.")
D <- read.table("tuno.txt", header=TRUE, sep=" ", 
                as.is=TRUE)
## Dimensions of D (number of rows and columns)
dim(D)
##  Column/variable names
names(D)
## The first rows/observations
head(D)
## The last rows/observations
tail(D)
## Selected summary statistics
summary(D)
## Another type of summary of the dataset
str(D)

D$date <- as.Date("2003-01-01")-1+D$r.day
D$pow.obs.norm <- D$pow.obs/5000

meltD <- D %>%
  select(-r.day, -month, -day, -pow.obs) %>%
  melt(id.vars = "date")

ggplot(meltD)+
  geom_histogram(aes(x = value, fill = variable), colour = "white")+
  facet_wrap(~ variable, scales = "free")

### Heatmap ###
D %>%
  select(pow.obs.norm, wd30, ws30) %>%
  cor() %>% 
  pheatmap()


hist(D$pow.obs, xlab="Power production pr day", prob=TRUE)
plot(D$r.day, D$pow.obs, type = 'l', xlab="Date", ylab="Power production pr. day", col=1)

