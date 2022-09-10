rm(list = ls())
library(lubridate)
library(tidyverse)
library(reshape)
library(pheatmap)
library(openair)
setwd("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1")
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

par(mfrow=c(1,2))
plot(D$date, D$pow.obs, type = 'l', xlab="Date", ylab="Average daily power production [kW]", 
     main = 'Development in average daily power production over time', cex.main = 0.8, col=1)
hist(D$pow.obs, xlab="Power production [kW]", main='Distribution of average daily power production', cex.main=0.8)
### for outliers ###
outlierFUN <- function(data, quantiles){
  v <- quantile(x = data, probs = quantiles)
  IQR <- v[2] - v[1]
  outliers <- ( ( data < (v[1] - 1.5 * IQR) ) | ( data > (v[2] + 1.5 * IQR) ) ) * data
  return (outliers)
}
D$outlierws30 <- outlierFUN(data = D$ws30, quantiles = c(0.25,0.75))
###              ###
par(mfrow=c(1,2))
plot(D$date, D$ws30, type = 'l', xlab="Date", ylab="Wind speeds [m/s]", cex.main = 0.8, col=1,
     main='Development in wind speeds over time')
points(x = D$date ,y = D$outlierws30, type = 'p', pch = 19, col = "red", cex = 0.25)
hist(D$ws30, xlab="Wind speeds [m/s]", main='Distribution of wind speeds', cex.main = 0.8)

par(mfrow=c(1,2))
plot(D$date, D$wd30, type = 'l', xlab="Date", ylab=expression(paste("Wind direction. N = 0, E = ", pi/2)), col=1,
     main='Development in wind directions over time', cex.main = 0.8)
lines(D$date, replicate(length(D$date), 3*pi/2), type='l', col=2) #W
lines(D$date, replicate(length(D$date), pi), type='l', col=3) #S
lines(D$date, replicate(length(D$date), pi/2), type='l', col=4) #E
lines(D$date, replicate(length(D$date), 0), type='l', col=5) #N
legend('topleft', legend = c('wd30', 'W', 'S', 'E', 'N'), col = 1:5, lty = 1, cex = 0.5)
#
hist(D$wd30, xlab=expression(paste('Wind direction. N = 0, E = ', frac(pi,2))),
     main='Distribution of wind directions', cex.main = 0.8, freq = TRUE) #####hist to show that wind rose is fine
abline(v = 3*pi/2, col=2)
abline(v = pi, col=3)
abline(v = pi/2, col=4)
abline(v = 0, col=5)
legend('topleft', legend = c('wd30', 'W', 'S', 'E', 'N'), col = 1:5, lty = 1, cex = 0.5)

D$wd30dg <- D$wd30 * 180/pi
###wind d
intv <- 4
#Dwinter <- D[1:54,] #for testing the season plots above :) D$date[55] = 2003-03-01
#Dsummer <- D[140:230,] #for testing the season plots above :) D$date[140] = 2003-06-01
windRose(D, ws = "ws30", wd = "wd30dg", ws2 = NA, wd2 = NA,
         ws.int = intv, angle = 30, type = "default", bias.corr = TRUE, cols
         = "heat", grid.line = list(value=4, lty=4, col="lightgrey"), width = 1, seg = NULL, auto.text
         = TRUE, breaks = round(max(D$ws30)/intv), offset = 10, normalise = FALSE, max.freq =
           NULL, paddle = FALSE, key.header = "Wind speed at 30 m", key.footer = "(m/s)",
         key.position = "bottom", key = list(height=2), dig.lab = 3, statistic =
           "prop.count", pollutant = NULL, annotate = FALSE, angle.scale =
           45, border = "black", main="Wind directions distribution (at 30 m)",
         cex.main=0.75) #change TYPE = "season" to "default" for a single plot
#angle = 10 v 45, ws.int = 0.5, 1 v 2, bias.corr = try both, cols = increment, heat, jet, hue, u.d.
summary(D)

