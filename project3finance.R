rm(list = ls())
library(lubridate)
library(tidyverse)
library(reshape)
library(pheatmap)
library(openair)
library(PowerNormal)
library(sn)
library(gnorm)
library(emg)
#setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/")
setwd("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1")
D <- read.table("finance_data.csv", header=TRUE, sep=";", 
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

D$date <- as.Date(D$time)
D$year <- year(D$date)

#ggplot(D, aes(x = date, y = SLV, colour = factor(year))) + geom_point()
ggplot(D, aes(x = date, y = SLV)) + geom_point()
#plot(D$SLV)
ggplot(D, aes(x = SLV)) +
  geom_histogram(aes(y = ..density..), color = 'black')


#########################################################################################################
#########################################################################################################
source("testDistribution.R")
#########################################################################################################
#########################################################################################################
par <- nlminb(start = c(1,1), objective = testDistribution,
              distribution = "normal",
              x = D$SLV)
ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') +#color, fill
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = par$par[1], sd = par$par[2]), color='red')

plot(ecdf(D$SLV), verticals = T)
xseq <- seq(0.9*min(D$SLV), 1.1*max(D$SLV), length.out=100)
lines(xseq, pnorm(xseq, mean(D$SLV), sd(D$SLV)), col='red')
#plot(xseq, pnorm(xseq, mean(D$SLV), sd(D$SLV)), col='red')
qqnorm(D$SLV)
qqline(D$SLV)
#################################################
#################################################
lcauchyFUNC <- function(p, data){
  x0 <- p[1] #location R
  gam <- p[2] #scale R > 0
  return(-sum(dcauchy(x = data, location = x0, scale = gam, log = T)))
}
lpownormFUNC <- function(p, data){
  alpha <- p
  return(-sum(log(dpn(x = data, p))))
}
ltFUNC <- function(p, data){
  return(-sum(dt(x = data, df = p, log = T)))
}
lsnFUNC <- function(p, data){ #skewed normal dist
  return(-sum(dsn(x = data, xi = p[1], omega = p[2], alpha = p[3], log = T)))
}
lgnFUNC <- function(p, data){ #symmetric generalized normal dist
  return(-sum(dgnorm(x = data, mu = p[1], alpha = p[2], beta = p[3], log = T)))
}
# lgnFUNC <- function(p, data){
#   mu <- p[1]
#   alpha <- p[2]
#   beta <- p[3]
#   -sum( log( beta/(2 * alpha * gamma(1/beta)) * exp(-(abs(data - mu)/alpha)^beta) ) )
# }
lasgnFUNC <- function(p, data){ #asymmetric generalized normal dist, when K = 0 has already been checked
  epsilon <- p[1]
  alpha <- p[2]
  kappa <- p[3]
  return(-sum( log(dnorm(x = -1/kappa * log(1 - kappa * (data - epsilon) / alpha) ) /
              (alpha - kappa * (data - epsilon)) ) ) )
}

lemgFUNC <- function(p, data){ #exponential modified gaussian dist
  return(-sum(demg(x = data, mu = p[1], sigma = p[2], lambda = p[3], log = T)))
}

par.cauchy <- nlminb(start = c(0,1), objective = lcauchyFUNC, data = D$SLV)
par.pownorm <- nlminb(start = 1, objective = lpownormFUNC, data = D$SLV)
par.t <- nlminb(start = 1, objective = ltFUNC, data = D$SLV)
par.sn <- nlminb(start = c(1,1,1), objective = lsnFUNC, data = D$SLV)
par.gn <- nlminb(start = c(1,1,1), objective = lgnFUNC, data = D$SLV)
par.asgn <- nlminb(start = c(1,1,1), lower = c(-Inf, -Inf, 0), objective = lasgnFUNC, data = D$SLV)
par.emg <- nlminb(start = c(1,1,1), lower = c(-Inf, 1/1000, 1/1000), objective = lemgFUNC, data = D$SLV)

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = par$par[1], sd = par$par[2]), aes(colour = "norm")) +
  stat_function(fun = dcauchy, n = dim(D)[1], args = list(location = par.cauchy$par[1],
                                                          scale = par.cauchy$par[2]), aes(colour = "cauchy")) +
  stat_function(fun = dpn, n = dim(D)[1], args = list(alpha = par.pownorm$par), aes(colour = "powernorm")) +
  stat_function(fun = dt, n = dim(D)[1], args = list(df = par.t$par), aes(colour = "t")) +
  stat_function(fun = dsn, n = dim(D)[1], args = list(xi = par.sn$par[1], omega = par.sn$par[2],
                                                      alpha = par.sn$par[3]), aes(colour = "sn")) +
  stat_function(fun = dgnorm, n = dim(D)[1], args = list(mu = par.gn$par[1], alpha = par.gn$par[2],
                                                         beta = par.gn$par[3]), aes(colour = "gnorm")) +
  stat_function(fun = demg, n = dim(D)[1], args = list(mu = par.emg$par[1], sigma = par.emg$par[2],
                                                       lambda = par.emg$par[3]), aes(colour = "emg"))+
  scale_colour_manual(values = c("blue", "red", "yellow", "black", "pink", "grey", "purple"))+
  labs(colour = "Distribution")
#legend('topright', legend=c('normal', 'cauchy', 'power normal', 't'), col=c('red', 'blue', 'green', 'yellow'))

AIC.norm <- -2 * sum(dnorm(x = D$SLV, mean = par$par[1], sd = par$par[2], log = T))
+ 2 * length(par$par)
AIC.cauchy <- -2 * sum(dcauchy(x = D$SLV, location = par.cauchy$par[1],
                               scale = par.cauchy$par[2], log = T)) + 2 * length(par.cauchy$par)
AIC.pownorm <- -2 * sum(log(dpn(x = D$SLV, alpha = par.pownorm$par)))
+ 2 * length(par.pownorm$par)
AIC.t <- -2 * sum(dt(x=D$SLV, df = par.t$par, log = T)) + 2 * length(par.t$par)
AIC.sn <- -2 * sum(dsn(x=D$SLV, xi = par.sn$par[1], omega = par.sn$par[2], alpha = par.sn$par[3], 
                       log = T)) + 2 * length(par.sn$par)
AIC.gn <- -2 * sum(dgnorm(x=D$SLV, mu = par.gn$par[1], alpha = par.gn$par[2], beta = par.gn$par[3], 
                          log = T)) + 2 * length(par.gn$par)
AIC.asgn <- -2 * sum( log(dnorm(x = -1/par.asgn$par[3] * log(1 - par.asgn$par[3] * (D$SLV - par.asgn$par[1]) / par.asgn$par[2]) ) /
                            (par.asgn$par[2] - par.asgn$par[3] * (D$SLV - par.asgn$par[1])) ) ) + 2 * length(par.asgn$par)
AIC.emg <- -2 * sum(demg(x=D$SLV, mu = par.emg$par[1], sigma = par.emg$par[2], lambda = par.emg$par[3], 
                         log = T)) + 2 * length(par.emg$par)


round(rbind(AIC.norm, AIC.cauchy, AIC.pownorm, AIC.t, AIC.sn, AIC.gn, AIC.asgn, AIC.emg), digits=5)

# boxplot(D$SLV)
n <- 100000
par(mfrow=c(2,3))
hist(rnorm(n, mean = par$par[1], sd = par$par[2]))
hist(D$SLV)
hist(rsn(n, xi = par.sn$par[1], omega = par.sn$par[2], alpha = par.sn$par[3]))
hist(rgnorm(n, mu = par.gn$par[1], alpha = par.gn$par[2], beta = par.gn$par[3]))
hist(rcauchy(n, location = par.cauchy$par[1], scale = par.cauchy$par[2]))
hist(remg(n, mu = par.emg$par[1], sigma = par.emg$par[2], lambda = par.emg$par[3]))

# ggplot(D)+
#   #geom_histogram(aes(x = rgnorm(dim(D)[1], mu = par.gn$par[1], alpha = par.gn$par[2], beta = par.gn$par[3]), y=..density..), color='black') +
#   geom_histogram(aes(x = SLV, y= ..density..,), color='red') + #color, fill
#   geom_histogram(aes(x = rgnorm(dim(D)[1], mu = par.gn$par[1], alpha = par.gn$par[2], beta = par.gn$par[3]), y=..density..), color='black')
lgamFUNC <- function(p, norm_data){
  k <- p[1] #shape
  beta <- p[2] # rate
  return(-sum(dgamma(x = norm_data, shape = k, rate = beta, log = T)))
}
lbetaFUNC <- function(p, norm_data){
  alpha <- p[1] #shape
  beta <- p[2] #¿shape?
  -sum(dbeta(x = norm_data, shape1 = alpha, shape2 = beta, log = T))
}

D$SLV.norm <- ( D$SLV - min(D$SLV) ) / (max(D$SLV) - min(D$SLV))
par.gam <- nlminb(start = c(0.5, 0.5), objective = lgamFUNC, norm_data = D$SLV.norm)
par.beta <- nlminb(start = c(0.5, 0.5), objective = lbetaFUNC, norm_data = D$SLV.norm)
par.gn.norm <- nlminb(start = c(1,1,1), objective = lgnFUNC, data = D$SLV.norm)

ggplot(D)+
  geom_histogram(aes(x = SLV.norm, y= ..density..,), color='black') + 
  stat_function(fun = dgamma, n = dim(D)[1], args = list(shape = par.gam$par[1], scale = par.gam$par[2]), aes(colour = 'gamma')) +
  stat_function(fun = dbeta, n = dim(D)[1], args = list(shape1 = par.beta$par[1], shape2 = par.beta$par[2]), aes(colour = 'beta')) +
  stat_function(fun = dgnorm, n = dim(D)[1], args = list(mu = par.gn.norm$par[1],
      alpha = par.gn.norm$par[2], beta = par.gn.norm$par[3]), aes(colour='gnorm')) +
  labs(colour = 'Distributions')

-2 * sum(dgamma(x=D$SLV.norm, shape = par.gam$par[1], rate = par.gam$par[2], log = T)) + 2 * length(par.gam$par)
-2 * sum(dbeta(x=D$SLV.norm, shape1 = par.beta$par[1], shape2 = par.beta$par[2], log = T)) + 2 * length(par.beta$par)
-2 * sum(dgnorm(x = D$SLV.norm, mu = par.gn.norm$par[1], alpha = par.gn.norm$par[2],
                beta = par.gn.norm$par[3], log = T)) + 2 * length(par.gn.norm$par)

par(mfrow=c(1,3))
plot(ecdf(D$SLV), verticals = T)
xseq <- seq(0.9*min(D$SLV), 1.1*max(D$SLV), length.out=100)
#lines(xseq, pnorm(xseq, mean(D$SLV), sd(D$SLV)), col='red')
lines(xseq, pnorm(xseq, mean = par$par[1], sd = par$par[2]), col='blue')
plot(ecdf(D$SLV), verticals = T)
lines(xseq, pgnorm(xseq, mu = par.gn$par[1], alpha = par.gn$par[2], beta = par.gn$par[3]), col='green')
plot(ecdf(D$SLV), verticals = T)
lines(xseq, psn(xseq, xi = par.sn$par[1], omega = par.sn$par[2], alpha = par.sn$par[3]), col='red')

