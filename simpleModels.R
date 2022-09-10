rm(list = ls())
library(lubridate)
library(tidyverse)
library(reshape)
library(stringr)
library(pheatmap)

#Retrieve data from the descriptive statistics script
##### TILPAS DEN HER TIL DEN COMPUTER DEER BRUGES
if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/Statistical Modelling/Project-1/")
} else {
  print("Fejl: Viktor, Indsæt lokationen på din egen folder her")
}

source("descriptiveStatistics.R")

testDistribution <- function(p, x, distribution = "Normal", giveDistributions = F){
  if (giveDistributions == T){
    NLL <- c("normal", "gamma", "beta", "negative binomial")
    distribution = "none"
  }
  if (str_to_lower(distribution) == "normal"){
    mu <- p[1]
    sigma <- p[2]
    NLL <- -sum(dnorm(x, mean = mu, sd = sigma, log = T))
  }
  if (str_to_lower(distribution) == "gamma"){
    shape <- p[1]
    rate <- p[2]
    NLL <- -sum(dgamma(x, shape = shape, rate = rate, log = T))
  }
  if (str_to_lower(distribution) == "beta"){
    #The beta distribution only ranges from [0;1] and thus it is
    #exclusively relevant to the normalized wind power statistic
    #and not the fitting of the other two parameters.
    shape1 <- p[1]
    shape2 <- p[2]
    NLL <- -sum(dbeta(x, shape = shape1, shape2 = shape2, log = T))
  }
  if (str_to_lower(distribution) == "negative binomial"){
    alpha <- p[1] #target number of succesfull trials
    pi <- p[2] #probability of succes in each trial
    NLL <- -sum(dnbinom(alpha, size = x, prob = pi, log = T))
  }
  return(NLL)
}

#### WINDPOWER ####
#Wind power have been normalized to lie between 0-1, so we use the second transformation in the
#assigment to see if we can get a normal distribution:
lambda <- 0.5
D$pow.obs.norm.trans <- 2*log(D$pow.obs.norm^lambda/(1-D$pow.obs.norm)^(1-lambda))
qqplot(rnorm(length(D$pow.obs.norm.trans), mean = mean(D$pow.obs.norm.trans), sd = sd(D$pow.obs.norm.trans)), D$pow.obs.norm.trans)

#Fitting models to wind power:
dists <- testDistribution(c(NA,NA), NA, giveDistributions = T)
opt.parms <- matrix(NA, nrow = length(dists), ncol = 2)

for (i in 1:length(dists)){
  opt.temp <- nlminb(start = c(1,0.5), objective = testDistribution,
         x = D$pow.obs.norm,
         distribution = dists[i],
         lower = c(0,0))
  opt.parms[i,] <- opt.temp$par
}

ggplot(D)+
  geom_histogram(aes(x = pow.obs.norm)
                 , bins = 30)+
  theme_bw()+
  stat_function(fun = dnorm,   n = length(D$pow.obs.norm), args = list(mean = opt.parms[1,1], sd = opt.parms[1,2]))+
  stat_function(fun = dgamma,  n = length(D$pow.obs.norm), args = list(shape = opt.parms[2,1], rate = opt.parms[2,2]))+
  stat_function(fun = dbeta,   n = length(D$pow.obs.norm), args = list(shape1 = opt.parms[3,1], shape2 = opt.parms[3,2]))+
  stat_function(fun = dnbinom, n = length(D$pow.obs.norm), args = list(size = opt.parms[4,1], prob = opt.parms[4,2]))
  

g.pow.obs
#### WIND SPEED ####
ws30_N <- nlminb(start = c(0,1), objective = testDistribution
            , x = D$ws30
            , lower = c(0,0))

ws30_G <- nlminb(start = c(1,1), objective = testDistribution
                 , x = D$ws30
                 , distribution = "gamma"
                 , lower = c(0,0))


ggplot(D)+
  geom_histogram(aes(x = ws30, y = ..count../sum(..count..))
                 , colour = "white"
                 , bins = 30)+
  stat_function(fun = dnorm, n = 101, args = list(mean = ws30_N$par[1], sd = ws30_N$par[2]))+
  stat_function(fun = dgamma
                , n = 101
                , args = list(shape = ws30_G$par[1], rate = ws30_G$par[2])
                , colour = "blue")

