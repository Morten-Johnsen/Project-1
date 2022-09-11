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
  setwd("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1")
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
  if (str_to_lower(distribution) == "exponential"){
    lambda = p
    NLL <- -sum(dexp(x, rate = lambda, log = T))
  }
  
  if (str_to_lower(distribution) == "weibull"){
    shape = p[1]
    scale = p[2]
    NLL <- -sum(dweibull(x, shape = shape, scale = scale, log = T))
  }
  if (str_to_lower(distribution) == "negative binomial"){
    alpha <- p[1] #target number of succesfull trials
    probs <- p[2] #probability of succes in each trial
    NLL <- -sum(dnbinom(x = alpha, size = x, prob = probs, log = T))
  }
  return(NLL)
}

#### WINDPOWER ####
#Fitting models to wind power:
par <- nlminb(start = 0.2, objective = testDistribution,
              distribution = "exponential",
              x = D$pow.obs.norm)

ggplot(D)+
  geom_histogram(aes(x = pow.obs.norm, y = ..density..), bins = 20)+
  theme_bw()+
  stat_function(fun = dexp, n = length(D$pow.obs.norm), args = list(rate = par$par))

g.pow.obs
#### WIND SPEED ####
par.ws30 <- nlminb(start = c(1,1), objective = testDistribution
            , x = D$ws30
            , distribution = "weibull"
            , lower = c(0,0))


ggplot(D)+
  geom_histogram(aes(x = ws30, y = ..count../sum(..count..))
                 , colour = "white"
                 , bins = 30)+
  theme_bw()+
  stat_function(fun = dweibull, n = dim(D)[1], args = list(shape = par.ws30$par[1], scale = par.ws30$par[2]))

#### WIND DIRECTION ####
#centrerer fordelingen omkring 3/2pi
D$wd30.centered <- D$wd30 - pi/2; D$wd30.centered[D$wd30.centered < 0] = D$wd30.centered[D$wd30.centered < 0] + 2*pi
par.wd30 <- nlminb(start = c(4,4),
                   objective = testDistribution,
                   x = D$wd30.centered,
                   distribution = "normal")

ggplot(D)+
  theme_bw()+
  geom_density(aes(x = wd30.centered, y = ..density..), alpha = .8, colour = "white", fill = "red", colour = "white")+
  geom_density(aes(x = wd30, y = ..density..), colour = "white", alpha = .2, fill = "blue")+
  scale_x_continuous(breaks = c(0,pi/2,pi,3/2*pi,2*pi)
                     , labels =c("0", "pi/2", "pi", "3/2pi", "2pi"))+
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = par.wd30$par[1], sd = par.wd30$par[2]))
