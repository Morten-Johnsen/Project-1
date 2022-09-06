rm(list = ls())
library(lubridate)
library(tidyverse)
library(reshape)
library(stringr)

#Retrieve data from the descriptive statistics script
##### TILPAS DEN HER TIL DEN COMPUTER DEER BRUGES
if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  setwd("OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/Statistical Modelling/Project-1/")
} else {
  print("Fejl: Viktor, Indsæt lokationen på din egen folder her")
}

source("descriptiveStatistics.R")

testDistribution <- function(p, x, distribution = "Normal"){
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
  return(NLL)
}


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

