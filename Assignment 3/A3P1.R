rm(list = ls())
library(lubridate)
library(latex2exp)
library(circular)
library(tidyverse)
library(reshape)
library(pheatmap)
library(openair)
library(stringr)
library(numDeriv)
library(gridExtra)
#library(openair)
#library(PowerNormal)
#library(sn)
#library(gnorm)
#library(emg)
#library(survival)
#library(survminer)

if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/")
} else {
  setwd("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1")
  D <- read.table("tuno.txt", header=TRUE, sep=" ", 
                  as.is=TRUE)
  setwd("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1/Assignment 3")
}
source("A1.R")
source("A2.R")

D$date <- as.Date("2003-01-01")-1+D$r.day
D$pow.obs.norm <- D$pow.obs/5000

  #1.
Trans.eq1 <- function(lambda, y){ #se tidl. projekter
  y_lambda <- 1/lambda * log(y^lambda/(1-y^lambda))#, lambda > 0
  return(y_lambda)
}
lambda_NLL <- function(lambda, x = D$pow.obs.norm){
  y <- Trans.eq1(lambda, x)
  NLL <- -as.numeric(shapiro.test(y)$statistic)
  return(NLL)
}
lambda.hat <- nlminb(start = 0.2, objective = lambda_NLL)
lambda <- round(lambda.hat$par, 2)
D$transformed.pow.obs.norm <- Trans.eq1(lambda, D$pow.obs.norm)

acf(D$pow.obs.norm)
acf(D$transformed.pow.obs.norm)

AR1m <- arima(D$transformed.pow.obs.norm, order=c(1,0,0))
e <- D$transformed.pow.obs.norm[-1] - AR1m$coef[1] * D$transformed.pow.obs.norm[-dim(D)[1]]
acf(e)
acf(AR1m$residuals)
( e - AR1m$residuals[-1])





