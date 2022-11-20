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
library(mvtnorm)
#library(openair)
#library(PowerNormal)
#library(sn)
#library(gnorm)
#library(emg)
#library(survival)
#library(survminer)

if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/")
  D <- read.table("tuno.txt", header=TRUE, sep=" ", 
                  as.is=TRUE)
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/Assignment 3")
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

ws.and.ws.squared <- lm( (D$transformed.pow.obs.norm ~ ws30 + I(ws30^2)), data = D)
D$ws.ws2.pred <- ws.and.ws.squared$coefficients[1] + ws.and.ws.squared$coefficients[2] * D$ws30 + ws.and.ws.squared$coefficients[3] * D$ws30^2

ggplot(data = D)+
  geom_point(aes(x=ws30, y=transformed.pow.obs.norm, colour="Data"), size = 4, alpha = 0.6)+
  geom_line(aes(x=ws30, y=ws.ws2.pred, colour="ws + ws^2"))+
  scale_colour_manual(values = c("red", "blue"))+
  #geom_line(aes(x=ws30, y=ws.ws2.wd.pred, colour="ws + ws^2 + cos(wd)"))+
  labs(x = "Wind speeds", y="Transformed, norm power obs", colour = "")+
  theme_bw()+
  ggtitle("Fitted model in Transformed Domain")

y_p <- D$ws.ws2.pred
D$y_inv_trans <- 1/( exp(y_p*lambda)+1 )^(1/lambda) * exp(y_p)

ggplot(data = D)+
  geom_point(aes(x=ws30, y=pow.obs.norm, colour = "Observed power production (Transformed)"), size = 4, alpha = 0.6)+
  geom_line(aes(x=ws30, y=y_inv_trans, colour="Inverse Transformation of the Normal model"), size = 2)+
  #geom_line(aes(x=ws30, y=betaModel.pow, colour="The directly fitted beta model"), size = 2)+
  scale_colour_manual(values = c("blue","red"))+
  labs(x = "Wind speeds", y="Norm power obs", colour ="")+
  theme_bw()+
  labs(colour = "")+
  ggtitle("Model fits")+
  theme(legend.background = element_rect(colour = "black"), legend.position = "top")

################################################################################################
###Analysis of auto-correlation
##TASK 1 of 7...
#Residuals are:
e1 <- ws.and.ws.squared$residuals[1:( length(ws.and.ws.squared$residuals)-1 )] #for e_1 - e_n-1
e2 <- ws.and.ws.squared$residuals[2:length(ws.and.ws.squared$residuals)] #for e_2 - e_n

m_res <- matrix(c(e1,e2), ncol=2, byrow=F)
#

acf(D$pow.obs.norm)
acf(D$transformed.pow.obs.norm)

AR1m <- arima(D$transformed.pow.obs.norm, order=c(1,0,0))
e <- D$transformed.pow.obs.norm[-1] - AR1m$coef[1] * D$transformed.pow.obs.norm[-dim(D)[1]]
acf(e)
acf(AR1m$residuals)
( e - AR1m$residuals[-1])





