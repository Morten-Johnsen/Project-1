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
library(extraDistr)

if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/")
  D <- read.table("finance_data.csv", header=TRUE, sep=";", 
                  as.is=TRUE)
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/Assignment 3")
} else {
  setwd("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1")
  D <- read.table("finance_data.csv", header=TRUE, sep=";", 
                  as.is=TRUE)
  setwd("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1/Assignment 3")
}
source("A1.R")
source("A2.R")

D$date <- as.Date(D$time)

llstFUNC <- function(p, data){ #location-scale t-distribution
  return(-sum(dlst(x = data, df = p[1], mu = p[2], sigma = p[3], log = T)))
}
par.lst <- nlminb(start = c(1,1,1), objective = llstFUNC, data = D$SLV)

mle.lst <- par.lst$par

temp1 <- paste("df == ", round(mle.lst[1], digits=2))
temp2 <- paste("mu == ", round(mle.lst[2], digits=5))
temp3 <- paste("sigma == ", round(mle.lst[3], digits=4))

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(fun = dlst, n = dim(D)[1], args = list(df = par.lst$par[1], mu = par.lst$par[2],
                                                       sigma = par.lst$par[3]), color = 'red') +
  annotate( "text", x = 2.7/5*max(D$SLV), y = c(10.5, 10.0, 9.5), label = c(temp1,temp2,temp3), parse = T  ) +
  ggtitle("Location-scale t-distribution and distribution of the weekly returns")
####################################################################################
####Task 1 of 2
###Subtask a of e
#These are the functions that we have to modify: so that they apply for N instead of pois:
#The first one we don't really need.
# pois.mix.pn2pw <- function(m,lambda,delta){
#   if(sum(delta) >= 1){return("sum(delta) should be < 1")}
#   if(length(lambda) != m){return("length(lambda) should be m")}
#   if(length(delta) != (m-1)){return("length(delta) should be m-1")}
#   eta <- log(lambda)
#   tau <- log(delta/(1-sum(delta)))
#   return(list(eta=eta,tau=tau))
# }
# pois.mix.pw2pn <- function(m,eta,tau){
#   if(m==1){return(exp(eta))}
#   if(length(eta) != m){return("length(lambda) should be m")}
#   if(length(tau) != (m-1)){return("length(delta) should be m-1")}
#   lambda <- exp(eta) #they do this because lambda > 0
#   delta <- exp(tau)/(1+sum(exp(tau))) #they do this because the delta's are constrained in that they sum to 1
#   delta <- c(1-sum(delta),delta) #this is to get the last delta. Remember, only m-1 delta's because they sum to 1.
#   return(list(lambda=lambda,delta=delta))
# }

# nll <- function(theta,m,y){
#   if(m==1){
#     return(-sum(dpois(y,lambda=exp(theta),log=TRUE)))
#   }
#   eta <- theta[1:m]
#   tau <- theta[(m+1):(2*m-1)]
#   n.pars <- pois.mix.pw2pn(m,eta,tau)
#   n <- length(y)
#   nll <- 0
#   for(i in 1:n){
#     nll <- nll - log(sum(n.pars$delta * dpois(y[i],lambda=n.pars$lambda)))
#   }
#   return(nll)
# }

pois.mix.pn2pw <- function(m,lambda,delta){
  if(sum(delta) >= 1){return("sum(delta) should be < 1")}
  if(length(lambda) != m){return("length(lambda) should be m")}
  if(length(delta) != (m-1)){return("length(delta) should be m-1")}
  eta <- log(lambda)
  tau <- log(delta/(1-sum(delta)))
  return(list(eta=eta,tau=tau))
}
N.mix.pn2pw <- function(m,mu,sd,delta){
  if(sum(delta) >= 1){return("sum(delta) should be < 1")}
  if(length(lambda) != m){return("length(lambda) should be m")}
  if(length(sd) != m){return("length(sd) should be m")}
  if(length(delta) != (m-1)){return("length(delta) should be m-1")}
}

N.mix.pw2pn <- function(m,theta, eta, tau){#theta = mu, eta = log(sigma), tau = log(delta / (1-sumdelta_other))
  if(m==1){return(theta, exp(eta))}
  if(length(theta) != m){return("length(mu) should be m")}
  if(length(eta) != m){return("length(sigma) should be m")}
  if(length(tau) != (m-1)){return("length(delta) should be m-1")}
  mu <- theta
  sigma <- exp(eta)
  delta <- exp(tau)/(1+sum(exp(tau)))
  delta <- c(1-sum(delta),delta)
  return(list(mu=mu,sigma=sigma,delta=delta))
}

N.nll <- function(theta,m,y){
  if(m==1){return( -sum(dnorm(y, mean=theta[1], sd=theta[2], log=T)))
  }
  theta_mu <- theta[1:m]
  eta <- theta[(m+1):2*m]
  tau <- theta[(2*m+1):(3*m-1)]
  n.pars <- N.mix.pw2pn(m,theta_mu,eta,tau)
  n <- length(y)
  nll <- 0
  for(i in 1:n){
    nll <- nll - log(sum(n.pars$delta * dnorm(y[i], mean=n.pars$mu, sd=n.pars$sigma)))
  }
  return(nll)
}
m <- 2 #2 components
mu <- c(1/2,3/2)*mean(D$SLV) #init vals for the optim
sd <- c(1/2,3/2)*sd(D$SLV) #init vals for the optim
delta <- c(0.5) #fifty-fifty for the two dists
wpars <- N.mix.pn2pw(m=m, mu=mu, sd=sd, delta=delta)
theta <- c(wpars$theta, wpars$eta, wpars$tau)
opt2 <- nlminb(theta=theta, N.nll, m=m, y=D$SLV)




