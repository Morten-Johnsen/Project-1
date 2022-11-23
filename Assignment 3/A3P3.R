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

#####
N.mix.pn2pw <- function(n.dist,mu,sigma,delta){ #n.dist is m
  if(sum(delta) >= 1){return("sum(delta) should be < 1")}
  if(length(mu) != n.dist){return("length(mu) should be m")}
  if(length(sigma) != n.dist){return("length(sd) should be m")}
  if(length(delta) != (n.dist-1)){return("length(delta) should be m-1")}
  theta <- mu
  eta <- log(sigma)
  tau = log(delta / (1-sum(delta)))
  return(list(theta=theta, eta=eta, tau=tau))
}

N.mix.pw2pn <- function(n.dist, theta, eta, tau){#theta = mu, eta = log(sigma), tau = log(delta / (1-sumdelta_other))
  if(n.dist==1){return(list(theta=theta, sigma=exp(eta)))}
  if(length(theta) != n.dist){return("length(mu) should be m")}
  if(length(eta) != n.dist){return("length(sigma) should be m")}
  if(length(tau) != (n.dist-1)){return("length(delta) should be m-1")}
  mu <- theta
  sigma <- exp(eta)
  delta <- exp(tau)/(1+sum(exp(tau)))
  delta <- c(1-sum(delta),delta)
  return(list(mu=mu,sigma=sigma,delta=delta))
}

N.nll <- function(p,n.dist,y){
  if(n.dist==1){return( -sum(dnorm(y, mean=p[1], sd=exp(p[2]), log=T)))
  }
  theta <- p[1:n.dist]
  eta <- p[(n.dist+1):(2*n.dist)]#important to remember these parentheses
  tau <- p[(2*n.dist+1):(3*n.dist-1)]
  n.pars <- N.mix.pw2pn(n.dist=n.dist, theta=theta, eta=eta, tau=tau)
  n <- length(y)
  nll <- 0
  for(i in 1:n){
    nll <- nll - log(sum(n.pars$delta * dnorm(y[i], mean=n.pars$mu, sd=n.pars$sigma)))
  }
  return(nll)
}
##################################################
m <- 1 # 1 components
mu <- mean(D$SLV)
sigma <- sd(D$SLV)
delta <- c()
wpars <- N.mix.pn2pw(n.dist=m,mu=mu,sigma=sigma,delta=delta)
N.nll(p=c(wpars$theta,wpars$eta),n.dist=m, y=D$SLV) #-sum(dnorm(D$SLV,mean=mu,sd=sigma,log=T)), check: so far so good
p <- c(wpars$theta, wpars$eta)
opt1 <- nlminb(start=p, N.nll, n.dist=m, y=D$SLV);opt1$objective #it gives the same as above and for the check: nice
(pars1 <- N.mix.pw2pn(n.dist=m,opt1$par[1:m],opt1$par[(m+1):(2*m)]))
ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = pars1$theta,
                                                       sd = pars1$sigma), color = 'red') +
  ggtitle("Normal distribution fitted to SLV")
##################################################
m <- 2 #2 components
mu <- c(1/2,3/2)*mean(D$SLV) #init vals for the optim
sigma <- c(1/2,3/2)*sd(D$SLV) #init vals for the optim
delta <- c(0.5) #fifty-fifty for the two dists
wpars <- N.mix.pn2pw(n.dist=m, mu=mu, sigma=sigma, delta=delta)
( N.mix.pw2pn(n.dist=m, theta=wpars$theta, eta=wpars$eta, tau=wpars$tau) ) #checking that pn2pw -> pw2pn works
p <- c(wpars$theta, wpars$eta, wpars$tau)
N.nll(p=p, n.dist=m, y=D$SLV)
opt2 <- nlminb(start=p, N.nll, n.dist=m, y=D$SLV)
(pars2 <- N.mix.pw2pn(n.dist=m,opt2$par[1:m],opt2$par[(m+1):(2*m)],opt2$par[(2*m+1):(3*m-1)]))
ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = pars2$mu[1],
                                                        sd = pars2$sigma[1]), color = 'red') +
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = pars2$mu[2],
                                                        sd = pars2$sigma[2]), color = 'orange') +
  ggtitle("Normal distribution fitted to SLV")
#with delta = 0.91, 0.09 for red, orange
##################################################
m <- 3 #3 components
mu <- c(1/2,1,3/2)*mean(D$SLV) #init vals for the optim
sigma <- c(1/2,1,3/2)*sd(D$SLV) #init vals for the optim
delta <- c(0.34,0.33) #fifty-fifty for the two dists
wpars <- N.mix.pn2pw(n.dist=m, mu=mu, sigma=sigma, delta=delta)
( N.mix.pw2pn(n.dist=m, theta=wpars$theta, eta=wpars$eta, tau=wpars$tau) ) #checking that pn2pw -> pw2pn works
p <- c(wpars$theta, wpars$eta, wpars$tau)
N.nll(p=p, n.dist=m, y=D$SLV)
opt3 <- nlminb(start=p, N.nll, n.dist=m, y=D$SLV)
(pars3 <- N.mix.pw2pn(n.dist=m,opt3$par[1:m],opt3$par[(m+1):(2*m)],opt3$par[(2*m+1):(3*m-1)]))
ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = pars3$mu[1],
                                                        sd = pars3$sigma[1]), color = 'red') +
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = pars3$mu[2],
                                                        sd = pars3$sigma[2]), color = 'orange') +
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = pars3$mu[3],
                                                        sd = pars3$sigma[3]), color = 'blue') +
  ggtitle("Normal distribution fitted to SLV")
#with delta = 0.83, 0.06, 0.10, for red, orange, blue
##################################################
#comparing the nlls of the lts and the three normal dists:
round(rbind(par.lst$objective,opt1$objective,opt2$objective,opt3$objective),digits=3) #NLLs
round(rbind(2*par.lst$objective+2*length(par.lst$par), #AICs
            2*opt1$objective+2*length(opt1$par),
            2*opt2$objective+2*length(opt2$par),
            2*opt3$objective+2*length(opt3$par)),digits=3)
round(rbind(2*par.lst$objective + length(par.lst$par)*log(length(D$SLV)), #BICs
            2*opt1$objective + length(opt1$par)*log(length(D$SLV)),
            2*opt2$objective + length(opt2$par)*log(length(D$SLV)),
            2*opt3$objective + length(opt3$par)*log(length(D$SLV))),digits=3)
#Den med 3 komponenter bliver straffet for hårdt pga. de mange parametre. Den med 2 kan næsten hamle op med lst.
####################################################################################################
###Subtask b of e
##CIs of the parameters of the best mixture model. This is the two-modal model.
opt2$par #we need 5 CIs
m <- 2
V.w2 <- diag(solve( hessian(func=N.nll, x=opt2$par, n.dist=m, y=D$SLV) ))
CIs.w2 <- matrix(c(opt2$par-sqrt(V.w2), opt2$par, opt2$par+sqrt(V.w2)), nrow=3, byrow=T) #1st: lower, 2nd: MLE, 3rd: upper
CIs.n2 <- matrix(rep(0,( dim(CIs.w2)[1] * ( dim(CIs.w2)[2]+1 ) ) ), nrow=3, byrow=T)
for (i in 1:3){
  n2_temp <- N.mix.pw2pn(n.dist=m, theta=CIs.w2[i,(1:m)], eta=CIs.w2[i,(m+1):(2*m)], tau=CIs.w2[i,(2*m+1):(3*m-1)])
  CIs.n2[i,1:dim(CIs.n2)[2]] <- c(n2_temp$mu, n2_temp$sigma, n2_temp$delta) #1st: lower, 2nd: MLE, 3rd: upper
}
####################################################################################################
###Subtask c of e
##PL of one of the variance parameters of the two component model. Driller lidt :/
#Modifying function N.nll from earlier to accept single value for on of the variances whilst optimising the others:
N.pl <- function(sd0, n.dist, y){
  #if(n.dist==1){return( -sum(dnorm(y, mean=p[1], sd=exp(p[2]), log=T)))
  #}
  nll_temp <- function(p, n.dist, y, sd0){
    theta <- p[1:n.dist]
    p[(n.dist+1)] <- sd0
    eta <- p[(n.dist+1):(2*n.dist)]#important to remember these parentheses
    tau <- p[(2*n.dist+1):(3*n.dist-1)]
    n.pars <- N.mix.pw2pn(n.dist=n.dist, theta=theta, eta=eta, tau=tau)
    n <- length(y)
    nll <- 0
    for(i in 1:n){
      nll <- nll - log(sum(n.pars$delta * dnorm(y[i], mean=n.pars$mu, sd=n.pars$sigma)))
    }
    return(nll)
  }
  m <- 2
  mu <- c(1/2,3/2)*mean(D$SLV) #init vals for the optim
  sigma <- c(1/2,3/2)*sd(D$SLV) #init vals for the optim
  delta <- c(0.5) #fifty-fifty for the two dists
  wpars <- N.mix.pn2pw(n.dist=m, mu=mu, sigma=sigma, delta=delta)
  opt2_temp <- nlminb(start=c(wpars$theta, wpars$eta, wpars$tau), N.nll, n.dist=m, y=D$SLV)
  n.pars2_temp <- N.mix.pw2pn(n.dist=m,opt2_temp$par[1:m],opt2_temp$par[(m+1):(2*m)],opt2_temp$par[(2*m+1):(3*m-1)])
  n <- length(y)
  ll <- 0
  for (i in 1:n){
    ll <- ll + log(sum(n.pars2_temp$delta * dnorm(y[i], mean=n.pars2_temp$mu, sd=n.pars2_temp$sigma)))
  }
  return(ll)
}
sd2.n <- seq(max(opt2$par[m+1]-2*sqrt(V.w2)[m+1],-1),
            min(opt2$par[m+1]+2*sqrt(V.w2)[m+1],1),
            length=10)
N.llp <- sapply(X = sd2.n, FUN=N.pl, n.dist=m, y=D$SLV)

plot(sd2.n,exp(N.llp-max(N.llp)),type="l")
lines(range(V2.n),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)


