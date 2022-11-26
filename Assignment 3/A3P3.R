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
alpha <- 0.05
V.w2 <- diag(solve( hessian(func=N.nll, x=opt2$par, n.dist=m, y=D$SLV) ))
CIs.w2 <- matrix(c(opt2$par-qnorm(1-alpha/2)*sqrt(V.w2), opt2$par, opt2$par+qnorm(1-alpha/2)*sqrt(V.w2)), nrow=3, byrow=T) #1st: lower, 2nd: MLE, 3rd: upper
CIs.n2 <- matrix(rep(0,( dim(CIs.w2)[1] * ( dim(CIs.w2)[2]+1 ) ) ), nrow=3, byrow=T)
for (i in 1:3){
  n2_temp <- N.mix.pw2pn(n.dist=m, theta=CIs.w2[i,(1:m)], eta=CIs.w2[i,(m+1):(2*m)], tau=CIs.w2[i,(2*m+1):(3*m-1)])
  CIs.n2[i,1:dim(CIs.n2)[2]] <- c(n2_temp$mu, n2_temp$sigma, n2_temp$delta) #1st: lower, 2nd: MLE, 3rd: upper
}
####################################################################################################
###Subtask c of e
##PL of one of the variance parameters of the two component model. Driller lidt :/
#Modifying function N.nll from earlier to accept single value for on of the variances whilst optimising the others:
N.pl <- function(eta0, n.dist, y, first=T, robust=F){ #to be used for PL of the 2nd variance parameter (working) of the 2 components normal mixture
  #if(n.dist==1){return( -sum(dnorm(y, mean=p[1], sd=exp(p[2]), log=T)))
  #}
  nll_temp <- function(p, n.dist=n.dist, y=y, eta0=eta0, bool=first){ #only works m = 2
    theta <- p[1:n.dist]
    eta <- p[(n.dist+1):(2*n.dist)]#important to remember these parentheses
    if(bool){
      eta[1] <- eta0
    } else {
      eta[n.dist] <- eta0
    }
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
  mu <- c(1/2,3/2)*mean(y) #init vals for the optim
  sigma <- c(1/4,8)*sd(y) #init vals for the optim
  delta <- c(0.5) #fifty-fifty for the two dists
  wpars <- N.mix.pn2pw(n.dist=m, mu=mu, sigma=sigma, delta=delta)
  opt2_temp <- nlminb(start=c(wpars$theta, wpars$eta, wpars$tau), objective=nll_temp, n.dist=m, y=y, eta0=eta0)
  #n.pars2_temp <- N.mix.pw2pn(n.dist=m,opt2_temp$par[1:m],opt2_temp$par[(m+1):(2*m)],opt2_temp$par[(2*m+1):(3*m-1)])
  if(first){
    n.pars2_temp <- N.mix.pw2pn(n.dist=m,opt2_temp$par[1:m],c(eta0, opt2_temp$par[(2*m)]),opt2_temp$par[(2*m+1):(3*m-1)])
  } else{
    n.pars2_temp <- N.mix.pw2pn(n.dist=m,opt2_temp$par[1:m],c(opt2_temp$par[(m+1)], eta0),opt2_temp$par[(2*m+1):(3*m-1)])
  }
  if(robust){
    n.pars2_temp$sigma = sort(n.pars2_temp$sigma)
  }
  n <- length(y)
  ll <- 0
  for (i in 1:n){
    ll <- ll + log(sum(n.pars2_temp$delta * dnorm(y[i], mean=n.pars2_temp$mu, sd=n.pars2_temp$sigma)))
  }
  return(ll)
}

eta1.w <- seq(opt2$par[m+1]-6*sqrt(V.w2)[m+1],
              opt2$par[m+1]+30*sqrt(V.w2)[m+1],
              length=100)
eta2.w <- seq(opt2$par[2*m]-6*sqrt(V.w2)[2*m],
             opt2$par[2*m]+6*sqrt(V.w2)[2*m],
            length=100)
m <- 2
N.llp <- sapply(X = eta1.w, FUN=N.pl, n.dist=m, y=D$SLV)
N.llp2 <- sapply(X = eta2.w, FUN=N.pl, n.dist=m, y=D$SLV, first=F)

plot(eta1.w,exp(N.llp-(max(N.llp))),type="l",
     main = "Profile likelihood of FIRST variance param of 2 component norm mixture dist") #smart Jan-plot
lines(range(eta1.w),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)
abline(v=opt2$par[m+1], col="red") #first MLE var
abline(v=opt2$par[2*m], col="blue") # second MLE var

plot(eta2.w,exp(N.llp2-(max(N.llp2))),type="l",
     main = "Profile likelihood of SECOND variance param of 2 component norm mixture dist") #smart Jan-plot
lines(range(eta2.w),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)
abline(v=opt2$par[2*m], col="blue")
abline(v=opt2$par[m+1], col="red")

# plot(eta2.w,N.llp, type="l") #grimt plot.
# abline(v = opt2$par[2*m])
# abline(h=-opt2$objective)

####################################################################################################
###Task d of e
N.llp.rep <- sapply(X = eta1.w, FUN=N.pl, n.dist=m, y=D$SLV, robust=T)
N.llp2.rep <- sapply(X = eta2.w, FUN=N.pl, n.dist=m, y=D$SLV, first=F, robust=T)

plot(eta1.w,exp(N.llp.rep-(max(N.llp.rep))),type="l",
     main = "Profile likelihood of FIRST variance param of 2 component norm mixture dist") #smart Jan-plot
lines(range(eta1.w),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)
abline(v=opt2$par[m+1], col="red") #first MLE var
abline(v=opt2$par[2*m], col="blue") # second MLE var

plot(eta2.w,exp(N.llp2.rep-(max(N.llp2.rep))),type="l",
     main = "Profile likelihood of SECOND variance param of 2 component norm mixture dist") #smart Jan-plot
lines(range(eta2.w),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)
abline(v=opt2$par[2*m], col="blue")
abline(v=opt2$par[m+1], col="red")
####################################################################################################
####################################################################################################
####Task 2 of 2
###Subtask a of e
pois.HMM.pn2pw <- function(m,lambda,gamma)
{                                              
  tlambda <- log(lambda)                         
  tgamma  <- NULL                              
  if(m>1)                                        
  {                                            
    foo   <- log(gamma/diag(gamma))           
    tgamma<- as.vector(foo[!diag(m)])             
  }                                             
  parvect <- c(tlambda,tgamma)                    
  parvect                                         
}  

m <- 2
lambda <- c(1/2,3/2) * 1 / mean(D$SLV)
gamma <- matrix(c(0.95,0.05,0.05,0.95), ncol=2,byrow=T)
wpars2 <- pois.HMM.pn2pw(m=m, lambda=lambda, gamma=gamma) #seeing how it works

pois.HMM.pw2pn <- function(m,parvect)                 
{                                                     
  epar   <- exp(parvect)                              
  lambda <- epar[1:m]                                   
  gamma  <- diag(m)                                    
  if(m>1)                                               
  {                                                  
    gamma[!gamma] <- epar[(m+1):(m*m)]                  
    gamma         <- gamma/apply(gamma,1,sum)          
  }                                                   
  delta  <- solve(t(diag(m)-gamma+1),rep(1,m))          
  list(lambda=lambda,gamma=gamma,delta=delta)           
}  

npars2 <- pois.HMM.pw2pn(m=m, parvect=wpars2)

pois.HMM.mllk <- function(parvect,x,m,...)       
{
  #    print(parvect)
  if(m==1) return(-sum(dpois(x,exp(parvect),log=TRUE))) 
  n          <- length(x)                            
  pn         <- pois.HMM.pw2pn(m,parvect)            
  allprobs   <- outer(x,pn$lambda,dpois)             
  allprobs   <- ifelse(!is.na(allprobs),allprobs,1)  
  lscale     <- 0                                    
  foo        <- pn$delta                             
  for (i in 1:n)                                    
  {                                                
    foo    <- foo%*%pn$gamma*allprobs[i,]            
    sumfoo <- sum(foo)                               
    lscale <- lscale+log(sumfoo)                    
    foo    <- foo/sumfoo                            
  }                                               
  mllk       <- -lscale                            
  mllk                                              
}    

mllk2 <- pois.HMM.mllk(parvect=wpars2,x=D$SLV,m=m)

pois.HMM.mle.nlminb <- function(x,m,lambda0,gamma0,...)
{                                                      
  parvect0 <- pois.HMM.pn2pw(m,lambda0,gamma0)
  np        <- length(parvect0)                          
  lower    <- rep(-10,np)
  upper    <- c(rep(max(x),m),rep(10,np-m))
  mod      <- nlminb(parvect0,pois.HMM.mllk,x=x,m=m,
                     lower=lower,upper=upper)
  if(mod$convergence!=0){
    print(mod)
  }
  pn       <- pois.HMM.pw2pn(m,mod$par)
  mllk     <- mod$objective
  AIC       <- 2*(mllk+np)
  n         <- sum(!is.na(x))
  BIC       <- 2*mllk+np*log(n)
  list(lambda=pn$lambda,gamma=pn$gamma,delta=pn$delta,   
       code=mod$convergence,mllk=mllk,AIC=AIC,BIC=BIC)   
}  

pois.HMM.mle.nlminb(x=D$SLV,m=m,lambda0=lambda,gamma0=gamma)


