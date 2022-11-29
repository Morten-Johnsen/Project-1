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
##PL of one of the variance parameters of the two component model.
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
              length=25)
eta2.w <- seq(opt2$par[2*m]-6*sqrt(V.w2)[2*m],
             opt2$par[2*m]+6*sqrt(V.w2)[2*m],
            length=25)
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
###Subtask d of e
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
###Subtask e of e


####################################################################################################
####################################################################################################
####Task 2 of 2 HMMs
###Subtask a of e
N.HMM.pn2pw <- function(m,mu,sigma,gamma)
{                                              
  theta <- mu
  eta <- log(sigma)
  tgamma  <- NULL                              
  if(m>1)                                        
  {                                            
    foo   <- log(gamma/diag(gamma))           
    tgamma<- as.vector(foo[!diag(m)])             
  }                                             
  parvect <- c(theta,eta,tgamma)                    
  parvect                                         
}  

#dat <- read.table("earthquakes.txt",header=FALSE)
m <- 2
mu <- c(-1,1) * mean(D$SLV) #arbitrary
sigma <- c(1/4,8) * sd(D$SLV) #arbitrary
gamma <- matrix(c(0.95,0.05,0.05,0.95), ncol=2,byrow=T) #arbitrary
wpars2.HMM <- N.HMM.pn2pw(m=m, mu=mu, sigma=sigma, gamma=gamma) #seeing how it works

N.HMM.pw2pn <- function(m,parvect)                 
{                                                     
  epar   <- exp(parvect)   
  mu <- parvect[1:m]
  sigma <- epar[(m+1):(2*m)]
  gamma  <- diag(m)                                    
  if(m>1)                                               
  {                                                  
    gamma[!gamma] <- epar[(2*m+1):(m*m+m)]                  
    gamma         <- gamma/apply(gamma,1,sum)          
  }                                                   
  delta  <- solve(t(diag(m)-gamma+1),rep(1,m)) #p. 19 in HMM, solving system of equations
  list(mu=mu,sigma=sigma,gamma=gamma,delta=delta)           
}  

npars2.HMM <- N.HMM.pw2pn(m=m, parvect=wpars2.HMM)

N.HMM.mllk <- function(parvect,x,m,...)       
{
  findNorm <- function(theta){
    mu<-theta[1]
    sd<-theta[2]
    dnorm(x,mean = mu, sd = sd)
    }
  #    print(parvect)
  if(m==1) return(-sum(dnorm(x,mean=parvect[1],sd=exp(parvect[2]),log=TRUE)))
  n          <- length(x)    
  pn         <- N.HMM.pw2pn(m,parvect)
  allprobs <- apply(cbind(pn$mu,pn$sigma), MARGIN = 1, FUN = findNorm)
  #allprobs   <- outer(x,pn$mu,pn$sigma,FUN=dnorm)
  allprobs   <- ifelse(!is.na(allprobs),allprobs,1)
  lscale     <- 0                                    
  foo        <- pn$delta                             
  for (i in 1:n)                                    
  {                                                
    foo    <- foo%*%pn$gamma*allprobs[i,] #pn$delta = pn$delta %*% pn$gamma for i = 1. Beregning er jf. slide 14, lecture 11
    sumfoo <- sum(foo)                               
    lscale <- lscale+log(sumfoo)                    
    foo    <- foo/sumfoo
  }                                               
  mllk       <- -lscale                            
  mllk                                              
}    

mllk2 <- N.HMM.mllk(parvect=wpars2.HMM,x=D$SLV,m=m)


N.HMM.mle.nlminb <- function(x,m,mu0,sigma0,gamma0,...)
{                                                      
  parvect0 <- N.HMM.pn2pw(m,mu0,sigma0,gamma0)
  np        <- length(parvect0)                          
  lower    <- rep(-10,np)
  upper    <- c(rep(max(x),m),rep(10,np-m))
  mod      <- nlminb(parvect0,N.HMM.mllk,x=x,m=m,
                     lower=lower,upper=upper)
  if(mod$convergence!=0){
    print(mod)
  }
  pn       <- N.HMM.pw2pn(m,mod$par)
  mllk     <- mod$objective
  AIC       <- 2*(mllk+np)
  n         <- sum(!is.na(x))
  BIC       <- 2*mllk+np*log(n)
  list(mu=pn$mu,sigma=pn$sigma,gamma=pn$gamma,delta=pn$delta,   
       code=mod$convergence,mllk=mllk,AIC=AIC,BIC=BIC)   
}
#following five lines are repeated from above
m <- 2
mu <- c(-2,2) * mean(D$SLV) #arbitrary
sigma <- c(1/2,4) * sd(D$SLV) #arbitrary
gamma <- matrix(c(0.95,0.05,0.05,0.95), ncol=2,byrow=T) #arbitrary
wpars2.HMM <- N.HMM.pn2pw(m=m, mu=mu, sigma=sigma, gamma=gamma) #seeing how it works
mle2.HMM <- N.HMM.mle.nlminb(x=D$SLV,m=m,mu0=mu,sigma0=sigma,gamma0=gamma)

#mllk1 <- N.HMM.mllk(parvect=c(mean(D$SLV),log(sd(D$SLV))),x=D$SLV,m=1)
#mle1.HMM <- N.HMM.mle.nlminb(x=D$SLV,m=1,mu0=mean(D$SLV),sigma0=sd(D$SLV),gamma0=c())

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = mle2.HMM$mu[1],
                                                        sd = mle2.HMM$sigma[1]), color = 'red') +
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = mle2.HMM$mu[2],
                                                        sd = mle2.HMM$sigma[2]), color = 'orange') +
  ggtitle("HMM with 2 normal distributions fitted to SLV")

m3 <- 3
mu3 <- c(0,0,0) * mean(D$SLV) #arbitrary -1,1/2,1
sigma3 <- c(1, 1, 1) * sd(D$SLV) #arbitrary 1/4,1,4
gamma3 <- matrix(rep(1/3,9), ncol=3,byrow=T) #arbitrary 0.95, 0.25, 0.25
#wpars3.HMM <- N.HMM.pn2pw(m=m0, mu=mu0, sigma=sigma0, gamma=gamma0)
#mllk3 <- N.HMM.mllk(parvect=wpars3.HMM,x=D$SLV,m=m0)

mle3.HMM <- N.HMM.mle.nlminb(x=D$SLV,m=m3,mu0=mu3,sigma0=sigma3,gamma0=gamma3)

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = mle3.HMM$mu[1],
                                                        sd = mle3.HMM$sigma[1]), color = 'red') +
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = mle3.HMM$mu[2],
                                                        sd = mle3.HMM$sigma[2]), color = 'orange') +
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = mle3.HMM$mu[3],
                                                        sd = mle3.HMM$sigma[3]), color = 'blue') +
  #xlim(c(-0.985,-0.90))+
  ggtitle("HMM with 3 normal distributions fitted to SLV")

mle2.HMM$mllk;mle3.HMM$mllk #2 comp normal HMM has 6 params, 3 comp normal HMM has 12... #params=m^2+m
mle2.HMM$AIC;mle3.HMM$AIC
mle2.HMM$BIC;mle3.HMM$BIC #so we use 2 components normal HMM
####################################################################################################
###Subtask b of e

N.HMM.mllk.mod <- function(parvect,y,n.dist,...)       
{
  m <- n.dist
  findNorm <- function(theta){
    mu<-theta[1]
    sd<-theta[2]
    dnorm(y,mean = mu, sd = sd)
  }
  #    print(parvect)
  if(m==1) return(-sum(dnorm(y,mean=parvect[1],sd=exp(parvect[2]),log=TRUE)))
  n          <- length(y)    
  pn         <- N.HMM.pw2pn(m,parvect)
  allprobs <- apply(cbind(pn$mu,pn$sigma), MARGIN = 1, FUN = findNorm)
  #allprobs   <- outer(x,pn$mu,pn$sigma,FUN=dnorm)
  allprobs   <- ifelse(!is.na(allprobs),allprobs,1)
  lscale     <- 0                                    
  foo        <- pn$delta                             
  for (i in 1:n)                                    
  {                                                
    foo    <- foo%*%pn$gamma*allprobs[i,] #pn$delta = pn$delta %*% pn$gamma for i = 1. Beregning er jf. slide 14, lecture 11
    sumfoo <- sum(foo)                               
    lscale <- lscale+log(sumfoo)                    
    foo    <- foo/sumfoo
  }                                               
  mllk       <- -lscale                            
  mllk                                              
}    
wpars2.HMM.opt <- N.HMM.pn2pw(m=m, mu=mle2.HMM$mu, sigma=mle2.HMM$sigma, gamma=mle2.HMM$gamma) #seeing how it works
H.w2.opt.HMM <- hessian(func=N.HMM.mllk.mod, x=wpars2.HMM.opt, n.dist=m, y=D$SLV)
sd.w2.opt.HMM <- sqrt(diag(solve(H.w2.opt.HMM)))
alpha <- 0.05
#CI.w2.HMM.l <- wpars2.HMM.opt -1 * qnorm(1-alpha/2) * sd.w2.opt.HMM
#CI.w2.HMM.u <- wpars2.HMM.opt + 1 * qnorm(1-alpha/2) * sd.w2.opt.HMM
#CI.w2.HMM <- matrix(cbind(CI.w2.HMM.l,wpars2.HMM.opt,CI.w2.HMM.u),nrow=3,byrow=T)
CI.w2.HMM <- replicate(3, wpars2.HMM.opt) + qnorm(0.975)*replicate(3, sd.w2.opt.HMM)*cbind(0,rep(-1,length(wpars2.HMM.opt)), 1)
CI.w2.HMM <- data.frame(CI.w2.HMM, row.names=c("theta1","theta2","eta1","eta2","gamma1","gamma2"))
colnames(CI.w2.HMM) <- c("mle","lower","upper")
CI.w2.HMM
####################################################################################################
###Subtask c of e
mle2.HMM$mu #mean of the two normal HMMs
mle2.HMM$sigma #sd of the two normal HMMs
mle2.HMM$gamma #t.p.m. for the transitioning between the two normal HMMs. Sølje: til, række: fra 
mle2.HMM$delta #Long term probabilities for ending up in each group.
               #delta 2 > delta 1, so the second model is most likely at t -> inf.
#solve(t(diag(m)-mle2.HMM$gamma+1),rep(1,m)) #this the way delta is calculated. System of two equations with two unknowns.

####################################################################################################
###Subtask d of e
#summing the two dists for the long term
N.pdf <- function(x,mu,sigma){
  return(1/(sigma*sqrt(2*pi)) * exp(-1/2 * (x-mu)^2/sigma^2))
}

final.Dist <- function(x){
  return(mle2.HMM$delta[1]*N.pdf(x, mle2.HMM$mu[1], mle2.HMM$sigma[1]) +
           mle2.HMM$delta[2]*N.pdf(x, mle2.HMM$mu[2], mle2.HMM$sigma[2]))
}

par(mfrow = c(1,1))
interval <- seq(min(D$SLV), max(D$SLV), length.out = length(D$SLV))
ggplot(D)+
  geom_histogram(aes(x = SLV, y=..density..))+
  geom_line(aes(x = interval, y = final.Dist(interval)))
plot(interval, final.Dist(interval), "l")

givenState1 <- function(x){
  return(mle2.HMM$gamma[1,1]*N.pdf(x, mle2.HMM$mu[1], mle2.HMM$sigma[1]) +
           mle2.HMM$gamma[1,2]*N.pdf(x, mle2.HMM$mu[2], mle2.HMM$sigma[2]))
}

givenState2 <- function(x){
  return(mle2.HMM$gamma[2,1]*N.pdf(x, mle2.HMM$mu[1], mle2.HMM$sigma[1]) +
           mle2.HMM$gamma[2,2]*N.pdf(x, mle2.HMM$mu[2], mle2.HMM$sigma[2]))
}

ggplot(D)+
  geom_histogram(aes(x = SLV, y=..density..), colour = "white")+
  geom_line(aes(x = interval, y = givenState1(interval), colour = "Given state 1"))+
  geom_line(aes(x = interval, y = givenState2(interval), colour = "Given state 2"))+
  #geom_line(aes(x = interval, y = givenState2(interval)*mle2.HMM$delta[2] + givenState1(interval)*mle2.HMM$delta[1], colour = "Testt stonks"), size = 2)+
  geom_line(aes(x = interval, y = final.Dist(interval), colour = "Long term returns"))+
  scale_colour_manual(values = c("darkgreen", "red", "orange", "purple"))+
  labs(colour = "", x = "SLV")+
  ggtitle("Long term and 1-step ahead return distributions")+
  theme_bw()

####################################################################################################
###Subtask e of f
N.HMM.generate_sample <- function(n,m,mu,sigma,gamma,delta=NULL)
{
  if(is.null(delta))delta<-solve(t(diag(m)-gamma+1),rep(1,m))
  mvect <- 1:m
  state <- numeric(n)
  state[1] <- sample(mvect,1,prob=delta)
  for (i in 2:n)
    state[i]<-sample(mvect,1,prob=gamma[state[i-1],])
  x <- rnorm(n,mean = mu[state], sd = sigma[state])
  return(x)
}

HMM_sample <- N.HMM.generate_sample(length(D$SLV), m = 2, mu = mle2.HMM$mu, sigma = mle2.HMM$sigma
                                    ,gamma = mle2.HMM$gamma, delta = mle2.HMM$delta)
ggplot(data = D)+
  geom_histogram(aes(x=SLV,colour="data"),alpha=0.3)+
  geom_histogram(aes(x=HMM_sample,colour="sample"), alpha=0.3)

######################### SORTING SO PARAMS ARE N~##############################
N.HMM.mllk.sort <- function(parvect,x,m,...)       
{
  findNorm <- function(theta){
    mu<-theta[1]
    sd<-theta[2]
    dnorm(x,mean = mu, sd = sd)
  }
  #    print(parvect)
  if(m==1) return(-sum(dnorm(x,mean=parvect[1],sd=exp(parvect[2]),log=TRUE)))
  n          <- length(x)    
  pn         <- N.HMM.pw2pn(m,parvect)
  allprobs <- apply(cbind(pn$mu,pn$sigma), MARGIN = 1, FUN = findNorm)
  #allprobs   <- outer(x,pn$mu,pn$sigma,FUN=dnorm)
  allprobs   <- ifelse(!is.na(allprobs),allprobs,1)
  lscale     <- 0                                    
  foo        <- pn$delta                           
  for (i in 1:n)                                    
  {                                                
    foo    <- foo%*%pn$gamma*allprobs[i,] #pn$delta = pn$delta %*% pn$gamma for i = 1. Beregning er jf. slide 14, lecture 11
    sumfoo <- sum(foo)                               
    lscale <- lscale+log(sumfoo)                    
    foo    <- foo/sumfoo
  }                                               
  mllk       <- -lscale                            
  mllk                                              
}    

N.HMM.mle.nlminb.sort <- function(x,m,mu0,sigma0,gamma0,...)
{                                                      
  parvect0 <- N.HMM.pn2pw(m,mu0,sigma0,gamma0)
  np        <- length(parvect0)                          
  lower    <- rep(-10,np)
  upper    <- c(rep(max(x),m),rep(10,np-m))
  mod      <- nlminb(parvect0,N.HMM.mllk.sort,x=x,m=m,
                     lower=lower,upper=upper)
  if(mod$convergence!=0){
    print(mod)
  }
  pn       <- N.HMM.pw2pn(m,mod$par)
  mllk     <- mod$objective
  AIC       <- 2*(mllk+np)
  n         <- sum(!is.na(x))
  BIC       <- 2*mllk+np*log(n)
  list(mu=pn$mu,sigma=pn$sigma,gamma=pn$gamma,delta=pn$delta,   
       code=mod$convergence,mllk=mllk,AIC=AIC,BIC=BIC)   
}


mu <- c(0,0) * mean(D$SLV) #arbitrary
sigma <- c(1,1) * sd(D$SLV) #arbitrary
gamma <- matrix(c(0.5,0.5,0.5,0.5), ncol=2,byrow=T) #arbitrary

k <- 200
GAMMA <- matrix(ncol=m*m,nrow=k)
Mu <- matrix(ncol=m,nrow=k)
Sigma <- matrix(ncol=m,nrow=k)
Delta <- matrix(ncol=m,nrow=k)
Code <- numeric(k)
for(i in 1:k){
  set.seed(i)
  ## generate sample
  y.sim <- N.HMM.generate_sample(length(D$SLV), m = 2, mu = mle2.HMM$mu, sigma = mle2.HMM$sigma
                                 ,gamma = mle2.HMM$gamma, delta = mle2.HMM$delta)
  ## fit model to sample
  mle.tmp <- N.HMM.mle.nlminb.sort(x=y.sim, m=m, mu0=mu, sigma0=sigma, gamma0=gamma)
  ## Store result
  index <- order(mle.tmp$sigma)#TILFØJET FOR AT SORTERE!
  mle.sort <- mle.tmp#
  mle.sort$mu <- mle.tmp$mu[index]#
  mle.sort$sigma <- mle.tmp$sigma[index]#
  mle.sort$delta <- mle.tmp$delta[index]#
  mle.tmp <- mle.sort#
  GAMMA[i, ] <- c(mle.tmp$gamma[index[1], ],#
                  mle.tmp$gamma[index[2], ])#
  Mu[i, ] <-  mle.tmp$mu
  Sigma[i, ] <- mle.tmp$sigma
  Delta[i, ] <-   mle.tmp$delta
  Code[i] <-      mle.tmp$code
  print(c(i,Code[i]))
}
sum(Code!=1)

par(mfrow=c(1,3))
hist(Mu)
hist(Sigma)
hist(Delta)

par(mfrow=c(1,2))
hist(Mu[,1])
abline(v=mle2.HMM$mu[1],col=2,lwd=3)
hist(Mu[,2])
abline(v=mle2.HMM$mu[2],col=2,lwd=3)

hist(Sigma[,1])
abline(v=mle2.HMM$sigma[1],col=2,lwd=3)
hist(Sigma[,2])
abline(v=mle2.HMM$sigma[2],col=2,lwd=3)

hist(Delta[,1])
abline(v=mle2.HMM$delta[1],col=2,lwd=3)
hist(Delta[,2])
abline(v=mle2.HMM$delta[2],col=2,lwd=3)

alpha <- 0.05
quantile(Mu[,1],c(alpha/2, 1-alpha/2)) #[-0.00628  0.0104  ] mle: 0.00109
quantile(Mu[,2],c(alpha/2, 1-alpha/2)) #[-0.00440  0.00948 ] mle: 0.00183
quantile(Sigma[,1],c(alpha/2, 1-alpha/2))#[0.0248 0.0505] mle: 0.0317
quantile(Sigma[,2],c(alpha/2, 1-alpha/2)) #[0.0494 0.0647] mle: 0.0596
quantile(Delta[,1],c(alpha/2, 1-alpha/2)) #[0.3018 0.500 ] mle: 0.353
quantile(Delta[,2],c(alpha/2, 1-alpha/2)) #[0.500 0.698] mle: 0.647

M <- diag(m^2+m) #We expect all zeros apart from the diagonal.
diagonal <- c(1,1,exp(wpars2.HMM.opt[3]),exp(wpars2.HMM.opt[4]),exp(wpars2.HMM.opt[5]),wpars2.HMM.opt[6])
M <- diagonal * M
#M <- mu1/theta1,  mu2/theta1, sigma1/theta1, sigma2/theta1, gamma1/theta1, gamma2/theta1
  #   mu1/theta2,  mu2/theta2, sigma1/theta2, sigma2/theta2, gamma1/theta2, gamma2/theta2
  #   mu1/eta1,    mu2/eta1,   sigma1/eta1,   sigma2/eta1,   gamma1/eta1,   gamma2/eta1
  #   mu1/eta2,    mu2/eta2,   sigma1/eta2,   sigma2/eta2,   gamma1/eta2,   gamma2/eta2
  #   mu1/tgamma1, mu2/tgamma, sigma1/tgamme, sigma2/tgamme, gamma1/tgamma, gamma2/tgamma

invG <- t(M) %*% solve(H.w2.opt.HMM) %*% M #p. 54 in HMM
sd.n2.opt.HMM <- sqrt(diag(invG))
CI.n2.HMM <- replicate(3, c(mle2.HMM$mu, mle2.HMM$sigma, mle2.HMM$delta)) + qnorm(0.975)*replicate(3, sd.n2.opt.HMM)*cbind(0,rep(-1,length(c(mle2.HMM$mu, mle2.HMM$sigma, mle2.HMM$delta))), 1)
CI.n2.HMM <- data.frame(CI.n2.HMM, row.names=c("mu1","mu2","sigma1","sigma2","gamma1","gamma2")) #"gamma/delta" er diag eller off-diag
colnames(CI.n2.HMM) <- c("mle","lower","upper")
CI.n2.HMM #OBS! på DELTA!, compare to:
#CI.w2.HMM
mle2.HMM$delta
mle2.HMM$gamma
#mle for gamma1 og gamma2 skal forstås som hhv. sandsynlighederne for at være i state 1 og blive dér og være i state 2 og blive dér.

