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
  sigma <- cumsum(exp(eta))
  delta <- exp(tau)/(1+sum(exp(tau)))
  delta <- c(1-sum(delta),delta)
  return(list(mu=mu,sigma=sigma,delta=delta))
}

N.nll <- function(p,n.dist,y){
  if(n.dist==1){return( -sum(dnorm(y, mean=p[1], sd=exp(p[2]), log=T)))
  }
  theta <- p[1:n.dist]
  eta <- p[(n.dist+1):(2*n.dist)]
  tau <- p[(2*n.dist+1):(3*n.dist-1)]
  n.pars <- N.mix.pw2pn(n.dist=n.dist, theta=theta, eta=eta, tau=tau)
  n <- length(y)
  nll <- 0
  for(i in 1:n){
    nll <- nll - log(sum(n.pars$delta * dnorm(y[i], mean=n.pars$mu, sd=n.pars$sigma)))
  }
  return(nll)
}
#####
N.pdf <- function(x,mu,sigma){
  return(1/(sigma*sqrt(2*pi)) * exp(-1/2 * (x-mu)^2/sigma^2))
}

final.Dist <- function(x, delta, mu, sigma){
  L <- 0
  for (i in 1:length(delta)){
    L <- L + delta[i] * N.pdf(x, mu[i], sigma[i])
  }
  return(L)
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
#opt2$par[(m+1+1):(2*m)] = opt2$par[(m+1):(2*m-1)] + exp(opt2$par[(m+1+1):(2*m)])
(pars2 <- N.mix.pw2pn(n.dist=m,opt2$par[1:m],opt2$par[(m+1):(2*m)],opt2$par[(2*m+1):(3*m-1)]))

D$interval <- seq(min(D$SLV), max(D$SLV), length.out = length(D$SLV))
D$dist2 <- sapply(D$interval, final.Dist, delta = pars2$delta, mu = pars2$mu, sigma = pars2$sigma)

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(aes(colour='1st'), fun = dnorm, n = dim(D)[1], args = list(mean = pars2$mu[1],
                                                        sd = pars2$sigma[1])) +
  stat_function(aes(colour='2nd'), fun = dnorm, n = dim(D)[1], args = list(mean = pars2$mu[2],
                                                        sd = pars2$sigma[2])) +
  geom_line(aes(x = interval, y = dist2, colour = "Combined"), size = 1.5)+
  scale_colour_manual(values = c("blue", "orange", "springgreen4"))+
  labs(colour = "", x = "SLV")+
  theme_bw()+
  theme(legend.position = "top")+
  ggtitle("(2) Normal mixture model fitted to SLV")


ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  geom_line(aes(x = interval, y = dist2), colour = "blue")+
  ggtitle("(2) Normal distribution fitted to SLV")
#with delta = 0.91, 0.09 for red, orange
##################################################
m <- 3 #3 components
mu <- c(-1/2,1,3/2)*mean(D$SLV) #init vals for the optim
sigma <- c(1/2,1,3/2)*sd(D$SLV) #init vals for the optim
delta <- c(0.34,0.33) #fifty-fifty for the two dists
wpars <- N.mix.pn2pw(n.dist=m, mu=mu, sigma=sigma, delta=delta)
( N.mix.pw2pn(n.dist=m, theta=wpars$theta, eta=wpars$eta, tau=wpars$tau) ) #checking that pn2pw -> pw2pn works
p <- c(wpars$theta, wpars$eta, wpars$tau)
N.nll(p=p, n.dist=m, y=D$SLV)
opt3 <- nlminb(start=c(wpars$theta, wpars$eta, wpars$tau), N.nll, n.dist=m, y=D$SLV)
opt3$par
#opt3$par[(m+1+1):(2*m)] = opt3$par[(m+1):(2*m-1)] + cumsum(exp(opt3$par[(m+1+1):(2*m)]))
#opt3$par
(pars3 <- N.mix.pw2pn(n.dist=m,opt3$par[1:m],opt3$par[(m+1):(2*m)],opt3$par[(2*m+1):(3*m-1)]))

D$interval <- seq(min(D$SLV), max(D$SLV), length.out = length(D$SLV))
D$dist3 <- sapply(D$interval, final.Dist, delta = pars3$delta, mu = pars3$mu, sigma = pars3$sigma)

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(aes(colour='1st'), fun = dnorm, n = dim(D)[1], args = list(mean = pars3$mu[1],
                                                        sd = pars3$sigma[1]))+
  stat_function(aes(colour='2nd'), fun = dnorm, n = dim(D)[1], args = list(mean = pars3$mu[2],
                                                        sd = pars3$sigma[2])) +
  stat_function(aes(colour='3rd'), fun = dnorm, n = dim(D)[1], args = list(mean = pars3$mu[3],
                                                        sd = pars3$sigma[3])) +
  geom_line(aes(x = interval, y = dist3, colour = "Combined"), size = 1.5)+
  scale_colour_manual(values = c("blue", "orange","turquoise3", "springgreen4"))+
  labs(colour = "", x = "SLV")+
  theme_bw()+
  theme(legend.position = "top")+
  ggtitle("(3) Normal mixture model fitted to SLV")

D$interval <- seq(min(D$SLV), max(D$SLV), length.out = length(D$SLV))
D$dist3 <- sapply(D$interval, final.Dist, delta = pars3$delta, mu = pars3$mu, sigma = pars3$sigma)
ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  geom_line(aes(x = interval, y = dist3), colour = "blue")+
  ggtitle("(3) Normal distribution fitted to SLV") 
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
CIs.w2 <- data.frame(CIs.w2, row.names=c("lower","mle","upper"))
colnames(CIs.w2) = c("theta1", "theta2", "eta1", "eta2", "tau")
#CIs.n2 <- matrix(rep(0,( dim(CIs.w2)[1] * ( dim(CIs.w2)[2]+1 ) ) ), nrow=3, byrow=T)
#for (i in 1:3){
#  n2_temp <- N.mix.pw2pn(n.dist=m, theta=CIs.w2[i,(1:m)], eta=CIs.w2[i,(m+1):(2*m)], tau=CIs.w2[i,(2*m+1):(3*m-1)])
#  CIs.n2[i,1:dim(CIs.n2)[2]] <- c(n2_temp$mu, n2_temp$sigma, n2_temp$delta) #1st: lower, 2nd: MLE, 3rd: upper
#}

## parametric bootstrap ##
rmix.normal <- function(k,delta,mu,sigma){
  states <- sample(1:length(delta),size=k,replace=TRUE,
                   prob=delta)
  return(rnorm(k, mean = mu[states], sd = sigma[states]))
}

#### Run for new mixture model parametric bootstraps
# k <- length(D$SLV)
# m <- 3
# mu0 <- c(1/2,0,-1/2)*mean(D$SLV)
# sigma0 <- c(1,2,5)*sd(D$SLV)
# delta0 <- c(0.34,0.33) #fifty-fifty for the two dists
# wpars0 <- N.mix.pn2pw(n.dist=m, mu=mu0, sigma=sigma0, delta=delta0)
# 
# runs_for_parametric_bootstrap <- 2000
# 
# Mu <- matrix(ncol=m,nrow=runs_for_parametric_bootstrap)
# Sigma <- matrix(ncol=m,nrow=runs_for_parametric_bootstrap)
# Delta <- matrix(ncol=m,nrow=runs_for_parametric_bootstrap)
# Code <- numeric(runs_for_parametric_bootstrap)
# 
# for(i in 1:runs_for_parametric_bootstrap){
#   set.seed(i)
#   ## generate sample
#   y.sim <- rmix.normal(k, delta = pars2$delta, mu = pars2$mu, sigma = pars2$sigma)
#   ## fit model to sample
#   p <- c(wpars0$theta, wpars0$eta, wpars0$tau)
#   opt2.tmp <- nlminb(start=p, N.nll, n.dist=m, y=y.sim)
#   pars2.tmp <- N.mix.pw2pn(n.dist=m, opt2.tmp$par[1:m], opt2.tmp$par[(m+1):(2*m)], opt2.tmp$par[(2*m+1):(3*m-1)])
# 
#   ## Store result
#   Mu[i, ] <-  pars2.tmp$mu
#   Sigma[i, ] <- pars2.tmp$sigma
#   Delta[i, ] <-   pars2.tmp$delta
#   Code[i] <-      opt2.tmp$convergence
#   print(c(i,Code[i]))
# }
# sum(Code!=1)
#save(Mu, Sigma, Delta, Code, file = "./3MixtureModelParametricBootstrap2000.Rdata")
load("./MixtureModelParametricBootstrap2000.Rdata")
melt(data.frame("Mu1" = Mu[Code!=1,1], "Mu2" = Mu[Code!=1,2])) %>%
  union(melt(data.frame("Sigma1" = Sigma[Code!=1,1], "Sigma2" = Sigma[Code!=1,2]))) %>%
  union(melt(data.frame("Delta1" = Delta[Code!=1,1], "Delta2" = Delta[Code!=1,2]))) %>%
  ggplot()+
  geom_histogram(aes(x = value), colour = "white")+
  facet_wrap(~variable, scales =  "free")

# load("./3MixtureModelParametricBootstrap2000.Rdata")
# melt(data.frame("Mu1" = Mu[Code!=1,1], "Mu2" = Mu[Code!=1,2], "Mu3" = Mu[Code!=1,3])) %>%
#   union(melt(data.frame("Sigma1" = Sigma[Code!=1,1], "Sigma2" = Sigma[Code!=1,2], "Sigma3" = Sigma[Code!=1,3]))) %>%
#   union(melt(data.frame("Delta1" = Delta[Code!=1,1], "Delta2" = Delta[Code!=1,2], "Delta3" = Delta[Code!=1,3]))) %>%
#   ggplot()+
#   geom_histogram(aes(x = value), colour = "white")+
#   facet_wrap(~variable, scales =  "free")

ParametricBootCI <- data.frame(round(cbind(c(pars2$mu, pars2$sigma^2, pars2$delta),
            rbind(quantile(Mu[ ,1],prob=c(0.025,0.975)),
            quantile(Mu[ ,2],prob=c(0.025,0.975)),
            quantile(Sigma[ ,1],prob=c(0.025,0.975))^2,
            quantile(Sigma[ ,2],prob=c(0.025,0.975))^2,
            quantile(Delta[ ,1],prob=c(0.025,0.975)),
            quantile(Delta[ ,2],prob=c(0.025,0.975)))),
      digits = 5)); colnames(ParametricBootCI) <- c("MLE", "2.5%", "97.5%"); rownames(ParametricBootCI) <- c("Mu1", "Mu2", "Sigma1^2", "Sigma2^2", "Delta1", "Delta2")
ParametricBootCI[,c(2,1,3)]

#Interpretations:
#Der er grundlæggende sæt to fordelinger som dataene kommer fra
#én med nær-positiv gennemsnit of smal varians, og én med primært negativ gennemsnit
#og bred varians. Ca. 90 af observationerne kommer fra den positive, smalle fordeling.
#Dvs. vi bruger 90% af tiden i den mindre volatile fordeling, hvor der er svag evidens
#for at udviklingen er positiv - det er dog ikke muligt på at 95% confidens-niveau at sige at gennemsnittet for fordeling
#1 er positiv og heller ikke om gennemsnittet af den anden (mere volatile) fordeling
#er negativ.

####################################################################################################
###Subtask c of e
##PL of one of the variance parameters of the two component model.
#Modifying function N.nll from earlier to accept single value for on of the variances whilst optimising the others:
N.pl <- function(eta0, n.dist, y, first=T){ #to be used for PL of the 2nd variance parameter (working) of the 2 components normal mixture
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
  n <- length(y)
  ll <- 0
  for (i in 1:n){
    ll <- ll + log(sum(n.pars2_temp$delta * dnorm(y[i], mean=n.pars2_temp$mu, sd=n.pars2_temp$sigma)))
  }
  return(ll)
}

eta1.w <- seq(opt2$par[m+1]-6*sqrt(V.w2)[m+1],
              opt2$par[m+1]+6*sqrt(V.w2)[m+1],
              length=50)
eta2.w <- seq(opt2$par[2*m]-6*sqrt(V.w2)[2*m],
             opt2$par[2*m]+6*sqrt(V.w2)[2*m],
            length=50)
m <- 2
N.llp <- sapply(X = eta1.w, FUN=N.pl, n.dist=m, y=D$SLV)
N.llp2 <- sapply(X = eta2.w, FUN=N.pl, n.dist=m, y=D$SLV, first=F)

plot(eta1.w,exp(N.llp-(max(N.llp))),type="l",
     ylab = "Normalized likelihood",
     xlab = expression(paste(eta[1])),
     main = "First variance param of 2 NMM") #smart Jan-plot
grid()
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
#see above


####################################################################################################
###Subtask e of e
#Model interpretation:
#See above interpreations from b of e regarding two states: growth, less-volatile and volatile (and possibly decreasing).

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
  sigma <- cumsum(epar[(m+1):(2*m)])
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
#following five lines are repeated from above?
m <- 2
mu <- c(-2,2) * mean(D$SLV) #arbitrary
sigma <- c(1/2,4) * sd(D$SLV) #arbitrary
gamma <- matrix(c(0.95,0.05,0.05,0.95), ncol=2,byrow=T) #arbitrary
wpars2.HMM <- N.HMM.pn2pw(m=m, mu=mu, sigma=sigma, gamma=gamma) #seeing how it works
mle2.HMM <- N.HMM.mle.nlminb(x=D$SLV,m=m,mu0=mu,sigma0=sigma,gamma0=gamma)

#mllk1 <- N.HMM.mllk(parvect=c(mean(D$SLV),log(sd(D$SLV))),x=D$SLV,m=1)
#mle1.HMM <- N.HMM.mle.nlminb(x=D$SLV,m=1,mu0=mean(D$SLV),sigma0=sd(D$SLV),gamma0=c())

D$interval <- seq(min(D$SLV), max(D$SLV), length.out = length(D$SLV))
D$distHMM2 <- sapply(D$interval, final.Dist, delta = mle2.HMM$delta, mu = mle2.HMM$mu,
                     sigma = mle2.HMM$sigma)

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(aes(colour='1st'), fun = dnorm, n = dim(D)[1], args = list(mean = mle2.HMM$mu[1],
                                                        sd = mle2.HMM$sigma[1])) +
  stat_function(aes(colour='2nd'), fun = dnorm, n = dim(D)[1], args = list(mean = mle2.HMM$mu[2],
                                                        sd = mle2.HMM$sigma[2])) +
  geom_line(aes(x = interval, y = distHMM2, colour = "Combined"), size = 1.5)+
  scale_colour_manual(values = c("blue", "orange","springgreen4"))+
  labs(colour = "", x = "SLV")+
  theme_bw()+
  theme(legend.position = "top")+
  ggtitle("HMM with 2 normal distributions fitted to SLV")

m3 <- 3
mu3 <- c(0,0,0) * mean(D$SLV) #arbitrary -1,1/2,1
sigma3 <- c(1, 1, 1) * sd(D$SLV) #arbitrary 1/4,1,4
gamma3 <- matrix(rep(1/3,9), ncol=3,byrow=T) #arbitrary 0.95, 0.25, 0.25
#wpars3.HMM <- N.HMM.pn2pw(m=m0, mu=mu0, sigma=sigma0, gamma=gamma0)
#mllk3 <- N.HMM.mllk(parvect=wpars3.HMM,x=D$SLV,m=m0)

mle3.HMM <- N.HMM.mle.nlminb(x=D$SLV,m=m3,mu0=mu3,sigma0=sigma3,gamma0=gamma3)
D$distHMM3 <- sapply(D$interval, final.Dist, delta = mle3.HMM$delta, mu = mle3.HMM$mu,
                     sigma = mle3.HMM$sigma)

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(aes(colour='1st'), fun = dnorm, n = dim(D)[1], args = list(mean = mle3.HMM$mu[1],
                                                        sd = mle3.HMM$sigma[1])) +
  stat_function(aes(colour='2nd'), fun = dnorm, n = dim(D)[1], args = list(mean = mle3.HMM$mu[2],
                                                        sd = mle3.HMM$sigma[2])) +
  stat_function(aes(colour='3rd'), fun = dnorm, n = dim(D)[1], args = list(mean = mle3.HMM$mu[3],
                                                        sd = mle3.HMM$sigma[3])) +
  geom_line(aes(x = interval, y = distHMM3, colour = "Combined"), size = 1.5)+
  scale_colour_manual(values = c("blue", "orange","turquoise3", "springgreen4"))+
  labs(colour = "", x = "SLV")+
  theme_bw()+
  theme(legend.position = "top")+
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
  geom_histogram(aes(x = SLV, y=..density..), colour = "white")+
  geom_line(aes(x = interval, y = final.Dist(interval)))+
  theme_bw()+
  ggtitle("(2) HMM")
#plot(interval, final.Dist(interval), "l")

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
  geom_line(aes(x = interval, y = givenState1(interval), colour = "One-step ahead Given state 1"),size = 1)+
  geom_line(aes(x = interval, y = givenState2(interval), colour = "One-step ahead Given state 2"),size = 1)+
  #geom_line(aes(x = interval, y = givenState2(interval)*mle2.HMM$delta[2] + givenState1(interval)*mle2.HMM$delta[1], colour = "Testt stonks"), size = 2)+
  geom_line(aes(x = interval, y = final.Dist(interval), colour = "Long term returns"), size = 1)+
  scale_colour_manual(values = c("blue", "red", "orange", "darkgreen"))+
  labs(colour = "", x = "SLV")+
  ggtitle("Long term and 1-step ahead return distributions")+
  theme_bw()+
  theme(legend.position = "top")

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
  geom_line(aes(x = date, y=SLV,colour="data"),alpha=0.8, size = 1)+
  geom_line(aes(x = date, y=HMM_sample,colour="sample"), alpha=0.8, size = 1)

#PARAMETRIC BOOTSTRAP
mu <- c(0,0) * mean(D$SLV) #arbitrary
sigma <- c(1,1) * sd(D$SLV) #arbitrary
gamma <- matrix(c(0.95,0.05,0.05,0.95), ncol=2,byrow=T) #arbitrary

k <- 1000
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
  mle.tmp <- N.HMM.mle.nlminb(x=y.sim, m=m, mu0=mu, sigma0=sigma, gamma0=gamma)
  ## Store resultater
  GAMMA[i, ] <- c(mle.tmp$gamma[1, ],#
                  mle.tmp$gamma[2, ])#
  Mu[i, ] <-  mle.tmp$mu
  Sigma[i, ] <- mle.tmp$sigma
  Delta[i, ] <-   mle.tmp$delta
  Code[i] <-      mle.tmp$code
  print(c(i,Code[i]))
}
sum(Code!=1)
#save(Mu, Sigma, Delta, Code, GAMMA, file = "./HMMBOot2States.Rdata")
load("./HMMBOot2States.Rdata")

melt(data.frame("Mu1" = Mu[Code!=1,1], "Mu2" = Mu[Code!=1,2],
                "Sigma1" = Sigma[Code!=1,1], "Sigma2" = Sigma[Code!=1,2],
                "Delta1" = Delta[Code!=1,1], "Delta2" = Delta[Code!=1,2],
                "Gamma11" = GAMMA[Code!=1,1], "Gamma12" = GAMMA[Code!=1,2],
                "Gamma21" = GAMMA[Code!=1,3], "Gamma22" = GAMMA[Code!=1,4])) %>%
  ggplot()+
  geom_histogram(aes(x = value), colour = "white")+
  facet_wrap(~variable, scales =  "free")+
  theme_bw()


ParametricBootCI <- data.frame(round(cbind(c(mle2.HMM$mu, mle2.HMM$sigma^2, mle2.HMM$delta, c(mle2.HMM$gamma[1,],mle2.HMM$gamma[2,])),
                                           rbind(quantile(Mu[Code!=1 ,1],prob=c(0.025,0.975)),
                                                 quantile(Mu[Code!=1 ,2],prob=c(0.025,0.975)),
                                                 quantile(Sigma[Code!=1 ,1],prob=c(0.025,0.975))^2,
                                                 quantile(Sigma[Code!=1 ,2],prob=c(0.025,0.975))^2,
                                                 quantile(Delta[Code!=1 ,1],prob=c(0.025,0.975)),
                                                 quantile(Delta[Code!=1 ,2],prob=c(0.025,0.975)),
                                                 quantile(GAMMA[Code!=1 ,1],prob=c(0.025,0.975)),
                                                 quantile(GAMMA[Code!=1 ,2],prob=c(0.025,0.975)),
                                                 quantile(GAMMA[Code!=1 ,3],prob=c(0.025,0.975)),
                                                 quantile(GAMMA[Code!=1 ,4],prob=c(0.025,0.975)))),
                                     digits = 5)); colnames(ParametricBootCI) <- c("MLE", "2.5%", "97.5%"); rownames(ParametricBootCI) <- c("Mu1", "Mu2",
                                                                                                                                            "Sigma1^2", "Sigma2^2",
                                                                                                                                            "Delta1", "Delta2",
                                                                                                                                            "GAMMA11", "GAMMA12",
                                                                                                                                            "GAMMA21", "GAMMA22")
ParametricBootCI[,c(2,1,3)]

#HMM f - Short term predictions
y <- D$SLV
x_first <- 1
x <- c()
x[1] <- x_first
x_test <- sample(x = c(1,2), size = length(D$SLV), replace <- TRUE, prob = mle2.HMM$delta)

L.HMM <- function(theta){
  for (i in 2:length(x)){
    P_Xi_given_Ximinus1 <- mle2.HMM$gamma[x[i-1],x[i]]
    P_Yi_given_Xi <- dnorm(y[i], mean = mle2.HMM$mu[x[i]], sd = mle2.HMM$sigma[x[i]])
    L <- L + log(P_Xi_given_Ximinus1) + log(P_Yi_given_Xi)
  }
  return(L)
}


expand.grid(rep(1,length(D$SLV)), rep(2,length(D$SLV)))[,1]

theta <- sample(c(4,0.6), size = length(D$SLV), replace = T)

x <- nlminb(theta, L.HMM, lower = 0)

L.HMM(x_test)

x <- c()
x[1] <- 1
for (i in 2:length(D$SLV)){
  ll_x_given_x_minus1 <- log(mle2.HMM$gamma[x[i-1],])
  ll_y_given_x <- log(dnorm(D$SLV[i], mean = mle2.HMM$mu, sd = mle2.HMM$sigma))
  ll_x_given_y <- ll_y_given_x + ll_x_given_x_minus1
  cat(ll_x_given_y, "SLV=",round(D$SLV[i],4),", y given x =", round(ll_y_given_x,4),"\n")
  x[i] <- order(ll_x_given_y)[2]
}
x
