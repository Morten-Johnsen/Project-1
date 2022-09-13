rm(list = ls())
library(lubridate)
library(tidyverse)
library(reshape)
library(stringr)
library(pheatmap)
library(numDeriv)

#Retrieve data from the descriptive statistics script
##### TILPAS DEN HER TIL DEN COMPUTER DEER BRUGES
if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/Statistical Modelling/Project-1/")
} else {
  setwd("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1")
}

load("dataWindPower.Rdata")
source("testDistribution.R")

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

#### WIND POWER ####
#Fitting models to wind power:
par <- nlminb(start = 0.2, objective = testDistribution,
              distribution = "exponential",
              x = D$pow.obs.norm)

ggplot(D)+
  geom_histogram(aes(x = pow.obs.norm, y = ..density..), bins = 20)+
  theme_bw()+
  stat_function(fun = dexp, n = length(D$pow.obs.norm), args = list(rate = par$par))

#g.pow.obs
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

## CI ## WIND POWER
par(mfrow=c(1,1))
alpha <- 0.05
c <- exp(-0.5 * qchisq(1-alpha, df = 1))
#likelihood-based
mle.pow.exp <- par$par

pow.fun <- function(lambda, data){
  prod(dexp(x = data, rate = lambda, log = F))
}

l.pow.fun <- function(lambda, data){
  sum(dexp(x = data, rate = lambda, log = T))
}

CIfun.pow <- function(y){
  sum(dexp(x = D$pow.obs.norm, rate = mle.pow.exp, log = T)) -
    sum(dexp(x = D$pow.obs.norm, rate = y, log = T)) -
    0.5 * qchisq(1-alpha, df = 1)
}
lambdas <- seq(0.2,7, by = 0.1)
pow <- sapply(X = lambdas, FUN = pow.fun, data = D$pow.obs.norm)
plot(lambdas, pow/max(pow), col = 1, type = "l", xlab = expression(paste(lambda)),
     main = "Parameter values for exponential model of power production")
CI.pow <- c(uniroot(f = CIfun.pow, interval = c(0.2, mle.pow.exp))$root,
            uniroot(f = CIfun.pow, interval = c(mle.pow.exp, 7))$root)
lines(range(lambdas), c*c(1,1), col = 2)

#wald
H.pow <- hessian(l.pow.fun, mle.pow.exp, data = D$pow.obs.norm)
V.pow <- as.numeric(-1/H.pow)
wald.pow <- mle.pow.exp + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.pow)

## CI ## WIND SPEED
par(mfrow=c(1,2))
#likelihood-based
mle.ws30.weib <- par.ws30$par

ws30.fun <- function(shape, scale, data){#####
  prod(dweibull(x = data, shape = shape, scale = scale, log = F)*2)#to not get full zeros
}

l.ws30.fun <- function(shape, scale, data){#####
  sum(dweibull(x = data, shape = shape, scale = scale, log = T))
}

CIfun.ws30 <- function(y, shape = T){##### T for shape, F for scale
  if(shape){
    sum(dweibull(x = D$ws30, shape = mle.ws30.weib[1], scale = mle.ws30.weib[2], log = T)) -
      sum(dweibull(x = D$ws30, shape = y, scale = mle.ws30.weib[2], log = T)) - 
      0.5 * qchisq(1-alpha, df = 1)
  } else {
    sum(dweibull(x = D$ws30, shape = mle.ws30.weib[1], scale = mle.ws30.weib[2], log = T)) -
      sum(dweibull(x = D$ws30, shape = mle.ws30.weib[1], scale = y, log = T)) - 
      0.5 * qchisq(1-alpha, df = 1) 
  }
}
shapes <- seq(1, 3.5, by = 0.1)
ws30.shape <- sapply(X = shapes, FUN = ws30.fun, scale = mle.ws30.weib[2], data = D$ws30)
plot(shapes, ws30.shape/max(ws30.shape), col = 1, type = "l", xlab = "shape",
     main = "Parameter value for shape for weibull model of wind speed")
CI.ws30.shape <- c(uniroot(f = CIfun.ws30, interval = c(1, mle.ws30.weib[1]), shape = T)$root,
                uniroot(f = CIfun.ws30, interval = c(mle.ws30.weib[1], 3.5), shape = T)$root)
lines(range(shapes), c*c(1,1), col = 2)

scales <- seq(7, 12, by = 0.1)
ws30.scale <- sapply(X = scales, FUN = ws30.fun, shape = mle.ws30.weib[1], data = D$ws30)
plot(scales, ws30.scale/max(ws30.scale), col = 1, type = "l", xlab = "scale",
     main = "Parameter value for scale for weibull model of wind speed")
CI.ws30.scale <- c(uniroot(f = CIfun.ws30, interval = c(7, mle.ws30.weib[2]), shape = F)$root,
                   uniroot(f = CIfun.ws30, interval = c(mle.ws30.weib[2], 12), shape = F)$root)
lines(range(scales), c*c(1,1), col = 2)

#wald
n <- dim(D)[1]
H.ws30.shape <- hessian(l.ws30.fun, mle.ws30.weib[1], scale = mle.ws30.weib[2], data = D$ws30)
V.ws30.shape <- as.numeric(-1/H.ws30.shape)
H.ws30.scale <- hessian(l.ws30.fun, mle.ws30.weib[2], shape = mle.ws30.weib[1], data = D$ws30)
V.ws30.scale <- as.numeric(-1/H.ws30.scale)
wald.ws30.shape <- mle.ws30.weib[1] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.ws30.shape)
wald.ws30.scale <- mle.ws30.weib[2] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.ws30.scale)


## CI ## WIND DIRECTION
par(mfrow=c(1,2))
#likelihood-based
mle.wd30.norm <- par.wd30$par

wd30.fun <- function(mu, sigma, data){#####
  prod(dnorm(x = data, mean = mu, sd = sigma, log = F))
}

l.wd30.fun <- function(mu, sigma, data){#####
  sum(dnorm(x = data, mean = mu, sd = sigma, log = T))
}

CIfun.wd30 <- function(y, mu = T){##### T from mean, F for sigma
  if(mu){
    sum(dnorm(x = D$wd30.centered, mean = mle.wd30.norm[1], sd = mle.wd30.norm[2], log = T)) -
      sum(dnorm(x = D$wd30.centered, mean = y, sd = mle.wd30.norm[2], log = T)) - 
      0.5 * qchisq(1-alpha, df = 1)
  } else {
    sum(dnorm(x = D$wd30.centered, mean = mle.wd30.norm[1], sd = mle.wd30.norm[2], log = T)) -
      sum(dnorm(x = D$wd30.centered, mean = mle.wd30.norm[1], sd = y, log = T)) - 
      0.5 * qchisq(1-alpha, df = 1) 
  }
}
mus <- seq(0, 6, by = 0.1)
wd30.mu <- sapply(X = mus, FUN = wd30.fun, sigma = mle.wd30.norm[2], data = D$wd30.centered)
plot(mus, wd30.mu/max(wd30.mu), col = 1, type = "l", xlab = expression(paste(mu)),
     main = "Parameter value for mean for normal model of wind direction")
CI.wd30.mu <- c(uniroot(f = CIfun.wd30, interval = c(0, mle.wd30.norm[1]), mu = T)$root,
            uniroot(f = CIfun.wd30, interval = c(mle.wd30.norm[1], 6), mu = T)$root)
lines(range(mus), c*c(1,1), col = 2)

sigmas <- seq(0, 2.5, by = 0.1)
wd30.sigma <- sapply(X = sigmas, FUN = wd30.fun, mu = mle.wd30.norm[1], data = D$wd30.centered)
plot(sigmas^2, wd30.sigma/max(wd30.sigma), col = 1, type = "l", xlab = expression(paste(sigma^2)),
     main = "Parameter value for var for normal model of wind direction")
CI.wd30.sigma <- c(uniroot(f = CIfun.wd30, interval = c(0, mle.wd30.norm[2]), mu = F)$root,
             uniroot(f = CIfun.wd30, interval = c(mle.wd30.norm[2], 2.5), mu = F)$root)
CI.wd30.sigmasq <- CI.wd30.sigma^2
lines(range(sigmas^2), c*c(1,1), col = 2)

#wald
n <- dim(D)[1]
H.wd30.mu <- hessian(l.wd30.fun, mle.wd30.norm[1], sigma = mle.wd30.norm[2], data = D$wd30.centered)
V.wd30.mu <- as.numeric(-1/H.wd30.mu)
H.wd30.sigma <- hessian(l.wd30.fun, mle.wd30.norm[2], mu = mle.wd30.norm[1], data = D$wd30.centered)
V.wd30.sigma <- as.numeric(-1/H.wd30.sigma)
I11 <- n/mle.wd30.norm[2]^2;1/I11;V.wd30.mu# P. 59
I22 <- n/( 2 * mle.wd30.norm[2]^4 );1/I22;V.wd30.sigma#Det driller med denne :(
wald.wd30.mu <- mle.wd30.norm[1] + c(-1,1) * qnorm(1-alpha/2) * sqrt(1/I11)
wald.wd30.sigmasq <- mle.wd30.norm[2]^2 + c(-1,1) * qnorm(1-alpha/2) * sqrt(1/I22)

#All CIs of parameters
mle.wd30.norm.sq <- c(mle.wd30.norm[1], mle.wd30.norm[2]^2)
round(rbind(CI.pow, wald.pow, mle.pow.exp, CI.ws30.shape, wald.ws30.shape,
            CI.ws30.scale, wald.ws30.scale,mle.ws30.weib, CI.wd30.mu, wald.wd30.mu,
            CI.wd30.sigmasq, wald.wd30.sigmasq, mle.wd30.norm.sq), digits=5)


      