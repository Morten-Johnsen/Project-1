rm(list = ls())
library(lubridate)
library(tidyverse)
library(reshape)
library(stringr)
library(pheatmap)
library(numDeriv)

#Retrieve data from the descriptive statistics script
##### TILPAS DEN HER TIL DEN COMPUTER DER BRUGES
if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/")
} else {
  setwd("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1")
}

load("dataWindPower.Rdata")
source("testDistribution.R")

#### WIND POWER ####
#Fitting models to wind power:
par.exp <- nlminb(start = 0.2, objective = testDistribution,
                  distribution = "exponential",
                  x = D$pow.obs.norm)

par.exp$objective

par.beta <- nlminb(start = c(2,5)
                   , objective = testDistribution
                   , distribution = "beta"
                   , x = D$pow.obs.norm
                   , lower = c(0,0.8))
par.beta$objective

par.gamma <- nlminb(start = c(2,5)
                    ,objective = testDistribution
                    ,distribution = "gamma"
                    ,x = D$pow.obs.norm)
par.gamma$objective

#Box-Cox transformation of pow.obs.norm
#Examine different transformations and the achieved fit when fitting a normal
#distribution
# 
# lambda <- c(0.0, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
# pal <- palette.colors(length(lambda))
# BoxCoxPlot <- list()
# for (i in 1:length(lambda)){
#   xData <- 2*log(D$pow.obs.norm^lambda[i]/(1-D$pow.obs.norm)^(1-lambda[i]))
#   n <- nlminb(start = c(-1,1)
#               , objective = testDistribution
#               , x = xData
#               , distribution = "normal")
#   
#   simData <- rnorm(n = length(D$pow.obs.norm), mean = n$par[1], sd = n$par[2])
#   D$sim <- simData
#   
#   print(c(n$par, lambda[i], mean(D$sim), sd(D$sim)))
#   BoxCoxPlot[[paste0(lambda[i])]] <- ggplot(D)+
#     geom_histogram(aes(x = sim, y = ..density..)
#                    , colour = "white"
#                    , alpha = 0.5
#                    , fill = "black")+
#     geom_histogram(aes(x = 2*log(pow.obs.norm^lambda[i]/(1-pow.obs.norm)^(1-lambda[i])), y = ..density..)
#                    , colour = "white"
#                    , alpha = 0.6
#                    , fill = pal[[i]])+
#     labs(x = "BoxCox(pow.obs.norm)")+
#     ggtitle(paste0("BoxCox Transformation of Windpower, ", expression(lambda), " = ", lambda[i]))+
#     geom_text(x = -5, y = 0.23, label = paste0("NLL = ", round(n$objective,2)))+
#     geom_text(x = -5, y = 0.2, label = paste0("(mean, sd) = (", round(n$par[1],2), ",", round(n$par[2],2), ")"))+
#     stat_function(fun = dnorm, n = length(D$pow.obs.norm), args = list(mean = n$par[1], sd = n$par[2]), alpha = 0.6, colour = "black")
# }
# library(gridExtra)
# grid.arrange(grobs = BoxCoxPlot)


#Compare distributions
for (i in 1:1){
  D$simdata <- rbeta(length(D$pow.obs.norm), shape1 = par.beta$par[1]
        ,shape2 = par.beta$par[2])
  b <- ggplot(D)+
    geom_histogram(aes(x = pow.obs.norm, y =..density..), colour = "white", alpha = 0.6)+
    geom_histogram(aes(x = simdata, y =..density..), alpha = 0.2, fill = "blue")+
    theme_bw()+
    ylim(c(0,10))+
    stat_function(fun = dbeta, n = length(D$pow.obs.norm), args = list(shape1 = par.beta$par[1],shape2 = par.beta$par[2]))
  show(b)
  #Sys.sleep(2)
}

for (i in 1:1){
  D$simdata <- rexp(length(D$pow.obs.norm), rate = par.exp$par)
  g <- ggplot(D)+
    geom_histogram(aes(x = pow.obs.norm, y = ..density..))+
    geom_histogram(aes(x = simdata, y = ..density..)
                   , alpha = 0.2, fill = "blue")+
    theme_bw()+
    stat_function(fun = dexp, n = length(D$pow.obs.norm), args = list(rate = par.exp$par))
  show(g)
  #Sys.sleep(2)
}

for (i in 1:1){
  D$simdata <- rgamma(length(D$pow.obs.norm), shape = par.gamma$par[1], rate = par.gamma$par[2])
  g <- ggplot(D)+
    geom_histogram(aes(x = pow.obs.norm, y = ..density..))+
    geom_histogram(aes(x = simdata, y = ..density..)
                   , alpha = 0.2, fill = "blue")+
    theme_bw()+
    ylim(c(0,10))+
    stat_function(fun = dgamma, n = length(D$pow.obs.norm), args = list(shape = par.gamma$par[1], rate = par.gamma$par[2]))
  show(g)
  #Sys.sleep(2)
}
D <- D %>%
  select(-simdata)
#Det er tydeligvist beta modellen der er bedst. Det er også det der fremgår af likelihoods,
#så alt er godt.

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

#Vi skal enten sammenligne med andre modeller eller referere til relevant litteratur
#som siger at man bruger weibull tætheden og sådan er det.

#### WIND DIRECTION ####
library(circular)

nll.wrappedNormal <- function(p,x){
  nll <- -sum(log(dwrappednormal(x, mu = circular(p[1]), rho = NULL, sd = p[2])))
  return(nll)
}

nll.wrappedCauchy <- function(p,x){
  nll <- -sum(log(dwrappedcauchy(x, mu = circular(p[1]), rho = p[2])))
  return(nll)
}

nll.vonMises <- function(p,x){
  nll <- -sum(dvonmises(x, mu = circular(p[1]), kappa = p[2], log = T))
  return(nll)
}

wrapped.par <- nlminb(start = c(2,1), objective = nll.wrappedNormal, x = D$wd30)
wrapped.cauc.par <- nlminb(start = c(1,1/10000), lower = c(-Inf, 1/10000), upper = c(Inf, 1),
                           objective = nll.wrappedCauchy, x = D$wd30)
wrapped.vonMises <- nlminb(start = c(0,1), objective = nll.vonMises, x = D$wd30, lower = c(-1000, 0))
#centrerer fordelingen omkring 3/2pi
#D$wd30.centered <- D$wd30 - pi/2; D$wd30.centered[D$wd30.centered < 0] = D$wd30.centered[D$wd30.centered < 0] + 2*pi
#par.wd30 <- nlminb(start = c(4,4),
#                   objective = testDistribution,
#                   x = D$wd30.centered,
#                   distribution = "normal")

ggplot(D)+
  theme_bw()+
  #geom_density(aes(x = wd30.centered, y = ..density..), alpha = .8, colour = "white", fill = "red", colour = "white")+
  geom_histogram(aes(x = wd30, y = ..density..), colour = "white", alpha = .4, fill = "blue", bins = 20)+
  scale_x_continuous(breaks = c(0,pi/2,pi,3/2*pi,2*pi)
                     , labels =c("0", "pi/2", "pi", "3/2pi", "2pi"))+
  #stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = par.wd30$par[1], sd = par.wd30$par[2]))+
  stat_function(fun = dwrappednormal, n = dim(D)[1], args = list(mu = wrapped.par$par[1], sd = wrapped.par$par[2]), aes(colour = "Wrapped Normal"))+
  stat_function(fun = dwrappedcauchy, n = dim(D)[1], args = list(mu = wrapped.cauc.par$par[1], rho = 0.36), aes(colour = "Wrapped Cauchy"))+
  stat_function(fun = dvonmises, n = dim(D)[1], args = list(mu = wrapped.vonMises$par[1], kappa = wrapped.vonMises$par[2]), aes(colour = "Von Mises"))+
  labs(x = "Wind Direction", colour = "")+
  scale_colour_manual(values = c("yellow", "red", "black"))

#Calculate AICs
print(paste0("AIC wrapped normal: ", round(-2*log(wrapped.par$objective)+2*2,4), "|"
            ,"AIC wrapped cauchy: ", round(-2*log(wrapped.cauc.par$objective)+2,4), "|"
            ,"AIC von Mises: "     , round(-2*log(wrapped.vonMises$objective)+2*2,4)))


## CI ## WIND POWER
par(mfrow=c(1,1))
alpha <- 0.05
c <- exp(-0.5 * qchisq(1-alpha, df = 1))
#likelihood-based
mle.pow <- par.beta$par

pow.fun <- function(shape1, shape2, data){
  return( prod( dbeta(x = data, shape1 = shape1, shape2 = shape2, log = F) ) )
}

l.pow.fun <- function(shape1, shape2, data){
  return( sum( dbeta(x = data, shape1 = shape1, shape2 = shape2, log = T) ) )
}

CIfun.pow <- function(y, first = T){##### T for shape, F for scale
  if(first){
    return( sum( dbeta(x = D$pow.obs.norm, shape1 = mle.pow[1], shape = mle.pow[2], log = T) ) -
      sum( dbeta(x = D$pow.obs.norm, shape1 = y, shape2 = mle.pow[2], log = T) ) - 
      0.5 * qchisq(1-alpha, df = 1) )
  } else {
    return( sum( dbeta(x = D$pow.obs.norm, shape1 = mle.pow[1], shape = mle.pow[2], log = T) ) -
      sum( dbeta(x = D$pow.obs.norm, shape1 = mle.pow[1], shape2 = y, log = T) ) - 
      0.5 * qchisq(1-alpha, df = 1) ) 
  }
}

par(mfrow=c(1,2))
shape1s <- seq(0, 1, by = 0.01)
pow1 <- sapply(X = shape1s, FUN = pow.fun, data = D$pow.obs.norm, shape2 = mle.pow[2])
plot(shape1s, pow1/max(pow1), col = 1, type = "l", xlab = expression(paste(alpha)),
     main = "Parameter value shape1 for beta model of power production")
CI.pow1 <- c(uniroot(f = CIfun.pow, interval = c(0, mle.pow[1]), first = T)$root,
            uniroot(f = CIfun.pow, interval = c(mle.pow[1], 1), first = T)$root)
lines(range(shape1s), c*c(1,1), col = 2)

shape2s <- seq(1, 2, by = 0.01)
pow2 <- sapply(X = shape2s, FUN = pow.fun, data = D$pow.obs.norm, shape1 = mle.pow[1])
plot(shape2s, pow2/max(pow2), col = 1, type = "l", xlab = expression(paste(beta)),
     main = "Parameter value shape2 for beta model of power production")
CI.pow2 <- c(uniroot(f = CIfun.pow, interval = c(1, mle.pow[2]), first = F)$root,
             uniroot(f = CIfun.pow, interval = c(mle.pow[2], 2), first = F)$root)
lines(range(shape2s), c*c(1,1), col = 2)

#wald
n <- dim(D)[1]
H.pow.shape1 <- hessian(l.pow.fun, mle.pow[1], shape2 = mle.pow[2], data = D$pow.obs.norm)
V.pow.shape1 <- as.numeric(-1/H.pow.shape1)
H.pow.shape2 <- hessian(l.pow.fun, mle.pow[2], shape1 = mle.pow[1], data = D$pow.obs.norm)
V.pow.shape2 <- as.numeric(-1/H.pow.shape2)
wald.pow.shape1 <- mle.pow[1] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.pow.shape1)
wald.pow.shape2 <- mle.pow[2] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.pow.shape2)

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
plot(shapes, ws30.shape/max(ws30.shape), col = 1, type = "l", xlab = "shape, k",
     main = "Parameter value for shape for weibull model of wind speed")
CI.ws30.shape <- c(uniroot(f = CIfun.ws30, interval = c(1, mle.ws30.weib[1]), shape = T)$root,
                   uniroot(f = CIfun.ws30, interval = c(mle.ws30.weib[1], 3.5), shape = T)$root)
lines(range(shapes), c*c(1,1), col = 2)

scales <- seq(7, 12, by = 0.1)
ws30.scale <- sapply(X = scales, FUN = ws30.fun, shape = mle.ws30.weib[1], data = D$ws30)
plot(scales, ws30.scale/max(ws30.scale), col = 1, type = "l", xlab = expression(paste("scale, ", lambda)),
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

####### E[x] ####### START 1
k <- par.ws30$par[1]
lambda <- par.ws30$par[2]

#E[X]
print("E[X]=")
lambda*gamma(1+1/k)
#Var[X]
print("Var[X]=")
lambda^2 * (gamma(1+2/k) - (gamma(1+1/k))^2)

#function for numerical estimate of the shape parameter based on mean and variance.
k.func <- function(k,mu,var){
  var.est <- mu^2 * (gamma(1+2/k) - gamma(1+1/k)^2)/gamma((k+1)/k)^2
  return((var - var.est)^2)
}

nll.wei <- function(theta, x){
  mu <- theta[1]
  var <- theta[2]
  k <- nlminb(start = 1, objective = k.func, mu = mu, var = var)$par
  nll <- -sum(dweibull(x, shape = k, scale = mu/gamma((k+1)/k), log = T))
  return(nll)
}

par.repar.wei <- nlminb(start = c(1,0.5), objective = nll.wei, x = D$ws30, lower = c(0,0))

#E[X]
par.repar.wei$par[1]
#Var[X]
par.repar.wei$par[2]

####### E[x] ####### END 1
####### E[x] ####### START 2
#We know that the expected value is given by: weibull_E[x] = lambda * gamma(1 + 1/k). So: lambda ) E[x]/gamma(1+1/k)
nll.repar.ws30 <- function(p, data){
  k <- p[1]
  theta <- p[2]
  return( -sum(dweibull(x = data, shape = k, scale = theta / gamma(1 + 1/k), log = T)) )
}
repar.E.ws30 <- nlminb(start = c(1,1), objective = nll.repar.ws30, data = D$ws30) #E[x] = par.E.ws30[2] 

### equiv ^ to the code above ###

repar.ws30.fun <- function(shape, theta, data){#####
  prod(dweibull(x = data, shape = shape, scale = theta / gamma(1 + 1/shape), log = F)*2)#to not get full zeros
}

repar.l.ws30.fun <- function(shape, theta, data){#####
  sum(dweibull(x = data, shape = shape, scale = theta / gamma(1 + 1/shape), log = T))
}

CIfun.repar.ws30 <- function(y, data, mles){##### T for shape, F for scale
  k <- mles[1]
  theta <- mles[2]
  sum(dweibull(x = data, shape = k, scale = theta / gamma(1 + 1/k), log = T)) -
    sum(dweibull(x = data, shape = k, scale = y / gamma(1 + 1/k), log = T)) - 
    0.5 * qchisq(1-alpha, df = 1) 
}

Es <- seq(8, 10.4, by = 0.001)
ws30.E <- sapply(X = Es, FUN = repar.ws30.fun, shape = repar.E.ws30$par[1], data = D$ws30)
plot(Es, ws30.E/max(ws30.E), col = 1, type = "l", xlab = expression(paste("Expected value E[x] or ", theta)),
     main = "Expected value for weibull model of wind speed")
CI.ws30.E <- c(uniroot(f = CIfun.repar.ws30, interval = c(min(Es), repar.E.ws30$par[2]), data = D$ws30, mles = repar.E.ws30$par)$root,
               uniroot(f = CIfun.repar.ws30, interval = c(repar.E.ws30$par[2], max(Es)), data = D$ws30, mles = repar.E.ws30$par)$root)
lines(range(Es), c*c(1,1), col = 2)

H.ws30.E <- hessian(repar.l.ws30.fun, repar.E.ws30$par[2], shape = repar.E.ws30$par[1], data = D$ws30)
V.ws30.E <- as.numeric(-1/H.ws30.E)
wald.ws30.E <- repar.E.ws30$par[2] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.ws30.E)

round( rbind(CI.ws30.E, wald.ws30.E), digits=3);round( rbind(repar.E.ws30$par[2]), digits=3)
####### E[x] ####### END 2

## CI ## WIND DIRECTION
par(mfrow=c(1,2))
#likelihood-based
mle.wd30 <- wrapped.cauc.par$par

wd30.fun <- function(mu, rho, data){#####
  prod(dwrappedcauchy(x = data, mu = mu, rho = rho))
}

l.wd30.fun <- function(mu, rho, data){#####
  sum( log( dwrappedcauchy(x = data, mu = mu, rho = rho) ) )
}

CIfun.wd30 <- function(y, mu = T){##### T from mean, F for sigma
  if(mu){
    return( sum( log( dwrappedcauchy(x = D$wd30, mu = mle.wd30[1], rho = mle.wd30[2]) ) ) -
      sum( log( dwrappedcauchy(x = D$wd30, mu = y, rho = mle.wd30[2]) ) ) -
      0.5 * qchisq(1-alpha, df = 1) )
  } else {
    return( sum( log( dwrappedcauchy(x = D$wd30, mu = mle.wd30[1], rho = mle.wd30[2]) ) ) - 
      sum( log( dwrappedcauchy(x = D$wd30, mu = mle.wd30[1], rho = y) ) ) -
      0.5 * qchisq(1-alpha, df = 1) )
  }
}

mus <- seq(-2.5, -1, by = 0.01)
wd30.mu <- sapply(X = mus, FUN = wd30.fun, rho = mle.wd30[2], data = D$wd30)
plot(mus, wd30.mu/max(wd30.mu), col = 1, type = "l", xlab = expression(paste("peak, ", mu)),
     main = "Parameter value for peak for wrapped cauchy model of wind direction")
CI.wd30.mu <- c(uniroot(f = CIfun.wd30, interval = c(-2.5, mle.wd30[1]), mu = T)$root,
                uniroot(f = CIfun.wd30, interval = c(mle.wd30[1], -1), mu = T)$root)
lines(range(mus), c*c(1,1), col = 2)

rhos <- seq(0, 0.5, by = 0.005)
wd30.rho <- sapply(X = rhos, FUN = wd30.fun, mu = mle.wd30[1], data = D$wd30)
plot(rhos, wd30.rho/max(wd30.rho), col = 1, type = "l", xlab = expression(paste("concentration, ", rho)),
     main = "Parameter value for concentration factor for wrapped cauchy model of wind direction")
CI.wd30.rho <- c(uniroot(f = CIfun.wd30, interval = c(0, mle.wd30[2]), mu = F)$root,
                   uniroot(f = CIfun.wd30, interval = c(mle.wd30[2], 0.5), mu = F)$root)
lines(range(rhos), c*c(1,1), col = 2)

#wald
n <- dim(D)[1]
H.wd30.mu <- hessian(l.wd30.fun, mle.wd30[1], rho = mle.wd30[2], data = D$wd30)
V.wd30.mu <- as.numeric(-1/H.wd30.mu)
H.wd30.rho <- hessian(l.wd30.fun, mle.wd30[2], mu = mle.wd30[1], data = D$wd30)
V.wd30.rho <- as.numeric(-1/H.wd30.rho)
wald.wd30.mu <- mle.wd30[1] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.wd30.mu)
wald.wd30.rho <- mle.wd30[2] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.wd30.rho)

#All CIs of parameters
round( rbind( CI.pow1, wald.pow.shape1, CI.pow2, wald.pow.shape2, mle.pow,
              CI.ws30.shape, wald.ws30.shape, CI.ws30.scale, wald.ws30.scale, mle.ws30.weib, 
              CI.wd30.mu, wald.wd30.mu, CI.wd30.rho, wald.wd30.rho, mle.wd30 ), digits = 3 )

#CI.E.pow.obs <- 1/par$par + c(-1,1) * qnorm(1-alpha/2) * sqrt(1/par$par^2) / dim(D)[1]
CI.E.pow.obs <- mean(D$pow.obs.norm) + c(-1,1) * qnorm(1-alpha/2) * sd(D$pow.obs.norm) / dim(D)[1]
#par.ws30$par[2]*gamma(1+1/par.ws30$par[1]) #mean = lambda * Gamma(1 + 1/k); lambda = scale, k = shape
E.ws30 <- par.ws30$par[2]*gamma(1+1/par.ws30$par[1])
V.ws30 <- par.ws30$par[2]^2*( gamma(1+2/par.ws30$par[1]) - (gamma(1+1/par.ws30$par[1]))^2)
#CI.E.ws30 <- E.ws30 + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.ws30) / dim(D)[1] #according to Central Limit Theorem
CI.E.ws30 <- mean(D$ws30) + c(-1,1) * qnorm(1-alpha/2) * sd(D$ws30) / dim(D)[1]
#CI.E.wd30 <- par.wd30$par[1] + c(-1,1) * qnorm(1-alpha/2) * par.wd30$par[2] / dim(D)[1] #according to Central Limit Theorem
CI.E.wd30 <- mle.wd30[1] + c(-1,1) * qnorm(1-alpha/2) * sd(D$wd30) / dim(D)[1] #mean(D$wd30) instead gives another result

round(rbind(c(CI.E.pow.obs[1], 1/par.exp$par, CI.E.pow.obs[2]) , c(CI.E.ws30[1], E.ws30, CI.E.ws30[2]) , c(CI.E.wd30[1], mle.wd30[1], CI.E.wd30[2])), digits=5)
##########

par(mfrow=c(1,3))

temp1 <- paste("alpha == ", mle.pow[1]) #par.beta$par[1]
temp2 <- paste("beta == ", mle.pow[2]) #par.beta$par[2]
temp <- c(temp1, temp2)

ggplot(D)+
  geom_histogram(aes(x = pow.obs.norm, y = ..density..), colour='white', alpha=0.6, bins=30)+
  theme_bw()+
  stat_function(fun = dbeta, n = dim(D)[1], args = list(shape1 = mle.pow[1], shape2 = mle.pow[2]))+
  ylim(c(0,10))+
  annotate( "text", x = 4/5*max(D$pow.obs.norm), y = c(9.5, 9.0), label = temp, parse = T  ) +
  ggtitle("Beta distribution and distribution of normalized power production")

temp1 <- paste("k == ", mle.ws30.weib[1]) #par.ws30$par[1]
temp2 <- paste("lambda == ", mle.ws30.weib[2]) #par.ws30$par[2]
temp <- c(temp1, temp2)

ggplot(D)+
  geom_histogram(aes(x = ws30, y = ..count../sum(..count..)) , colour = "white", alpha=0.6, bins = 30)+
  theme_bw()+
  stat_function(fun = dweibull, n = dim(D)[1], args = list(shape = par.ws30$par[1], scale = par.ws30$par[2])) +
  annotate( "text", x = 4/5*max(D$ws30), y = c(0.09,0.085), label = temp, parse = T  ) +
  ggtitle("Weibull distribution and distribution of the wind speed")

temp1 <- paste("mu ==", mle.wd30[1]) #wrapped.cauc.par$par[1]
temp2 <- paste("rho ==", mle.wd30[2]) #wrapped.cauc.par$par[2]
temp <- c(temp1, temp2)

ggplot(D)+
  theme_bw()+
  geom_histogram(aes(x = wd30, y = ..density..), colour = "white", alpha = 0.6, bins = 20)+
  scale_x_continuous(breaks = c(0,pi/2,pi,3/2*pi,2*pi)
                     , labels =c("0", "pi/2", "pi", "3/2pi", "2pi"))+
  stat_function(fun = dwrappedcauchy, n = dim(D)[1], args = list(mu = wrapped.cauc.par$par[1], rho = 0.36))+
  annotate( "text", x = 1/5*max(D$wd30), y = c(0.35, 0.325), label = c(temp1,temp2), parse = T  ) +
  ggtitle("Wrapped cauchy distribution and distribution of the wind direction")

