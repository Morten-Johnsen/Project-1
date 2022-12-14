rm(list = ls())
library(lubridate)
library(tidyverse)
library(reshape)
library(pheatmap)
library(openair)
library(PowerNormal)
library(sn)
library(gnorm)
library(emg)
library(numDeriv)
library(rlang)
library(extraDistr)
library(fitdistrplus)
library(SMPracticals)
#setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/")
setwd("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1")
D <- read.table("finance_data.csv", header=TRUE, sep=";", 
                as.is=TRUE)
## Dimensions of D (number of rows and columns)
dim(D)
##  Column/variable names
names(D)
## The first rows/observations
head(D)
## The last rows/observations
tail(D)
## Selected summary statistics
summary(D)
## Another type of summary of the dataset
str(D)

D$date <- as.Date(D$time)
D$year <- year(D$date)

#ggplot(D, aes(x = date, y = SLV, colour = factor(year))) + geom_point()
ggplot(D, aes(x = date, y = SLV)) + geom_point()
#plot(D$SLV)
ggplot(D, aes(x = SLV)) +
  geom_histogram(aes(y = ..density..), color = 'black')


#########################################################################################################
#########################################################################################################
source("testDistribution.R")
#########################################################################################################
#########################################################################################################
par <- nlminb(start = c(1,1), objective = testDistribution,
              distribution = "normal",
              x = D$SLV)
ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') +#color, fill
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = par$par[1], sd = par$par[2]), color='red')

plot(ecdf(D$SLV), verticals = T)
xseq <- seq(0.9*min(D$SLV), 1.1*max(D$SLV), length.out=100)
lines(xseq, pnorm(xseq, mean(D$SLV), sd(D$SLV)), col='red')
#plot(xseq, pnorm(xseq, mean(D$SLV), sd(D$SLV)), col='red')
qqnorm(D$SLV)
qqline(D$SLV)
#################################################
#################################################
lcauchyFUNC <- function(p, data){
  x0 <- p[1] #location R
  gam <- p[2] #scale R > 0
  return(-sum(dcauchy(x = data, location = x0, scale = gam, log = T)))
}
llstFUNC <- function(p, data){ #location-scale t-distribution
  return(-sum(dlst(x = data, df = p[1], mu = p[2], sigma = p[3], log = T)))
}

lsnFUNC <- function(p, data){ #skewed normal dist
  return(-sum(dsn(x = data, xi = p[1], omega = p[2], alpha = p[3], log = T)))
}
lgnFUNC <- function(p, data){ #symmetric generalized normal dist
  return(-sum(dgnorm(x = data, mu = p[1], alpha = p[2], beta = p[3], log = T)))
}
# lgnFUNC <- function(p, data){
#   mu <- p[1]
#   alpha <- p[2]
#   beta <- p[3]
#   -sum( log( beta/(2 * alpha * gamma(1/beta)) * exp(-(abs(data - mu)/alpha)^beta) ) )
# }
lasgnFUNC <- function(p, data){ #asymmetric generalized normal dist, when K = 0 has already been checked
  epsilon <- p[1]
  alpha <- p[2]
  kappa <- p[3]
  return(-sum( log(dnorm(x = -1/kappa * log(1 - kappa * (data - epsilon) / alpha) ) /
              (alpha - kappa * (data - epsilon)) ) ) )
}

lemgFUNC <- function(p, data){ #exponential modified gaussian dist
  return(-sum(demg(x = data, mu = p[1], sigma = p[2], lambda = p[3], log = T)))
}

par.cauchy <- nlminb(start = c(0,1), objective = lcauchyFUNC, data = D$SLV)
par.lst <- nlminb(start = c(1,1,1), objective = llstFUNC, data = D$SLV)
par.sn <- nlminb(start = c(1,1,1), objective = lsnFUNC, data = D$SLV)
par.gn <- nlminb(start = c(1,1,1), objective = lgnFUNC, data = D$SLV)
par.asgn <- nlminb(start = c(1,1,1), lower = c(-Inf, -Inf, 0), objective = lasgnFUNC, data = D$SLV)
par.emg <- nlminb(start = c(1,1,1), lower = c(-Inf, 1/1000, 1/1000), objective = lemgFUNC, data = D$SLV)

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(fun = dnorm, n = dim(D)[1], args = list(mean = par$par[1], sd = par$par[2]), aes(colour = "norm")) +
  stat_function(fun = dcauchy, n = dim(D)[1], args = list(location = par.cauchy$par[1],
                                                          scale = par.cauchy$par[2]), aes(colour = "cauchy")) +
  stat_function(fun = dlst, n = dim(D)[1], args = list(df = par.lst$par[1], mu = par.lst$par[2],
                                                       sigma = par.lst$par[3]), aes(colour = "lst")) +
  stat_function(fun = dsn, n = dim(D)[1], args = list(xi = par.sn$par[1], omega = par.sn$par[2],
                                                      alpha = par.sn$par[3]), aes(colour = "sn")) +
  stat_function(fun = dgnorm, n = dim(D)[1], args = list(mu = par.gn$par[1], alpha = par.gn$par[2],
                                                         beta = par.gn$par[3]), aes(colour = "gnorm")) +
  stat_function(fun = demg, n = dim(D)[1], args = list(mu = par.emg$par[1], sigma = par.emg$par[2],
                                                       lambda = par.emg$par[3]), aes(colour = "emg"))+
  scale_colour_manual(values = c("blue", "red", "yellow", "black", "grey", "purple"))+
  labs(colour = "Distribution")
#legend('topright', legend=c('normal', 'cauchy', 'power normal', 't'), col=c('red', 'blue', 'green', 'yellow'))

AIC.norm <- -2 * sum(dnorm(x = D$SLV, mean = par$par[1], sd = par$par[2], log = T))
+ 2 * length(par$par)
AIC.cauchy <- -2 * sum(dcauchy(x = D$SLV, location = par.cauchy$par[1],
                               scale = par.cauchy$par[2], log = T)) + 2 * length(par.cauchy$par)
AIC.lst <- -2 * sum(dlst(x=D$SLV, df = par.lst$par[1], mu = par.lst$par[2], sigma = par.lst$par[3], log = T)) + 2 * length(par.lst$par)
AIC.sn <- -2 * sum(dsn(x=D$SLV, xi = par.sn$par[1], omega = par.sn$par[2], alpha = par.sn$par[3], 
                       log = T)) + 2 * length(par.sn$par)
AIC.gn <- -2 * sum(dgnorm(x=D$SLV, mu = par.gn$par[1], alpha = par.gn$par[2], beta = par.gn$par[3], 
                          log = T)) + 2 * length(par.gn$par)
AIC.asgn <- -2 * sum( log(dnorm(x = -1/par.asgn$par[3] * log(1 - par.asgn$par[3] * (D$SLV - par.asgn$par[1]) / par.asgn$par[2]) ) /
                            (par.asgn$par[2] - par.asgn$par[3] * (D$SLV - par.asgn$par[1])) ) ) + 2 * length(par.asgn$par)
AIC.emg <- -2 * sum(demg(x=D$SLV, mu = par.emg$par[1], sigma = par.emg$par[2], lambda = par.emg$par[3], 
                         log = T)) + 2 * length(par.emg$par)


round(rbind(AIC.norm, AIC.cauchy, AIC.lst, AIC.sn, AIC.gn, AIC.asgn, AIC.emg), digits=5)

# boxplot(D$SLV)
n <- 100000
par(mfrow=c(2,3))
hist(rnorm(n, mean = par$par[1], sd = par$par[2]))
hist(D$SLV)
hist(rsn(n, xi = par.sn$par[1], omega = par.sn$par[2], alpha = par.sn$par[3]))
hist(rgnorm(n, mu = par.gn$par[1], alpha = par.gn$par[2], beta = par.gn$par[3]))
hist(rlst(n, df = par.lst$par[1], mu = par.lst$par[2], sigma = par.lst$par[3]))
hist(remg(n, mu = par.emg$par[1], sigma = par.emg$par[2], lambda = par.emg$par[3]))

# ggplot(D)+
#   #geom_histogram(aes(x = rgnorm(dim(D)[1], mu = par.gn$par[1], alpha = par.gn$par[2], beta = par.gn$par[3]), y=..density..), color='black') +
#   geom_histogram(aes(x = SLV, y= ..density..,), color='red') + #color, fill
#   geom_histogram(aes(x = rgnorm(dim(D)[1], mu = par.gn$par[1], alpha = par.gn$par[2], beta = par.gn$par[3]), y=..density..), color='black')
# lgamFUNC <- function(p, norm_data){
#   k <- p[1] #shape
#   beta <- p[2] # rate
#   return(-sum(dgamma(x = norm_data, shape = k, rate = beta, log = T)))
# }
# lbetaFUNC <- function(p, norm_data){
#   alpha <- p[1] #shape
#   beta <- p[2] #??shape?
#   -sum(dbeta(x = norm_data, shape1 = alpha, shape2 = beta, log = T))
# }
# 
# D$SLV.norm <- ( D$SLV - min(D$SLV) ) / (max(D$SLV) - min(D$SLV))
# par.gam <- nlminb(start = c(0.5, 0.5), objective = lgamFUNC, norm_data = D$SLV.norm)
# par.beta <- nlminb(start = c(0.5, 0.5), objective = lbetaFUNC, norm_data = D$SLV.norm)
# par.gn.norm <- nlminb(start = c(1,1,1), objective = lgnFUNC, data = D$SLV.norm)
# 
# ggplot(D)+
#   geom_histogram(aes(x = SLV.norm, y= ..density..,), color='black') + 
#   stat_function(fun = dgamma, n = dim(D)[1], args = list(shape = par.gam$par[1], scale = par.gam$par[2]), aes(colour = 'gamma')) +
#   stat_function(fun = dbeta, n = dim(D)[1], args = list(shape1 = par.beta$par[1], shape2 = par.beta$par[2]), aes(colour = 'beta')) +
#   stat_function(fun = dgnorm, n = dim(D)[1], args = list(mu = par.gn.norm$par[1],
#       alpha = par.gn.norm$par[2], beta = par.gn.norm$par[3]), aes(colour='gnorm')) +
#   labs(colour = 'Distributions')
# 
# -2 * sum(dgamma(x=D$SLV.norm, shape = par.gam$par[1], rate = par.gam$par[2], log = T)) + 2 * length(par.gam$par)
# -2 * sum(dbeta(x=D$SLV.norm, shape1 = par.beta$par[1], shape2 = par.beta$par[2], log = T)) + 2 * length(par.beta$par)
# -2 * sum(dgnorm(x = D$SLV.norm, mu = par.gn.norm$par[1], alpha = par.gn.norm$par[2],
#                 beta = par.gn.norm$par[3], log = T)) + 2 * length(par.gn.norm$par)

par(mfrow=c(1,3))
plot(ecdf(D$SLV), verticals = T, main = "Normal distribution")
xseq <- seq(0.9*min(D$SLV), 1.1*max(D$SLV), length.out=100)
#lines(xseq, pnorm(xseq, mean(D$SLV), sd(D$SLV)), col='red')
lines(xseq, pnorm(xseq, mean = par$par[1], sd = par$par[2]), col='pink')
plot(ecdf(D$SLV), verticals = T, main = "Generalized normal distribution")
lines(xseq, pgnorm(xseq, mu = par.gn$par[1], alpha = par.gn$par[2], beta = par.gn$par[3]), col='green')
plot(ecdf(D$SLV), verticals = T, main = "Location-scale t-distribution")
lines(xseq, plst(xseq, df = par.lst$par[1], mu = par.lst$par[2], sigma = par.lst$par[3]), col='red')

########################################################################################################
#Relevant keynumbers for the estimates:
#CIs of params for norm and gnorm
#CI for E[x] for norm and gnorm
alpha <- 0.05
c <- exp(-0.5 * qchisq(1-alpha, df = 1))
par(mfrow=c(1,2))
#likelihood-based CI for normal distribution
mle.norm <- par$par

fun.norm <- function(mu, sigma, data){#####
  prod(dnorm(x = data, mean = mu, sd = sigma, log = F) / 2) #to avoid inf-values
}

l.fun.norm <- function(mu, sigma, data){#####
  sum(dnorm(x = data, mean = mu, sd = sigma, log = T))
}

CIfun.norm <- function(y, data, mu = T){##### T from mean, F for sigma
  if(mu){
    sum(dnorm(x = data, mean = mle.norm[1], sd = mle.norm[2], log = T)) -
      sum(dnorm(x = data, mean = y, sd = mle.norm[2], log = T)) - 
      0.5 * qchisq(1-alpha, df = 1)
  } else {
    sum(dnorm(x = data, mean = mle.norm[1], sd = mle.norm[2], log = T)) -
      sum(dnorm(x = data, mean = mle.norm[1], sd = y, log = T)) - 
      0.5 * qchisq(1-alpha, df = 1) 
  }
}
mus <- seq(-0.01, 0.013, by = 0.00001)
mu.norm <- sapply(X = mus, FUN = fun.norm, sigma = mle.norm[2], data = D$SLV)
plot(mus, mu.norm/max(mu.norm), col = 1, type = "l", xlab = expression(paste(mu)),
     main = "Parameter value for mean for normal model of SLV")
CI.mu.norm <- c(uniroot(f = CIfun.norm, interval = c(-0.003, mle.norm[1]), data = D$SLV, mu = T)$root,
                uniroot(f = CIfun.norm, interval = c(mle.norm[1], 0.006), data = D$SLV, mu = T)$root)
lines(range(mus), c*c(1,1), col = 2)

sigmas <- seq(0.033, 0.060, by = 0.00001)
sigma.norm <- sapply(X = sigmas, FUN = fun.norm, mu = mle.norm[1], data = D$SLV)
plot(sigmas^2, sigma.norm/max(sigma.norm), col = 1, type = "l", xlab = expression(paste(sigma^2)),
     main = "Parameter value for var for normal model of SLV")
CI.sigma.norm <- c(uniroot(f = CIfun.norm, interval = c(0.033, mle.norm[2]), data = D$SLV, mu = F)$root,
                   uniroot(f = CIfun.norm, interval = c(mle.norm[2], 0.063), data = D$SLV, mu = F)$root)
CI.sigmasq.norm <- CI.sigma.norm^2
lines(range(sigmas^2), c*c(1,1), col = 2)

#wald
n <- dim(D)[1]
H.mu.norm <- hessian(l.fun.norm, mle.norm[1], sigma = mle.norm[2], data = D$SLV)
V.mu.norm <- as.numeric(-1/H.mu.norm)
H.sigma.norm <- hessian(l.fun.norm, mle.norm[2], mu = mle.norm[1], data = D$SLV)
V.sigma.norm <- as.numeric(-1/H.sigma.norm)
wald.mu.norm <- mle.norm[1] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.mu.norm)
wald.sigmasq.norm <- (mle.norm[2] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.sigma.norm))^2

mle.norm.sq <- c(mle.norm[1], mle.norm[2]^2)

#####Likelihood based CI for GENERALIZED normal distribution
mle.gn <- par.gn$par

lgnFUNC <- function(p, data){ #symmetric generalized normal dist
  return(-sum(dgnorm(x = data, mu = p[1], alpha = p[2], beta = p[3], log = T)))
}

fun.Gnorm <- function(mu, alpha, beta, data){#####
  prod(dgnorm(x = data, mu = mu, alpha = alpha, beta = beta, log = F) / 2)#to avoid inf-values
}

l.fun.Gnorm <- function(mu, alpha, beta, data){#####
  sum(dgnorm(x = data, mu = mu, alpha = alpha, beta = beta, log = T))
}

CIfun.Gnorm <- function(y, data, p = "mu"){##### T from mean, F for sigma
  if(p == "mu"){
    sum(dgnorm(x = data, mu = mle.gn[1], alpha = mle.gn[2], beta = mle.gn[3], log = T)) -
      sum(dgnorm(x = data, mu = y, alpha = mle.gn[2], beta = mle.gn[3], log = T)) - 
      0.5 * qchisq(1-alpha, df = 1)
  } else if(p == "alpha") {
    sum(dgnorm(x = data, mu = mle.gn[1], alpha = mle.gn[2], beta = mle.gn[3], log = T)) -
      sum(dgnorm(x = data, mu = mle.gn[1], alpha = y, beta = mle.gn[3], log = T)) - 
      0.5 * qchisq(1-alpha, df = 1) 
  } else {
    sum(dgnorm(x = data, mu = mle.gn[1], alpha = mle.gn[2], beta = mle.gn[3], log = T)) -
      sum(dgnorm(x = data, mu = mle.gn[1], alpha = mle.gn[2], beta = y, log = T)) - 
      0.5 * qchisq(1-alpha, df = 1) 
  }
}
###PROFILE likelihoods
par(mfrow = c(1,3))
#mu
mus.gn <- seq(-0.0055, 0.011, by = 0.00001)
mu.gn <- sapply(X = mus.gn, FUN = fun.Gnorm, alpha = mle.gn[2], beta = mle.gn[3], data = D$SLV)
plot(mus.gn, mu.gn/max(mu.gn), col = 1, type = "l", xlab = expression(paste(mu)),
     main = "Parameter value for location of generalized normal model of SLV")
CI.mu.gn <- c(uniroot(f = CIfun.Gnorm, interval = c(min(mus.gn), mle.gn[1]), data = D$SLV, p = "mu")$root,
                uniroot(f = CIfun.Gnorm, interval = c(mle.gn[1], max(mus.gn)), data = D$SLV, p = "mu")$root)
lines(range(mus.gn), c*c(1,1), col = 2)
#alpha
alphas <- seq(0.042, 0.062, by = 0.0001)
alpha.gn <- sapply(X = alphas, FUN = fun.Gnorm, mu = mle.gn[1], beta = mle.gn[3], data = D$SLV)
plot(alphas, alpha.gn/max(alpha.gn), col = 1, type = "l", xlab = expression(paste(alpha)),
     main = "Parameter value for scale for generalized normal model of SLV")
CI.alpha.gn <- c(uniroot(f = CIfun.Gnorm, interval = c(min(alphas), mle.gn[2]), data = D$SLV, p = "alpha")$root,
                   uniroot(f = CIfun.Gnorm, interval = c(mle.gn[2], max(alphas)), data = D$SLV, p = "alpha")$root)
lines(range(alphas), c*c(1,1), col = 2)
#beta
betas <- seq(1.19, 1.59, by = 0.001)
beta.gn <- sapply(X = betas, FUN = fun.Gnorm, mu = mle.gn[1], alpha = mle.gn[2], data = D$SLV)
plot(betas, beta.gn/max(beta.gn), col = 1, type = "l", xlab = expression(paste(beta)),
     main = "Parameter value for shape for generalized normal model of SLV")
CI.beta.gn <- c(uniroot(f = CIfun.Gnorm, interval = c(min(betas), mle.gn[3]), data = D$SLV, p = "beta")$root,
                uniroot(f = CIfun.Gnorm, interval = c(mle.gn[3], max(betas)), data = D$SLV, p = "beta")$root)
lines(range(betas), c*c(1,1), col = 2)

#Wald CIs
n <- dim(D)[1]
H.mu.gn <- hessian(l.fun.Gnorm, mle.gn[1], alpha = mle.gn[2], beta = mle.gn[3], data = D$SLV)
V.mu.gn <- as.numeric(-1/H.mu.gn)
H.alpha.gn <- hessian(l.fun.Gnorm, mle.gn[2], mu = mle.gn[1], beta = mle.gn[3], data = D$SLV)
V.alpha.gn <- as.numeric(-1/H.alpha.gn)
H.beta.gn <- hessian(l.fun.Gnorm, mle.gn[3], mu = mle.gn[1], alpha = mle.gn[3], data = D$SLV)
V.beta.gn <- as.numeric(-1/H.beta.gn)
wald.mu.gn <- mle.gn[1] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.mu.gn)
wald.alpha.gn <- mle.gn[2] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.alpha.gn)
wald.beta.gn <- mle.gn[3] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.beta.gn)

round( rbind(CI.mu.norm, wald.mu.norm, CI.sigmasq.norm, wald.sigmasq.norm, mle.norm.sq), digits = 5)
round( rbind(CI.mu.gn, wald.mu.gn, CI.alpha.gn, wald.alpha.gn, CI.beta.gn, wald.beta.gn), digits=5 );round(rbind (mle.gn), digits = 5)

#E[x] for gnorm, x_bar +- 1.96 * sd/sqrt(n)
E.gn <- mle.gn[1]
V.gn <- mle.gn[2]^2 * gamma(3/mle.gn[3]) / gamma(1/mle.gn[3])
CI.E.gn <- E.gn + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.gn / n)
rbind(E.gn);rbind(CI.E.gn)

temp1 <- paste("mu == ", mle.gn[1])
temp2 <- paste("alpha == ", mle.gn[2])
temp3 <- paste("beta == ", mle.gn[3])

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(fun = dgnorm, n = dim(D)[1], args = list(mu = par.gn$par[1], alpha = par.gn$par[2],
                                                         beta = par.gn$par[3]), color = 'red') +
  annotate( "text", x = 2.7/5*max(D$SLV), y = c(10.5, 10.0, 9.4), label = c(temp1,temp2,temp3), parse = T  ) +
  ggtitle("Generalized normal distribution and distribution of the weekly returns")

#####Likelihood based CI for location-scale t-distribution
mle.lst <- par.lst$par

#lstFUNC w/ NLL = T, negative log-likelihood using p's

#lstFUNC w/ NLL = F, log-likelihood using p's

lstFUNC <- function(df, mu, sigma, data, log = F){
  if(!log){
    return(prod(dlst(x = data, df = df, mu = mu, sigma = sigma, log = F) / 2)) #to avoid inf values
  } else {
    return(sum(dlst(x = data, df = df, mu = mu, sigma = sigma, log = T)))
  }
}  


CIfun.lst <- function(y, data, p = "df"){##### T from mean, F for sigma
  if(p == "df"){
    sum(dlst(x = data, df = mle.lst[1], mu = mle.lst[2], sigma = mle.lst[3], log = T)) -
      sum(dlst(x = data, df = y, mu = mle.lst[2], sigma = mle.lst[3], log = T)) - 
      0.5 * qchisq(1-alpha, df = 1)
  } else if(p == "mu") {
    sum(dlst(x = data, df = mle.lst[1], mu = mle.lst[2], sigma = mle.lst[3], log = T)) -
      sum(dlst(x = data, df = mle.lst[1], mu = y, sigma = mle.lst[3], log = T)) - 
      0.5 * qchisq(1-alpha, df = 1)
  } else { #p == "sigma"
    sum(dlst(x = data, df = mle.lst[1], mu = mle.lst[2], sigma = mle.lst[3], log = T)) -
      sum(dlst(x = data, df = mle.lst[1], mu = mle.lst[2], sigma = y, log = T)) - 
      0.5 * qchisq(1-alpha, df = 1)
  }
}
###PROFILE likelihoods
par(mfrow = c(1,3))
#df
dfs.lst <- seq(mle.lst[1]-2.5, mle.lst[1]+4.25, by = 0.0001)
df.lst <- sapply(X = dfs.lst, FUN = lstFUNC, mu = mle.lst[2], sigma = mle.lst[3], data = D$SLV, log = F)
plot(dfs.lst, df.lst/max(df.lst), col = 1, type = "l", xlab = "df",
     main = "Parameter value for df for location-scale t-distribution model of SLV")
CI.df.lst <- c(uniroot(f=CIfun.lst, interval = c(min(dfs.lst), mle.lst[1]), data = D$SLV, p = "df")$root,
               uniroot(f=CIfun.lst, interval = c(mle.lst[1], max(dfs.lst)), data = D$SLV, p = "df")$root)
lines(range(dfs.lst), c*c(1,1), col = 2)
#mu
mus.lst <- seq(mle.lst[2]-0.01, mle.lst[2]+0.01, by = 0.00001)
mu.lst <- sapply(X = mus.lst, FUN = lstFUNC, df = mle.lst[1], sigma = mle.lst[3], data = D$SLV, log = F)
plot(mus.lst, mu.lst/max(mu.lst), col = 1, type = "l", xlab = expression(paste(mu)),
     main = "Parameter value for location of location-scale t-distribution model of SLV")
CI.mu.lst <- c(uniroot(f = CIfun.lst, interval = c(min(mus.lst), mle.lst[2]), data = D$SLV, p = "mu")$root,
              uniroot(f = CIfun.lst, interval = c(mle.lst[2], max(mus.lst)), data = D$SLV, p = "mu")$root)
lines(range(mus.lst), c*c(1,1), col = 2)
#sigma
sigmas.lst <- seq(mle.lst[3]-0.01,mle.lst[3]+0.01, by = 0.0001)
sigma.lst <- sapply(X = sigmas.lst, FUN = lstFUNC, df = mle.lst[1], mu = mle.lst[2], data = D$SLV, log = F)
plot(sigmas.lst, sigma.lst/max(sigma.lst), col = 1, type = "l", xlab = expression(paste(sigma)),
     main = "Parameter value for scale for generalized normal model of SLV")
CI.sigma.lst <- c(uniroot(f = CIfun.lst, interval = c(min(sigmas.lst), mle.lst[3]), data = D$SLV, p = "sigma")$root,
                 uniroot(f = CIfun.lst, interval = c(mle.lst[3], max(sigmas.lst)), data = D$SLV, p = "sigma")$root)
lines(range(sigmas.lst), c*c(1,1), col = 2)

#Wald CIs
n <- dim(D)[1]
H.df.lst <- hessian(lstFUNC, mle.lst[1], mu = mle.lst[2], sigma = mle.lst[3], data = D$SLV, log = T)
V.df.lst <- as.numeric(-1/H.df.lst)
H.mu.lst <- hessian(lstFUNC, mle.lst[2], df = mle.lst[1], sigma = mle.lst[3], data = D$SLV, log = T)
V.mu.lst <- as.numeric(-1/H.mu.lst)
H.sigma.lst <- hessian(lstFUNC, mle.lst[3], df = mle.lst[1], mu = mle.lst[2], data = D$SLV, log = T)
V.sigma.lst <- as.numeric(-1/H.sigma.lst)
wald.df.lst <- mle.lst[1] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.df.lst)
wald.mu.lst <- mle.lst[2] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.mu.lst)
wald.sigma.lst <- mle.lst[3] + c(-1,1) * qnorm(1-alpha/2) * sqrt(V.sigma.lst)

round( rbind(CI.mu.norm, wald.mu.norm, CI.sigmasq.norm, wald.sigmasq.norm, mle.norm.sq), digits = 5)
round( rbind(CI.mu.gn, wald.mu.gn, CI.alpha.gn, wald.alpha.gn, CI.beta.gn, wald.beta.gn), digits=5 );round(rbind (mle.gn), digits = 5)
round( rbind(CI.df.lst, wald.df.lst, CI.mu.lst, wald.mu.lst, CI.sigma.lst, wald.sigma.lst), digits=5);round( rbind(mle.lst), digits=5)

temp1 <- paste("df == ", round(mle.lst[1], digits=2))
temp2 <- paste("mu == ", round(mle.lst[2], digits=5))
temp3 <- paste("sigma == ", round(mle.lst[3], digits=4))

ggplot(D)+
  geom_histogram(aes(x = SLV, y= ..density..,), color='black') + #color, fill
  stat_function(fun = dlst, n = dim(D)[1], args = list(df = par.lst$par[1], mu = par.lst$par[2],
                                                         sigma = par.lst$par[3]), color = 'red') +
  annotate( "text", x = 2.7/5*max(D$SLV), y = c(10.5, 10.0, 9.5), label = c(temp1,temp2,temp3), parse = T  ) +
  ggtitle("Location-scale t-distribution and distribution of the weekly returns")


qqnorm2 <- function (y, line = FALSE, ...) #function is taken from package SMPracticals and was originally for qqexp :)
{
  y <- y[!is.na(y)]
  n <- length(y)
  x <- qnorm(c(1:n)/(n + 1))
  m <- mean(y)
  ylim <- c(min(y), max(y))
  qqplot(x, y, xlab = "normal plotting position", ylim = ylim, 
         ylab = "Ordered sample", ...)
  if (line) 
    abline(0, m, lty = 2)
  invisible()
}
qqlst <- function (y, line = FALSE, ...)
{
  y <- y[!is.na(y)]
  n <- length(y)
  x <- qlst(c(1:n)/(n + 1), df=mle.lst[1])
  m <- mean(y)
  ylim <- c(min(y), max(y))
  qqplot(x, y, xlab = "location-scale t-distribution plotting position", ylim = ylim, 
         ylab = "Ordered sample", ...)
  if (line) 
    abline(0, m, lty = 2)
  invisible()
}
par(mfrow=c(1,2))
qqnorm(D$SLV)
qqline(D$SLV)
qqlst(D$SLV)
qqline(D$SLV)
