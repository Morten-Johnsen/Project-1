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
library(numDeriv)
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

#Our model from assignment 2
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
D$x_pred <- 1/( exp(y_p*lambda)+1 )^(1/lambda) * exp(y_p)

ggplot(data = D)+
  geom_point(aes(x=ws30, y=pow.obs.norm, colour = "Observed power production (Transformed)"), size = 4, alpha = 0.6)+
  geom_line(aes(x=ws30, y=x_pred, colour="Inverse Transformation of the Normal model"), size = 2)+
  #geom_line(aes(x=ws30, y=betaModel.pow, colour="The directly fitted beta model"), size = 2)+
  scale_colour_manual(values = c("blue","red"))+
  labs(x = "Wind speeds", y="Norm power obs", colour ="")+
  theme_bw()+
  labs(colour = "")+
  ggtitle("Model fits")+
  theme(legend.background = element_rect(colour = "black"), legend.position = "top")

## Calculate the model residuals in original pow.obs/5000 domain.
e_vec <- D$pow.obs.norm - D$x_pred

# Show fit vs residuals
fit <- ggplot(data = D)+
  theme_bw()+ #THEME
  geom_point(aes(x = ws30, y = pow.obs.norm))+ #OBSERVED DATA
  geom_line(aes(x = ws30, y = x_pred), colour = "red", size = 1)+ #MODEL FIT (LINE)
  labs(x = "Wind speed [m/s]", y = "Generated Power [1/5000 kW]")+ #AXIS LABELS
  ggtitle("Fit") #TITLE
res <- ggplot(data = D)+
  theme_bw()+ #THEME
  geom_point(aes(x = ws30, y = e_vec))+ #ACTUAL PLOT
  labs(x = "Wind speed [m/s]", y = "Residuals [1/5000 kW]")+ #AXIS LABELS
  geom_line(aes(x = ws30, y = 0), size = 1, colour = "red")+ #LINE GOING THROUGH ZERO
  ggtitle("Residuals") #TITLE
grid.arrange(fit,res,nrow = 1)

################################################################################################
###Analysis of auto-correlation
##TASK 1 of 7...
#Residuals in the transformed domain are:
e1 <- ws.and.ws.squared$residuals[1:( length(ws.and.ws.squared$residuals)-1 )] #for e_1 - e_n-1
e2 <- ws.and.ws.squared$residuals[2:length(ws.and.ws.squared$residuals)] #for e_2 - e_n

m_res <- matrix(c(e1,e2), ncol=2, byrow=F)

#1.1
##acf plot of model residuals in the x domain
acf(e_vec)
##arrange model residuals in a matrix consisting of [e_1...e_n-1; e_2...e_n]
n <- length(e_vec)
e <- matrix(cbind("e_One_to_nMinusOne" = e_vec[-n], "e_Two_to_n" = e_vec[-1]), ncol = 2)
#1.2
##calculate rho and sigma
rho_matrix <- cor(e)#?
(phi <- sum(e[,2]*e[,1])/sum(e[,1]^2))#?
#1/(n-1) * sum (e_{i+1} - phi * e_{i})
(sigma.sq <- 1/(length(e_vec)-1) * sum((e[,2] - phi * e[,1])^2))#?


e_corrected <- e[,2]-phi*e[,1] #?
acf(e_corrected) #?

AR1m <- arima(e_vec, order=c(1,0,0))
acf(AR1m$residuals)
#### FROM LECT 10 ####
## Moment estimates
(phi.hat <- acf(e_vec,plot=FALSE)$acf[2]) ## ACF in lag 1 (index 1 is itself)
cor(e_vec[-1], e_vec[-n]) #same calculation but by hand.
(sigma <- var(e_vec)*(1-phi.hat^2))


## Likelihood estimation
nll <- function(theta,y){
  n <- length(y) - 1
  sigma <- theta[1]
  phi <- theta[2]
  return(n/2 * log(sigma) + 1/(2*sigma) * sum((y[-1]-phi*y[-(n+1)])^2))
}

## MLE
(opt <- nlminb(c("sigma"=1,"phi"=1/2),nll,y=e_vec,lower=c(0.001,-0.999),upper=c(Inf,0.999)))

## Check MLEs
(phi <- sum(e_vec[-1]*e_vec[-n])/sum(e_vec[-n]^2))
(sigma.sq <- 1/(n-1) * sum((e_vec[-1]-phi * e_vec[-n])^2))

## Did it help?
e_corrected <- e_vec[-1]-opt$par[2]*e_vec[-n]
acf(e_corrected)

## standard errors
V <- solve(hessian(nll,opt$par,y=e_vec))
(se <- sqrt(diag(V)))

## Profile likelihood for phi
llp.phi <- function(phi,x,x0){
  n <- length(x) - 1
  return(-n/2 * log(sum((x[-1]-phi*x[-(n+1)])^2)))
}

## Plot profile likelihood
phi_interval <- seq(max(opt$par[2]-3*se[2],-1),
                    min(opt$par[2]+3*se[2],1),
                    length=200)
sigma.sq_interval <- seq(max(opt$par[1]-3*se[1],-1),
                         min(opt$par[1]+3*se[1],1),
                         length=200)

llp <- sapply(phi_interval,llp.phi,x=e_vec[-1],x0)

plot(phi_interval,exp(llp-max(llp)),type="l")
lines(range(phi_interval),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)

#### 2.1 - MLEs and confidence intervals ####
#phi confidence interval
#phi MLE
phi.hat.opt
## pf based
phi.lower.pf <- min(phi_interval[exp(llp-max(llp)) >= exp(-qchisq(0.95,df=1)/2)])
phi.upper.pf <- max(phi_interval[exp(llp-max(llp)) >= exp(-qchisq(0.95,df=1)/2)])

## wald based
sd.phi <- se[2]
phi.hat.opt <- opt$par[2]

phi.lower.wald <- phi.hat.opt - sd.phi*qnorm(0.975)
phi.upper.wald <- phi.hat.opt + sd.phi*qnorm(0.975)

#sigma.sq confidence interval
llp.sigma.sq <- function(sigma,phi,y){
  #funktionen er taget fra lect10.R scriptet, og ved ikke helt hvor han har den fra,
  #men det ligner meget log(pdf) for en normal fordeling, hvor han har tilføjet leddet
  #n/2 * log(sigma) og fjernet log(1/(sigma*sqrt(2*pi)))
  n <- length(y) - 1
  return(-n/2 * log(sigma) - 1/(2*sigma) * sum((y[-1]-phi*y[-(n+1)])^2))
}

ll.pf.sigma.sq <- sapply(sigma.sq_interval,llp.sigma.sq,y=e_vec[-1], phi = phi.hat)

plot(sigma.sq_interval,exp(ll.pf.sigma.sq-max(ll.pf.sigma.sq)),type="l")
lines(range(sigma.sq_interval),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)

##profile likelihood based CI for sigma
sigma.sq.lower.pf <- min(sigma.sq_interval[exp(ll.pf.sigma.sq-max(ll.pf.sigma.sq)) >= exp(-qchisq(0.95,df=1)/2)])
sigma.sq.upper.pf <- max(sigma.sq_interval[exp(ll.pf.sigma.sq-max(ll.pf.sigma.sq)) >= exp(-qchisq(0.95,df=1)/2)])

## wald based CI for sigma
sd.sigma.sq <- se[1]
sigma.sq.hat.opt <- opt$par[1]

sigma.sq.lower.wald <- sigma.sq.hat.opt - sd.sigma.sq*qnorm(0.975)
sigma.sq.upper.wald <- sigma.sq.hat.opt + sd.sigma.sq*qnorm(0.975)

## Directly in R
#arima(e_vec,order=c(1,0,0))
#opt$par ## rather close. 

#### 2.2 - Contour plot ####
library(plotly)
nll.contour <- function(theta){
  n <- length(e_vec) - 1
  sigma <- theta[1]
  phi <- theta[2]
  return(n/2 * log(sigma) + 1/(2*sigma) * sum((e_vec[-1]-phi*e_vec[-(n+1)])^2))
}

intervals <- expand.grid(sigma.sq = sigma.sq_interval, phi = phi_interval)
z <- apply(intervals, MARGIN = 1, FUN = nll.contour)
intervals$likelihood <- exp(-z-max(-z))
fig <- plot_ly(intervals, x = ~sigma.sq, y = ~phi, z = ~likelihood
               , type = "contour")
c <- round(exp(-qchisq(0.95,df=1)/2), 4)
fig2 <- add_trace(p = fig, data = intervals, x = ~sigma.sq, y = ~phi, z = ~likelihood
                , type = "contour", 
                contours = list(showlabels = T
                                ,start = c
                                ,end = c
                                ,coloring = "lines"))
fig2 #den gule linje angiver contour-linjen for 95% confidence regions for de to
#parametre

#### 2.3 - p-value for LRT and Wald test for the null hypothesis: H_0: rho = 0 ####
nll.null.sigma.sq <- function(sigma,phi,y){
  #funktionen er taget fra lect10.R scriptet, og ved ikke helt hvor han har den fra,
  #men det ligner meget log(pdf) for en normal fordeling, hvor han har tilføjet leddet
  #n/2 * log(sigma) og fjernet log(1/(sigma*sqrt(2*pi)))
  n <- length(y) - 1
  return(n/2 * log(sigma) + 1/(2*sigma) * sum((y[-1]-phi*y[-(n+1)])^2))
}
null_model <- nlminb(c("sigma" = 1), nll.null.sigma.sq, y=e_vec[-1], phi = 0)
null_model$par

# LRT
chi.squared <- - 2 * (opt$objective - null_model$objective)
p.value.LRT <- 1 - pchisq(chi.squared, df = 1)

# WALD TEST
waldTestStatistic <- phi.hat.opt/sd.phi
wald.p.value <- 2*dnorm(waldTestStatistic)

#### 3.1 - Compare the found numerical information matrix with the algebraic form ####
I_theory <- function(sigma, phi, y){
  n <- length(y)
  I_theoretical <- matrix(c(sum(-1/(2*sigma^2) + (-phi*y[-(n)] + y[-1])^2/sigma^3)
                            ,sum((-phi*y[-(n)]+y[-1])*y[-(n)]/sigma^2)
                            ,sum((-phi*y[-(n)]+y[-1])*y[-(n)]/sigma^2)
                            ,sum(y[-(n+1)]^2)/sigma), nrow = 2)
  return(solve(I_theoretical))
}
I_theory(sigma.sq.hat.opt, phi.hat.opt, y = e_vec)
V
#MAGNIFIQUE!!!

#### 3.2 - two tasks (see below) ####
  #1) Plot of profile likelihood of phi and compare it two the quadratic approximation
  #2) Use z-transform (exercise 3.23) for phi and log-transform for sigma.sq
      #plot the profile likelihood for z and compare with the quadratic approx

# 3.2.1
par(mfrow = c(1,1))
plot(phi_interval,exp(llp-max(llp)),type="l", main = TeX("Profile likelihood for $\\rho$"))
grid()
lines(range(phi_interval),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)
