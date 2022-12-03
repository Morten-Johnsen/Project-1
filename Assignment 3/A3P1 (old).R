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

## Calculate the model residuals in original pow.obs/5000 domain.
D$residuals <- D$transformed.pow.obs.norm - D$ws.ws2.pred

# Show fit vs residuals
fit <- ggplot(data = D)+
  theme_bw()+ #THEME
  geom_point(aes(x = ws30, y = transformed.pow.obs.norm))+ #OBSERVED DATA
  geom_line(aes(x = ws30, y = ws.ws2.pred), colour = "red", size = 1)+ #MODEL FIT (LINE)
  labs(x = "Wind speeds", y="Transformed, norm power obs", colour = "")+
  theme_bw()+
  ggtitle("Fitted model in Transformed Domain") #TITLE
res <- ggplot(data = D)+
  theme_bw()+ #THEME
  geom_point(aes(x = ws30, y = residuals))+ #ACTUAL PLOT
  labs(x = "Wind speed [m/s]", y = "Residuals [1/5000 kW]")+ #AXIS LABELS
  geom_line(aes(x = ws30, y = 0), size = 1, colour = "red")+ #LINE GOING THROUGH ZERO
  ggtitle("Residuals")+ #TITLE
  theme_bw()
D$residuals2 <- ws.and.ws.squared$residuals
res2 <- ggplot(data = D)+
  theme_bw()+ #THEME
  geom_point(aes(x = ws30, y = residuals2))+ #ACTUAL PLOT
  labs(x = "Wind speed [m/s]", y = "Residuals [1/5000 kW]")+ #AXIS LABELS
  geom_line(aes(x = ws30, y = 0), size = 1, colour = "red")+ #LINE GOING THROUGH ZERO
  ggtitle("Residuals (Tester/direct with ws.and.ws.squared$residuals)")+ #TITLE
  theme_bw()
grid.arrange(fit,res,res2, nrow = 1)

################################################################################################
###Analysis of auto-correlation
##TASK 1 of 7...
#Residuals in the transformed domain are:
n <- dim(D)[1]
e1 <- D$residuals[-n] #for e_1 - e_n-1
e2 <- D$residuals[-1] #for e_2 - e_n
##arrange model residuals in a matrix consisting of [e_1...e_n-1; e_2...e_n]
e <- matrix(c(e1,e2), ncol=2, byrow=F)

#1.1
##acf plot of model residuals in the x domain
par(mfrow = c(1,1))
e_vec <- D$residuals
acf(e_vec)

#1.2
##calculate rho and sigma
rho_matrix <- cor(e)#?
(phi <- sum(e[,2]*e[,1])/sum(e[,1]^2))# direct calculation of correlation
#1/(n-1) * sum (e_{i+1} - phi * e_{i})
(sigma.sq <- 1/(length(e_vec)-1) * sum((e[,2] - phi * e[,1])^2))#?

#acf(e) plot of the found model based on direct MLE estimates
e_corrected <- e[,2]-phi*e[,1] #?
acf(e_corrected) #?

#directly with arima:
AR1m <- arima(e_vec, order=c(1,0,0))
acf(AR1m$residuals)

#### FROM LECT 10 ####
## Moment estimates
(phi.hat <- acf(e_vec,plot=FALSE)$acf[2]) ## ACF in lag 1 (index 1 is itself)
phi.hat <- cor(e1, e2) #same calculation but by hand.
(sigma <- var(e_vec)*(1-phi.hat^2))


## Likelihood estimation
nll <- function(theta,y){
  n <- length(y) - 1
  sigma <- theta[1]
  phi <- theta[2]
  return(n/2 * log(sigma) + 1/(2*sigma) * sum((y[-1]-phi*y[-(n+1)])^2)) #multivariate normal dist?
}

## MLE
(opt <- nlminb(c("sigma"=1,"phi"=1/2),nll,y=e_vec,lower=c(0.001,-0.999),upper=c(Inf,0.999)))

## Check MLEs
(phi <- sum(e2*e1)/sum(e1^2))
(sigma.sq <- 1/(n-1) * sum((e2 - phi * e1)^2))

#acf(e) plot of the found model based on OPTIMIZED MLE estimates
e_corrected <- e2-opt$par[2]*e1
acf(e_corrected)

## standard errors
V <- solve(hessian(nll,opt$par,y=e_vec))
(se <- sqrt(diag(V)))

## Profile likelihood for phi
llp.phi <- function(phi,x){
  n <- length(x) - 1
  return(-n/2 * log(sum((x[-1]-phi*x[-(n+1)])^2)))
}

## Plot profile likelihood
phi_interval <- seq(opt$par[2]-3*se[2],
                    opt$par[2]+3*se[2],
                    length=200)
sigma.sq_interval <- seq(opt$par[1]-3*se[1],
                         opt$par[1]+3*se[1],
                         length=200)

llp <- sapply(phi_interval,llp.phi,x=e_vec) #Har ændret x = e_vec[-1] til x = e_vec. (I jans kode var der anvendt e_vec[-1], men dette giver en pf, som ikke passer med den fundne MLE + det er jo tydeligt at man laver e_vec[-1], e_vec[-n] stratificering inde i funktionen så ved ikke hvor meget mening de giver at inputte e_vec[-1])

plot(phi_interval,exp(llp-max(llp)),type="l", main = "phi profile likelihood")
grid()
lines(range(phi_interval),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)

#### 2.1 - MLEs and confidence intervals ####
#phi confidence interval
#phi MLE
phi.hat.opt <- opt$par[2]
phi.hat.opt
abline(v = phi.hat.opt, lty = 2, col = 4)
## pf based
(phi.lower.pf <- min(phi_interval[exp(llp-max(llp)) >= exp(-qchisq(0.95,df=1)/2)]))
(phi.upper.pf <- max(phi_interval[exp(llp-max(llp)) >= exp(-qchisq(0.95,df=1)/2)]))

## wald based
sd.phi <- se[2]

(phi.lower.wald <- phi.hat.opt - sd.phi*qnorm(0.975))
(phi.upper.wald <- phi.hat.opt + sd.phi*qnorm(0.975))

#sigma.sq confidence interval
llp.sigma.sq <- function(sigma,phi,y){
  #funktionen er taget fra lect10.R scriptet, og ved ikke helt hvor han har den fra,
  #men det ligner meget log(pdf) for en normal fordeling, hvor han har tilføjet leddet
  #n/2 * log(sigma) og fjernet log(1/(sigma*sqrt(2*pi)))
  n <- length(y) - 1
  return(-n/2 * log(sigma) - 1/(2*sigma) * sum((y[-1]-phi*y[-(n+1)])^2))
}

ll.pf.sigma.sq <- sapply(sigma.sq_interval,llp.sigma.sq, y=e_vec, phi = phi.hat)

plot(sigma.sq_interval,exp(ll.pf.sigma.sq-max(ll.pf.sigma.sq)),type="l"
     ,main = "sigma.squared profile likelihood")
grid()
lines(range(sigma.sq_interval),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)

#sigma MLE
sigma.sq.hat.opt <- opt$par[1]
sigma.sq.hat.opt
abline(v = sigma.sq.hat.opt, lty = 2, col = 4)

##profile likelihood based CI for sigma
(sigma.sq.lower.pf <- min(sigma.sq_interval[exp(ll.pf.sigma.sq-max(ll.pf.sigma.sq)) >= exp(-qchisq(0.95,df=1)/2)]))
(sigma.sq.upper.pf <- max(sigma.sq_interval[exp(ll.pf.sigma.sq-max(ll.pf.sigma.sq)) >= exp(-qchisq(0.95,df=1)/2)]))

## wald based CI for sigma
sd.sigma.sq <- se[1]

(sigma.sq.lower.wald <- sigma.sq.hat.opt - sd.sigma.sq*qnorm(0.975))
(sigma.sq.upper.wald <- sigma.sq.hat.opt + sd.sigma.sq*qnorm(0.975))

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
               , type = "contour", contours = list(showlines = FALSE, size = 0.05), colorscale='Viridis')
c <- round(exp(-qchisq(0.95,df=1)/2), 4) #95% CI

fig2 <- add_trace(p = fig, data = intervals, x = ~sigma.sq, y = ~phi, z = ~likelihood
                  , type = "contour", showlegend = F,
                  contours = list(showlabels = T
                                  ,start = c
                                  ,end = c
                                  ,coloring = "lines"
                                  ,showlegend = F))
fig3 <- add_trace(p = fig2
                  , type = "scatter"
                  , x = sigma.sq.hat.opt
                  , y = phi.hat.opt)

fig3 #den gule linjer angiver contour-linjen for 95% confidence region for de to
  #parametre

#### 2.3 - p-value for LRT and Wald test for the null hypothesis: H_0: rho = 0 ####
nll.null.sigma.sq <- function(sigma,phi,y){
  #funktionen er taget fra lect10.R scriptet, og ved ikke helt hvor han har den fra,
  #men det ligner meget log(pdf) for en normal fordeling, hvor han har tilføjet leddet
  #n/2 * log(sigma) og fjernet log(1/(sigma*sqrt(2*pi)))
  n <- length(y) - 1
  return(n/2 * log(sigma) + 1/(2*sigma) * sum((y[-1]-phi*y[-(n+1)])^2))
}
null_model <- nlminb(c("sigma" = 1), nll.null.sigma.sq, y=e_vec, phi = 0)
null_model$par

# LRT
chi.squared <- - 2 * (opt$objective - null_model$objective)
p.value.LRT <- 1 - pchisq(chi.squared, df = 1)

# WALD TEST
waldTestStatistic <- phi.hat.opt/sd.phi
wald.p.value <- 2 * (1 - pnorm(waldTestStatistic))#fejl

#### 3.1 - Compare the found numerical information matrix with the algebraic form ####
#Differentierede udtrykket fra vores nll i maple for at lave en hessian.
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

#### 4.2 - two tasks (see below) ####
#1) Plot of profile likelihood of phi and compare it two the quadratic approximation
#2) Use z-transform (exercise 3.23) for phi and log-transform for sigma.sq
#plot the profile likelihood for z and compare with the quadratic approx

# 4.2.1
par(mfrow=c(1,2))
Hessian <- solve(I_theory(sigma.sq.hat.opt, phi.hat.opt, y = e_vec))

#til hvis man vil vise at den kvadratiske approksimation afviger på virkelig stor skala (uden for det interval det er relevant at kigge på)
# phi_interval <- seq(opt$par[2]-20*se[2],
#                     opt$par[2]+20*se[2],
#                     length=2000)
# llp <- sapply(phi_interval,llp.phi,x=e_vec) #Har ændret x = e_vec[-1] til x = e_vec. (I jans kode var der anvendt e_vec[-1], men dette giver en pf, som ikke passer med den fundne MLE + det er jo tydeligt at man laver e_vec[-1], e_vec[-n] stratificering inde i funktionen så ved ikke hvor meget mening de giver at inputte e_vec[-1])
#

plot(phi_interval,exp(llp-max(llp)),type="l", main = TeX("Profile likelihood for $\\rho$"))
grid()
lines(range(phi_interval),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)
abline(v = phi.hat.opt)
lines(phi_interval
      #quadratic approximation expression found in lect01.R og sol2.pdf.
      #Here 
      , exp(-0.5 * Hessian[2,2] * (phi_interval - phi.hat.opt)^2)
      , col = 4
      , lwd = 2
      , lty = 3)
text(x = 0.45, y = 0.8, "Quadratic \napproximation", col = 4)
plot(phi_interval,llp-max(llp),type="l", main = TeX("Profile log likelihood for $\\rho$"))
grid()
lines(range(phi_interval),
      -c(1,1) * qchisq(0.95,df=1)/2,
      col=2,lty=2)
abline(v = phi.hat.opt)
lines(phi_interval
      #quadratic approximation expression found in lect01.R og sol2.pdf.
      #Here 
      , -0.5 * Hessian[2,2] * (phi_interval - phi.hat.opt)^2
      , col = 4
      , lwd = 2
      , lty = 3)
text(x = 0.45, y = 0.8, "Quadratic \napproximation", col = 4)

#Perfekt kvadratisk approximation - måske for perfekt?

# 4.2.2
#z-transformation af rho og log-transformation af sigma
#Sol 2 under exercise 3.23 indeholder en god beskrivelse af de korrelations beregninger som er lavet
#for at finde sigma og rho.

#Føler der er noget i vejen med det her afsnit. z-transformationen gør den jo ringere? -.-

transformation_nll <- function(rho_sigmasq_vector, y){
  rho <- rho_sigmasq_vector[2]
  sigma.sq <- rho_sigmasq_vector[1]
  
  sigma.sq.log <- exp(sigma.sq)
  rho.z <- (exp(2*rho)-1)/(1+exp(2*rho))
  Sigma <- diag(2)
  Sigma[!Sigma] <- rho.z
  Sigma <- sigma.sq.log * Sigma
  return(-sum(dmvnorm(y, sigma = Sigma, log = TRUE)))
}
opt.z <- nlminb(c(1,0), transformation_nll, y = e)

H <- hessian(transformation_nll, opt.z$par ,y=e)
sd <- sqrt(solve(H))

z_interval <- seq(max(opt.z$par[2]-5*sd[2],-1),
                    min(opt.z$par[2]+5*sd[2],1),
                    length=200)
log.sigma.sq_interval <- seq(max(opt.z$par[1]-5*sd[1],-1),
                         min(opt.z$par[1]+5*sd[1],1),
                         length=200)

#Ved ikke om man bør optimere sigma i profile likelihooden for phi (rho) og vice-versa.
ll_z <- function(z, log.sigma.sq.hat, y){
  rho <- z
  sigma.sq <- log.sigma.sq.hat
  
  sigma.sq.log <- exp(sigma.sq)
  rho.z <- (exp(2*rho)-1)/(1+exp(2*rho))
  Sigma <- diag(2)
  Sigma[!Sigma] <- rho.z
  Sigma <- sigma.sq.log * Sigma
  return(sum(dmvnorm(y, sigma = Sigma, log = TRUE)))
}

llp.z <- sapply(z_interval, FUN = ll_z, log.sigma.sq.hat = opt.z$par[1], y = e)

plot(z_interval,exp(llp.z-max(llp.z)),type="l", main = TeX("Profile likelihood for z"))
grid()
lines(range(z_interval),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)
abline(v = opt.z$par[2])
lines(z_interval
      #quadratic approximation expression found in lect01.R og sol2.pdf.
      #Here 
      , exp(-0.5 * H[2,2] * (z_interval - opt.z$par[2])^2)
      , col = 2
      , lwd = 2
      , lty = 2)
text(x = 0.45, y = 0.8, "Quadratic \napproximation", col = 2)
plot(z_interval,llp.z-max(llp.z),type="l", main = TeX("Profile log-likelihood for z"))
grid()
lines(range(z_interval),
      -c(1,1) * qchisq(0.95,df=1)/2,
      col=2,lty=2)
abline(v = opt.z$par[2])
lines(z_interval
      #quadratic approximation expression found in lect01.R og sol2.pdf.
      #Here 
      , -0.5 * H[2,2] * (z_interval - opt.z$par[2])^2
      , col = 2
      , lwd = 2
      , lty = 2)
text(x = 0.45, y = 0.8, "Quadratic \napproximation", col = 2)

# 5 - Estimate the parameters of the AR model (Example 11.1) when conditioning on e1 and then full estimation
#compare with the estimation above

# p. 298 conditioning on the first observation or not. If not then it is the "full estimation"
#this can also be done by using the the arima model and setting the "method" input to different settings.
#in arima one can say "include mean = false" as we know the mean of the residuals should be zero
#Example 11.1: L2(theta) is a conditional likelihood THIS CONDITIONAL LIKELIHOOD IS COMMONLY ASSUMED IN ROUTINE DATAANALYSIS
#MÅSKE SKAL DETTE GØRES I OVENSTÅENDE MANUELLE BEREGNINGER.

arima(e_vec, order = c(1,0,0), method = 'CSS', n.cond = 1)
arima(e_vec, order = c(1,0,0))
#?



# 6 - simultaneous parameter estimation of our wind speed model and AR(1) model 
#  + likelihood comparison between the combined model and the purely linear model.
#can maybe be done in arima by specifying xreg.

pdf_normal <- function(x, mu = 0, sigma = 1, log_out = F){
  if (log_out == T){
    return(log(1/(sigma*sqrt(2*pi))) - 1/2*((x-mu)/sigma)^2)
  } else {
    return(1/(sigma*sqrt(2*pi)) * exp(-1/2 * ((x-mu)/sigma)^2))
  }
}

combined_nll <- function(theta, pow.obs, ws){
  beta0 <- theta[1]
  beta1 <- theta[2]
  beta2 <- theta[3]
  phi <- theta[4]
  sigma <- theta[5]
  
  e <- pow.obs - (beta0 + beta1*ws + beta2*ws^2)
  px <- -sum(pdf_normal(e, log_out = T))
  
  n <- length(e)
  px_given_xminus1 <- sum(n/2 * log(sigma) + 1/(2*sigma) * sum((e[-1]-phi*e[-n])^2))
  return(px + px_given_xminus1)
}

combined_model <- nlminb(start = c(1,0,0,0,1)
                       , objective = combined_nll
                       , pow.obs = D$transformed.pow.obs.norm
                       , ws = D$ws30)

#Linear model
linear_model_nll = function(theta){
  pred <-  theta[1]  + theta[2]*D$ws30 + theta[3]*D$ws30^2
  e <- D$transformed.pow.obs.norm - pred
  return(-sum(pdf_normal(e, log_out = T)))
}
linear_model <- nlminb(c(1,0,0), linear_model_nll)

#comparing likelihoods through LRT:
chi.squared <- - 2 * (combined_model$objective - linear_model$objective)
p.value.LRT <- 1 - pchisq(chi.squared, df = 1)

opt_lect10 <- nlminb(c(1,0), nll, y = e_vec)
opt_lect10$par

p298_full <- function(theta, e, e1){
  sigma.sq <- theta[1]
  phi <- theta[2]
  
  n <- length(e)
  L1 <- sqrt(2*pi*sigma.sq)^(-1) * sqrt((1 - phi^2)) * exp(-(1 - phi^2)/(2*sigma.sq) * ((e1 - 0)/(1 - phi))^2)
  L2 <- sqrt(2*pi*sigma.sq)^(-1) * exp(-1/(2*sigma.sq) * (e[-1] - 0 - phi*e[-n])^2)
  
  return(-log(L1) - sum(log(L2)))
}

p298_cond <- function(theta, e){
  sigma.sq <- theta[1]
  phi <- theta[2]
  n <- length(e)
  
  L2 <- sqrt(2*pi*sigma.sq)^(-1) * exp(-1/(2*sigma.sq) * (e[-1] - 0 - phi*e[-n])^2)
  return(-sum(log(L2)))
}

bivariate_normal <- function(x, mu = c(0,0), Sigma = diag(2)){
  return(1/(2*pi) * sqrt(solve(det(Sigma))) * exp(-1/2 * t(x) %*% solve(Sigma) %*% x))
}

mn_nll <- function(theta, e){
  sigma.sq <- theta[1]
  phi <- theta[2]
  Sigma <- diag(2)
  Sigma[!Sigma] <- phi
  Sigma <- Sigma * sigma.sq
  
  L <- apply(e, MARGIN = 1, FUN = bivariate_normal, Sigma = Sigma)
  return(-sum(log(L)))
}

opt_p298_full <- nlminb(c("sigma.sq" = 1, "phi" = 0), p298_full, e = e_vec, e1 = e_vec[1])
opt_p298_full$par

opt_p298_cond <- nlminb(c("sigma.sq" = 1, "phi" = 0), p298_cond, e = e_vec)
opt_p298_cond$par

opt_bivariate_normal <- nlminb(c("sigma.sq" = 1, "phi" = 0), mn_nll, e = e)
opt_bivariate_normal$par

opt_lect10$par

p298_cond_2 <- function(theta, e){
  sigma.sq <- theta[1]
  phi <- theta[2]
  n <- length(e)
  
  L2 <- sqrt(2*pi*sigma.sq)^(-1) * exp(-1/(2*sigma.sq) * (e[-1] - phi*e[-n])^2)
  return(-sum(log(L2)))
}

opt_p298_cond_2 <- nlminb(c("sigma.sq" = 1, "phi" = 0), p298_cond_2, e = e_vec)
opt_p298_cond_2$par

#discussion:

##husk husk: det er PREDICTIONS derfor er det kun tidsskridtet igår der er kendt
#og derfor vil effekten af AR(1) modellen aftage på long-term predictions.

#på short term gør AR(1) helt klar modellen bedre, da den inkluderer
#x_t-1 korrelationen i modellen og dermed opnår bedre prædiktioner, og modellens
#residualer bliver nærmere i.d.d.