rm(list = ls())
library(lubridate)
library(latex2exp)
library(circular)
library(tidyverse)
library(reshape)
library(pheatmap)
library(openair)
library(stringr)
library(gridExtra)
library(mvtnorm)
library(numDeriv)
library(plotly)
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
D$residuals <- ws.and.ws.squared$residuals
D$residuals2 <- D$transformed.pow.obs.norm - D$ws.ws2.pred

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
  ggtitle("Residuals (direct with ws.and.ws.squared$residuals)")+ #TITLE
  theme_bw()
res2 <- ggplot(data = D)+
  theme_bw()+ #THEME
  geom_point(aes(x = ws30, y = residuals2))+ #ACTUAL PLOT
  labs(x = "Wind speed [m/s]", y = "Residuals [1/5000 kW]")+ #AXIS LABELS
  geom_line(aes(x = ws30, y = 0), size = 1, colour = "red")+ #LINE GOING THROUGH ZERO
  ggtitle("Residuals")+ #TITLE
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
par(mfrow = c(1,1))
e_vec <- D$residuals

plot(e1, e2, xlab = "e1", ylab = "e2", main="Correlation between error and the error at the previous time step")
grid()

#1.2
##calculate rho and sigma
bivariate_normal <- function(x, mu = c(0,0), Sigma = diag(2)){
  return(1/(2*pi) * sqrt(solve(det(Sigma))) * exp(-1/2 * t(x) %*% solve(Sigma) %*% x))
}

bn_nll <- function(theta, e){
  sigma.sq <- theta[1]
  phi <- theta[2]
  Sigma <- diag(2)
  Sigma[!Sigma] <- phi
  Sigma <- Sigma * sigma.sq
  
  L <- apply(e, MARGIN = 1, FUN = bivariate_normal, Sigma = Sigma)
  return(-sum(log(L)))
}

opt_bivariate_normal <- nlminb(c("sigma.sq" = 1, "phi" = 0), bn_nll, e = e)
opt_bivariate_normal$par
sigma.sq.hat <- opt_bivariate_normal$par[1]
phi.hat <- opt_bivariate_normal$par[2]

(Sigma <- opt_bivariate_normal$par[1] * matrix(c(1,opt_bivariate_normal$par[2],opt_bivariate_normal$par[2],1),2))
#this is equivalent to:
var(e)

## standard errors
V <- solve(hessian(bn_nll, opt_bivariate_normal$par, e = e))
(se <- sqrt(diag(V)))

#### 2.1 - MLEs and confidence intervals ####
#phi confidence interval
#phi MLE
#abline(v = phi.hat, lty = 2, col = 4)

## wald based
phi.hat
sd.phi <- se[2]

(phi.lower.wald <- phi.hat - sd.phi*qnorm(0.975))
(phi.upper.wald <- phi.hat + sd.phi*qnorm(0.975))


#sigma MLE
sigma.sq.hat

## wald based CI for sigma
sd.sigma.sq <- se[1]

(sigma.sq.lower.wald <- sigma.sq.hat - sd.sigma.sq*qnorm(0.975))
(sigma.sq.upper.wald <- sigma.sq.hat + sd.sigma.sq*qnorm(0.975))

str <-expression(paste(sigma^2))

cat(" = ", round(sigma.sq.hat,2), "95% WALD CI [",  round(sigma.sq.lower.wald,2),",",round(sigma.sq.upper.wald,2),"]"
     ,"\n = ", round(phi.hat,3), "95% WALD CI [",round(phi.lower.wald,3),",",round(phi.upper.wald,3),"]"
)#first is variance, second is rho

phi_interval <- seq(phi.hat-4*se[2],
                    phi.hat+4*se[2],
                    length=200)
sigma.sq_interval <- seq(sigma.sq.hat-4*se[1],
                         sigma.sq.hat+4*se[1],
                         length=200)

## Directly in R
#### 2.2 - Contour plot ####
nll.contour <- function(theta){
  sigma.sq <- theta[1]
  phi <- theta[2]
  Sigma <- diag(2)
  Sigma[!Sigma] <- phi
  Sigma <- Sigma * sigma.sq
  
  L <- sum(dmvnorm(e, sigma = Sigma, log = T))
  return(L)
}

intervals <- expand.grid(sigma.sq = sigma.sq_interval, phi = phi_interval)
z_temp <- apply(intervals, MARGIN = 1, FUN = nll.contour)
intervals$likelihood <- exp(z_temp-max(z_temp))
z <- matrix(intervals$likelihood, nrow = length(sigma.sq_interval))

alpha <- c(0.01, .05, .1, .5, 1)
contour(x = sigma.sq_interval, y = phi_interval,
        z = z,
        levels = exp(-qchisq(1 - alpha, df = 2) / 2),
        labels = alpha, xlab=expression(paste(sigma^2)),ylab=expression(paste(rho))
        ,main = "Contour plot with confidence regions")
grid()

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
                  , x = sigma.sq.hat
                  , y = phi.hat
                  , xlab=expression(paste(sigma^2)), ylab=expression(paste(rho))
                  )

fig3


#### 2.3 - p-value for LRT and Wald test for the null hypothesis: H_0: phi = 0 ####
nll.null <- function(sigma.sq, e){
  Sigma <- diag(2)
  Sigma <- Sigma * sigma.sq
  
  nll <- -sum(dmvnorm(e, sigma = Sigma, log = T))
  return(nll)
}

null_model <- nlminb(c("sigma" = 1), nll.null, e=e)
null_model$par

# LRT
chi.squared <- - 2 * (opt_bivariate_normal$objective - null_model$objective)
p.value.LRT <- 1 - pchisq(chi.squared, df = 1)

# WALD TEST
waldTestStatistic <- phi.hat/sd.phi
wald.p.value <- 2 * (1 - pnorm(waldTestStatistic)) 


#### 3.1 - Compare the found numerical information matrix with the algebraic form ####
#FROM PAGE 72 in the course textbook
V_theory <- function(sigma.sq, phi, y){
  n <- length(y)
  I_theoretical <- matrix(c(n/sigma.sq^2
                            ,-n*phi/(sigma.sq*(1 - phi^2))
                            ,-n*phi/(sigma.sq*(1 - phi^2))
                            ,n*(1+phi^2)/(1-phi^2)^2), nrow = 2)
  return(solve(I_theoretical))
}
V_theory(sigma.sq.hat, phi.hat, y = ws.and.ws.squared$residuals)
V
#The information:
hessian(bn_nll, opt_bivariate_normal$par, e = e)
solve(V)
solve(V_theory(sigma.sq.hat, phi.hat, y = ws.and.ws.squared$residuals))
#MAGNIFIQUE!!!

#Relatively large diagnonal elements in V as compared to the diagonal elements (meaning that
#the two parameters are somewhat correlated). Therefore when calculating the profile
#likelihoods, it is necessary to account for this correlation, by reoptimizing 
#the other parameter (for phi pf we reoptimize sigma.sq etc). This would not be necessary
#if the off-diagonal elements had been very small.

#### 4.2 - two tasks (see below) ####
#1) Plot of profile likelihood of phi and compare it two the quadratic approximation
#2) Use z-transform (exercise 3.23) for phi and log-transform for sigma.sq
#plot the profile likelihood for z and compare with the quadratic approx

# 4.2.1
## Profile likelihood for phi
llp.phi <- function(phi, e){
  fun_tmp <- function(sigma.sq, phi){
    Sigma <- diag(2)
    Sigma[!Sigma] <- phi
    Sigma <- Sigma * sigma.sq
    
    return(-sum(dmvnorm(e, sigma = Sigma, log = T)))
  }
  
  L <- nlminb(start = 1, objective = fun_tmp, phi = phi)
  return(L$objective)
}

## Plot profile likelihood
nllp <- sapply(X = phi_interval, FUN = llp.phi, e=e) #Har ændret x = e_vec[-1] til x = e_vec. (I jans kode var der anvendt e_vec[-1], men dette giver en pf, som ikke passer med den fundne MLE + det er jo tydeligt at man laver e_vec[-1], e_vec[-n] stratificering inde i funktionen så ved ikke hvor meget mening de giver at inputte e_vec[-1])
par(mfrow = c(2,2))

Information_matrix <- solve(V_theory(sigma.sq.hat, phi.hat, y = ws.and.ws.squared$residuals))

plot(phi_interval,exp(-(nllp-min(nllp))),type="l", main = expression(paste(rho, " profile likelihood")),
     xlab=expression(paste(rho)))
grid()
lines(range(phi_interval),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)
abline(v = phi.hat, col = 4, lty = 3)
lines(phi_interval
      #quadratic approximation p. 33 in course textbook
      #expression also found in lect01.R og sol2.pdf.
      , exp(-0.5 * Information_matrix[2,2] * (phi_interval - phi.hat)^2)
      , col = 4
      , lwd = 2
      , lty = 3)
text(x = 0.45, y = 0.8, "Quadratic \napproximation", col = 4)

plot(phi_interval,-(nllp-min(nllp)),type="l", main = expression(paste(rho," log-profile likelihood")),
     xlab=expression(paste(rho)))
grid()
lines(range(phi_interval),
      -c(1,1) * qchisq(0.95,df=1)/2,
      col=2,lty=2)
abline(v = phi.hat, col = 4, lty = 3)
lines(phi_interval
      #quadratic approximation p. 33 in course textbook
      #expression also found in lect01.R og sol2.pdf.
      , -0.5 * Information_matrix[2,2] * (phi_interval - phi.hat)^2
      , col = 4
      , lwd = 2
      , lty = 3)
text(x = 0.45, y = -1.3, "Quadratic \napproximation", col = 4)

## pf based CI
(phi.lower.pf <- min(phi_interval[exp(-(nllp-min(nllp))) >= exp(-qchisq(0.95,df=1)/2)]))
(phi.upper.pf <- max(phi_interval[exp(-(nllp-min(nllp))) >= exp(-qchisq(0.95,df=1)/2)]))
#God kvadratisk approximation


# 4.2.2
#z-transformation af rho og log-transformation af sigma
#Sol 2 under exercise 3.23 indeholder en god beskrivelse af de korrelations beregninger som er lavet
#for at finde sigma og rho.

transformed_nll <- function(theta, e){
  rho <- theta[2]
  sigma.sq <- theta[1]
  z <- (exp(2*rho) - 1) / (1 + exp(2*rho))
  sigma.sq_log <- exp(sigma.sq)
  Sigma <- diag(2)
  Sigma[!Sigma] <- z
  Sigma <- Sigma * sigma.sq_log
  
  L <- -sum(dmvnorm(e, sigma = Sigma, log = T))
  return(L)
}

transformed_nll_pf <- function(rho, e){
  z <- (exp(2*rho) - 1) / (1 + exp(2*rho))
  fun_tmp <- function(sigma.sq, z){
    sigma.sq_log <- exp(sigma.sq)
    Sigma <- diag(2)
    Sigma[!Sigma] <- z
    Sigma <- Sigma * sigma.sq_log
    
    return(-sum(dmvnorm(e, sigma = Sigma, log = T)))
  }
  
  L <- nlminb(start = 1, objective = fun_tmp, z = z)
  return(L$objective)
}

opt.z <- nlminb(start = c(1,0), objective = transformed_nll, e = e)
Information_matrix_transformed <- hessian(transformed_nll, opt.z$par ,e=e)
sd <- diag(sqrt(solve(Information_matrix_transformed)))

z_interval <- seq(opt.z$par[2]-5*sd[2],
                  opt.z$par[2]+5*sd[2],
                  length=200)

nllp.z <- sapply(z_interval, FUN = transformed_nll_pf, e = e)

plot(z_interval,exp(-(nllp.z-min(nllp.z))), type="l", main = expression(paste("Profile likelihood for z")),
     xlab="z")
grid()
lines(range(z_interval),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
      col=2,lty=2)
abline(v = opt.z$par[2], col = 4, lty = 3)
lines(z_interval
      , exp(-0.5 * Information_matrix_transformed[2,2] * (z_interval - opt.z$par[2])^2)
      , col = 4
      , lwd = 2
      , lty = 3)
text(x = 0.45, y = 0.8, "Quadratic \napproximation", col = 4)
plot(z_interval,-(nllp.z-min(nllp.z)),type="l", main = expression(paste("Profile log-likelihood for z")),
     xlab="z")
grid()
lines(range(z_interval),
      -c(1,1) * qchisq(0.95,df=1)/2,
      col=2,lty=2)
abline(v = opt.z$par[2], col = 4, lty = 3)
lines(z_interval
      #quadratic approximation expression found in lect01.R og sol2.pdf.
      #Here 
      , -0.5 * Information_matrix_transformed[2,2] * (z_interval - opt.z$par[2])^2
      , col = 4
      , lwd = 2
      , lty = 3)
text(x = 0.45, y = -1.5, "Quadratic \napproximation", col = 4)

# 5 - Estimate the parameters of the AR model (Example 11.1) when conditioning on e1 and then full estimation
#compare with the estimation above

# p. 298 conditioning on the first observation or not. If not then it is the "full estimation"
#this can also be done by using the the arima model and setting the "method" input to different settings.
#in arima one can say "include mean = false" as we know the mean of the residuals should be zero
#Example 11.1: L2(theta) is a conditional likelihood THIS CONDITIONAL LIKELIHOOD IS COMMONLY ASSUMED IN ROUTINE DATAANALYSIS
#MÅSKE SKAL DETTE GØRES I OVENSTÅENDE MANUELLE BEREGNINGER.

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

opt_p298_full <- nlminb(c("sigma.sq" = 1, "phi" = 0), p298_full, e = e_vec, e1 = e_vec[1])
opt_p298_full$par

opt_p298_cond <- nlminb(c("sigma.sq" = 1, "phi" = 0), p298_cond, e = e_vec)
opt_p298_cond$par

opt_bivariate_normal$par
cat("Conditional = ", round(opt_p298_full$par[1],2),",", round(opt_p298_full$par[2],3)
    ,"\nFull = ", round(opt_p298_cond$par[1],2),",", round(opt_p298_cond$par[2],3)
    ,"\nMVN = ", round(opt_bivariate_normal$par[1],2),",", round(opt_bivariate_normal$par[2],3)
)
# giver god mening at variansen bliver mindre ift. det der ses i den direkte covarians
#beregning, da vi nu har anvendt en AR(1) til at reducere spredningen i residualerne
#dermed er estimatet på≈8 et udtryk for spredningen sigma.sq_u (se egne noter på bagerste
#side i bogen) og det originale estimat på ≈ 9 er sigma.sq.eps. TL:DR AR(1) modellen
#reducerer spredningen i modellens fejl (læs den gør modellen bedre).

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
  rho <- theta[4]
  sigma.sq <- theta[5]
  
  e <- pow.obs - (beta0 + beta1*ws + beta2*ws^2)
  e1 <- e[1]
  
  n <- length(e)
  pe_L1 <- -sum(log(sqrt(2*pi*sigma.sq)^(-1) * sqrt((1 - rho^2)) * exp(-(1 - rho^2)/(2*sigma.sq) * (e1/(1 - rho))^2)))
  pe_L2 <- -sum(log(1/sqrt(2*pi*sigma.sq) * exp(-1/(2*sigma.sq) * (e[-1] - rho*e[-n])^2)))
  return(pe_L1 + pe_L2)
}

combined_model <- nlminb(start = c(1,0,0,0,1)
                         , objective = combined_nll
                         , pow.obs = D$transformed.pow.obs.norm
                         , ws = D$ws30)
#Linear model
linear_model_nll = function(theta){
  pred <-  theta[1]  + theta[2]*D$ws30 + theta[3]*D$ws30^2
  e <- D$transformed.pow.obs.norm - pred
  return(-sum(pdf_normal(e, sigma = sqrt(var(e)), log_out = T)))
}
linear_model <- nlminb(c(1,0,0), linear_model_nll)

#comparing likelihoods through LRT:
chi.squared <- - 2 * (combined_model$objective - linear_model$objective)
p.value.LRT <- 1 - pchisq(chi.squared, df = 1)

(AIC_combined <- 2*combined_model$objective + 2*5)#mangler change of variable, men sammenligningen holder vel stadig
(AIC_linear <- 2*linear_model$objective + 2*3)

logLik.normTrans <- -linear_model$objective+sum(log(-1/(D$pow.obs.norm*(-1+D$pow.obs.norm^0.26))))#change of variable
(AIC.normTrans <- -2*(logLik.normTrans) + 2*length(linear_model$par))

logLik.normTrans.comb <- -combined_model$objective+sum(log(-1/(D$pow.obs.norm*(-1+D$pow.obs.norm^0.26))))
(AIC.normTrans.comb <- -2*(logLik.normTrans.comb) + 2*length(combined_model$par))

#til beregning af sigma:
( var.lin.mod <- 1/length(D$transformed.pow.obs.norm) * sum( (D$transformed.pow.obs.norm - ( linear_model$par[1] + linear_model$par[2]*D$ws30 + linear_model$par[3]*D$ws30^2 ))^2  ) )
( var.comb.mod <- combined_model$par[5] )

#discussion:
#(se bagerste side i bogen med egne noter)
##husk husk: det er PREDICTIONS derfor er det kun tidsskridtet igår der er kendt
#og derfor vil effekten af AR(1) modellen aftage på long-term predictions.

#på short term gør AR(1) helt klar modellen bedre, da den inkluderer
#x_t-1 korrelationen i modellen og dermed opnår bedre prædiktioner, og modellens
#residualer bliver nærmere i.d.d.