rm(list = ls())
library(tidyverse)
library(gridExtra)

if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/")
} else {
  print("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1")
}
source("testDistribution.R")
#### Analysis of the Binary Data ####
log.data <- read.table("Logistic.txt", header=TRUE, sep="", 
                       as.is=TRUE)

#Data visualization:


#all data from one population:
bin.par <- nlminb(start = 0.1, objective = testDistribution
                  , x = c(sum(log.data$AIDS_yes), sum(log.data$n))
                  , distribution = "binomial")

#separately for the groups
x.AZT <- log.data %>%
  filter(AZT == "Yes") %>%
  select(AIDS_yes, n) %>%
  as.numeric()

AZT.par <- nlminb(start = 0.1, objective = testDistribution
                  , x = c(x.AZT[1], x.AZT[2])
                  , distribution = "binomial")

x.no.AZT <- log.data %>%
  filter(AZT == "No") %>%
  select(AIDS_yes, n) %>%
  as.numeric()

no.AZT.par <- nlminb(start = 0.1, objective = testDistribution
                     , x = c(x.no.AZT[1], x.no.AZT[2])
                     , distribution = "binomial")

#Compare proportions:
#There's generally two approaches to doing this: 1) The large sample hypothesis test
#2) Likelihood analysis


## Large sample hypothesis test:
#Null-hypothesis: There's no difference between the groups and our best estimate of p
# is the combined p (probability for developing AIDS)
p.hat <- sum(log.data$AIDS_yes)/sum(log.data$n)#bin.par$par


#Calculate expected values for this group based on each group size:
e.A.AZT <- log.data$n[log.data$AZT == "Yes"]*p.hat
e.A.no_AZT <- log.data$n[log.data$AZT == "No"]*p.hat

e.nA.AZT <- log.data$n[log.data$AZT == "Yes"]*(1-p.hat)
e.nA.no_AZT <- log.data$n[log.data$AZT == "No"]*(1-p.hat)

e <- c(e.A.AZT, e.A.no_AZT, e.nA.AZT, e.nA.no_AZT)
#by hand
chi_squared <- sum((c(log.data$AIDS_yes,log.data$n-log.data$AIDS_yes)-e)^2/e)
(chi_squared)
#probability of observing this chi-squared test statistic given that the null-hypothesis is true
rows <- dim(log.data)[1]
columns <- dim(log.data)[2]-1 #-1 because of the AZT column
pchisq(chi_squared,df=(rows-1)*(columns-1),lower.tail=FALSE)

#WITH CONTINUITY CORRECTION:
#https://en.wikipedia.org/wiki/Yates%27s_correction_for_continuity
chi_squared_yates <- sum((abs(c(log.data$AIDS_yes,log.data$n-log.data$AIDS_yes)-e)-0.5)^2/e)
(chi_squared_yates)
#probability of observing this chi-squared test statistic given that the null-hypothesis is true
rows <- dim(log.data)[1]
columns <- dim(log.data)[2]-1 #-1 because of the AZT column
pchisq(chi_squared_yates,df=(rows-1)*(columns-1),lower.tail=FALSE)

#direct using R:
log.data.for.chi <- log.data; log.data.for.chi$f <- log.data.for.chi$n - log.data.for.chi$AIDS_yes
prop.test(log.data$AIDS_yes, log.data$n)
#or
chisq.test(as.matrix(log.data.for.chi[,c(2,4)]))
### Result: There's a difference between the two groups.
print(paste0("Mean p for group with AZT treatment: ", round(AZT.par$par,3)))
print(paste0("Mean p for group with no AZT treatment: ", round(no.AZT.par$par,3)))


## Likelihood analysis (Den her del er formodentligt ikke færdig)
x1 <- log.data$AIDS_yes[1]
x2 <- log.data$AIDS_yes[2]
n1 <- log.data$n[2]
n2 <- log.data$n[2]

(theta.hat <- log(x1/(n1-x1)*(n2-x2)/x2))

eta.start <- log(x2/(n2-x2))

L_theta <- function(eta){ 
  return(-log(exp(theta.hat*x1+eta*(x1+x2))/((1+exp(theta.hat+eta))^n1*(1+exp(eta))^n2)))
}
(eta.hat <- nlminb(eta.start, objective = L_theta)$par)

#Now we can construct CIs based on the Likelihood function
L <- function(theta){
  return(-log(exp(theta*x1+eta.hat*(x1+x2))/((1+exp(theta+eta.hat))^n1*(1+exp(eta.hat))^n2)))
}

library(numDeriv)
var <- solve(hessian(func <- L, x = theta.hat))
sd <- as.numeric(sqrt(var))

CI <- p.hat + c(-1,1)*qnorm(c(0.975))*sd/sqrt(sum(log.data$n))


#### Bullet point 4
#Estimate parameters in the model and report a confidence interval for the parameter 
#describing the difference, compare with the result above.
#p_0: Probability of aids in control group
#p_1: Probability of aids in treatment group


#calculate likelihood
nll.p_0 <- function(beta, x = log.data$AIDS_yes[2], n = log.data$n[2]){
  p <- exp(beta)/(1+exp(beta))
  nll <- -sum(dbinom(x, size = n, prob = p, log = T))
  return(nll)
}
opt.p_0 <- nlminb(start = 1, objective = nll.p_0, x = log.data$AIDS_yes[2], n = log.data$n[2])
beta_0 <- opt.p_0$par

nll.p_1 <- function(beta_1, beta_0, x = log.data$AIDS_yes[1], n = log.data$n[1]){
  p <- exp(beta_0+beta_1)/(1+exp(beta_0+beta_1))
  nll <- -sum(dbinom(x, size = n, prob = p, log = T))
}
opt.p_1 <- nlminb(start = 1
                  , objective = nll.p_1
                  , beta_0 = beta_0
                  , x = log.data$AIDS_yes[1]
                  , n = log.data$n[1])
beta_1 <- opt.p_1$par

(p_0 <- exp(beta_0)/(1 + exp(beta_0)))
(p_1 <- exp(beta_0 + beta_1) / (1 + exp(beta_0 + beta_1)))

log.data
logistic <- data.frame("AZT" = c(rep(1,170), rep(0,168))
                       ,"AIDS_yes" = c(rep(c(1,0),c(25,170-25)), rep(c(1,0), c(44, 168-44))))

fit.glm <- glm(AIDS_yes ~ AZT, data = logistic, family = binomial)
print(paste0("with glm model: ", coef(fit.glm)))
print(paste0("By hand (according to slide 19 lect 4): "))
print(paste0("beta_0 = ", beta_0, ", beta_1 = ", beta_1))
summary(fit.glm)

#results show: -0.72 logits(?) for developing AIDS when using the treatment

# Confidence interval for the two beta parameters. 
# Bør det her være for beta eller p????
confint(fit.glm)
#calculate profile likelihoods
prof.b0 <- function(beta0, x, n){
  p <- exp(beta0)/(1+exp(beta0))
  return(-dbinom(x, size = n, prob = p, log = T))
}

prof.b1 <- function(beta1, beta0, x, n){
  p <- exp(beta0+beta1)/(1+exp(beta0+beta1))
  return(-dbinom(x, size = n, prob = p, log = T))
}
beta.zero.sims <- seq(-1.5,-0.6,0.01)
beta.one.sims <- seq(-1.3,-0.2,0.01)
pL.b0 <- prof.b0(beta.zero.sims
                 , x = log.data$AIDS_yes[2]
                 , n = log.data$n[2])
pL.b1 <- prof.b1(beta.one.sims
                 , x = log.data$AIDS_yes[1]
                 , n = log.data$n[1]
                 , beta0 = beta_0)
par(mfrow=c(1,2))
plot(beta.zero.sims
     , -(pL.b0+max(-pL.b0))
     , "l"
     ,main = "Profile likelihood for Beta_0")
abline(h = -qchisq(0.95, df = 1)/2, lty = "dashed")
plot(beta.one.sims
     , -(pL.b1+max(-pL.b1))
     , "l"
     ,main = "Profile likelihood for Beta_1")
abline(h = -qchisq(0.95, df = 1)/2, lty = "dashed")

#From these figures it can be concluded that the quadratic approximation
#of the CI through use of fischers information matrix, is a 
#good approksimation.
#redefine because x is used

sd_0 <- as.numeric(sqrt(solve(hessian(beta_0, func = nll.p_0))))
sd_1 <- as.numeric(sqrt(solve(hessian(beta_1, func = nll.p_1, beta_0 = beta_0))))

#Wald 95% CIs and profile-likelihoods with approx 95% CI
(W.CI.0 <- beta_0 + c(-1,1)*qnorm(0.975)*sd_0)
(W.CI.1 <- beta_1 + c(-1,1)*qnorm(0.975)*sd_1)

#Direkte numerisk approksimation:
(CI.0 <- c(min(beta.zero.sims[-(pL.b0+max(-pL.b0)) > -qchisq(0.95, df = 1)/2])
           ,max(beta.zero.sims[-(pL.b0+max(-pL.b0)) > -qchisq(0.95, df = 1)/2])))
(CI.1 <- c(min(beta.one.sims[-(pL.b1+max(-pL.b1)) > -qchisq(0.95, df = 1)/2])
           ,max(beta.one.sims[-(pL.b1+max(-pL.b1)) > -qchisq(0.95, df = 1)/2])))

plot(beta.zero.sims
     , -(pL.b0+max(-pL.b0))
     , "l"
     ,main = "Profile likelihood for Beta 0"
     ,xlab = expression(beta[0])
     ,ylab = "lp - max(lp)")
abline(h = -qchisq(0.95, df = 1)/2, col = 2)
abline(v = c(W.CI.0), col = 6)
text(x = W.CI.0[1]+0.2, y = -3, "Wald CI", col = 6)
text(x = CI.0[1]+0.1, y = -2.5, "CI", col = 2)
abline(v = c(CI.0), lty = "dashed", col = 2)
plot(beta.one.sims
     , -(pL.b1+max(-pL.b1))
     , "l"
     ,main = "Profile likelihood for beta 1"
     ,xlab = expression(beta[1])
     ,ylab = "lp - max(lp)")
abline(h = -qchisq(0.95, df = 1)/2, col = 2)
abline(v = c(W.CI.1), col = 6)
text(x = W.CI.1[1]+0.2, y = -3, "Wald CI", col = 6)
text(x = CI.1[1]+0.1, y = -2.5, "CI", col = 2)
abline(v = c(CI.1),lty = "dashed", col = 2)

#Hvordan kan man sige at de her værdier beskriver forskellen?
#her er det tydeligt at beta_1 confint ikke krydser 0 og dermed kan det siges
#at der er forskel på sandsynligheden for at få AIDS p_0 og p_1 blandt de to
#grupper.

#### Analysis of the Survival Data ####
#tx: Treatment indicator. 1 = New treatment, 0 = Control treatment
#event: Indicator for AIDS or death. 1 = AIDS diagnosis or death, 0 = Otherwise
#time: Time to AIDS diagnosis or death. Days
#så tiden for event = 0 må angive at personen har været med i studiet time[X] dage uden at være enten død eller fået AIDS.
actg320 <- read.table("actg320.txt", header=TRUE, sep="", 
                      as.is=TRUE)

#select time, event and tx as they are the only relevant variables in this project
actg <- actg320 %>%
  select(time, event, tx)

#### Bulletpoint 2) #### in the project description (MISSING: "Other relevant number that could be calculated?")
actg %>%
  group_by(tx) %>%
  summarise("Got AIDS or DIED" = sum(event),
            "Proportion" = sum(event)/n(),
            "Participants Given the Treatment" = n())

#### Bulletpoint 3) ####
#Fitting an exponential model to time for both and for each treatment
#only use times for event = 1, to filter out all the time of event indices with are longer than the reported time
#given the fact that the participants in the event = 0 group, has not 'experienced' the event yet.
#Ved sgu ikke om ovenstående er en passende antagelse....

actg_event <- actg %>%
  filter(event == 1)

both <- nlminb(start = 2
               , objective = testDistribution
               , x = actg_event$time
               , distribution = "exponential")

#separate exponential models
t1 <- nlminb(start = 2
             , objective = testDistribution
             , x = filter(actg_event, tx == 1)$time
             , distribution = "exponential")

t0 <- nlminb(start = 2
             , objective = testDistribution
             , x = filter(actg_event, tx == 0)$time
             , distribution = "exponential")
#Potato plots:
p.both <- ggplot(actg_event)+
  geom_histogram(aes(x = time, y = ..density.., fill = "Data"), alpha = 0.5)+
  stat_function(aes(colour = "Exp. Model"), fun = dexp, n = dim(actg_event)[1], args = list(rate = both$par))+
  ggtitle("Ignoring Treatment Effect")+
  theme(legend.position = "top")+
  lims(x = c(0,max(actg_event$time)+10), y = c(0,0.012))+
  labs(fill = "", colour = "", x = "Time to Event")+
  scale_colour_manual(values = "purple")+
  scale_fill_manual(values = "purple")

p.t1 <- ggplot(actg_event[actg_event$tx == 1,])+
  geom_histogram(aes(x = time, y = ..density.., fill = "Data"), alpha = 0.5)+
  stat_function(aes(colour = "Exp. Model"), fun = dexp, n = dim(actg_event)[1], args = list(rate = t1$par))+
  ggtitle("Treatment")+
  theme(legend.position = "top")+
  lims(x = c(0,max(actg_event$time)+10), y = c(0,0.012))+
  labs(fill = "", colour = "", x = "Time to Event")+
  scale_colour_manual(values = "blue")+
  scale_fill_manual(values = "blue")

p.t2 <- ggplot(actg_event[actg_event$tx == 0,])+
  geom_histogram(aes(x = time, y = ..density.., fill = "Data"), alpha = 0.5)+
  stat_function(aes(colour = "Exp. Model"), fun = dexp, n = dim(actg_event)[1], args = list(rate = t0$par))+
  ggtitle("No Treatment")+
  theme(legend.position = "top")+
  lims(x = c(0,max(actg_event$time)+10), y = c(0,0.012))+
  scale_colour_manual(values = "red")+
  labs(fill = "", colour = "", x = "Time to Event")+
  scale_fill_manual(values = "red")

grid.arrange(p.both, p.t1, p.t2, nrow = 1)

#### bullet point 4: Compare likelihoods ####
#Likelihood Ratio Test (LRT) comparison
#one model:
chi_squared <- - 2 * ((t1$objective + t0$objective) - both$objective)
(p_value <- 1 - pchisq(chi_squared, df = 1))
#no difference as p_value ≈ 0.46 > 0.05.

#Overvej også at lave et residual plot for at vise at det også er et elendigt fit.
#Overvej at beregne time ratio og hazard ratio


#### bullet point 5: Model with a parameter indicating the treatment effect. ####
#library(survival)
#library(survminer)
kaplan.meier <- survfit(Surv(time) ~ tx, data = actg)
ggsurvplot_add_all(kaplan.meier, data = actg, conf.int = T
                   ,risk.table = T, pval = T)

fit <- survreg(Surv(time, event = event) ~ tx, data = actg,
               dist = "exponential")
summary(fit)
confint(fit)

#Overvej residual plot

#Ifølge ovenstående:
#beta0 = 7.62 95% CI [7.38; 7.87]
#beta1 = 0.699 85% CI [0.28; 1.12]
# => Significant difference.

#ifølge ovenstående er der statistisk signifikant forskel. Herunder regnes i hånden i stedet, så
#vi ved hvad der foregår.

#I hånden (jvf. slides fra uge 7):
#model: T = exp(B0 + B1*tx)*epsilon, epsilon ~ exp(1)
#Der kan opstilles to forskellige modeller afhængigt af tx = 0 eller tx = 1.
#tx = 0: E[T] = exp(b0)*epsilon
#tx = 1: E[T] = exp(b0 + b1)*epsilon

#Likelihood
nll.exp <- function(beta, time = actg$time, event = actg$event, treatment = actg$tx){
  beta0 <- beta[1]
  #dont want to make two functions so let beta1 = 0 if no treatment is not considered/used:
  if (max(treatment) == 0){
    beta1 <- 0
  } else {
    beta1 <- beta[2]
  }
  
  h <- exp(- beta0 - beta1 * treatment)
  H <- time/exp(beta0 + beta1*treatment)
  nll <- -sum(event*log(h) - H)
  return(nll)
}

beta_hat <- nlminb(start = c(1,1)
                   , objective = nll.exp
                   , time = actg$time
                   , event = actg$event
                   , treatment = actg$tx)
beta_hat$par


#Comparing likelihoods with the result from bullet-point 4
beta_hat$objective #Ved ikke lige om der skal sammenlignes med de to modeller eller den ene? ahh
#måske skal man undersøge om begge værdier er statistisk signifikante således at vi kan argumentere for
#at der er tale om at behandlingen virker og sammenligne dette resultat med bullet-point 4.

#Calculate LRT
#optimise model without beta1 (no treatment):
beta_no_treatment_effect <- nlminb(start = 1
                                   , objective = nll.exp
                                   , time = actg$time
                                   , event = actg$event
                                   , treatment = rep(0, length(actg$tx)))

beta_no_treatment_effect$par

#LRT:
chi_squared <- - 2 * (beta_hat$objective - beta_no_treatment_effect$objective)
(p_value <- 1 - pchisq(chi_squared, df = 1))
#Here, we see that the treatment effect is statistically significant.

#### Bullet point 6 - Wald CI for the treatment parameters beta0 and beta1 ####
#Calculate profile likelihoods to ensure that the quadratic approximation by using Fischers Information matrix
#is acceptable.

beta.zero.sims <- seq(7.35,7.9,0.01)
beta.one.sims <- seq(0.3,1.2,0.01)
pL.beta0 <- apply(X = data.frame(beta.zero.sims,beta_hat$par[2]), MARGIN = 1 , FUN = nll.exp, time = actg$time, event = actg$event, treatment = actg$tx)
pL.beta1 <- apply(X = data.frame(beta_hat$par[1],beta.one.sims), MARGIN = 1 , FUN = nll.exp, time = actg$time, event = actg$event, treatment = actg$tx)

par(mfrow=c(1,2))
plot(beta.zero.sims
     , -(pL.beta0+max(-pL.beta0))
     , "l"
     ,main = "Profile likelihood for Beta 0"
     ,xlab = expression(beta[0])
     ,ylab = "lp - max(lp)"
     ,ylim = c(-3,0))
abline(h = -qchisq(0.95, df = 1)/2, col = 2)
plot(beta.one.sims
     , -(pL.beta1+max(-pL.beta1))
     , "l"
     ,main = "Profile likelihood for beta 1"
     ,xlab = expression(beta[1])
     ,ylab = "lp - max(lp)"
     ,ylim = c(-3,0))
abline(h = -qchisq(0.95, df = 1)/2, col = 2)

#CI:
sd <- as.numeric(sqrt(diag(solve(hessian(beta_hat$par, func = nll.exp)))))

#Wald 95% CIs and profile-likelihoods with approx 95% CI
#Måske er der et eller andet i vejen med de her WALD CIs
(Wald.CI <- beta_hat$par + matrix(c(-1,1), 2,2, byrow = T) * matrix(qnorm(0.975)*sd, 2,2, byrow = F))

#Direkte numerisk approksimation:
(CI.0 <- c(min(beta.zero.sims[-(pL.beta0+max(-pL.beta0)) > -qchisq(0.95, df = 1)/2])
           ,max(beta.zero.sims[-(pL.beta0+max(-pL.beta0)) > -qchisq(0.95, df = 1)/2])))
(CI.1 <- c(min(beta.one.sims[-(pL.beta1+max(-pL.beta1)) > -qchisq(0.95, df = 1)/2])
           ,max(beta.one.sims[-(pL.beta1+max(-pL.beta1)) > -qchisq(0.95, df = 1)/2])))

#### Bullet Point 6 - (and partially also 5) ####
plot(beta.zero.sims
     , -(pL.beta0+max(-pL.beta0))
     , "l"
     ,main = "Profile likelihood for Beta 0"
     ,xlab = expression(beta[0])
     ,ylab = "lp - max(lp)"
     ,ylim = c(-3,0))
abline(h = -qchisq(0.95, df = 1)/2, col = 2)
abline(v = Wald.CI[1,], col = 6)
text(x = Wald.CI[1,1]+.15, y = -2.4, "Wald CI", col = 6)
text(x = CI.0[1]+.05, y = -2.5, "CI", col = 2)
abline(v = c(CI.0), lty = "dashed", col = 2)
plot(beta.one.sims
     , -(pL.beta1+max(-pL.beta1))
     , "l"
     ,main = "Profile likelihood for beta 1"
     ,xlab = expression(beta[1])
     ,ylab = "lp - max(lp)"
     ,ylim = c(-3,0))
abline(h = -qchisq(0.95, df = 1)/2, col = 2)
abline(v = Wald.CI[2,], col = 6)
text(x = Wald.CI[2,1]+0.25, y = -1.6, "Wald CI", col = 6)
text(x = CI.1[1]+0.09, y = -1.7, "CI", col = 2)
abline(v = c(CI.1),lty = "dashed", col = 2)

#Consider residual plot

#### Try Bullet Point 4-6 again but with a weibull distribution ####

nll.wei <- function(theta, time = actg$time, event = actg$event, treatment = actg$tx){
  sigma <- theta[1]
  beta0 <- theta[2]
  #dont want to make two functions so let beta1 = 0 if no treatment is not considered/used:
  if (max(treatment) == 0){
    beta1 <- 0
  } else {
    beta1 <- theta[3]
  }
  
  h <-1/sigma * time^(1/sigma - 1) * exp(-1/sigma * (beta0 + beta1 * treatment))
  H <- time^(1/sigma) * exp(-1/sigma * (beta0 + beta1*treatment))
  nll <- -sum(event*log(h) - H)
  return(nll)
}

theta_hat <- nlminb(start = c(1,1,0.5)
                   , objective = nll.wei)

theta_hat$par


#Comparing likelihoods with the result from bullet-point 4
theta_hat$objective 

#Calculate LRT
#optimise model without beta1 (no treatment):
theta_no_treatment_effect <- nlminb(start = c(1,1)
                                   , objective = nll.wei
                                   , treatment = rep(0, length(actg$tx)))

theta_no_treatment_effect$par

#LRT:
chi_squared <- - 2 * (theta_hat$objective - theta_no_treatment_effect$objective)
(p_value <- 1 - pchisq(chi_squared, df = 1))
#Here, we see that the treatment effect is statistically significant.

#Calculate profile likelihoods to ensure that the quadratic approximation by using Fischers Information matrix
#is acceptable.

beta.zero.sims <- seq(7.8,8.7,0.01)
beta.one.sims <- seq(0.2,1.6,0.01)
pL.beta0 <- apply(X = data.frame(theta_hat$par[1], beta.zero.sims,theta_hat$par[3]), MARGIN = 1 , FUN = nll.wei, time = actg$time, event = actg$event, treatment = actg$tx)
pL.beta1 <- apply(X = data.frame(theta_hat$par[1], theta_hat$par[2],beta.one.sims), MARGIN = 1 , FUN = nll.wei, time = actg$time, event = actg$event, treatment = actg$tx)

par(mfrow=c(1,2))
plot(beta.zero.sims
     , -(pL.beta0+max(-pL.beta0))
     , "l"
     ,main = "Profile likelihood for Beta 0"
     ,xlab = expression(beta[0])
     ,ylab = "lp - max(lp)"
     ,ylim = c(-3,0))
abline(h = -qchisq(0.95, df = 1)/2, col = 2)
plot(beta.one.sims
     , -(pL.beta1+max(-pL.beta1))
     , "l"
     ,main = "Profile likelihood for beta 1"
     ,xlab = expression(beta[1])
     ,ylab = "lp - max(lp)"
     ,ylim = c(-3,0))
abline(h = -qchisq(0.95, df = 1)/2, col = 2)

#CI:
sd <- as.numeric(sqrt(diag(solve(hessian(theta_hat$par, func = nll.wei)))))

#Wald 95% CIs and profile-likelihoods with approx 95% CI
#Der er et eller andet ved de her wald CI som ikke fungerer?
(Wald.CI <- theta_hat$par + matrix(c(-1,1), 3,2, byrow = T) * matrix(qnorm(0.975)*sd, 3,2, byrow = F))

#Direkte numerisk approksimation:
(CI.0 <- c(min(beta.zero.sims[-(pL.beta0+max(-pL.beta0)) > -qchisq(0.95, df = 1)/2])
           ,max(beta.zero.sims[-(pL.beta0+max(-pL.beta0)) > -qchisq(0.95, df = 1)/2])))
(CI.1 <- c(min(beta.one.sims[-(pL.beta1+max(-pL.beta1)) > -qchisq(0.95, df = 1)/2])
           ,max(beta.one.sims[-(pL.beta1+max(-pL.beta1)) > -qchisq(0.95, df = 1)/2])))

plot(beta.zero.sims
     , -(pL.beta0+max(-pL.beta0))
     , "l"
     ,main = "Profile likelihood for Beta 0"
     ,xlab = expression(beta[0])
     ,ylab = "lp - max(lp)"
     ,ylim = c(-3,0))
abline(h = -qchisq(0.95, df = 1)/2, col = 2)
abline(v = Wald.CI[2,], col = 6)
text(x = Wald.CI[2,1]+.15, y = -2.4, "Wald CI", col = 6)
text(x = CI.0[1]+.05, y = -2.5, "CI", col = 2)
abline(v = c(CI.0), lty = "dashed", col = 2)
plot(beta.one.sims
     , -(pL.beta1+max(-pL.beta1))
     , "l"
     ,main = "Profile likelihood for beta 1"
     ,xlab = expression(beta[1])
     ,ylab = "lp - max(lp)"
     ,ylim = c(-3,0))
abline(h = -qchisq(0.95, df = 1)/2, col = 2)
abline(v = Wald.CI[3,], col = 6)
text(x = Wald.CI[3,1]+0.25, y = -1.6, "Wald CI", col = 6)
text(x = CI.1[1]+0.09, y = -1.7, "CI", col = 2)
abline(v = c(CI.1),lty = "dashed", col = 2)

#### Er rimelig sikker på at alt herunder ikke nødvendigvis er korrekt, og at man i stedet bør anvende noget ####
#kaplan-meier.
#BØR MAN HER ANVENDE DET FULDE DATASÆT?#
#BØR MAN HER ANVENDE DET FULDE DATASÆT?#
#BØR MAN HER ANVENDE DET FULDE DATASÆT?#

#set theta_hat_0 = exp(beta_0)
#    theta_hat_1 = exp(beta_0 + beta_1)

#calculate likelihood
#define theta as E[time_tx]. Thus, theta = 1/lambda
nll.theta_0 <- function(beta0, x = actg_event$time[actg_event$tx == 0]){
  theta = exp(beta0)
  nll <- -sum(dexp(x, rate =  1/theta, log = T))
  return(nll)
}
opt.theta_0 <- nlminb(start = 1, objective = nll.theta_0)
(beta_0_hat <- opt.theta_0$par)

nll.theta_1 <- function(beta1, beta0 = beta_0_hat, x = actg_event$time[actg_event$tx == 1]){
  theta = exp(beta0 + beta1)
  nll <- -sum(dexp(x, rate =  1/theta, log = T))
  return(nll)
}
opt.theta_1 <- nlminb(start = 1, objective = nll.theta_1)
(beta_1_hat <- opt.theta_1$par)

#calculate estimated thetas
(theta_0_hat <- exp(beta_0_hat))
(theta_1_hat <- exp(beta_0_hat + beta_1_hat))

#### Bullet point 6: WALD CI ####
#Profile likelihoods:
beta.zero.sims <- seq(4.43,5,0.01)
beta.one.sims <- seq(-0.5,0.23,0.01)
pL.beta0 <- sapply(X = beta.zero.sims, FUN = nll.theta_0)
pL.beta1 <- sapply(X = beta.one.sims, FUN = nll.theta_1, beta0 = beta_0_hat)
par(mfrow=c(1,2))
plot(beta.zero.sims
     , -(pL.beta0+max(-pL.beta0))
     , "l"
     ,main = "Profile likelihood for Beta_0")
abline(h = -qchisq(0.95, df = 1)/2, lty = "dashed")
plot(beta.one.sims
     , -(pL.beta1+max(-pL.beta1))
     , "l"
     ,main = "Profile likelihood for Beta_1")
abline(h = -qchisq(0.95, df = 1)/2, lty = "dashed")

#From these figures it can be concluded that the quadratic approximation
#of the CI through use of fischers information matrix, is a 
#good approksimation.

sd_0 <- as.numeric(sqrt(diag(solve(hessian(beta_0_hat, func = nll.theta_0)))))
sd_1 <- as.numeric(sqrt(diag(solve(hessian(beta_1_hat, func = nll.theta_1, beta0 = beta_0_hat)))))

#Wald 95% CIs and profile-likelihoods with approx 95% CI
(W.CI.0 <- beta_0_hat + c(-1,1)*qnorm(0.975)*sd_0)
(W.CI.1 <- beta_1_hat + c(-1,1)*qnorm(0.975)*sd_1)

#Direkte numerisk approksimation:
(CI.0 <- c(min(beta.zero.sims[-(pL.beta0+max(-pL.beta0)) > -qchisq(0.95, df = 1)/2])
           ,max(beta.zero.sims[-(pL.beta0+max(-pL.beta0)) > -qchisq(0.95, df = 1)/2])))
(CI.1 <- c(min(beta.one.sims[-(pL.beta1+max(-pL.beta1)) > -qchisq(0.95, df = 1)/2])
           ,max(beta.one.sims[-(pL.beta1+max(-pL.beta1)) > -qchisq(0.95, df = 1)/2])))

#### Bullet point 7: Theoretical results, SDs and Profile likelihood ####
plot(beta.zero.sims
     , -(pL.beta0+max(-pL.beta0))
     , "l"
     ,main = "Profile likelihood for Beta 0"
     ,xlab = expression(beta[0])
     ,ylab = "lp - max(lp)")
abline(h = -qchisq(0.95, df = 1)/2, col = 2)
abline(v = c(W.CI.0), col = 6)
text(x = W.CI.0[1]+.1, y = -2.3, "Wald CI", col = 6)
text(x = CI.0[1]+.05, y = -2.5, "CI", col = 2)
abline(v = c(CI.0), lty = "dashed", col = 2)
plot(beta.one.sims
     , -(pL.beta1+max(-pL.beta1))
     , "l"
     ,main = "Profile likelihood for beta 1"
     ,xlab = expression(beta[1])
     ,ylab = "lp - max(lp)")
abline(h = -qchisq(0.95, df = 1)/2, col = 2)
abline(v = c(W.CI.1), col = 6)
text(x = W.CI.1[1]+0.15, y = -1.5, "Wald CI", col = 6)
text(x = CI.1[1]+0.07, y = -1.7, "CI", col = 2)
abline(v = c(CI.1),lty = "dashed", col = 2)

#Theoretical results: The effect of the treatment is not statistically significant based on above results.