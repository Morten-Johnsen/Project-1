rm(list = ls())
library(tidyverse)


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

#results show: -0.72 logits? for developing AIDS when using the treatment

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

# Hvordan kan man sige at de her værdier beskriver forskellen?
#medmindre det på en eller anden måde skal være for p?

#### Analysis of the Survival Data ####
actg320 <- read.table("actg320.txt", header=TRUE, sep="", 
                as.is=TRUE)

#select time, event and tx as they are the only relevant variables in this project
actg <- actg320 %>%
  select(time, event, tx)

#Bulletpoint 2) in the project description (MISSING: "Other relevant number that could be calculated?")
actg %>%
  group_by(tx) %>%
  summarise("Got AIDS or DIED" = sum(event),
            "Proportion" = sum(event)/n(),
            "Participants Given the Treatment" = n())

#Bulletpoint 3)
#Fitting an exponential model to time for both and for each treatment
#
library(survival)
library(survminer)
kaplan.meier <- survfit(Surv(time) ~ tx, data = filter(actg, event == 1))
ggsurvplot_add_all(kaplan.meier, data = filter(actg, event == 1), conf.int = T
                   ,risk.table = T, pval = T)

both <- nlminb(start = 2
       , objective = testDistribution
       , x = actg$time[actg$event == 1]
       , distribution = "exponential")

#separate exponential models
t1 <- nlminb(start = 2
               , objective = testDistribution
               , x = filter(actg, tx == 1)$time
               , distribution = "exponential")

t0 <- nlminb(start = 2
             , objective = testDistribution
             , x = filter(actg, tx == 0)$time
             , distribution = "exponential")

ggplot(actg)+
  geom_histogram(aes(x = time, y = ..density..))+
  stat_function(fun = dexp, n = dim(actg)[1], args = list(rate = both$par))

ggplot(actg)+
  geom_histogram(aes(x = time, fill = factor(tx)), colour = "white")

#### bullet point 4: Compare likelihoods
#Anvende metode fra bogen som vi snakkede om sidste gang.