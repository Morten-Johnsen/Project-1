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

x1 <- log.data$AIDS_yes[1]
x2 <- log.data$AIDS_yes[2]
n1 <- log.data$n[2]
n2 <- log.data$n[2]

## Likelihood analysis (Den her del er formodentligt ikke færdig)
(theta.hat <- log(x1/(n1-x1)*(n2-x2)/x2))

eta.start <- log(x2/(n2-x2))

L_theta <- function(eta){ 
  return(-log(exp(theta.hat*x1+eta*(x1+x2))/((1+exp(theta.hat+eta))^n1*(1+exp(eta))^n2)))
}
(eta.hat <- nlminb(eta.start, objective = L_theta)$par)

#Now we can construct CIs the 'normal' way based on the Likelihood function
L <- function(theta){
  return(-log(exp(theta*x1+eta.hat*(x1+x2))/((1+exp(theta+eta.hat))^n1*(1+exp(eta.hat))^n2)))
}

library(numDeriv)
var <- solve(hessian(func <- L, x = theta.hat))
sd <- as.numeric(sqrt(var))

CI <- p.hat + c(-1,1)*qnorm(c(0.975))*sd/sqrt(sum(log.data$n))

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
#Fitting an exponential model to time for both and for each treatments
both <- nlminb(start = 1
       , objective = testDistribution
       , x = actg$time
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
