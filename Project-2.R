rm(list = ls())

if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/Statistical Modelling/Project-1/")
} else {
  print("Fejl: Viktor, Indsæt lokationen på din egen folder her")
}
source("testDistribution.R")
#### Analysis of the Binary Data ####
log.data <- read.table("Logistic.txt", header=TRUE, sep="", 
                      as.is=TRUE)

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
