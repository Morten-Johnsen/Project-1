rm(list = ls())

if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/Statistical Modelling/Project-1/")
} else {
  print("Fejl: Viktor, Indsæt lokationen på din egen folder her")
}

#### Analysis of Logistic Data ####

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
