#kaplan.meier survival curves function
library(tidyverse)

actg320 <- read.table("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/actg320.txt", header=TRUE, sep="", 
                      as.is=TRUE)

#select time, event and tx as they are the only relevant variables in this project
actg <- actg320 %>%
  dplyr::select(time, event, tx, cd4)
#actg$id <- 1:dim(actg)[1]

kaplan.meier <- function(data = actg, event = "event", time = 'time'){
  
  data %>% 
    select(time, event) -> data_selected
  
  data_selected %>%
    dplyr::rename( t = time
                   ,d = event) %>%
    arrange(t) %>%
    mutate(R = n():1
           ,censored = 1-d) %>%
    group_by(t) %>%
    summarise(d = sum(d), R = max(R), censored = sum(censored)) %>%
    mutate(EventsPerAvailable = 1 - d/R
           ,S = cumprod(EventsPerAvailable)
           ,logS = log(S)
           ,Varloglog = 1/logS^2 * cumsum(d/(R*(R-d)))
           ,loglog_minus = log(-logS) - qnorm(0.975) * sqrt(Varloglog)
           ,loglog_plus = log(-logS) + qnorm(0.975) * sqrt(Varloglog)
           ,lower = exp(-exp(loglog_plus))
           ,upper = exp(-exp(loglog_minus))) %>%
    select(-loglog_minus, -loglog_plus, -logS, -Varloglog, -EventsPerAvailable)-> init.table
  
  return(init.table)
}

#edit room


#edit room end



kaplan.meier()
