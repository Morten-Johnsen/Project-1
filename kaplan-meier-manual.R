#kaplan.meier survival curves function
library(tidyverse)

WHAS <- read.delim(file = "Downloads/WHAS.txt") 


kaplan.meier <- function(data = WHAS, event = 'fstat', time = 'lenfol'){
  data %>% 
    select(id, time, event) -> data_selected
  
  data_selected %>%
    dplyr::rename( t = time
                   ,d = event) %>%
    arrange(t) %>%
    mutate(R = n():1
          ,censored = 1-d
          ,EventsPerAvailable = 1 - d/R
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

# Old: 
# WHAS %>% 
#   select(lenfol, age, bmi, fstat) -> WHAS_selected
# 
# WHAS_selected %>%
#   dplyr::rename( time = lenfol
#                 ,d = fstat) %>%
#   arrange(time) %>%
#   mutate(R = n():1
#          ,censored = 1-d
#          ,EventsPerAvailable = 1 - d/R
#          ,S = cumprod(EventsPerAvailable)
#          ,logS = log(S)
#          ,Varloglog = 1/logS^2 * cumsum(d/(R*(R-d)))
#          ,loglog_minus = log(-logS) - qnorm(0.975) * sqrt(Varloglog)
#          ,loglog_plus = log(-logS) + qnorm(0.975) * sqrt(Varloglog)
#          ,lower = exp(-exp(loglog_plus))
#          ,upper = exp(-exp(loglog_minus))) %>%
#   select(-loglog_minus, -loglog_plus, -logS, -Varloglog, -EventsPerAvailable)-> init.table
# 
# init.table