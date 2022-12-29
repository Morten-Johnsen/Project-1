#Comparing Survival Functions
if (Sys.getenv("LOGNAME") == "mortenjohnsen"){
  source("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/kaplan-meier-manual.R")
} else {
  source("~/Documents/02418 Statistical Modelling/Assignments/Assignment 1/Project-1/kaplan-meier-manual.R")
}

logRank <- function(data = actg, group = "tx", event = "event", time = "time"){
  
  groups <- unique(data[[group]])
  
  outs <- list()
  e <- list()
  
  outs[["none"]] <- kaplan.meier(data = data, event = event, time = time)
  
  #Iterate over groups and merge on timestep
  for (i in groups){
    outs[[paste(i)]] <- kaplan.meier(data = filter(data, !!sym(group) == i), event = event, time = time)
    
    tmp <- outs$none %>%
      left_join(outs[[paste(i)]], by = 't')
    
    #Udfylder alle NA værdier i grupperne (0 eller 1) med værdien fra det forudgående tidsskridt, da der i tidsskridtene med
    #NA værdier ikke sker nogen ændring.
    for (j in dim(tmp)[1]:1){
      if(is.na(tmp[j,8])){
        tmp[j,-(1:7)] <- tmp[j+1,-(1:7)]
        tmp$d.y[j] <- 0
        tmp$censored.y[j] <- 0
      }
    }
    
    tmp %>%
      mutate(e = R.y * d.x/R.x) -> e[[paste(i)]] #slide 32 uge 36
    #mutate(e = R.y * d.x / R.x)-> e[[paste(i)]]
  }
  
  w <- 1 
  #}
  
  e$`0` %>%
    mutate(v = e*((R.x - d.x)/R.x)*((R.x-R.y)/(R.x - 1))) -> e$`0`
  Q0 <- sum(w*(e$`0`$d.y - e$`0`$e), na.rm = T)^2 / sum(w^2 * e$`0`$v, na.rm = T)
  
  e$`1` %>%
    mutate(v = e*((R.x - d.x)/R.x)*((R.x-R.y)/(R.x - 1))) -> e$`1`
  Q1 <- sum(w*(e$`1`$d.y - e$`1`$e))^2 / sum(w^2 * e$`1`$v, na.rm = T)
  
  cat("",group,": ",0, " | "
      ,"n = ", max(e$`0`$R.y)
      , "Observed: ", sum(e$`0`$d.y, na.rm = T)
      , ", Expected: ", sum(e$`0`$e, na.rm = T)
      , ", (O-E)^2/V =", (sum(e$`0`$d.y) - sum(e$`0`$e))^2/sum(e$`0`$v)
      , ", Q = ", Q0
      ,"\n",group,": ",1, " | "
      ,"n = ", max(e$`1`$R.y)
      , "Observed: ", sum(e$`1`$d.y)
      , ", Expected: ", sum(e$`1`$e)
      , ", (O-E)^2/V =", (sum(e$`1`$d.y) - sum(e$`1`$e))^2/sum(e$`1`$v)
      , ", Q = ", Q1
      , "\n\n\np-value for chi-squared test for q = ",Q0," ", "with df = ", 1," : ", 1-pchisq(q = Q0, df = 1))
  
  
}

# 
# library(survival)
# fit <- survfit(Surv(time, event) ~ tx, data = actg)
# fit$surv
# fit$time
# 
# plot(fit$time[1:228], fit$surv[1:228], "l", lwd = 5)
# grid()
# lines(fit$time[229:(228+205)], fit$surv[229:(228+205)], "l", lwd = 5)
# lines(fit$time[229:(228+205)], fit$lower[229:(228+205)], "l", lwd = 5)
# lines(fit$time[229:(228+205)], fit$upper[229:(228+205)], "l", lwd = 5)
# lines(fit$time[1:(228)], fit$lower[1:(228)], "l", lwd = 5)
# lines(fit$time[1:(228)], fit$upper[1:(228)], "l", lwd = 5)
# lines(outs$`0`$t, outs$`0`$S, col = "yellow")
# lines(outs$`1`$t, outs$`1`$S, col = "yellow")
# lines(outs$`1`$t, outs$`1`$lower, lty = "dashed")
# lines(outs$`1`$t, outs$`1`$upper, lty = "dashed")
# lines(outs$`0`$t, outs$`0`$lower, lty = "dashed")
# lines(outs$`0`$t, outs$`0`$upper, lty = "dashed")

