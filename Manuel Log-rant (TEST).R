#Comparing Survival Functions
source("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/kaplan-meier-manual.R")

#group
group = "tx"
data = actg
event = "event"
time = "time"
head(data)

#data$id <- 1:dim(data)[1]
groups <- unique(data[[group]])

outs <- list()
e <- list()

outs[["none"]] <- kaplan.meier(data = data, event = event, time = time)

for (i in groups){
  outs[[paste(i)]] <- kaplan.meier(data = filter(data, !!sym(group) == i), event = event, time = time)
  
  outs$none %>%
    left_join(outs[[paste(i)]], by = 't') %>%
    group_by('t') %>%
    summarise(e = sum(R.y * d.x/R.x, na.rm = T)) -> e[[paste(i)]]
    #mutate(e = R.y * d.x / R.x)-> e[[paste(i)]]
}

#test for groups == 1
w <- 1

e$`0` %>%
  mutate(v = e*((R.x - d.x)/R.x)*((R.x-R.y)/(R.x - 1))) -> e$`0`
Q0 <- sum(w*(e$`0`$d.y - e$`0`$e), na.rm = T)^2 / sum(w^2 * e$`0`$v, na.rm = T)
cat("Observed: ", sum(e$`0`$d.y, na.rm = T), ", Expected: ", sum(e$`0`$e, na.rm = T), ", (O-E)^2/E =", sum((e$`0`$d.y - e$`0`$e)^2 / e$`0`$e, na.rm = T))

 e$`1` %>%
  mutate(v = e*((R.x - d.x)/R.x)*((R.x-R.y)/(R.x - 1))) -> e$`1`
Q1 <- sum(w*(e$`1`$d.y - e$`1`$e))^2 / sum(w^2 * e$`1`$v, na.rm = T)
cat("Observed: ", sum(e$`1`$d.y), ", Expected: ", sum(e$`1`$e), ", (O-E)^2/E =", sum((e$`1`$d.y - e$`1`$e)^2 / e$`1`$e, na.rm = T))

library(survival)
fit <- survfit(Surv(time, event) ~ tx, data = actg)
fit$surv
fit$time

plot(fit$time[1:228], fit$surv[1:228], "l", lwd = 5)
grid()
lines(fit$time[229:(228+205)], fit$surv[229:(228+205)], "l", lwd = 5)
lines(fit$time[229:(228+205)], fit$lower[229:(228+205)], "l", lwd = 5)
lines(fit$time[229:(228+205)], fit$upper[229:(228+205)], "l", lwd = 5)
lines(fit$time[1:(228)], fit$lower[1:(228)], "l", lwd = 5)
lines(fit$time[1:(228)], fit$upper[1:(228)], "l", lwd = 5)
lines(outs$`0`$t, outs$`0`$S, col = "yellow")
lines(outs$`1`$t, outs$`1`$S, col = "yellow")
lines(outs$`1`$t, outs$`1`$lower, lty = "dashed")
lines(outs$`1`$t, outs$`1`$upper, lty = "dashed")
lines(outs$`0`$t, outs$`0`$lower, lty = "dashed")
lines(outs$`0`$t, outs$`0`$upper, lty = "dashed")

