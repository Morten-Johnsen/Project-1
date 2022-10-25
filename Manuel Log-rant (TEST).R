#Comparing Survival Functions
source("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/9. Semester/02418 - Statistical Modelling/Project-1/kaplan-meier-manual.R")

#group
group = "tx"
data = actg
event = "event"
time = "time"
head(data)

data$id <- 1:dim(data)[1]
groups <- unique(data[[group]])

outs <- list()
e <- list()

outs[["none"]] <- kaplan.meier(data = data, event = event, time = time)

for (i in groups){
  outs[[paste(i)]] <- kaplan.meier(data = filter(data, !!sym(group) == i), event = event, time = time)
  
  outs$none %>%
    inner_join(outs[[paste(i)]], by = 'id') %>%
    mutate(e = R.y * d.x / R.x)-> e[[paste(i)]]
}

#test for groups == 1
w <- 1

e$`0` %>%
  mutate(v = e*((R.x - d.x)/R.x)*((R.x-R.y)/(R.x - 1))) -> e$`0`
Q0 <- sum(w*(e$`0`$d.x - e$`0`$e))^2 / sum(w^2 * e$`0`$v, na.rm = T)
cat("Observed: ", sum(e$`0`$d.x), ", Expected: ", sum(e$`0`$e))

e$`1` %>%
  mutate(v = e*((R.x - d.x)/R.x)*((R.x-R.y)/(R.x - 1))) -> e$`1`
Q1 <- sum(w*(e$`1`$d.x - e$`1`$e))^2 / sum(w^2 * e$`1`$v, na.rm = T)
cat("Observed: ", sum(e$`1`$d.x), ", Expected: ", sum(e$`1`$e))

plot(outs$`0`$t, outs$`0`$S, "l", ylim = c(0.8,1))
grid()
lines(outs$`1`$t, outs$`1`$S)
lines(outs$`1`$t, outs$`1`$lower, lty = "dashed")
lines(outs$`1`$t, outs$`1`$upper, lty = "dashed")
lines(outs$`0`$t, outs$`0`$lower, lty = "dashed")
lines(outs$`0`$t, outs$`0`$upper, lty = "dashed")
