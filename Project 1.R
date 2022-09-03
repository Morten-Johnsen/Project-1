rm(list = ls())
library(tidyverse)
library(reshape)


test.data <- data.frame("id" = sample(c("1","2","3","4"), size = 100, replace = T, prob = c(0.1,0.4,0.2,0.3))
                        ,"N" = rnorm(100))
ggplot(test.data)+
  geom_histogram(aes(x = N, fill = id), colour = "white")+
  theme_linedraw()+
  labs(x = "", y = "", fill = "ID")+
  ggtitle("Test hist")
