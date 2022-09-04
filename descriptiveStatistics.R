#setwd("Replace with path to directory containing project files.")
D <- read.table("tuno.txt", header=TRUE, sep=" ", 
                as.is=TRUE)
## Dimensions of D (number of rows and columns)
dim(D)
##  Column/variable names
names(D)
## The first rows/observations
head(D)
## The last rows/observations
tail(D)
## Selected summary statistics
summary(D)
## Another type of summary of the dataset
str(D)

hist(D$pow.obs, xlab="Power production pr day", prob=TRUE)
plot(D$r.day, D$pow.obs, type = 'l', xlab="Date", ylab="Power production pr. day", col=1)

plot(D$r.day, D$pow.obs, type="l", xlim=as.Date(c("2003-01-01","2003-10-31")), 
     ylim=c(0,9), xlab="Date", ylab="Power production pr. day", col=2)


