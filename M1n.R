require(rjags)
require(runjags)
require(coda)
library(tidyverse)
library(haven)
library(psych) 
library (GGally)
library(dplyr)
library(foreign)



## reading in the data

FRAMINGHAM <- read_dta('FRAMINGHAM.dta ')

final <- FRAMINGHAM %>% 
  select(pp, age, sbp, chol)

## Structure and summary of the data
str(final)
summary(final)
describe(final)

## looking at the correlations and scatter plot
ggpairs(final)

# variables
y = final$pp
x1 = final$sbp
x2 = final$age
x3 = final$chol
nn = length(y)

dataList = list(  # Put the information into a list.
  y = y, x1 = x1, x2 = x2, x3 = x3, Ntotal = nn
)

# Define the model with no individual differences
modelString = "
data {
  ymean <- mean(y); ysd <- sd(y)
  x1mean <- mean(x1); x1sd <- sd(x1)
  x2mean <- mean(x2); x2sd <- sd(x2)
  x3mean <- mean(x3); x3sd <- sd(x3)
  for (i in 1:Ntotal){
    zy[i] <- (y[i] - ymean)/ysd
    z1[i] <- (x1[i] - x1mean)/x1sd
    z2[i] <- (x2[i] - x2mean)/x2sd
    z3[i] <- (x3[i] - x3mean)/x3sd
  }
}
model{
  for (i in 1:Ntotal){
    zy[i] ~ dnorm(zbeta0  + zbeta1*z1[i] + zbeta2*z2[i] + zbeta3*z3[i], 1/zsigma^2)
  }
  # Vague priors
  zbeta0 ~ dnorm(0,1/2^2)  # mean 0, sd = 2
  zbeta1 ~ dnorm(0,1/2^2) 
  zbeta2 ~ dnorm(0,1/2^2) 
  zbeta3 ~ dnorm(0,1/2^2) 
  zsigma ~ dunif(0.00001,.99999)
  #nu ~ dexp(.0333) # mean 30, sd = 30
  #sigmaBeta ~ dgamma(2.618, 1.618) # mode 1, sd 1.0
  
  # Transform to original scale
  beta0 <- zbeta0*ysd + ymean - zbeta1*ysd*x1mean/x1sd - zbeta2*ysd*x2mean/x2sd - zbeta3*ysd*x3mean/x3sd
  betaSbp <- zbeta1*ysd/x1sd 
  betaAge <- zbeta2*ysd/x2sd 
  betaChol <- zbeta3*ysd/x3sd 
  sigma <- zsigma*ysd
}
" # close quote for modelString
writeLines(modelString, con="MultipleRegression-M1h.txt" )

# MCMC run
myinits <-list(list(zbeta0=rnorm(1,0,1), zbeta1=rnorm(1,0,1), zbeta2=rnorm(1,0,1),zbeta3=rnorm(1,0,1),
                    zsigma=runif(1)),
               list(zbeta0=rnorm(1,0,1), zbeta1=rnorm(1,0,1), zbeta2=rnorm(1,0,1),zbeta3=rnorm(1,0,1),
                    zsigma=runif(1)),
               list(zbeta0=rnorm(1,0,1), zbeta1=rnorm(1,0,1), zbeta2=rnorm(1,0,1),zbeta3=rnorm(1,0,1),
                    zsigma=runif(1)))
out <- run.jags(model="MultipleRegression-M1h.txt",data=dataList, inits=myinits, n.chains=3,
                adapt=500,burnin=2000,sample=7000, 
                monitor=c("zbeta0","zbeta1","zbeta2","zbeta3","beta0","betaSbp","betaAge","betaChol","zsigma","DIC"))
print(out)
plot(out)

# Extract MCMC samples
aa0 <- as.mcmc.list(out, vars="beta0")
b0 <- as.numeric(unlist(aa0[,2]))
aa1 <- as.mcmc.list(out, vars="betaSbp")
b1 <- as.numeric(unlist(aa1[,1]))
aa2 <- as.mcmc.list(out, vars="betaAge")
b2 <- as.numeric(unlist(aa2[,1]))
aa3 <- as.mcmc.list(out, vars="betaChol")
b3 <- as.numeric(unlist(aa3[,1]))

# Posterior predictive
xp1 <- seq(92,240,1) # 92-240
xp2 <- seq(20,66,1) # 20 - 66
xp3 <- seq (100,386,1) # 100 - 386

mm <- 20
plot(x1, y,xlim=range(c(92,240)),ylim=range(c(15,150)),xlab="x1 - Sbp", ylab="y - PP", 
     main = "Posterior Predictive Check",col="black", pch=20, cex = 1)
for (i in 1:mm){
  y.prd1 <- b0[i*50] + b1[i*50]*xp1 + b2[i*50]*mean(x2) + b3[i*50]*mean(x3)
  lines(xp1, y.prd1, lwd=1, col='blue')
}

plot(x2, y,ylim=range(c(15,150)),xlab="x2 - Age", ylab="y - PP", 
     main = "Posterior Predictive Check",col="black", pch=20, cex = 1)
for (i in 1:mm){
  y.prd2 <- b0[i*50] + b1[i*50]*mean(x1) + b2[i*50]*xp2 + b3[i*50]*mean(x3)
  lines(xp2, y.prd2, lwd=1, col='blue')
}


plot(x3, y,ylim=range(c(15,150)),xlab="x3 - Chol", ylab="y - PP", 
     main = "Posterior Predictive Check",col="black", pch=20, cex = 1)
for (i in 1:mm){
  y.prd3 <- b0[i*50] + b1[i*50]*mean(x1) + b2[i*50]*mean(x2) + b3[i*50]*xp3
  lines(xp3, y.prd3, lwd=1, col='blue')
}


# Posterior predictive
mmm <- 50
pvaf <- rep(0,mmm)
for (i in 1:mmm){
  y.prd <- b0[i*50] + b1[i*50]*x1 + b2[i*100]*x2 + b3[i*50]*x3
  pvaf[i] <- 1-sum((y - y.prd)^2)/sum((y - mean(y))^2)
}
cat('Mean percent variance accounted for (PVAF): ', mean(pvaf))
# Normality check for residuals
residu <- mean(b0)+mean(b1)*x1+mean(b2)*x2 + mean(b3)*x3 - y
qqnorm(residu);qqline(residu)
hist(residu, breaks=11)








