#########################################################################
## This R scrip contains all nessassery code for 
## creating trace plot & density plots, sa well as 
## Geweke test and Gelman test

## created 2023-3-31 by Yuxin Huang 
##########################################################################

## load package

library(coda)
library(mcmcplots)
library(bayesplot)
library(ggplot2)
library(gridExtra)

## load .rds file save from mcmc_winbugs.R

MCMC_g5 <- readRDS("mcmc_gamma05_nodflat.rds")

## convert to mcmc object for diagonosis test

post_Gamma5 <- as.mcmc.bugs(MCMC_g5)

## trace plot & density plot

plot(post_Gamma5[,c("gamma[1]","gamma[2]","gamma[3]")],
     trace = TRUE, density = TRUE, smooth = TRUE, auto.layout = TRUE, 
     type = "l",
     col=c(rgb(red = 0, green = 0, blue = 1, alpha = 0.5),rgb(red = 1, green = 0.2, blue = 0.4, alpha = 0.5)),
     xlab = "Iterations", ylab = "value")

#it can also be done by traceplot() and densplot()

## to save plots in local computer:
png("trace_gamma.png")

#plot()

dev.off

## Geweke test
geweke.diag(post_Gamma5[[1]][,"ranuv[4]"])

## Gelman test
gelman.diag(post_Gamma5[,"ranuv[58]"],
            confidence = 0.95)

## autocorrelation plot (if needed)
autocorr.plot(post_Gamma5[,"gamma[1]"],lag.max = 4000,
              main="gamma1")