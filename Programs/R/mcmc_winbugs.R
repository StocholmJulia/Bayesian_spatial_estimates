#########################################################################
## This R scrip contains all nessassery code for 
## calling Winbugs from R and do MCMC estimation 


## created 2023-3-31 by Yuxin Huang 
##########################################################################

## load package

library(coda)
library(haven)
library(dplyr)
library(R2WinBUGS)
## load dataset save previous estimation (Stata: lle_est.do and cif_est.do)

cancer_spatial <- read_dta("C:/Users/n11117761/Work/Data/Analysis/mergeData_winbugs.dta")
nblist <- read_dta("C:/Users/n11117761/Work/StataCode_Yuliya/adj_matrix.dta")

## data preperation
# rename the variable in dataset as Winbugs don't like variable start with "_"
cancer_spatial <- cancer_spatial %>%
  rename("rcs1" = "_rcs1",
         "rcs2" = "_rcs2",
         "rcs3" = "_rcs3",
         "rcs4" = "_rcs4",
         "rcs5" = "_rcs5",
         "drcs1" = "_d_rcs1",
         "drcs2" = "_d_rcs2",
         "drcs3" = "_d_rcs3",
         "drcs4" = "_d_rcs4",
         "drcs5" = "_d_rcs5",
         "t" = "_t"
  )

############################# Data ########################################
## create dataset for Winbugs 
##(acually at the end I used dataset obtained externally: Data.txt)
## but it always more convinent to create data in R


bugs.dat <- list(
  
  # fix values for loops
  datarows = length(cancer_spatial$id),   # number of observations
  Nsla = 478,                             # number of grids
  sumNumNeigh = 2724,              # total neighbors (2724)
  
  # values for spatial data
  adj = as.numeric(nblist$adj_grid),                  # the index of neighbor
  num = as.numeric(nblist$num),
  grid = as.numeric(cancer_spatial$grid),# number of neighbors
  #cum = c(cumsum(nblist$num) - nblist$num, sum(nblist$num)), # cumulative sum
  
  # variables
  d = as.numeric(cancer_spatial$death),  # exist indicator (1,0)
  rcs1 = as.numeric(cancer_spatial$rcs1),
  rcs2 = as.numeric(cancer_spatial$rcs2),
  rcs3 = as.numeric(cancer_spatial$rcs3),
  rcs4 = as.numeric(cancer_spatial$rcs4),
  rcs5 = as.numeric(cancer_spatial$rcs5),
  drcs1 = as.numeric(cancer_spatial$drcs1),
  drcs2 = as.numeric(cancer_spatial$drcs2),
  drcs3 = as.numeric(cancer_spatial$drcs3),
  drcs4 = as.numeric(cancer_spatial$drcs4),
  drcs5 = as.numeric(cancer_spatial$drcs5),
  cage1 = as.numeric(cancer_spatial$cage1),
  cage2 = as.numeric(cancer_spatial$cage2),
  
  rate = as.numeric(cancer_spatial$rate),
  t = as.numeric(cancer_spatial$t)
)

########################### Inits #########################################

inits_fix_gamma05 <- function() {
  list(
    tauu=450,
    tauv=450,
    gamma=c(-1.305, 1.034, 0.169, -0.031, -0.046, -0.032), 
    beta=c( 0.359, -0.059),
    u = rep(0,478),
    v = rnorm(478,mean = 0,sd = 1)
  )
  
}


########################### MCMC  ################################

M.burnin <- 50000       # Number of burn-in iterations (discarded)
M <- 20000             # Number of iterations retained (after thinning)
n.thin <- 5           # Thinning factor


MCMCGamma_final <- bugs(
  data = "C:/Users/n11117761/Work/Data/WinBugs/Data.txt",
  inits = inits_fix_gamma05,
  parameters.to.save = c("gamma", "beta","u","v","ranuv","fracspatial"),
  model.file = "C:/Users/n11117761/Work/Data/WinBugs/Model_YH_gamma0.5.txt",
  n.chains = 2,
  n.burnin = M.burnin,
  n.iter = M.burnin + (M * n.thin), # Total iterations
  n.thin = n.thin,
  bugs.directory = "C:/Users/n11117761/WinBUGS14",
  DIC = FALSE,
  debug = TRUE
)
# Store the list of parameters for convenience (ignore the other stuff in MCMC)
pars <- MCMC5$sims.list
# save bug object
saveRDS(MCMCGamma_final, file = "mcmc_g5.rds")