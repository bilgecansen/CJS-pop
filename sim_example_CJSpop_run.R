
# Example script to run CJS-pop on simulated data

## See sim_models_notebook for extended notes on model structure.
## sim_models_notebook also creates additional directories.

rm(list = ls())

# Set up initial parameters -----------------------------------------------

## Code uses parallel computing. 
## set cores to 1 if not sure about the number of cores 
cores <- 7

## Iterations of MCMC
ni <- 50000 

## burn-in
nb <- 20000

## Number of chains 
## Jags.parallel function is used so each chain is run on a seperate core
## total core useage is cores*nc
nc <- 4

## Thinning
nt <- 20


# Load packages -----------------------------------------------------------

library(rjags)
library(R2jags)
library(MCMCvis)
library(magrittr)


# Load data and set up jags functions -------------------------------------

chdata <- readRDS("chdata_m_p10_k150.rds")

# ch.init function initilizes values for the latent state z, using the latent_state matrix (see
# data_generation notebook for the creation of this matrix). NA values in the latent state is given
# an initial value of 1 (Years before first or after the last capture of an individual). 
# Observations between the first and last year of captures are given NA as initial values 
# because they are provided as data below. init.func initilizes other model paramters 
# and objects usually byt random draws that correspond to their priors.

ch.init <- function(latent_state, first, last_pop, pop) {
  ch <- ifelse(is.na(latent_state),1,NA)
  for (i in 1:nrow(ch)) {
    ch[i,1:first[i]] <- NA
    if(last_pop[pop[i]]<ncol(ch)) ch[i,(last_pop[pop[i]]+1):ncol(ch)] <- NA
  }
  return(ch)    
}

jags.func <- function(chdata, ni, nb, nc, nt) {
  
  init.func <- function() {
    list(survival_ad = runif(1,0,1),
         beta = rnorm(1,0,10),
         gamma = rnorm(2,0,10),
         fecundity = runif(1,1,10),
         zeta = rnorm(1,0,3),
         sigma_s = runif(1,0,5),
         sigma_f = runif(1,0,5),
         z = ch.init(chdata$latent_state, chdata$first, 
                     chdata$last_pop, chdata$pop))
  }
  
  params <- c("survival_ad", "survival_juv", "beta", "gamma", "fecundity", "zeta", 
              "sigma_f", "sigma_s","pvalue.ch", "pvalue.fec")
  
  data <- list(y = chdata$ch,
               z = chdata$latent_state,
               Nobs_ad = chdata$Nobs_ad,
               Nobs_juv1 = chdata$Nobs_juv,
               Nobs_juv2 = chdata$fec_data,
               pop = chdata$pop,
               npop = chdata$npop,
               pop.index = chdata$pop_index,
               year.index = chdata$year_index,
               stage = chdata$stage,
               first = chdata$first,
               first_pop = chdata$first_pop,
               first_sub = chdata$first_sub,
               nind = chdata$nind,
               nyear = chdata$nyear,
               nsub = 4,
               nfec = chdata$nfec)
  
  model_type <- c("dd_model_sim.jags")
  
  results <- jags.parallel(model.file = model_type,
                           working.directory = getwd(),
                           data = data,
                           inits = init.func,
                           parameters.to.save = params,
                           n.iter = ni,
                           n.chains = nc,
                           n.burnin = nb,
                           n.thin = nt,
                           DIC = F)
  
  return(results)
}


# Run the jags model ------------------------------------------------------

cl <- makeCluster(cores, types = "SOCK")
registerDoSNOW(cl)

system.time({
  results <- foreach(i=1:length(chdata), .packages = c("rjags", "R2jags")) %dopar% {
    setwd(code_directory)
    jags.func(chdata[[i]], ni = ni, nb = nb, nc = nc, nt = nt)
  }
})

stopCluster(cl)

if (!"results" %in% list.files()) dir.create("results")
saveRDS(results, "ddresults_h_p10_k150.rds")
