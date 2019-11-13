
# CJS-pop for Brown Creeper MAPS data

## See brcr_models_notebook for extended notes on model structure.

rm(list = ls())

# Set up initial parameters -----------------------------------------------

# Models used in this script

models <- c("cjspop.jags", "cjspop_nores.jags", "cjspop_di.jags", "cjs_r.jags")

## Code uses parallel computing. 
## set cores to 1 if not sure about the number of cores 
cores <- length(models)

## Iterations of MCMC
ni <- 100000 

## burn-in
nb <- 50000

## Number of chains 
## Jags.parallel function is used so each chain is run on a seperate core
## total core useage is cores*nc
nc <- 4

## Thinning
nt <- 20


# Load packages -----------------------------------------------------------

if ("pacman" %in% rownames(installed.packages())==F) 
  install.packages("pacman", repos = "http://cran.rstudio.com/")
pacman::p_load(snow, doSNOW, foreach, rjags, R2jags)


# Load data and set up jags functions -------------------------------------

chdata <- readRDS("BRCR_chdata.rds")

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

jags.func <- function(model, chdata, ni, nb, nc, nt) {
  
  init.func <- function() {
    list(survival_ad = runif(1,0,1),
         beta = rnorm(1,0,3),
         gamma = rnorm(2,0,10),
         delta = rnorm(1,0,2.5),
         fecundity = runif(1,1,10),
         zeta = rnorm(1,0,3),
         pi = runif(2,0,1),
         rho = runif(2,0,1),
         sigma_s = runif(1,0,3),
         sigma_f = runif(1,0,3),
         psi = runif(1,0,1),
         R = rep(1, nrow(chdata$stage)),
         z = ch.init(chdata$latent_state, chdata$first, 
                     chdata$last_pop, chdata$pop))
  }
  
  params <- c("survival_ad", "survival_juv", "beta", "gamma", "delta",
              "fecundity", "zeta", "pi", "rho", "sigma_f", "psi",
              "sigma_s", "pvalue_ch", "pvalue_fec", "Njuv2")
  
  data <- list(y = chdata$ch_robust,
               z = chdata$latent_state,
               r = chdata$residents,
               effort = chdata$effort_cent,
               effort_year = chdata$effort_year,
               Nobs_ad = chdata$Nobs_ad,
               Nobs_juv1 = chdata$Nobs_juv,
               Nobs_juv2 = chdata$fec_data2,
               pop = chdata$pop,
               stage = chdata$stage,
               first = chdata$first,
               first_pop = chdata$first_pop,
               last_pop = chdata$last_pop,
               first_sub = chdata$first_sub,
               pop_index = chdata$pop_index2,
               year_index = chdata$year_index2,
               nind = chdata$nind,
               nyear = chdata$nyear,
               nsub = 4,
               npop = chdata$npop,
               nfec = chdata$nfec2)
  
  results <- jags.parallel(model.file = model,
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
  results <- foreach(i=1:length(models), .packages = c("rjags", "R2jags")) %dopar% {
    jags.func(models[i], chdata, ni = ni, nb = nb, nc = nc, nt = nt)
  }
})

stopCluster(cl)

saveRDS(results, file = "BRCR_results.rds")
