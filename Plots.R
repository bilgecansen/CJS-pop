
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(MCMCvis)
library(cowplot)
library(boot)
library(foreach)
library(grid)
library(gridExtra)
library(cowplot)
library(ggthemes)

# Change the main directory as necessary
main_directory <- paste(getwd(), "paper_results", sep = "/")
setwd(main_directory)

# Data Input --------------------------------------------------------------

## Simulation Results

filenames <- list.files() %>%
  .[grepl("dd", .) | grepl("di", .)]

modelnames <- gsub("ddresults_", "dd_", filenames) %>%
  gsub(".rds", "", .) %>%
  gsub("diresults_", "di_", .) %>%
  gsub(".rds", "", .) %>%
  .[grepl("dd", .) | grepl("di", .)]

model_sum <- foreach(i=1:length(filenames)) %do% {
  model_res <- readRDS(filenames[i])
  lapply(model_res, function(x) MCMCsummary(x, n.eff = T, digits = 3))
}
names(model_sum) <- modelnames

# remove iterations with no convergence
detect_con <- function(x) {
  y <- any(x[,"Rhat"]>1.05) | any(x[,"n.eff"]<100)
  return(!y)
}

index <- map(model_sum, function(x) map_lgl(x, detect_con))
for (i in 1:length(model_sum)) {
  model_sum[[i]] <- model_sum[[i]][index[[i]]]
}

model_sum_dd <- model_sum[1:27]
model_sum_di <- model_sum[28:30]

par_names <- c("survival_ad", "survival_juv","beta", "sigma_s", 
                      "fecundity", "zeta", "sigma_f")

## MAPS Results
maps_res <- readRDS("BRCR_results.rds")

maps_sum <- map(maps_res, function(x) MCMCsummary(x, n.eff = T, digits = 4))
names(maps_sum) <- c("dd_r", "dd", "di_r", "cjs_r")

maps_chains <- map(maps_res, MCMCchains)
names(maps_chains) <- c("dd_r", "dd", "di_r", "cjs_r")

# Distributions of posterior means ----------------------------------------

tidy_sum <- function(x, par_names) {
  map(x, function(x) as.data.frame(x[par_names, c(1,3,5)])) %>%
    map(function(x) mutate(x, par = row.names(x))) %>%
    map(function(x) gather(x, stat, est, -par)) %>%
    do.call(rbind, .)
}

sum_data_dd <- map(model_sum_dd, tidy_sum, par_names = par_names)
#sum_data_di <- map(model_sum_di, tidy_sum, par_names = par_names)

pknames <- c("p1_k50", "p1_k100", "p1_k150",
             "p5_k50", "p5_k100", "p5_k150",
             "p10_k50","p10_k100", "p10_k150")

sum_data_dd2 <- foreach (i = 1:length(pknames), .combine = "rbind") %do% {
  index <- grep(pknames[i], names(sum_data_dd))
  data_length <- map_dbl(sum_data_dd[index], nrow)
  do.call(rbind, sum_data_dd[index]) %>%
    mutate(dd = factor(c(rep("Strong", data_length[1]),
                         rep("Weak", data_length[2]),
                         rep("Moderate", data_length[3])), 
                       levels = c("Weak", "Moderate", "Strong"), ordered = T)) %>%
    mutate(pop = strsplit(pknames[i], "_")[[1]][1]) %>%
    mutate(K = strsplit(pknames[i], "_")[[1]][2])
}

sum_data_dd2$pop <- factor(sum_data_dd2$pop, levels = c("p1", "p5", "p10"), ordered = T)
levels(sum_data_dd2$pop) <- c("Pop = 1", "Pop = 5", "Pop = 10")

sum_data_dd2$K <- factor(sum_data_dd2$K, levels = c("k50", "k100", "k150"), ordered = T)
levels(sum_data_dd2$K) <- c("K = 50", "K = 100", "K = 150")

sum_data_dd2$par <- factor(sum_data_dd2$par, 
                        levels = c("survival_ad", "survival_juv", "beta", "sigma_s", 
                                   "fecundity", "zeta", "sigma_f"), 
                        ordered = T)
levels(sum_data_dd2$par) <- c("S(ad)", "S(juv)", "DD(lgt)", "SD(lgt)", 
                           "Fec", "DD(log)", "SD(log)")


# Simulation Parameters ---------------------------------------------------

state_l <- readRDS("state_l_p1_k50.rds")
state_m <- readRDS("state_m_p1_k50.rds")
state_h <- readRDS("state_h_p1_k50.rds")

sim_params_l <- unlist(state_l[[1]]$parameters)
sim_params_l[1] <- inv.logit(sim_params_l["alpha1"] + sim_params_l["beta"])
sim_params_l[2] <- inv.logit(sim_params_l["alpha2"] + sim_params_l["beta"])
sim_params_l[3] <- exp(sim_params_l["theta"] + sim_params_l["zeta"])
sim_params_l <- sim_params_l[c("alpha1","alpha2", "beta", "sigma_s", "theta", "zeta", "sigma_f")]

sim_params_l <- as.data.frame(sim_params_l) %>%
  mutate(par = par_names) %>%
  select(est = sim_params_l,
         par = par) %>%
  mutate(dd = rep("L", 7))

sim_params_m <- unlist(state_m[[1]]$parameters)
sim_params_m[1] <- inv.logit(sim_params_m["alpha1"] + sim_params_m["beta"])
sim_params_m[2] <- inv.logit(sim_params_m["alpha2"] + sim_params_m["beta"])
sim_params_m[3] <- exp(sim_params_m["theta"] + sim_params_m["zeta"])
sim_params_m <- sim_params_m[c("alpha1","alpha2", "beta", "sigma_s", "theta", "zeta", "sigma_f")]

sim_params_m <- as.data.frame(sim_params_m) %>%
  mutate(par = par_names) %>%
  select(est = sim_params_m,
         par = par) %>%
  mutate(dd = rep("M", 7))

sim_params_h <- unlist(state_h[[1]]$parameters)
sim_params_h[1] <- inv.logit(sim_params_h["alpha1"] + sim_params_h["beta"])
sim_params_h[2] <- inv.logit(sim_params_h["alpha2"] + sim_params_h["beta"])
sim_params_h[3] <- exp(sim_params_h["theta"] + sim_params_h["zeta"])
sim_params_h <- sim_params_h[c("alpha1","alpha2", "beta", "sigma_s", "theta", "zeta", "sigma_f")]

sim_params_h <- as.data.frame(sim_params_h) %>%
  mutate(par = par_names) %>%
  select(est = sim_params_h,
         par = par) %>%
  mutate(dd = rep("H", 7))

sim_params <- rbind(sim_params_l, sim_params_m, sim_params_h)
sim_params$dd <- factor(sim_params$dd, levels = c("L", "M", "H"), ordered = T)

sim_params$par <- factor(sim_params$par, 
                        levels = c("survival_ad", "survival_juv", "beta", "sigma_s", 
                                   "fecundity", "zeta", "sigma_f"), 
                        ordered = T)
levels(sim_params$par) <- c("S(ad)", "S(juv)", "DD(lgt)", "SD(lgt)", 
                           "Fec", "DD(log)", "SD(log)")

sim_params_dd <- data.frame(est = c(-0.05, -0.5, -1), 
                            dd = c("Weak", "Moderate", "Strong"))


# base Plots -------------------------------------------------------------

theme_set(theme_bw())
jitter <- position_jitter(width = 0.1, height = 0.1)

# Figure 3
dd_data <- filter(sum_data_dd2, (par == "DD(lgt)" | par == "DD(log)") 
                  & pop == "Pop = 10" & K == "K = 150")
dd_data$par <- recode_factor(dd_data$par, "DD(lgt)" = "Survival", "DD(log)" = "Fecundity")

p_dd <- ggplot(dd_data, aes(dd, est, color = stat)) +
  geom_jitter(position = jitter, size = 1, shape = 1) +
  geom_point(data = sim_params_dd, size = 8, col = "black", shape = "+")

p_dd + facet_wrap(~ par, ncol = 3, nrow = 2, scales = "free_y") +
  labs(y = "RD-pop DD Estimates", x = "Simulated DD Strength") +
  scale_color_manual(name ="Posterior Statistics",
                     breaks = c("2.5%", "97.5%","mean"),
                     labels = c("2.5%", "97.5%","Mean"),
                     values = c("#0072B2","#009E73","#E69F00")) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom")

ggsave("fig3.jpeg", width = 6, height = 6, units = "in")

# Figure S1
sum_data_dd3 <- sum_data_dd2
sum_data_dd3$dd <- as.factor(sum_data_dd3$dd)
levels(sum_data_dd3$dd) <- c("Weak DD", "Moderate DD", "Strong DD")

sim_params2 <- sim_params
sim_params2$dd <- as.factor(sim_params2$dd)
levels(sim_params2$dd) <- c("Weak DD", "Moderate DD", "Strong DD")

p_sur <- ggplot(filter(sum_data_dd3, 
                       pop == "Pop = 10" & K == "K = 150" & stat == "mean" & 
                       par %in% c("S(ad)", "S(juv)")), 
                aes(par, est)) +
  geom_boxplot() +
  geom_point(data = filter(sim_params2, par %in% c("S(ad)", "S(juv)")), size = 2, col = "#D55E00")

p_sur + facet_wrap(~ dd, nrow = 3, scales = "fixed") +
  labs(y = "Posterior Means") + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave("figS1.tiff", width = 6, height = 8, units = "in")

# Figure not in manuscript
p_dd <- ggplot(filter(sum_data_dd3, 
                       pop == "Pop = 10" & K == "K = 150" & stat == "mean" & 
                         par %in% c("DD(lgt)", "DD(log)")), 
                aes(par, est)) +
  geom_boxplot() +
  geom_point(data = filter(sim_params2, par %in% c("DD(lgt)", "DD(log)")), size = 2, col = "#D55E00")

p_dd + facet_wrap(~ dd, nrow = 3, scales = "fixed") +
  labs(y = "Posterior Means") + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave("figS2.tiff", width = 6, height = 8, units = "in")

# Figure S2
p_sig <- ggplot(filter(sum_data_dd3, 
                      pop == "Pop = 10" & K == "K = 150" & stat == "mean" & 
                        par %in% c("SD(lgt)", "SD(log)")), 
               aes(par, est)) +
  geom_boxplot() +
  geom_point(data = filter(sim_params2, par %in% c("SD(lgt)", "SD(log)")), size = 2, col = "#D55E00")

p_sig + facet_wrap(~ dd, nrow = 3, scales = "fixed") +
  labs(y = "Posterior Means") + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

ggsave("figS3.tiff", width = 6, height = 8, units = "in")

# Figure 2
p_fec <- ggplot(filter(sum_data_dd3, 
                      pop == "Pop = 10" & K == "K = 150" & stat == "mean" & 
                        par == "Fec"), 
               aes(par, est)) +
  geom_boxplot() +
  geom_point(data = filter(sim_params2, par == "Fec"), size = 2, col = "#D55E00")

p_fec + facet_wrap(~ dd, ncol = 3, scales = "fixed") +
  labs(y = "Posterior Means", x = "Fecundity") + 
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size = 14))

ggsave("fig2.jpeg", width = 8, height = 6, units = "in")

# Figure S3
p_beta <- ggplot(filter(sum_data_dd2, par == "DD(lgt)"), aes(dd, est, color = stat)) +
  geom_jitter(position = jitter, size = 1, shape = 1) +
  geom_point(data = sim_params_dd, size = 8, col = "black", shape = "+")

p_beta + facet_wrap(~ pop + K, ncol = 3, nrow = 3, scales = "free_y") +
  labs(y = "CJS-pop DD Estimates", x = "Simulated DD Strength",
       title = "Density Dependence on Survival") +
  scale_color_manual(name ="Posterior Statistics",
                       breaks = c("2.5%", "97.5%","mean"),
                       labels = c("2.5%", "97.5%","Mean"),
                       values = c("#0072B2","#009E73","#E69F00")) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.position = "bottom") 

ggsave("figS2.tiff", width = 8, height = 8, units = "in")

# Figure S4
p_zeta <- ggplot(filter(sum_data_dd2, par == "DD(log)"), aes(dd, est, color = stat)) +
  geom_jitter(position = jitter, size = 1, shape = 1) +
  geom_point(data = sim_params_dd, size = 8, col = "black", shape = "+")

p_zeta <- ggplot(filter(sum_data_dd2, par == "DD(log)"), aes(dd, est, color = stat)) +
  geom_jitter(position = jitter, size = 1, shape = 1) +
  geom_point(data = sim_params_dd, size = 8, col = "black", shape = "+")

p_zeta + facet_wrap(~ pop + K, ncol = 3, nrow = 3, scales = "free_y") +
  labs(y = "DD Estimates", x = "DD Strength",
       title = "Density Dependence on Fecundity") +
  scale_color_manual(name ="Posterior Statistics",
                      breaks = c("2.5%", "97.5%","mean"),
                      labels = c("2.5%", "97.5%","Mean"),
                     values = c("#0072B2","#009E73","#E69F00")) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.position = "bottom")

ggsave("figS3.tiff", width = 8, height = 8, units = "in")


# Survival/Fecundity - Density Graphs -------------------------------------

dens <- seq(0,2,0.01)
beta <- maps_sum[[1]]["beta",1]
alpha1 <- logit(maps_sum[[1]]["survival_ad",1])
alpha2 <- logit(maps_sum[[1]]["survival_juv",1])
theta <- log(maps_sum[[1]]["fecundity",1])
zeta <- maps_sum[[1]]["zeta",1]

beta_chain <- maps_chains[[1]][,"beta"]
alpha1_chain <- logit(maps_chains[[1]][,"survival_ad"])
alpha2_chain <- logit(maps_chains[[1]][,"survival_juv"])
theta_chain <- log(maps_chains[[1]][,"fecundity"])
zeta_chain <- maps_chains[[1]][,"zeta"]

survival_ad <- inv.logit(alpha1 + beta*(dens-1))
survival_juv <- inv.logit(alpha2 + beta*(dens-1))
fecundity <- exp(theta + zeta*(dens-1))

survival_ad_h <- c()
survival_ad_l <- c()
survival_juv_h <- c()
survival_juv_l <- c()
fecundity_h <- c()
fecundity_l <- c()

for (i in 1:length(dens)) {
  survival_ad_h[i] <- inv.logit(alpha1_chain + beta_chain*(dens[i]-1)) %>%
    quantile(.,0.975)
  survival_ad_l[i] <- inv.logit(alpha1_chain + beta_chain*(dens[i]-1)) %>%
    quantile(.,0.025)
  
  survival_juv_h[i] <- inv.logit(alpha2_chain + beta_chain*(dens[i]-1)) %>%
    quantile(.,0.975)
  survival_juv_l[i] <- inv.logit(alpha2_chain + beta_chain*(dens[i]-1)) %>%
    quantile(.,0.025)
  
  fecundity_h[i] <- exp(theta_chain + zeta_chain*(dens[i]-1)) %>%
    quantile(.,0.975)
  fecundity_l[i] <- exp(theta_chain + zeta_chain*(dens[i]-1)) %>%
    quantile(.,0.025)
}

dens_data <- data.frame(s_ad = survival_ad,
                        s_ad_h = survival_ad_h,
                        s_ad_l = survival_ad_l,
                        s_juv = survival_juv,
                        s_juv_h = survival_juv_h,
                        s_juv_l = survival_juv_l,
                        fec  = fecundity,
                        fec_h = fecundity_h,
                        fec_l = fecundity_l,
                        dens = dens)

g_sad <- ggplot(dens_data, aes(x=dens)) +
  geom_ribbon(aes(ymin = s_ad_l, ymax = s_ad_h), fill  = "grey", alpha = 0.7) +
  geom_line(aes(y=s_ad), col = "black") +
  labs(x = "Relative Population Density", y = "Adult Survival") +
  scale_y_continuous(breaks=seq(0, 0.7, 0.1), limits = c(0,0.7))

g_sjuv <- ggplot(dens_data, aes(x=dens)) +
  geom_ribbon(aes(ymin = s_juv_l, ymax = s_juv_h), fill  = "grey", alpha = 0.7) +
  geom_line(aes(y=s_juv)) +
  labs(x = "Relative Population Density", y = "Juvenile Survival")

g_fec <- ggplot(dens_data, aes(x=dens)) +
  geom_ribbon(aes(ymin = fec_l, ymax = fec_h), fill  = "grey", alpha = 0.7) +
  geom_line(aes(y=fec), col = "black") +
  labs(x = "Relative Population Density", y = "Fecundity")+
  scale_y_continuous(breaks=seq(0, 4, 1), limits = c(0,4))

# Projection Plots --------------------------------------------------------

iterations <- 1000
nt <- 20
K <- 1000

## Projections with density dependence

pop_fun_dens <- function(data) {
  Nad <- c()
  Njuv <- c()
  Ntotal<- c()
  Nad[1] <- round(K/2)
  Njuv[1] <- K - Nad[1]
  Ntotal[1] <- Nad[1] + Njuv[1]
  
  index1 <- sample(1:length(data), 1)
  index2 <- sample(c(1,3,5), 1)
  alpha1 <- logit(data[[index1]]["survival_ad",index2])
  alpha2 <- logit(data[[index1]]["survival_juv",index2])
  beta <- data[[index1]]["beta",index2]
  theta <- log(data[[index1]]["fecundity",index2])
  zeta <- data[[index1]]["zeta",index2]
  sigma_s <- data[[index1]]["sigma_s",index2]
  sigma_f <- data[[index1]]["sigma_f",index2]
  
  for (t in 1:(nt-1)) {
    epsilon <- rnorm(1, 0, sigma_s)
    s.ad <- inv.logit(alpha1 + beta*(Ntotal[t]/K-1) + epsilon)
    s.juv <- inv.logit(alpha2 + beta*(Ntotal[t]/K-1) + epsilon)
    
    Nad1 <- rbinom(1, Nad[t], s.ad)
    Nad2 <- rbinom(1, Njuv[t], s.juv)
    Nad[t+1] <- Nad1 + Nad2
    
    omega <- rnorm(1, 0, sigma_f)
    f <- exp(theta + zeta*(Ntotal[t]/K-1) + omega)
    Njuv[t+1] <- rpois(1,f*Nad[t+1])
    Ntotal[t+1] <- Nad[t+1] + Njuv[t+1]
  }
  res <- list(Ntotal = Ntotal,
              Nad = Nad,
              Njuv = Njuv)
  
  return(res)
}

Ntotal_dd_h <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_dd_h[i,] <- pop_fun_dens(data = model_sum_dd$dd_h_p10_k150)$Ntotal
}
ema_dd_h <- apply(Ntotal_dd_h, 1, min, na.rm = T)

Ntotal_dd_m <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_dd_m[i,] <- pop_fun_dens(model_sum_dd$dd_m_p10_k150)$Ntotal
}
ema_dd_m <- apply(Ntotal_dd_m, 1, min, na.rm = T)

Ntotal_dd_l <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_dd_l[i,] <- pop_fun_dens(model_sum_dd$dd_l_p10_k150)$Ntotal
}
ema_dd_l <- apply(Ntotal_dd_l, 1, min, na.rm = T)


## Projections with density independence

pop_fun_avg <- function(data) {
  Nad <- c()
  Njuv <- c()
  Ntotal<- c()
  Nad[1] <- round(K/2)
  Njuv[1] <- K - Nad[1]
  Ntotal[1] <- Nad[1] + Njuv[1]
  
  index1 <- sample(1:length(data), 1)
  index2 <- sample(c(1,3,5), 1)
  alpha1 <- logit(data[[index1]]["survival_ad",index2])
  alpha2 <- logit(data[[index1]]["survival_juv",index2])
  theta <- log(data[[index1]]["fecundity",index2])
  sigma_s <- data[[index1]]["sigma_s",index2]
  sigma_f <- data[[index1]]["sigma_f",index2]
  
  for (t in 1:(nt-1)) {
    epsilon <- rnorm(1, 0, sigma_s)
    s.ad <- inv.logit(alpha1 + epsilon)
    s.juv <- inv.logit(alpha2 + epsilon)
    
    Nad1 <- rbinom(1, Nad[t], s.ad)
    Nad2 <- rbinom(1, Njuv[t], s.juv)
    Nad[t+1] <- Nad1 + Nad2
    
    omega <- rnorm(1, 0, sigma_f)
    f <- exp(theta + omega)
    Njuv[t+1] <- rpois(1, f*Nad[t+1])
    Ntotal[t+1] <- Nad[t+1] + Njuv[t+1]
  }
  res <- list(Ntotal = Ntotal,
              Nad = Nad,
              Njuv = Njuv)
  
  return(res)
}


Ntotal_di_h <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_di_h[i,] <- pop_fun_avg(model_sum_di$di_h_p10_k150)$Ntotal
}
ema_di_h <- apply(Ntotal_di_h, 1, min, na.rm = T)

Ntotal_di_m <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_di_m[i,] <- pop_fun_avg(model_sum_di$di_m_p10_k150)$Ntotal
}
ema_di_m <- apply(Ntotal_di_m, 1, min, na.rm = T)

Ntotal_di_l <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_di_l[i,] <- pop_fun_avg(model_sum_di$di_l_p10_k150)$Ntotal
}
ema_di_l <- apply(Ntotal_di_l, 1, min, na.rm = T)

## Projections with true simulation parameters

pop_fun_true <- function(data) {
  Nad <- c()
  Njuv <- c()
  Ntotal<- c()
  Nad[1] <- round(K/2)
  Njuv[1] <- K - Nad[1]
  Ntotal[1] <- Nad[1] + Njuv[1]
  
  alpha1 <- logit(filter(data, par == "survival_ad")$est)
  alpha2 <- logit(filter(data, par == "survival_juv")$est)
  beta <- filter(data, par == "beta")$est
  theta <- log(filter(data, par == "fecundity")$est)
  zeta <- filter(data, par == "zeta")$est
  sigma_s <- filter(data, par == "sigma_s")$est
  sigma_f <- filter(data, par == "sigma_f")$est
  
  for (t in 1:(nt-1)) {
    epsilon <- rnorm(1, 0, sigma_s)
    s.ad <- inv.logit(alpha1 + beta*(Ntotal[t]/K-1) + epsilon)
    s.juv <- inv.logit(alpha2 + beta*(Ntotal[t]/K-1) + epsilon)
    
    Nad1 <- rbinom(1, Nad[t], s.ad)
    Nad2 <- rbinom(1, Njuv[t], s.juv)
    Nad[t+1] <- Nad1 + Nad2
    
    omega <- rnorm(1, 0, sigma_f)
    f <- exp(theta + zeta*(Ntotal[t]/K-1) + omega)
    Njuv[t+1] <- rpois(1,f*Nad[t+1])
    Ntotal[t+1] <- Nad[t+1] + Njuv[t+1]
  }
  res <- list(Ntotal = Ntotal,
              Nad = Nad,
              Njuv = Njuv)
  
  return(res)
}

Ntotal_tr_h <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_tr_h[i,] <- pop_fun_true(sim_params_h)$Ntotal
}
ema_tr_h <- apply(Ntotal_tr_h, 1, min)

Ntotal_tr_m <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_tr_m[i,] <- pop_fun_true(sim_params_m)$Ntotal
}
ema_tr_m <- apply(Ntotal_tr_m, 1, min)

Ntotal_tr_l <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_tr_l[i,] <- pop_fun_true(sim_params_l)$Ntotal
}
ema_tr_l <- apply(Ntotal_tr_l, 1, min)

## Projections with density dependent maps model parameters

pop_fun_mapsdd <- function(data, index) {
  Nad <- c()
  Njuv <- c()
  Ntotal<- c()
  Nad[1] <- round(K/2)
  Njuv[1] <- K - Nad[1]
  Ntotal[1] <- Nad[1] + Njuv[1]
  
  index <- sample(1:nrow(data), 1)
  alpha1 <- logit(data[index, "survival_ad"])
  alpha2 <- logit(data[index,"survival_juv"])
  beta <- data[index, "beta"]
  theta <- log(data[index, "fecundity"])
  zeta <- data[index, "zeta"]
  sigma_s <- data[index, "sigma_s"]
  sigma_f <- data[index, "sigma_f"]
  
  for (t in 1:(nt-1)) {
    epsilon <- rnorm(1, 0, sigma_s)
    s.ad <- inv.logit(alpha1 + beta*(Ntotal[t]/K-1) + epsilon)
    s.juv <- inv.logit(alpha2 + beta*(Ntotal[t]/K-1) + epsilon)
    
    Nad1 <- rbinom(1, Nad[t], s.ad)
    Nad2 <- rbinom(1, Njuv[t], s.juv)
    Nad[t+1] <- Nad1 + Nad2
    
    omega <- rnorm(1, 0, sigma_f)
    f <- exp(theta + zeta*(Ntotal[t]/K-1) + omega)
    Njuv[t+1] <- rpois(1,f*Nad[t+1])
    Ntotal[t+1] <- Nad[t+1] + Njuv[t+1]
  }
  res <- list(Ntotal = Ntotal,
              Nad = Nad,
              Njuv = Njuv)
  
  return(res)
}

Ntotal_mapsdd <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_mapsdd[i,] <- pop_fun_mapsdd(maps_chains$dd_r, i)$Ntotal
}

med_dd <- apply(Ntotal_mapsdd, 2, median, na.rm = T)
h_dd <- apply(Ntotal_mapsdd, 2, function(x) quantile(x, 0.775, na.rm = T))
l_dd <- apply(Ntotal_mapsdd, 2, function(x) quantile(x, 0.225, na.rm = T))
ema_mapsdd <- apply(Ntotal_mapsdd, 1, min, na.rm = T)

## Projections with density independent maps model parameters

pop_fun_mapsdi <- function(data, index) {
  Nad <- c()
  Njuv <- c()
  Ntotal<- c()
  Nad[1] <- round(K/2)
  Njuv[1] <- K - Nad[1]
  Ntotal[1] <- Nad[1] + Njuv[1]
  
  index <- sample(1:nrow(data), 1)
  alpha1 <- logit(data[index, "survival_ad"])
  alpha2 <- logit(data[index,"survival_juv"])
  theta <- log(data[index, "fecundity"])
  sigma_s <- data[index, "sigma_s"]
  sigma_f <- data[index, "sigma_f"]
  
  for (t in 1:(nt-1)) {
    epsilon <- rnorm(1, 0, sigma_s)
    s.ad <- inv.logit(alpha1 + epsilon)
    s.juv <- inv.logit(alpha2 + epsilon)
    
    Nad1 <- rbinom(1, Nad[t], s.ad)
    Nad2 <- rbinom(1, Njuv[t], s.juv)
    Nad[t+1] <- Nad1 + Nad2
    
    omega <- rnorm(1, 0, sigma_f)
    f <- exp(theta + omega)
    Njuv[t+1] <- rpois(1,f*Nad[t+1])
    Ntotal[t+1] <- Nad[t+1] + Njuv[t+1]
  }
  res <- list(Ntotal = Ntotal,
              Nad = Nad,
              Njuv = Njuv)
  
  return(res)
}

Ntotal_mapsdi <- matrix(nrow = iterations, ncol = nt)
for (i in 1:iterations) {
  Ntotal_mapsdi[i,] <- pop_fun_mapsdi(maps_chains$di_r,i)$Ntotal
}

med_di <- apply(Ntotal_mapsdi, 2, median, na.rm = T)
h_di <- apply(Ntotal_mapsdi, 2, function(x) quantile(x, 0.775, na.rm = T))
l_di <- apply(Ntotal_mapsdi, 2, function(x) quantile(x, 0.225, na.rm = T))
ema_mapsdi <- apply(Ntotal_mapsdi, 1, min, na.rm = T)

# Fig S5
med_data <- data.frame(med = c(med_dd, med_di),
                       ll = c(l_dd, l_di),
                       hl= c(h_dd, h_di),
                       type = rep(c("DD", "DI"), each = 20),
                       time = c(1:20,1:20))

med_plot <- ggplot(med_data, aes(x = time)) +
  geom_ribbon(aes(ymin = ll, ymax = hl), fill  = "grey", alpha = 0.7) +
  geom_line(aes(y = med)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  labs(x = "Time", y = "Population Size") +
  scale_y_continuous(limits = c(0,NA))

med_plot2 <- med_plot + facet_wrap(~type, nrow = 2, scales = "free_y")

g_s <- plot_grid(g_sad, g_fec, ncol = 2, labels = c("a", "b"))
plot_grid(g_s, med_plot2, nrow = 2, rel_heights = c(0.6, 1, 1), 
          labels = c("", "c"))

ggsave("figS5.tiff", width = 8, height = 8, units = "in")

# Fig 4
ema_data1 <- data.frame(ema = c(ema_dd_l,ema_dd_m,ema_dd_h,
                                ema_di_l,ema_di_m,ema_di_h,
                                ema_tr_l,ema_tr_m,ema_tr_h),
                        dd_strength = rep(c("Weak", "Moderate", "Strong"), each = 1000, times = 3),
                        type = rep(c("DD", "DI", "TRUE"), each = 3000))

ema_data2 <- data.frame(ema = c(ema_mapsdd, ema_mapsdi),
                        dd_strength = "BRCR",
                        type = rep(c("DD", "DI"), each = 1000))
  

ema_data <- rbind(ema_data1,ema_data2)
ema_data$dd_strength <- factor(ema_data$dd_strength,
                               levels = c("BRCR","Weak", "Moderate", "Strong"))
levels(ema_data$dd_strength) <-  c("A) MAPS: BRCR",
                                   "B) Simulations: Weak DD", 
                                   "C) Simulations: Moderate DD", 
                                   "D) Simulations: Strong DD")

ema_plot <- ggplot(ema_data, aes(x = ema)) + 
  geom_density(aes(fill = type), alpha = 0.5)

ema_plot + facet_wrap(~dd_strength, nrow = 4, scales = "free_y") +
  labs(x = "Minimum Abundance", y = "Probability Density") +
  scale_fill_manual(name = "Projection Type",
                    values = c("#D55E00","#56B4E9","#009E73")) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 12))

ggsave("fig4.jpeg", width = 8, height = 8, units = "in")

# Tables for BRCR ---------------------------------------------------------

maps_tab <- maps_sum[[1]][c(-8,-9),c(1,3,5,6,7)]
maps_tab["gamma[1]", 1:3] <- inv.logit(maps_tab["gamma[1]", 1:3])
maps_tab["gamma[2]", 1:3] <- inv.logit(maps_tab["gamma[2]", 1:3])
maps_tab <- round(maps_tab,3)
maps_tab <- maps_tab[c("survival_ad", "survival_juv", "beta", "sigma_s", "fecundity", 
                       "zeta", "sigma_f", "gamma[1]", "gamma[2]", "delta", "pi[1]",
                       "pi[2]", "rho[1]", "rho[2]"),]

write.csv(maps_tab, "table1.csv")

tabS1 <- maps_sum$cjs_r
tabS1["gamma[1]", 1:3] <- inv.logit(tabS1["gamma[1]", 1:3])
tabS1["gamma[2]", 1:3] <- inv.logit(tabS1["gamma[2]", 1:3])
tabS1 <- round(tabS1,3)
tabS1 <- tabS1[c("survival_ad", "survival_juv", "sigma_s", 
                       "gamma[1]", "gamma[2]", "delta", "pi[1]",
                       "pi[2]", "rho[1]", "rho[2]"),]

write.csv(tabS1, "tableS1.csv")

