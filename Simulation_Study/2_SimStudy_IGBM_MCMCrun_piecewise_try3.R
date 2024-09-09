###########################################################################
#
# Setup
#
###########################################################################

#---- Session info: 
sessionInfo()

#---- Libraries:
lib <- c("ggplot2",
         "tidyverse",
         "rvest",
         "vroom",
         "rstan",
         "bayesplot",
         "gridExtra",
         "mgcv",
         "Bernadette",
         "readxl",
         "compiler",
         "extraDistr"
)
lapply(lib, require, character.only = TRUE)

.libPaths()

###########################################################################
#
# Paths
#
###########################################################################
model_id        <- "MCMCrun_IGBM_run_piecewise"
output_id       <- model_id
data_path       <- paste0("...//", "Datageneration_IGBM_v3_piecewise_try3", ".RData")
stan_file_path  <- paste0("...//", model_id, ".stan")
store_posterior_estimates <- paste0("...//", output_id, "_try3.RData")

#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")

# #---- Data engineering; daily time series of new cases and deaths for a given period:
source(paste0(path_to_source, "2_Data_Engineering", ".R"))

###########################################################################
#
# Load the mortality data:
#
###########################################################################
load(data_path)
rm(simul_v2, standata_sim)

###########################################################################
#
# Age distribution and contact matrix:
#
###########################################################################
age_mapping <- c(rep("0-39",   8),
                 rep("40-59",  4),
                 rep("60+",    4))

age_mapping_contacts <- age_mapping

#---- Age distribution and contact matrix:
age_distr <- age_distribution(country = "United Kingdom", year = 2020)
cm        <- contact_matrix(country = "GBR")

#---- Aggregate the age distribution table:
lookup_table <- data.frame(Initial = age_distr$AgeGrp,
                           Mapping = age_mapping
)
aggr_age  <- aggregate_age_distribution(age_distr, lookup_table)
age_bands <- length( unique(aggr_age$AgeGrp) )

aggr_age$fraction <- aggr_age$PopTotal/sum(aggr_age$PopTotal)
n_pop                     <- sum(aggr_age$PopTotal)
aggr_age_sim              <- aggr_age

#---- Lookup table:
lookup_table_cm <- data.frame(Initial = rownames(cm),
                              Mapping = age_mapping_contacts)

#---- Aggregate the contact matrix (average value):
# one small change to the function aggregate_contact_matrix from bernadette package that returned the names of the groups wrong
# This has been pushed to the Github repository on 25 June 2024.
aggregate_contact_matrix <- function(object, lookup_table, age_distr) 
{
  if (!identical(unique(lookup_table$Mapping), age_distr$AgeGrp)) {
    stop("The mapped age group labels do not correspond to the age group labels of the aggregated age distribution matrix.\n")
  }
  indiv_age_df <- base::data.frame(indiv_age = lookup_table$Initial)
  object_df <- cbind(object, indiv_age_df)
  long_dt <- stats::reshape(object_df, varying = list(base::names(object_df)[-ncol(object_df)]), 
                            direction = "long", sep = "_", timevar = "contact_age", 
                            times = base::names(object_df)[-ncol(object_df)], v.names = "contact")
  long_dt$id <- NULL
  long_dt$indiv_age <- as.factor(long_dt$indiv_age)
  long_dt$contact_age <- as.factor(long_dt$contact_age)
  names(lookup_table) <- c("indiv_age", "indiv_age_agr")
  long_dt <- base::merge(long_dt, lookup_table, by = "indiv_age", 
                         all.x = TRUE)
  names(lookup_table) <- c("contact_age", "contact_age_agr")
  long_dt <- base::merge(long_dt, lookup_table, by = "contact_age", 
                         all.x = TRUE)
  long_dt$diag_element <- ifelse(long_dt$indiv_age == long_dt$contact_age, 
                                 TRUE, FALSE)
  dt_aggregate <- stats::aggregate(long_dt$contact, by = list(long_dt$indiv_age_agr, 
                                                              long_dt$contact_age_agr, long_dt$diag_element), mean)
  names(dt_aggregate) <- c("indiv_age_agr", "contact_age_agr", 
                           "diag_element", "mean_cm")
  dt_aggregate <- stats::aggregate(dt_aggregate$mean_cm, by = list(dt_aggregate$indiv_age_agr, 
                                                                   dt_aggregate$contact_age_agr), sum)
  names(dt_aggregate) <- c("indiv_age_agr", "contact_age_agr", 
                           "mean_cm")
  cm_to_rescale <- as.matrix(stats::reshape(dt_aggregate, 
                                            idvar = "indiv_age_agr", timevar = "contact_age_agr", 
                                            direction = "wide"))
  cm_to_rescale <- cm_to_rescale[, 2:ncol(cm_to_rescale)]
  class(cm_to_rescale) <- "numeric"
  age_bands <- length(unique(age_distr$AgeGrp))
  ret <- matrix(0, age_bands, age_bands)
  for (i in seq(1, age_bands, 1)) {
    for (j in seq(1, age_bands, 1)) {
      ret[i, j] <- 0.5 * (1/age_distr$PopTotal[i]) * (cm_to_rescale[i, 
                                                                    j] * age_distr$PopTotal[i] + cm_to_rescale[j, 
                                                                                                               i] * age_distr$PopTotal[j])
    }
  }
  rownames(ret) <- colnames(ret) <- age_distr$AgeGrp
  ret <- as.data.frame(ret)
  return(ret)
}

aggr_cm <- aggregate_contact_matrix(cm, lookup_table_cm, age_distr = aggr_age_sim)

Bernadette::plot_contact_matrix(aggr_cm)

#---- A symmetric transformation of the contact matrix, to work with a 
# 1. Feed the lower triangular matrix from the Cholesky decomposition in Stan:
aggr_cm_sym <- diag(aggr_age_sim$PopTotal) %*% as.matrix(aggr_cm)
L_cm        <- t( chol(aggr_cm_sym) )

###########################################################################
#
# Group-specific Infection-to-death distribution:
#
###########################################################################
ifr_react2 <- data.frame(AgeGrp = c("0-39", "40-59", "60+"),
                         IFR     = c(0.024, 0.22, 7.0)/100)

###########################################################################
#
# Infection-to-death distribution:
#
###########################################################################
# Infection-to-onset distribution: shape = 1.352082207, rate = 0.265114158
# Onset-to-death distribution: shape = 4.938271605, rate = 0.2626
n_obs <- 100 # Time horizon (days)
ditd  <- Bernadette::itd_distribution(ts_length  = n_obs,
                                      gamma_mean = 24.19231,
                                      gamma_cv   = 0.3987261)

###########################################################################
#
# Initial proportion of exposed (at time t_0). Riou et al suggests Beta(1,999)
#
###########################################################################
pi_perc         <- 0.10   # Assume that 10% of each age group are Exposed, rest 85% are Susceptible
pi_prior_params <- lapply(pi_perc, function(x) estBetaParams(x, (0.05*x)) ) # Variance = 5% around the mean
pi_prior_params <- data.frame(do.call(rbind, pi_prior_params))

E_deathsByAge_day1 <- c(2, 5, 8)
p_sigmaCM <- 0.05
dE        <- 3    # Mean Incubation period of 3 days
dI        <- 4    # Mean Infection period of 4 days

ecr_changes  <- 7
n_changes    <- ceiling(n_obs / ecr_changes)
n_remainder  <- (n_obs - (n_changes-1)*ecr_changes)

###########################################################################
#
# Stan model statement
#
###########################################################################
compilation_time_start <- Sys.time()
m1 <- rstan::stan_model(stan_file_path)
compilation_time_end <- Sys.time()
duration_compilation <- compilation_time_end - compilation_time_start

###########################################################################
#
# 8. HMC initialisation:
#
###########################################################################

#---- Number of states:
model_states   <- c("S", "E", "E", "I", "I", "C")

#---- Modify data into a form suitable for Stan:
cov_data <- list(A           = age_bands,
                 n_obs       = nrow(y_sim_draw5_v2),
                 n_pop       = n_pop,
                 ecr_changes = 7,
                 n_changes   = n_changes,
                 n_remainder = n_remainder,
                 age_dist    = aggr_age$fraction,
                 pop_diag    = 1/(aggr_age$PopTotal),
                 n_difeq     = length(model_states),
                 y_deaths    = y_sim_draw5_v2[,-1],
                 L_cm        = L_cm,
                 ifr_age     = ifr_react2[,-1],
                 # ODE solver input:
                 t0          = 0,
                 ts          = 1:n_obs,
                 left_t      = 1:n_obs,
                 right_t     = 2:(n_obs + 1),
                 # Add a small number to avoid zero expected deaths in the young age group:
                 E_deathsByAge_day1 = unlist(y_sim_draw5_v2[1,-1]) + 0.001,
                 # Fixed parameters:
                 dE        = dE,    # Mean Incubation period of 3 days
                 dI        = dI,    # Mean Infection period of 4 days
                 # Prior distributions:
                 eta0_sd   = 1, #5,
                 p_sigmaCM = p_sigmaCM,
                 p_sigmaBM = rep(1, age_bands),
                 p_phiD    = 1/5,
                 p_pi      = pi_prior_params,
                 # Infection-to-death distribution:
                 I_D       = ditd,
                 # Debugging:
                 inference = 1,
                 doprint   = 0
)

#---- Specify parameters to monitor:
parameters <- c("rho",
                "phiD",
                "sigmaBM",
                "cm_sample",
                "beta_trajectory",
                "E_casesByAge",
                "E_deathsByAge",
                "E_cases",
                "E_deaths",
                "Susceptibles",
                "log_lik",
                "deviance")

#---- Set initial values:
fit_optim <- rstan::optimizing(m1, data = cov_data, seed = 1, hessian = TRUE)

sampler_init <- function(){
  list(eta0      = fit_optim$par[names(fit_optim$par) %in% "eta0"],
       x_init    = fit_optim$par[grepl("x_init",   names(fit_optim$par), fixed = TRUE)],
       x_noise   = fit_optim$par[grepl("x_noise[", names(fit_optim$par), fixed = TRUE)],
       L_raw     = fit_optim$par[grepl("L_raw[",   names(fit_optim$par), fixed = TRUE)], 
       rho       = fit_optim$par[names(fit_optim$par) %in% "rho"],
       sigmaBM   = fit_optim$par[grepl("sigmaBM",  names(fit_optim$par), fixed = TRUE)],
       phiD      = fit_optim$par[names(fit_optim$par) %in% "phiD"]
  )
}# End function

###########################################################################
#
# 9. Fit and sample from the posterior using NUTS:
#
###########################################################################

parallel::detectCores()

#---- MCMC options:
nChains       <- 6
nBurn         <- 500 
nPost         <- 500
nThin         <- 1
adapt_delta   <- 0.85
max_treedepth <- 16

nIter   <- (nBurn + nPost) * nThin
nBurnin <- nBurn * nThin

#---- Stan options:
rstan_options(auto_write = TRUE)
options(mc.cores = nChains)

#---- Test execution:
test <- rstan::sampling(m1,
                        data    = cov_data,
                        init    = sampler_init,
                        pars    = parameters,
                        chains  = 1,
                        warmup  = 50,
                        iter    = 100,
                        seed    = 1,
                        control = list(max_treedepth = max_treedepth,
                                       adapt_delta   = adapt_delta
                        ),
                        show_messages   = FALSE)

rm(test)

time_start_nuts1 <- Sys.time()
nuts_fit_1 <- rstan::sampling(m1,
                              data    = cov_data,
                              init    = sampler_init,
                              pars    = parameters,
                              chains  = nChains,
                              warmup  = nBurnin,
                              iter    = nIter,
                              seed    = 1,
                              control = list(max_treedepth = max_treedepth,
                                             adapt_delta   = adapt_delta
                              ),
                              show_messages   = FALSE)
time_end_nuts1 <- Sys.time()
duration_nuts1 <- time_end_nuts1 - time_start_nuts1

# warmup and sampling times for each chain:
print( get_elapsed_time(nuts_fit_1) )

time_run <- get_elapsed_time(nuts_fit_1)
apply(time_run, 1, sum)

###########################################################################
#
# 8. Store HMC output and MCMC summaries
#
###########################################################################
save(cov_data,
     sampler_init,
     nuts_fit_1,
     duration_nuts1,
     file = store_posterior_estimates)

load(file = store_posterior_estimates)