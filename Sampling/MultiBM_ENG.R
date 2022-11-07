###########################################################################
#
# Setup:
#
###########################################################################

sessionInfo()

#---- Libraries:

# ACTION: Manually execute the code in the file /R/1_Libraries.R

#--- Paths:
stan_file_path  <- paste0("...", "/src/", "seeiir_model_validation", ".stan")
store_model_out <- paste0("...", ".RData")
rt_path         <- paste0("...", "model_name", "_Eff_Rt" ,".RData")

#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")

#---- Country and analysis period:
country    <- "England"
start_date <- "2020-03-02"
end_date   <- "2020-09-27" 

#---- Age mappings:
age_mapping_deaths <- c(rep("0-39", 2), 
                        "40-59", 
                        rep("60+",2))

age_mapping_cases <- c(rep("0-39",  8),
                       rep("40-59", 4),
                       rep("60+",   7))

# The contact matrix contains 16 age groups:
age_mapping_contacts <- c(rep("0-39",  8),
                          rep("40-59", 4),
                          rep("60+",   4))

age_mapping_vaccinations <- c(rep("0-39",  7),
                              rep("40-59", 4),
                              rep("60+",   7))

# Mapping of the age distribution:
age_mapping_UN <- c(rep("0-39",  8),
                    rep("40-59", 4),
                    rep("60+",   7))

###########################################################################
#
# Data import:
#
###########################################################################
data_all <- data_engineering_observations(country = country,
                                          start_date,
                                          end_date,
                                          data_cases_AT = NULL, 
                                          age_mapping_deaths,
                                          age_mapping_cases)
mortality_data <- data_all$Deaths

###########################################################################
#
# Age distribution and contact matrix:
#
###########################################################################
cm    <- Bernadette::contact_matrix(country = "GBR")

#---- Aggregate the age distribution table:
aggr_age  <- data_all$Age_distribution_aggr
age_bands <- length( unique(aggr_age$AgeGrp) )

#---- Lookup table:
lookup_table_cm <- data.frame(Initial = rownames(cm),
                              Mapping = age_mapping_contacts)

#---- Aggregate the contact matrix (average value):
aggr_cm <- Bernadette::aggregate_contact_matrix(cm, 
                                                lookup_table_cm, 
                                                age.distr = aggr_age)

#---- Feed the lower triangular matrix from the Cholesky decomposition in Stan:
aggr_cm_sym <- diag(aggr_age$PopTotal) %*% as.matrix(aggr_cm)
L_cm        <- t( chol(aggr_cm_sym) )

###########################################################################
#
# 3. Group-specific Infection-to-death distribution:
#
###########################################################################
ifr_table  <- Bernadette::aggregate_ifr_react(x           = data_all$Age_distribution, 
                                              user_AgeGrp = age_mapping_UN,
                                              data_cases  = data_all$Cum_Cases)
ifr        <- ifr_table[[2]]

#---- If the calculated IFRs are not appropriate to the given country, modify them:
ifr$AgrIFR <- c(0.024, 0.22, 7.0)/100

###########################################################################
#
# Infection-to-death distribution:
#
###########################################################################
ditd <- Bernadette::discretize_itd(ts_length  = nrow(mortality_data),
                                   gamma_mean = 24.19231,
                                   gamma_cv   = 0.3987261)

###########################################################################
#
# Initial proportion of exposed (at time t_0). Riou et al suggests Beta(1,999)
#
###########################################################################
pi_perc         <- 0.1   # Assume that 10% of each age group are Exposed, rest 90% are Susceptible
pi_prior_params <- lapply(pi_perc, function(x) estBetaParams(x, (0.05*x)) ) # Variance = 5% around the mean
pi_prior_params <- data.frame(do.call(rbind, pi_prior_params))

###########################################################################
#
# 7. Stan model statement
#
###########################################################################
m1 <- rstan::stan_model(stan_file_path)

###########################################################################
#
# Stan data and nitialisation:
#
###########################################################################

#---- Number of states:
model_states   <- c("S", "E", "E", "I", "I", "C")

#---- Modify data into a form suitable for Stan:
cov_data <- list(A           = age_bands,
                 n_obs       = nrow(mortality_data),
                 n_pop       = sum(aggr_age$PopTotal),
                 age_dist    = aggr_age$PopTotal/sum(aggr_age$PopTotal),
                 pop_diag    = 1/(aggr_age$PopTotal),
                 n_difeq     = length(model_states),
                 y_deaths    = mortality_data[,-c(1:5)],
                 L_cm        = L_cm,
                 ifr_age     = ifr[,-1],
                 # ODE solver input:
                 t0          = 0,
                 ts          = mortality_data$Index,
                 left_t      = mortality_data$Index,
                 right_t     = mortality_data$Right,
                 # Add a small number to avoid zero expected deaths in the young age group:
                 E_deathsByAge_day1 = unlist(mortality_data[1,-c(1:5)]) + 0.001,
                 # Fixed parameters:
                 dE        = 3,    # Mean Incubation period of 3 days
                 dI        = 4,    # Mean Infection period of 4 days
                 # Prior distributions:
                 eta0_sd   = 5,
                 p_sigmaCM = 0.05,
                 p_sigmaBM = rep(4, age_bands),
                 p_phiD    = 1/5,
                 p_pi      = pi_prior_params,
                 # Infection-to-death distribution:
                 I_D       = ditd,
                 # Debugging:
                 inference = 1,
                 doprint   = 0
)

#---- Specify parameters to monitor:
parameters <- c("pi",
                "phiD",
                "sigmaBM",
                "cm_sample",
                "beta_N",
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
  list(eta0        = fit_optim$par[names(fit_optim$par) %in% "eta0"],
       eta_init    = fit_optim$par[grepl("eta_init",   names(fit_optim$par), fixed = TRUE)],
       eta_noise   = fit_optim$par[grepl("eta_noise[", names(fit_optim$par), fixed = TRUE)],
       L_raw       = fit_optim$par[grepl("L_raw[",     names(fit_optim$par), fixed = TRUE)], 
       pi          = fit_optim$par[names(fit_optim$par) %in% "pi"],
       sigmaBM     = fit_optim$par[grepl("sigmaBM",  names(fit_optim$par), fixed = TRUE)],
       phiD        = fit_optim$par[names(fit_optim$par) %in% "phiD"]
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
nBurn         <- 500 ## Number of warm-up samples per chain after thinning
nPost         <- 500 ## Number of post-warm-up samples per chain after thinning
nThin         <- 1
adapt_delta   <- 0.80
max_treedepth <- 14

nIter   <- (nBurn + nPost) * nThin
nBurnin <- nBurn * nThin

#---- Stan options:
rstan_options(auto_write = TRUE)
options(mc.cores = nChains)

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

#---- Warmup and sampling times for each chain:
print( get_elapsed_time(nuts_fit_1) )

time_run <- get_elapsed_time(nuts_fit_1)
apply(time_run, 1, sum)

###########################################################################
#
# Store MCMC output
#
###########################################################################
save(cov_data,
     sampler_init,
     nuts_fit_1,
     duration_nuts1,
     file = store_model_out)

# load(file = store_model_out)
