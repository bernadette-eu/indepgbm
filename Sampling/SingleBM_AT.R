###########################################################################
#
# Setup:
#
###########################################################################

sessionInfo()

#---- Libraries:

# ACTION: Manually execute the code in the file /R/1_Libraries.R

#--- Paths:
store_model_out <- paste0("...", ".RData")
rt_path         <- paste0("...", "model_name", "_Eff_Rt" ,".RData")

#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")

#---- Country and analysis period:
country    <- "Austria"
start_date <- "2020-08-31"
end_date   <- "2021-03-28"

#---- Age mappings:
age_mapping_deaths <- c(rep("<65",7), rep("65+",3))

age_mapping_UN <- c(rep("<65", 13),
                    rep("65+",  3))

age_mapping_cases <- c(rep("<65", 13),
                       rep("65+", 8))

###########################################################################
#
# Data import:
#
###########################################################################

#---- Load the data | Follow the instructions in the file /R/0_Austria_Data_download.R

#---- Data processing:
data_all <- data_engineering_observations(country = country,
                                          start_date,
                                          end_date,
                                          data_cases_AT = Austria_5yo_bands_data, 
                                          age_mapping_deaths,
                                          age_mapping_cases)

mortality_data        <- data_all$Deaths
cusum_infections_data <- data_all$Cusum_Cases

###########################################################################
#
# Age distribution and contact matrix:
#
###########################################################################
age_distr <- Bernadette::age_distribution(country = country, year = 2020)
cm        <- Bernadette::contact_matrix(country = "AUT")

#---- Lookup table:
lookup_table <- data.frame(Initial = age_distr$AgeGrp,
                           Mapping = age_mapping_UN)

#---- Aggregate the age distribution table:
aggr_age  <- Bernadette::aggregate_age_distribution(age_distr, lookup_table)
age_bands <- length( unique(aggr_age$AgeGrp) )

#---- Aggregate the contact matrix (average value):
aggr_cm <- Bernadette::aggregate_contact_matrix(cm, lookup_table, age.distr = aggr_age)

#---- A symmetric transformation of the contact matrix, to work with a 
# 1. Feed the lower triangular matrix from the Cholesky decomposition in Stan:
aggr_cm_sym <- diag(aggr_age$PopTotal) %*% as.matrix(aggr_cm)
L_cm        <- t( chol(aggr_cm_sym) )

###########################################################################
#
# 3. Group-specific Infection-to-death distribution:
#
###########################################################################
ifr_table  <- Bernadette::aggregate_ifr_react(x           = age_distr, 
                                              user_AgeGrp = age_mapping_UN,
                                              data_cases  = cusum_infections_data)
ifr        <- ifr_table[[2]]
ifr$AgrIFR <- c(0.07, 6.8)/100

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
# Stan data and nitialisation:
#
###########################################################################

#---- Number of states:
model_states   <- c("S", "E", "E", "I", "I", "C")

#---- Modify data into a form suitable for Stan:
cov_data <- list(A           = age_bands,
                 n_obs       = nrow(mortality_data),
                 n_pop       = sum(aggr_age$PopTotal),
                 n_weeks     = max(mortality_data$Week_ID),
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
                 eta0_sd   = 1,
                 p_sigmaCM = 0.05,
                 p_sigmaBM = c(0, 1),
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
fit_optim <- rstan::optimizing(Bernadette:::stanmodels$seeiir_sbm, 
                               data    = cov_data,
                               seed    = 2, 
                               hessian = TRUE)

sampler_init <- function(){
  list(eta0        = fit_optim$par[names(fit_optim$par) %in% "eta0"],
       eta_init    = fit_optim$par[names(fit_optim$par) %in% "eta_init"],
       eta_noise   = fit_optim$par[grepl("eta_noise",  names(fit_optim$par), fixed = TRUE)],
       L_raw       = fit_optim$par[grepl("L_raw[",     names(fit_optim$par), fixed = TRUE)], 
       pi          = fit_optim$par[names(fit_optim$par) %in% "pi"],
       sigmaBM     = fit_optim$par[names(fit_optim$par) %in% "sigmaBM"],
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
nChains       <- 3
nBurn         <- 1500 ## Number of warm-up samples per chain after thinning
nPost         <- 1500 ## Number of post-warm-up samples per chain after thinning
nThin         <- 1
adapt_delta   <- 0.80
max_treedepth <- 14

nIter   <- (nBurn + nPost) * nThin
nBurnin <- nBurn * nThin

#---- Stan options:
rstan_options(auto_write = TRUE)
options(mc.cores = nChains)

time_start_nuts1 <- Sys.time()
nuts_fit_1 <- rstan::sampling(Bernadette:::stanmodels$seeiir_sbm,
                              data     = cov_data,
                              init     = sampler_init,
                              pars     = parameters,
                              chains   = nChains,
                              warmup   = nBurnin,
                              iter     = nIter,
                              seed     = 1,
                              control  = list(max_treedepth = max_treedepth,
                                              adapt_delta   = adapt_delta
                              ),
                              show_messages = FALSE)
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
