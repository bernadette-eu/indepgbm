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
         "cmdstanr",
         "bayesplot",
         "gridExtra",
         "mgcv",
         "Bernadette",
         "readxl",
         "compiler",
         "extraDistr"
)
lapply(lib, require, character.only = TRUE)

#---- Session info:
.libPaths()

###########################################################################
#
# Paths
#
###########################################################################
model_id        <- "Datageneration_IGBM_v3_piecewise"
output_id       <- model_id
stan_file_path  <- paste0("...//", model_id, ".stan")
store_simdata   <- paste0("...//", output_id, "_try3.RData")

#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")

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
n_pop             <- sum(aggr_age$PopTotal)
aggr_age_sim      <- aggr_age

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
                         IFR    = c(0.024, 0.22, 7.0)/100)

###########################################################################
#
# Infection-to-death distribution:
#
###########################################################################
# Infection-to-onset distribution: shape = 1.352082207, rate = 0.265114158
# Onset-to-death distribution: shape = 4.938271605, rate = 0.2626
#
n_obs <- 100 # Time horizon (days)
ditd  <- Bernadette::itd_distribution(ts_length  = n_obs,
                                      gamma_mean = 24.19231,
                                      gamma_cv   = 0.3987261)

###########################################################################
#
# Initial proportion of exposed (at time t_0). Riou et al suggests Beta(1,999)
#
###########################################################################
pi_perc            <- 0.10   # Assume that 10% of each age group are Exposed, rest 85% are Susceptible
E_deathsByAge_day1 <- c(2, 5, 8)

eta0     <- c(-0.06)
x_init   <- c(-0.34, -0.20, -1.02) 
phiD     <- 0.14
sigmaBM  <- c(0.70, 0.30, 0.60) 
p_sigmaCM<- 0.05
dE       <- 3    # Mean Incubation period of 3 days
dI       <- 4    # Mean Infection period of 4 days

ecr_changes  <- 7
n_changes    <- ceiling(n_obs / ecr_changes)
n_remainder  <- (n_obs - (n_changes-1)*ecr_changes)

###########################################################################
#
# Stan model statement (for simulation)
#
###########################################################################
compilation_time_start <- Sys.time()
m1 <- rstan::stan_model(stan_file_path)
compilation_time_end <- Sys.time()
duration_compilation <- compilation_time_end - compilation_time_start

#---- Number of states:
model_states   <- c("S", "E", "E", "I", "I", "C")

#---- Modify data into a form suitable for Stan:
standata_sim <- list(A           = dim(L_cm)[1],
                     n_obs       = n_obs,
                     n_pop       = n_pop,
                     ecr_changes = 7,
                     n_changes   = n_changes,
                     n_remainder = n_remainder,
                     age_dist    = aggr_age_sim$PopTotal/sum(aggr_age_sim$PopTotal),
                     pop_diag    = 1/(aggr_age_sim$PopTotal),
                     n_difeq     = length(model_states),
                     L_cm        = L_cm,
                     ifr_age     = ifr_react2[,-1],
                     # ODE solver input:
                     t0          = 0,
                     ts          = 1:n_obs,
                     left_t      = 1:n_obs,
                     right_t     = 2:(n_obs + 1),
                     # Infection-to-death distribution:
                     I_D       = ditd,
                     # Add a small number to avoid zero expected deaths in the young age group:
                     E_deathsByAge_day1 = E_deathsByAge_day1 + 0.001,
                     # Fixed parameters:
                     dE        = dE,    # Mean Incubation period of 3 days
                     dI        = dI,    # Mean Infection period of 4 days
                     # Prior distributions:
                     eta0      = eta0,
                     x_init    = x_init,
                     pi        = pi_perc,
                     sigmaBM   = sigmaBM,
                     phiD      = phiD,
                     p_sigmaCM = p_sigmaCM
                     )

#---- Specify parameters to monitor:
parameters_sim <- c("cm_sample",
                    "beta_trajectory",
                    "y_sim",
                    "Susceptibles")

max_treedepth <- 16
adapt_delta   <- 0.8

time_start_simul_v2 <- Sys.time()
simul_v2 <- rstan::sampling(m1,
                            data    = standata_sim,
                            pars    = parameters_sim,
                            chains  = 1,
                            warmup  = 500,
                            iter    = 510,
                            seed    = 1,
                            control = list(max_treedepth = max_treedepth,
                                            adapt_delta  = adapt_delta
                            ),
                            show_messages   = FALSE)
time_end_simul_v2 <- Sys.time()
duration_simul_v2 <- time_end_simul_v2 - time_start_simul_v2

# warmup and sampling times for each chain:
print( get_elapsed_time(simul_v2) )

simul_v2_time <- get_elapsed_time(simul_v2)
apply(simul_v2_time, 1, sum)

###########################################################################
#
# Store HMC output and MCMC summaries (simulation)
#
###########################################################################
posts_simul_v2 <- rstan::extract(simul_v2)

pb <- txtProgressBar(min = 0, max = dim(posts_simul_v2$beta_trajectory)[1], initial = 0, style = 3) 

for ( t in dim(posts_simul_v2$beta_trajectory)[1]:1){

  draw               <- t 
  data_mat_cols      <- c("Date", "Group", "Observed")
  data_mat           <- data.frame(matrix(ncol = length(data_mat_cols), nrow = 0))
  colnames(data_mat) <- data_mat_cols
  
  for (i in 1:standata_sim$A){
    
    fit_rt  <- posts_simul_v2$y_sim[draw,,i] #-c(1:rm_dates)
    
    dt_rt_group  <- data.frame(Date     = 1:standata_sim$n_obs,
                               Group    = rep(aggr_age_sim$AgeGrp[i], standata_sim$n_obs)
    )
    
    dt_rt_group$Observed <- fit_rt
    data_mat <- rbind(data_mat, dt_rt_group)
    
  }# End for
  
  
  eff_cont_rate_mat_cols      <- c("Date", "Group", "Observed")
  eff_cont_rate_mat           <- data.frame(matrix(ncol = length(eff_cont_rate_mat_cols), nrow = 0))
  colnames(eff_cont_rate_mat) <- eff_cont_rate_mat_cols
  
  for (i in 1:standata_sim$A){
    
    fit_effcr  <- posts_simul_v2$beta_trajectory[draw,,i] #-c(1:rm_dates)
    
    dt_effcr_group  <- data.frame(Date  = 1:standata_sim$n_obs,
                                  Group = rep(aggr_age_sim$AgeGrp[i], standata_sim$n_obs)
    )
    
    dt_effcr_group$Observed <- fit_effcr
    eff_cont_rate_mat       <- rbind(eff_cont_rate_mat, dt_effcr_group)
    
  }# End for

###########################################################################
#
# Time series - daily mortality data per age group (simulation)
#
###########################################################################
print(
    ggplot(data_mat,
         aes(Date       = Date,
             New_Deaths = Observed)) +
  geom_point(aes(x   = Date,
                 y   = Observed,
                 fill= "Reported deaths")) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right") +
  scale_fill_manual(values = c('Reported deaths' = 'black'), guide = "none") +
  labs(x = "Day",
       y = "Number of new deaths",
       title = paste0("Draw:", t)) +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(angle = 0, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())
)

###########################################################################
#
# Effective contact rate trajectories (simulation)
#
###########################################################################
print(
ggplot(eff_cont_rate_mat,
       aes(Date       = Date,
           New_Deaths = Observed)) +
  geom_point(aes(x   = Date,
                 y   = Observed,
                 fill= "Observed")) +
  geom_line(aes(x   = Date,
                 y   = Observed)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "right") +
  scale_fill_manual(values = c('Observed' = 'black'), guide = "none") +
  labs(x = "Day",
       y = "Effective contact rate",
       title = paste0("Draw:", t)) +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(angle = 0, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())
)
  
  setTxtProgressBar(pb, dim(posts_simul_v2$beta_trajectory)[1] - t + 1)
  
}# End for

#---- Outcome:
selected_draw <- 5

y_sim_draw5_v2 <- data.frame(Date = 1:standata_sim$n_obs)

for (i in 1:standata_sim$A) y_sim_draw5_v2 <- cbind(y_sim_draw5_v2, posts_simul_v2$y_sim[selected_draw,,i])
colnames(y_sim_draw5_v2) <- c("Date", ifr_react2$AgeGrp)

#---- Store the data:
save(standata_sim,
     simul_v2,
     duration_simul_v2,
     y_sim_draw5_v2,
     file = store_simdata)

#---- Remove redundant objects:
rm(m1, pb, simul_v2, posts_simul_v2, simul_v2_time, eff_cont_rate_mat, data_mat)


