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
         "readr",
         "dplyr",
         "tidyr",
         "rvest",
         "vroom",
         "rstan",
         "bayesplot",
         "gridExtra",
         "mgcv",
         "Bernadette",
         "readxl",
         "compiler",
         "extraDistr",
         "patchwork"
)
lapply(lib, require, character.only = TRUE)

.libPaths()

###########################################################################
#
# Paths
#
###########################################################################
model_id        <- "MCMCrun_IGBM_run_piecewise_try3"
output_id       <- model_id
store_posterior_estimates <- paste0("...//", output_id, ".RData")
store_simdata   <- paste0("...//", "Datageneration_IGBM_v3_piecewise_try3", ".RData")
rt_path         <- paste0("...//", output_id, "_Rt.RData")
  
#---- Set system locale to English:
Sys.setlocale("LC_ALL", "English")

#---- Data engineering; daily time series of new cases and deaths for a given period:
source(paste0(path_to_source, "2_Data_Engineering", ".R"))

###########################################################################
#
# Load libraries:
#
###########################################################################
load(file = store_simdata)
load(file = store_posterior_estimates)

posts_simul_v2 <- rstan::extract(simul_v2) 
posts_mcmc     <- rstan::extract(nuts_fit_1) # Posterior draws

#########################################################
#
# Age distribution and contact matrix:
#
#########################################################
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

#########################################################
# 
# Age-specific effective contact rate
# 
#########################################################
selected_draw <- 5

eff_cont_rate_mat_cols      <- c("Date", "Group", "Observed", "median", "low", "low25", "high75", "high")
eff_cont_rate_mat           <- data.frame(matrix(ncol = length(eff_cont_rate_mat_cols), nrow = 0))
colnames(eff_cont_rate_mat) <- eff_cont_rate_mat_cols

for (i in 1:standata_sim$A){
  
  observed_effcr  <- posts_simul_v2$beta_trajectory[selected_draw,,i]
  fit_rt          <- posts_mcmc$beta_trajectory[,,i] 
  
  dt_effcr_group  <- data.frame(Date     = 1:standata_sim$n_obs,
                                Group    = rep(aggr_age_sim$AgeGrp[i], standata_sim$n_obs),
                                Observed = observed_effcr
  )
  
  # Add quantiles from the model outputs:
  dt_effcr_group$median <- apply(fit_rt, 2, median)
  dt_effcr_group$low    <- apply(fit_rt, 2, quantile, probs = c(0.025))
  dt_effcr_group$low25  <- apply(fit_rt, 2, quantile, probs = c(0.25))
  dt_effcr_group$high75 <- apply(fit_rt, 2, quantile, probs = c(0.75))
  dt_effcr_group$high   <- apply(fit_rt, 2, quantile, probs = c(0.975))
  
  eff_cont_rate_mat  <- rbind(eff_cont_rate_mat, dt_effcr_group)
  
}# End for

eff_cont_rate_mat$Group_f <- factor(eff_cont_rate_mat$Group, 
                                    levels = aggr_age$AgeGrp)

colnames(eff_cont_rate_mat)[3] <- "True"

#---- Plot of posterior estimates:
posterior_eff_cont_rate <- 
ggplot(eff_cont_rate_mat,
       aes(Date = Date,
           true = True)) +
  facet_wrap(. ~ Group_f, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top") +
  geom_point(aes(x   = Date,
                 y   = True,
                 fill = "True value"
                 ),
             shape = 4
             ) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"
                ),
            size = 1.0
            ) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"
                  ),
              alpha = 0.5) +
  labs(x = "Day",
       y = "Contact rate",
       fill = "") +
  scale_fill_manual(values = c('True value' = "black")) +
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black", 
                                 "50% CrI" = "gray40"    
                      )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(size = 16, angle = 0, hjust = 1),
        axis.text.y      = element_text(size = 16),
        axis.title.x     = element_text(size = 18, face = "bold"),
        axis.title.y     = element_text(size = 18, face = "bold"),
        legend.text      = element_text(size = 16),
        strip.text.x     = element_text(size = 14),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())


posterior_eff_cont_rate_95CrI <- 
  ggplot(eff_cont_rate_mat,
         aes(Date = Date,
             true = True)) +
  facet_wrap(. ~ Group_f, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top") +
  geom_point(aes(x   = Date,
                 y   = True,
                 fill = "True value"
  ),
  shape = 4
  ) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"
  ),
  size = 1.0
  ) +
  geom_ribbon(aes(x    = Date,
                  ymin = low,
                  ymax = high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Day",
       y = "Contact rate",
       fill = "") +
  scale_fill_manual(values = c('True value' = "black")) +
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black", 
                                 "95% CrI" = "gray40"    
  
                      )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(size = 16, angle = 0, hjust = 1),
        axis.text.y      = element_text(size = 16),
        axis.title.x     = element_text(size = 18, face = "bold"),
        axis.title.y     = element_text(size = 18, face = "bold"),
        legend.text      = element_text(size = 16),
        strip.text.x     = element_text(size = 14),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())

#########################################################
# 
# Expected deaths - goodness of fit
# 
#########################################################

deaths_random_draws <- function(dispersion_type = "time-independent",
                                cov_data,
                                model_out
                                ){
  
  set.seed(1)
  
  #---- Checks:
  if(class(model_out)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(model_out)
  E_deathsByAge   <- posterior_draws[["E_deathsByAge"]]
  phi             <- posterior_draws[["phiD"]]
  
  mcmc_length     <- dim(E_deathsByAge)[1]
  ts_length       <- dim(E_deathsByAge)[2]
  age_grps        <- dim(E_deathsByAge)[3]
  dates           <- 1:nrow(cov_data$y_deaths)
  
  death_draws           <- array(NA, c(mcmc_length, ts_length, age_grps))
  data_deaths_cols      <- c("Date",
                             "Group",
                             'median',
                             "rng_low",
                             "rng_low25",
                             "rng_high75",
                             "rng_high")
  age_rng_draws           <- data.frame(matrix(ncol = length(data_deaths_cols),
                                               nrow = 0))
  colnames(age_rng_draws) <- data_deaths_cols
  aggregated_rng_draws    <- age_rng_draws
  
  message(" > Estimation for age-spefic deaths")
  
  for (k in 1:age_grps) {
    message(paste0(" > Estimation in group ", k))
    
    for (i in 1:mcmc_length) {
      for (j in 1:ts_length) {
        
        if(dispersion_type == "time-independent"){
          death_draws[i,j,k] <- rnbinom(1,
                                        mu   = E_deathsByAge[i,j,k],
                                        size = phi[i])
          
        } else {
          death_draws[i,j,k] <- rnbinom(1,
                                        mu   = E_deathsByAge[i,j,k],
                                        size = (E_deathsByAge[i,j,k] / phi[i]))
          
        }# End if
      }# End for
    }# End for
  }# End for
  
  for (k in 1:age_grps) {
    rng_draws_age_grp  <- data.frame(Date  = dates,
                                     Group = rep( colnames(cov_data$y_deaths)[k],
                                                  length(dates)
                                     ))
    rng_draws_age_grp$median     <- apply(death_draws[,,k], 2, quantile, probs = c(0.50))
    rng_draws_age_grp$rng_low    <- apply(death_draws[,,k], 2, quantile, probs = c(0.025))
    rng_draws_age_grp$rng_low25  <- apply(death_draws[,,k], 2, quantile, probs = c(0.25))
    rng_draws_age_grp$rng_high75 <- apply(death_draws[,,k], 2, quantile, probs = c(0.75))
    rng_draws_age_grp$rng_high   <- apply(death_draws[,,k], 2, quantile, probs = c(0.975))
    
    age_rng_draws <- rbind(age_rng_draws, rng_draws_age_grp)
  }# End for
  
  ###################################
  age_rng_draws <- age_rng_draws %>% distinct(Date, Group, .keep_all = TRUE)
  
  return(age_rng_draws)
  
}# End function
deaths_random_draws <- compiler::cmpfun(deaths_random_draws)

age_specific_deaths_simdata           <- posts_simul_v2$y_sim[selected_draw,,]
age_specific_deaths_simdata           <- data.frame(age_specific_deaths_simdata)
colnames(age_specific_deaths_simdata) <- colnames(cov_data$y_deaths)
age_specific_deaths_simdata$Date      <- 1:nrow(age_specific_deaths_simdata)
age_specific_deaths_simdata           <- age_specific_deaths_simdata[c("Date", colnames(cov_data$y_deaths))]
age_specific_deaths_simdata_long      <- tidyr::pivot_longer(age_specific_deaths_simdata, 
                                                             cols      = -c("Date"), 
                                                             values_to = "Deaths", 
                                                             names_to  = "Group") 

deaths_rng    <- deaths_random_draws(dispersion_type     = "time-dependent",
                                     cov_data,
                                     model_out           = nuts_fit_1
                                     )

age_specific_deaths_data <- age_specific_deaths_simdata_long %>% 
                            left_join(deaths_rng, by = c("Date", "Group"))

posterior_deaths <-
ggplot(age_specific_deaths_data,
       aes(Date       = Date,
           Deaths = Deaths)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top") +
  geom_point(aes(x   = Date, 
                 y   = Deaths,
                 fill= "True value"),
             shape = 4
             ) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"),
            size = 1.0) +
  geom_ribbon(aes(x    = Date,
                  ymin = rng_low,
                  ymax = rng_high,
                  fill = "95% CrI"),
              alpha = 0.5) +
  labs(x = "Day",
       y = "New daily deaths",
       fill = "") +
  scale_fill_manual(values = c('True value' = "black")) +
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black",   
                                 "95% CrI" = "gray60"
                      )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(size = 16, angle = 0, hjust = 1),
        axis.text.y      = element_text(size = 16),
        axis.title.x     = element_text(size = 18, face = "bold"),
        axis.title.y     = element_text(size = 18, face = "bold"),
        legend.text      = element_text(size = 16),
        strip.text.x     = element_text(size = 14),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())
  
#########################################################
# 
# Contact matrix
# 
#########################################################

niters       <- nrow(posts_mcmc[["cm_sample"]])
prior_cm     <- array(numeric(),c(niters, cov_data$A, cov_data$A)) 
plot_cm_list <- list()
plot_cm_indx <- 1

for (i in 1:cov_data$A){
  for (j in 1:cov_data$A){
    
    # Create the appropriate dataset in a format acceptable by ggplot:
    dt_cm_post  <- data.frame(Posterior = posts_mcmc[["cm_sample"]][,i, j])
    dt_cm_post2 <- reshape2::melt(dt_cm_post)
    
    p <- ggplot(dt_cm_post2, 
                aes(x = value)) +
      geom_density(alpha = 0.8) +
      scale_x_continuous(
        limits = c(min(dt_cm_post2$value)*0.95, max(dt_cm_post2$value)*1.05)
      ) +
      geom_vline(xintercept = aggr_cm[i,j], linetype = "dashed") +
      labs(x = "",
           y = "Density",
           title = paste0("C[",i,",",j,"]")) +
      theme(strip.placement  = "outside",
            strip.background = element_rect(fill = NA, colour="grey50"),
            panel.spacing    = unit(0,"cm"),
            legend.position  = "bottom",
            legend.title     = element_blank() )
    
    plot_cm_list[[plot_cm_indx]] <- p
    
    # Increment the index of stored graph:
    plot_cm_indx <- plot_cm_indx + 1
    
  }# End for
}# End for

# Remove redundant objects:
rm(i, dt_cm_post, dt_cm_post2)

# Plot the distributions for rows 1-2:
cm_posterior_store_graph <- do.call("grid.arrange", c(plot_cm_list, nrow = cov_data$A))

#########################################################
# 
# Transmission rate
# 
#########################################################
transm_rate_generation <- function(cov_data,
                                   model_out,
                                   posts_simul,
                                   aggr_cm,
                                   selected_draw
){
  
  `%nin%` <- Negate(`%in%`)
  
  #---- Checks:
  if(class(model_out)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(model_out)
  
  #---- Checks:
  if(ncol(posterior_draws$cm_sample) != cov_data$A) stop( paste0("The number of rows in the age distribution table must be equal to ", cov_data$A) )
  
  age_grps <- cov_data$A
  
  transm_rate_mat_cols      <- c("Date", "Group", "True", "median", "low", "low25", "high75", "high")
  transm_rate_mat           <- data.frame(matrix(ncol = length(transm_rate_mat_cols), nrow = 0))
  colnames(transm_rate_mat) <- transm_rate_mat_cols
  
  dates        <- 1:nrow(cov_data$y_deaths) 
  beta_draws   <- posterior_draws$beta_trajectory
  chain_length <- nrow(beta_draws)
  ts_length    <- length(dates)
  
  for (i in 1:age_grps){
    message(paste0(" > Estimation in group ", i))
    
    transm_rate_group_tmp <- data.frame(Date  = dates,
                                        Group = rep( colnames(cov_data$y_deaths)[i],
                                                     length(dates)
                                        ))
    trans_rate_temp       <- matrix(0L, nrow = chain_length, ncol = ts_length)
    
    #---- Observed:
    observed_effcr             <- posts_simul$beta_trajectory[selected_draw,,i]
    transm_rate_group_tmp$True <- observed_effcr * aggr_cm[age_grps,age_grps]
    
    #---- Posterior summaries:
    for (j in 1:ts_length) trans_rate_temp[,j] <- beta_draws[,j,i] * posterior_draws$cm_sample[,i,i]
    
    transm_rate_group_tmp$median  <- apply(trans_rate_temp, 2, median)
    transm_rate_group_tmp$low     <- apply(trans_rate_temp, 2, quantile, probs = c(0.025)) # c(0.025)
    transm_rate_group_tmp$low25   <- apply(trans_rate_temp, 2, quantile, probs = c(0.25))  # c(0.025)
    transm_rate_group_tmp$high75  <- apply(trans_rate_temp, 2, quantile, probs = c(0.75))  # c(0.975)
    transm_rate_group_tmp$high    <- apply(trans_rate_temp, 2, quantile, probs = c(0.975)) # c(0.975)
    
    transm_rate_mat <- rbind(transm_rate_mat,
                             transm_rate_group_tmp)
    
  }# End for
  
  #---- Output:
  output <- transm_rate_mat
  
  return(output)
  
}# End function

post_transm_rate <- transm_rate_generation(cov_data      = cov_data,
                                           model_out     = nuts_fit_1,
                                           posts_simul   = posts_simul_v2,
                                           aggr_cm       = aggr_cm,
                                           selected_draw = selected_draw
)


#---- Plot of posterior estimates:
posterior_transm_rate <- 
  ggplot(post_transm_rate,
         aes(Date = Date,
             true = True)) +
  facet_wrap(. ~ Group, 
             scales = "free_y", 
             ncol   = 1,
             strip.position = "top") +
  geom_point(aes(x   = Date,
                 y   = True,
                 fill = "True value"
  ),
  shape = 4
  ) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"
  ),
  size = 1.0
  ) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"
  ),
  alpha = 0.5) +
  labs(x = "Day",
       y = "Transmission rate",
       fill = "") +
  scale_fill_manual(values = c('True value' = "black")) +
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black", 
                                 "50% CrI" = "gray40"    
                      )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(size = 16, angle = 0, hjust = 1),
        axis.text.y      = element_text(size = 16),
        axis.title.x     = element_text(size = 18, face = "bold"),
        axis.title.y     = element_text(size = 18, face = "bold"),
        legend.text      = element_text(size = 16),
        strip.text.x     = element_text(size = 14),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())


###########################################################################
#
# Effective reproduction number 
#
###########################################################################
effective_rt_seeiir <- function(cov_data,
                                model_out,
                                posts_simul,
                                aggr_cm,
                                selected_draw,
                                progress_bar = FALSE){
  
  `%nin%` <- Negate(`%in%`)
  
  #---- Checks:
  if(class(model_out)[1] != "stanfit") stop("Provide an object of class 'stanfit' using rstan::sampling()")
  
  #---- Posterior draws:
  posterior_draws <- rstan::extract(model_out)
  
  #---- Checks:
  if(ncol(posterior_draws$cm_sample) != cov_data$A) stop( paste0("The number of rows in the age distribution table must be equal to ", cov_data$A) )
  
  #---- Functions:
  eigen_mat  <- function(mat) max( Re( eigen(mat)$values ) )
  
  #---- Assistant matrices:
  age_grps             <- cov_data$A
  zero_mat             <- matrix(0L, nrow = age_grps, ncol = age_grps)  
  identity_mat         <- diag(age_grps)
  reciprocal_age_distr <- matrix(rep(cov_data$pop_diag, age_grps),   ncol  = age_grps, nrow  = age_grps, byrow = TRUE)
  age_distr            <- matrix(rep(1/cov_data$pop_diag, age_grps), ncol  = age_grps, nrow  = age_grps, byrow = FALSE)
  
  # Step 1 - Calculate Q^{-1} at iteration i:  
  Q_inverse <- cov_data$dI * identity_mat
  
  #---- Output storage:
  beta_draws   <- posterior_draws$beta_trajectory
  chain_length <- nrow(beta_draws)
  ts_length    <- dim(beta_draws)[2]
  
  R_eff_mat     <- R_mat     <- matrix(0L, nrow = chain_length, ncol = ts_length)
  R_eff_t_obs <- R_t_obs <- rep(0L, lentgh = ts_length)
  
  if(progress_bar == TRUE) pb = txtProgressBar(min = 0, max = ts_length, initial = 0, style = 3) 
  
  #---- Calculation of the effective at iteration i and time point j:
  for (j in 1:ts_length){
      
    # Step 2.1: simulated quantities:
    observed_effcr <- posts_simul$beta_trajectory[selected_draw,,]
    
    B_tmp_obs <- 
      matrix( rep(observed_effcr[j,], age_grps),
                         ncol  = age_grps, 
                         nrow  = age_grps, 
                         byrow = FALSE) * 
      data.matrix(aggr_cm) *
     age_distr *
     reciprocal_age_distr
    
    B_eff_tmp_obs <- 
      matrix( rep(observed_effcr[j,], age_grps),
                             ncol  = age_grps, 
                             nrow  = age_grps, 
                             byrow = FALSE) *
      data.matrix(aggr_cm)  *
      matrix(rep( posts_simul$Susceptibles[selected_draw,j,], age_grps),
             ncol  = age_grps, 
             nrow  = age_grps, 
             byrow = FALSE) *
      reciprocal_age_distr
    
    BQinv_tmp_obs     <- B_tmp_obs     %*% Q_inverse
    BQinv_eff_tmp_obs <- B_eff_tmp_obs %*% Q_inverse
    
    R_t_obs[j]     <- eigen_mat(BQinv_tmp_obs)
    R_eff_t_obs[j] <- eigen_mat(BQinv_eff_tmp_obs)
    
    # Step 2.2: MCMC output:
    for (i in 1:chain_length) {
  
      B_tmp <- matrix( rep(beta_draws[i,,][j,], age_grps),
                       ncol  = age_grps, 
                       nrow  = age_grps, 
                       byrow = FALSE) *
        matrix( posterior_draws$cm_sample[i,,], 
                nrow = age_grps, 
                ncol = age_grps) *
        age_distr *
        reciprocal_age_distr
      
      B_eff_tmp <- matrix( rep(beta_draws[i,,][j,], age_grps),
                           ncol  = age_grps, 
                           nrow  = age_grps, 
                           byrow = FALSE) *
        matrix( posterior_draws$cm_sample[i,,], 
                nrow = age_grps, 
                ncol = age_grps) *
        matrix(rep( posterior_draws$Susceptibles[i,,][j,], age_grps),
               ncol  = age_grps, 
               nrow  = age_grps, 
               byrow = FALSE) *
        reciprocal_age_distr
      
      # Step 3 calculate the B_t %*% Q^{-1} matrix:
      BQinv_tmp     <- B_tmp     %*% Q_inverse
      BQinv_eff_tmp <- B_eff_tmp %*% Q_inverse
      
      # Step 4 - calculate the spectral radius of B_t %*% Q^{-1} at iteration i and time point j:
      R_mat[i,j]     <- eigen_mat(BQinv_tmp)
      R_eff_mat[i,j] <- eigen_mat(BQinv_eff_tmp)
      
      
    }# End for
    
    if(progress_bar == TRUE) setTxtProgressBar(pb,j)
    
  }# End for
  
  output <- list(R_t         = R_mat,
                 R_eff_t     = R_eff_mat,
                 R_t_obs     = R_t_obs,
                 R_eff_t_obs = R_eff_t_obs)
  
  #---- Export:
  return(output)
  
}# End function


#---- Compilation of the function:
effective_rt_seeiir <- compiler::cmpfun(effective_rt_seeiir)

#---- Estimation of the effective reproduction number:
time_start_Rt <- Sys.time()
post_Rt <- effective_rt_seeiir(cov_data     = cov_data,
                               model_out    = nuts_fit_1,
                               posts_simul  = posts_simul_v2,
                               aggr_cm      = aggr_cm,
                               selected_draw= selected_draw,
                               progress_bar = TRUE)
time_end_Rt <- Sys.time()
duration_Rt <- time_end_Rt - time_start_Rt

save(post_Rt,
     duration_Rt,
     file = rt_path)

data_R_eff_t         <- data.frame(Date = 1:cov_data$n_obs)
data_R_eff_t$True    <- post_Rt$R_eff_t_obs
data_R_eff_t$median  <- apply(post_Rt$R_eff_t, 2, median)
data_R_eff_t$low     <- apply(post_Rt$R_eff_t, 2, quantile, probs = c(0.025)) # c(0.025)
data_R_eff_t$low25   <- apply(post_Rt$R_eff_t, 2, quantile, probs = c(0.25))  # c(0.025)
data_R_eff_t$high75  <- apply(post_Rt$R_eff_t, 2, quantile, probs = c(0.75))  # c(0.975)
data_R_eff_t$high    <- apply(post_Rt$R_eff_t, 2, quantile, probs = c(0.975)) # c(0.975)

#---- Plot of posterior estimates:
posterior_eff_reprod_number <- 
  ggplot(data_R_eff_t,
         aes(Date = Date,
             true = True)) +
  geom_point(aes(x   = Date,
                 y   = True,
                 fill = "True value"
  ),
  shape = 4
  ) +
  geom_line(aes(x     = Date,
                y     = median,
                color = "Median"
  ),
  size = 1.0
  ) +
  geom_ribbon(aes(x    = Date,
                  ymin = low25,
                  ymax = high75,
                  fill = "50% CrI"
  ),
  alpha = 0.5) +
  geom_ribbon(aes(x    = Date,
                  ymin = low,
                  ymax = high,
                  fill = "95% CrI"
  ),
  alpha = 0.5) +
  labs(x    = "Day",
       y    = "Effective reproduction number",
       fill = "") +
  geom_hline(yintercept = 1, color = "black") +
  scale_y_continuous(breaks = c(0,1,2,4,6,8,10,12,14,16)) +
  scale_fill_manual(values = c('True value' = "black")) +
  guides(color = "none", fill = "none") + 
  scale_colour_manual(name = '', 
                      values = c('Median'  = "black",   
                                 "50% CrI" = "gray20",  
                                 "95% CrI" = "gray60"   
                      )
  ) +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(size = 16, angle = 0, hjust = 1),
        axis.text.y      = element_text(size = 16),
        axis.title.x     = element_text(size = 18, face = "bold"),
        axis.title.y     = element_text(size = 18, face = "bold"),
        legend.text      = element_text(size = 16),
        strip.text.x     = element_text(size = 14),
        legend.position  = "none",
        legend.title     = element_blank(),
        legend.box       = "horizontal", 
        legend.margin    = margin())

#########################################################
# 
# Percentage of exposed individuals per group
# 
#########################################################
pi_perc <- 0.10

# Create the appropriate dataset in a format acceptable by ggplot:
dt_pi_post  <- data.frame(Posterior = posts_mcmc[["rho"]])
dt_pi_post2 <- reshape2::melt(dt_pi_post)

# Create the actual col-row labels of the aggregated matrix, from the vectorized version:
# The aggregated contact matrix has been vectorized by major row order.
posterior_rho <- 
ggplot(dt_pi_post2, 
       aes(x    = value)) +
  geom_density(alpha = 0.8) +
  scale_x_continuous(
    limits = c(min(dt_pi_post2$value)*0.95, 
               max(dt_pi_post2$value)*1.05)
  ) +
  geom_vline(xintercept = standata_sim$pi, linetype = "dashed") +
  labs(x = expression(rho),
       y = "Density") +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(angle = 0, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())

#########################################################
#
# NegBin overdispersion parameter: 
#
#########################################################

# Create the appropriate dataset in a format acceptable by ggplot:
dt_phiD_post  <- data.frame(PPosterior = posts_mcmc[["phiD"]])
dt_phiD_post2 <- reshape2::melt(dt_phiD_post)

# Create the actual col-row labels of the aggregated matrix, from the vectorized version:
# The aggregated contact matrix has been vectorized by major row order.
posterior_phi <- 
ggplot(dt_phiD_post2, 
       aes(x    = value)) +
  geom_density(alpha = 0.8) +
  scale_x_continuous(
    limits = c(min(dt_phiD_post2$value)*0.95, 
               max(dt_phiD_post2$value)*1.05)
  ) +
  geom_vline(xintercept = standata_sim$phiD, linetype = "dashed") +
  labs(x = expression(phi),
       y = "Density") +
  theme_bw() +
  theme(panel.spacing    = unit(0.4,"cm"),
        axis.text.x      = element_text(angle = 0, hjust = 1),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.box       = "vertical", 
        legend.margin    = margin())

#########################################################
# 
# Age-specific volatility parameters
# 
#########################################################

niters            <- nrow(posts_mcmc[["sigmaBM"]])
prior_sigmaBM     <- array(numeric(),c(niters, cov_data$A)) 
plot_sigmaBM_list <- list()
plot_sigmaBM_indx <- 1

for (i in 1:cov_data$A){
    
  dt_sigmaBM_post <- data.frame(Posterior = posts_mcmc[["sigmaBM"]][,i])
  dt_sigmaBM_post2 <- reshape2::melt(dt_sigmaBM_post)
  
  p <- ggplot(dt_sigmaBM_post2, 
              aes(x    = value)) +
    geom_density(alpha = 0.8) +
    scale_x_continuous(
      limits = c(min(dt_sigmaBM_post2$value)*0.95, max(dt_sigmaBM_post2$value)*1.05)
    ) +
    geom_vline(xintercept = standata_sim$sigmaBM[i], linetype = "dashed") +
    labs(x = bquote(sigma[.(i)]), #bquote(sigma[x]^.(i))
         y = "Density",
         title = paste0("Age group: ", aggr_age$AgeGrp[i])) + 
    theme_bw() +
    theme(panel.spacing    = unit(0.4,"cm"),
          axis.text.x      = element_text(angle = 0, hjust = 1),
          axis.title.x     = element_text(size = 14, face = "bold"),
          axis.title.y     = element_text(size = 14, face = "bold"),
          legend.position  = "bottom",
          legend.title     = element_blank(),
          legend.box       = "vertical", 
          legend.margin    = margin())
  
  plot_sigmaBM_list[[plot_sigmaBM_indx]] <- p
  
  # Increment the index of stored graph:
  plot_sigmaBM_indx <- plot_sigmaBM_indx + 1
  
}# End for

# Remove redundant objects:
rm(i, dt_sigmaBM_post, dt_sigmaBM_post2)

# Plot the distributions for rows 1-2:
posterior_sigma <- do.call("grid.arrange", c(plot_sigmaBM_list, nrow = 1))

#########################################################
#
# Print the combined plot:
#
#########################################################
t1 <- patchwork::wrap_plots(plot_cm_list, 
                            nrow = cov_data$A, 
                            ncol = cov_data$A)

combined <- 
  posterior_eff_cont_rate_95CrI + 
  t1 +
  posterior_eff_reprod_number +
  posterior_deaths &
  theme(legend.position = "bottom") 
 
combined + 
  plot_annotation( tag_levels = list(c("A", "B", rep("",8), "C", "D", "E")) ) + 
  plot_layout(nrow = 2, guides = "collect")






