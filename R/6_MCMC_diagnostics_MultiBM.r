###########################################################################
#
# Path to posterior outputs:
#
###########################################################################
path_to_stan_output <- "..."

###########################################################################
#
# Load model output (stored in the object nuts_fit_1):
#
###########################################################################
load(file = path_to_stan_output)

#---- Session info:
sessionInfo()

mod1_diagnostics <- rstan::get_sampler_params(nuts_fit_1)
posts_1          <- rstan::extract(nuts_fit_1)
lp_cp            <- bayesplot::log_posterior(nuts_fit_1)
np_cp            <- bayesplot::nuts_params(nuts_fit_1)
posterior_1      <- as.array(nuts_fit_1)

###########################################################################
#
# Subset the parameters
#
###########################################################################

granularity <- "daily" # c("weekly", "daily")

model <- "multi-BM"    # vs "single-BM"

cm_params   <- paste0("cm_sample[", 
                      apply(expand.grid(1:cov_data$A, 1:cov_data$A), 1, paste, collapse = ","),
                      "]")

if (model == "single-BM") {
  sigmaBM_params <- "sigmaBM"
  beta_params    <- paste0("beta_N[", 1:cov_data$n_obs ,"]")
  
} else {
  sigmaBM_params <- paste0("sigmaBM[", 1:cov_data$A,"]")
  beta_params    <- paste0("beta_N[", apply(expand.grid(1:cov_data$n_obs, 1:cov_data$A), 1, paste, collapse = ","), "]")
}

rest_params <- c("pi",
                 sigmaBM_params,
                 "phiD")

##########################################################################
#
# MCMC summaries
#
###########################################################################
nuts_fit_1_summary <- summary(nuts_fit_1,
                              pars = c("lp__",
                                       rest_params,
                                       "cm_sample",
                                       beta_params
                                       ))$summary

options(scipen = 999)
round(nuts_fit_1_summary, 3)

###########################################################################
#
# Diagnostics for the No-U-Turn Sampler
#
###########################################################################

# Check for divergent transitions:
rstan::check_divergences(nuts_fit_1)
rstan::check_treedepth(nuts_fit_1)

check_div(nuts_fit_1)
check_treedepth(nuts_fit_1)

# Global patterns in the divergences:
color_scheme_set("darkgray")
bayesplot::mcmc_parcoord(posterior_1)

# Understand how the divergences interact with the model globally:
# Top: Distribution of the log-posterior when there was no divergence vs
# the distribution when there was a divergence
# Bottom: Same, for the NUTS acceptance statistic.
color_scheme_set("red")
bayesplot::mcmc_nuts_divergence(np_cp, lp_cp)

###########################################################################
#
# Energy and Bayesian fraction of missing information
#
###########################################################################

color_scheme_set("red")
bayesplot::mcmc_nuts_energy(np_cp)

check_energy(nuts_fit_1)

#---- Effective sample size per iteration (efficiency of the sampler):
check_n_eff(nuts_fit_1)

#---- Checks the potential scale reduction factors
check_rhat(nuts_fit_1)

###########################################################################
#
# Plotting MCMC output
#
###########################################################################
color_scheme_set("viridis")

#---- Marginal posterior distributions (combining all chains):
bayesplot::mcmc_dens(posterior_1,
                     pars = c(ifr_params))

bayesplot::mcmc_dens(posterior_1,
                     pars = c(rest_params))

bayesplot::mcmc_dens(posterior_1,
                     pars = cm_params)

#---- Prior of the NegBin overdispersion parameter: 
set.seed(1)
niters            <- nrow(posts_1[["phiD"]])
plot_phiD_indx <- 1

# Generate draws from the prior distribution:
prior_phiD   <- rexp(niters, cov_data$p_phiD) 

# Create the appropriate dataset in a format acceptable by ggplot:
dt_phiD_post <- data.frame(Prior     = prior_phiD,
                           Posterior = posts_1[["phiD"]])
dt_phiD_post2 <- reshape2::melt(dt_phiD_post)

# Create the actual col-row labels of the aggregated matrix, from the vectorized version:
# The aggregated contact matrix has been vectorized by major row order.
ggplot(dt_phiD_post2, 
       aes(x    = value, 
           fill = variable)) +
  geom_density(alpha = 0.8) +
  scale_x_continuous(
    limits = c(min(dt_phiD_post2$value)*0.95, 
               max(dt_phiD_post2$value)*1.05)
  ) +
  theme(strip.placement  = "outside",
        strip.background = element_rect(fill = NA, colour="grey50"),
        panel.spacing    = unit(0,"cm"),
        legend.position  = "bottom",
        legend.title     = element_blank() )

# Remove redundant objects:
rm(dt_phiD_post, dt_phiD_post2)

#---- Prior of the GBM volatility vs the posterior: 
set.seed(1)
niters            <- nrow(posts_1[["sigmaBM"]])
plot_sigmaBM_indx <- 1

# Generate draws from the prior distribution:
#prior_sigmaBM   <- extraDistr::rhcauchy(niters, cov_data$p_sigmaBM[2])
prior_sigmaBM <- extraDistr::rhnorm(niters, cov_data$p_sigmaBM[2])

# Create the appropriate dataset in a format acceptable by ggplot:
dt_sigmaBM_post <- data.frame(Prior     = prior_sigmaBM,
                              Posterior = posts_1[["sigmaBM"]]
                              )
dt_sigmaBM_post2 <- reshape2::melt(dt_sigmaBM_post)

# Create the actual col-row labels of the aggregated matrix, from the vectorized version:
# The aggregated contact matrix has been vectorized by major row order.
ggplot(dt_sigmaBM_post2, 
          aes(x    = value, 
              fill = variable)) +
geom_density(alpha = 0.8) +
scale_x_continuous(
        limits = c(min(dt_sigmaBM_post2$value)*0.95, 
                   max(dt_sigmaBM_post2$value)*1.05)
) +
theme(strip.placement  = "outside",
      strip.background = element_rect(fill = NA, colour="grey50"),
      panel.spacing    = unit(0,"cm"),
      legend.position  = "bottom",
      legend.title     = element_blank() )
  
# Remove redundant objects:
rm(dt_sigmaBM_post, dt_sigmaBM_post2)

#---- Pairs plots:
bayesplot::mcmc_pairs(posterior_1,
                      np = np_cp,
                      pars = c(cm_params),
                      off_diag_args = list(size = 1.5))

bayesplot::mcmc_pairs(posterior_1,
                      np = np_cp,
                      pars =  c(beta_params, rest_params[2]),
                      off_diag_args = list(size = 1.5))

bayesplot::mcmc_pairs(posterior_1,
                      np = np_cp,
                      pars =  rest_params,
                      off_diag_args = list(size = 1.5))

#---- Trace plot of MCMC draws:
bayesplot::mcmc_trace(posterior_1,
                      pars = c(beta_params, rest_params[2]),
                      facet_args = list(ncol = 5,
                                        strip.position = "left"))

bayesplot::mcmc_trace(posterior_1,
                      pars = cm_params,
                      facet_args = list(ncol = 2,
                                        strip.position = "left"))

bayesplot::mcmc_trace(posterior_1,
                      pars = rest_params,
                      facet_args = list(ncol = 4,
                                        strip.position = "left"))

#---- Autocorrelation:
bayesplot::mcmc_acf(posterior_1,
                    pars = cm_params,
                    lags = 100)  

bayesplot::mcmc_acf(posterior_1,
                    pars = rest_params,
                    lags = 100)

#---- Central posterior uncertainty intervals (50%, 90%):
bayesplot::mcmc_intervals(posterior_1,
                          pars = c(beta_params))

bayesplot::mcmc_intervals(posterior_1,
                          pars = c(cm_params))

bayesplot::mcmc_intervals(posterior_1,
                          pars = c(rest_params))

#---- Uncertainty intervals as shaded areas under the estimated posterior density curves:
bayesplot::mcmc_areas(posterior_1,
                      pars = c(cm_params),
                      prob = 0.8,        # 80% intervals
                      prob_outer = 0.99, # 99%
                      point_est = "mean")

bayesplot::mcmc_areas(posterior_1,
                      pars = c(rest_params),
                      prob = 0.8,        # 80% intervals
                      prob_outer = 0.99, # 99%
                      point_est = "mean")

###########################################################################
#
# DIC & Effective number of parameters
#
###########################################################################

model <- "Multi-BM, Austria"

#---- Dhat:
E_deathsByAge_params <- paste0(
  "E_deathsByAge[",
  apply(expand.grid(1:cov_data$n_obs, 1:cov_data$A), 1, paste, collapse = ","),
  "]")
E_deathsByAge_params_summary <- summary(nuts_fit_1, pars = E_deathsByAge_params)$summary

loglik <- 0

for (i in 1:cov_data$n_obs) {
  for (j in 1:cov_data$A) {
    if (model == "Multi-BM, Austria") {
      size <- nuts_fit_1_summary[rownames(nuts_fit_1_summary) %in% "phiD","mean"]
    } else {
      size <- E_deathsByAge_params_summary[
              grepl(paste0("E_deathsByAge[",i,",",j,"]"), 
                    rownames(E_deathsByAge_params_summary), fixed = TRUE), "mean"] /
              nuts_fit_1_summary[rownames(nuts_fit_1_summary) %in% "phiD","mean"]
    }
    
    loglik <- loglik + 
              dnbinom(
                cov_data$y_deaths[i,j], 
                mu   = E_deathsByAge_params_summary[
                  grepl(paste0("E_deathsByAge[",i,",",j,"]"), 
                        rownames(E_deathsByAge_params_summary), fixed = TRUE), "mean"],
                size = size,
                log = TRUE)
  }# End for
}# End for  

Dhat <- -2*loglik

#---- Dbar:
Dbar <- mean(posts_1$dev)

#---- Effective number of parameters:
pD  <- Dbar - Dhat
pV  <- var(posts_1$dev)/2
DIC <- Dbar + pD

#---- Gather the estimated values:
model_comparison_dic        <- round( c(Dbar, Dhat, pD, pV, DIC), 2)
names(model_comparison_dic) <- c("Dbar", "Dhat", "pD", "pV", "DIC")
print(model_comparison_dic)