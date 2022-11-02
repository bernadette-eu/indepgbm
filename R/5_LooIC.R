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

###########################################################################
#
# Calculate PSIS-LooIC and WAIC:
#
###########################################################################
parallel::detectCores()

nChains   <- 6

log_lik_1 <- loo::extract_log_lik(nuts_fit_1, merge_chains = FALSE)
r_eff_1   <- loo::relative_eff(exp(log_lik_1),    cores = nChains)
loo_1     <- loo::loo(log_lik_1, r_eff = r_eff_1, cores = nChains)
waic_1    <- loo::waic(log_lik_1)

print(loo_1)
print(waic_1)