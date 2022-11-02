###########################################################################
#
# Effective reproduction number - Next generation matrix approach 
# SEEIIR transmission model
#
###########################################################################
effective_rt_seeiir <- function(cov_data,
                                model_out,
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
  # The R matrix must have the same dimensions as the matrix posterior_draws$beta_N:
  if ( "beta_weekly" %in% names(posterior_draws) ){ 
    
    beta_draws   <- posterior_draws$beta_weekly 
    chain_length <- nrow(beta_draws)
    ts_length    <- dim(beta_draws)[2]
    
  } else if ( "beta_N" %in% names(posterior_draws) ) {
    
    beta_draws   <- posterior_draws$beta_N
    chain_length <- nrow(beta_draws)
    ts_length    <- dim(beta_draws)[2]
    
  }# End for
  
  R_eff_mat <- R_mat <- matrix(0L, nrow = chain_length, ncol = ts_length)
  
  if(progress_bar == TRUE) pb = txtProgressBar(min = 0, max = chain_length, initial = 0, style = 3) 
  
  #---- Calculation of the effective at iteration i and time point j:
  for (i in 1:chain_length) {
    for (j in 1:ts_length){
      
      # Step 2 - Create the B_t matrix:
      # NOTE: perform element-wise multiplications:
      
      # Common GBM across groups:
      if( is.na(dim(beta_draws)[3]) ){
        
        B_tmp <- beta_draws[i,j] *
                 matrix( posterior_draws$cm_sample[i,,], 
                         nrow = age_grps, 
                         ncol = age_grps) *
                 age_distr *
                 reciprocal_age_distr
        
        B_eff_tmp <- beta_draws[i,j] *
                     matrix( posterior_draws$cm_sample[i,,], 
                             nrow = age_grps, 
                             ncol = age_grps) *
                     matrix(rep( posterior_draws$Susceptibles[i,,][j,], age_grps),
                            ncol  = age_grps, 
                            nrow  = age_grps, 
                            byrow = FALSE) *
                    reciprocal_age_distr
        
      # Multiple GBMs: 
      } else if ( !is.na(dim(beta_draws)[3]) ){
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
        
      }# End if  
      
      # Step 3 calculate the B_t %*% Q^{-1} matrix:
      BQinv_tmp     <- B_tmp %*% Q_inverse
      BQinv_eff_tmp <- B_eff_tmp %*% Q_inverse
      
      # Step 4 - calculate the spectral radius of B_t %*% Q^{-1} at iteration i and time point j:
      R_mat[i,j]     <- eigen_mat(BQinv_tmp)
      R_eff_mat[i,j] <- eigen_mat(BQinv_eff_tmp)
      
    }# End for
    
    #---- Cleanup:
    #rm(B_inside_tmp, B_tmp, BQinv_tmp)
    
    if(progress_bar == TRUE) setTxtProgressBar(pb,i)
    
  }# End for
  
  # Keep every 7 days if we assume weekly GBMs:
  if("n_weeks" %in% names(cov_data) ){
    R_mat     <- t( t(R_mat)     %>% as.data.frame %>% slice(which(row_number() %% 7 == 1)))
    R_eff_mat <- t( t(R_eff_mat) %>% as.data.frame %>% slice(which(row_number() %% 7 == 1)))
  }# End if
  
  output <- list(R_t     = R_mat,
                 R_eff_t = R_eff_mat)
  
  #---- Export:
  return(output)
  
}# End function

#---- Compilation of the function:
effective_rt_seeiir <- compiler::cmpfun(effective_rt_seeiir)
