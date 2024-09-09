functions {
    // Source: https://discourse.mc-stan.org/t/lower-upper-triangular-matrix-data-types/4361/6
	//         https://www.mathworks.com/matlabcentral/answers/152818-how-to-covert-a-vector-of-lower-triangular-matrix-factors-back-to-the-original-matrix
	matrix to_triangular(vector x, int K) {
	
		matrix[K, K] y = rep_matrix(0.0, K, K); //Declare a matrix of zeros to avoid NaNs in the upper triangular part. 
		int pos = 1;
			
		for (col in 1:K) {
			for (row in col:K) {
			  y[row, col] = x[pos];
			  pos += 1;
			}
		}
		
		//for (row in 1:K) {
		//	for (col in 1:row) {
		//	  y[row, col] = x[pos];
		//	  pos += 1;
		//	}
		//}
		
		return y;
    }// End function
	
	// Source: https://stackoverflow.com/questions/60810640/what-is-the-equivalent-to-the-r-function-repx-each-n-in-stan
	real[] rep_each(real[] x, int K) {
		int N = size(x);
		real y[N  *  K];
		int pos = 1;
		
		for (n in 1:N) {
		  for (k in 1:K) {
			y[pos] = x[n];
			pos += 1;
		  }
		}
		return y;
	}// End function
	
	// Source: https://jrnold.github.io/ssmodels-in-stan/stan-functions.html
	real[] to_vector_rowwise(matrix x) {
	  
	  real res[num_elements(x)];
	  int n;
	  int m;
	  
	  n = rows(x);
	  m = cols(x);
	  
	  for (i in 1:n) for (j in 1:m) res[(i - 1) * m + j] = x[i, j];
  
	  return res;
	}// End function
	
	real[] to_vector_colwise(matrix x) {
		
	  real res[num_elements(x)];
	  int n;
	  int m;
	  
	  n = rows(x);
	  m = cols(x);
	  
	  for (i in 1:n) for (j in 1:m) res[n * (j - 1) + i] = x[i, j];

	  return res;
	}// End function
	
	matrix repeat_matrix(matrix input, int K) {
	  int N = rows(input);
	  int M = cols(input);
	  matrix[N * K, M] repmat; // stack N*M matrix K times
	  int pos = 1;

		for (n in 1:N) {
		  for (k in 1:K) {
		  repmat[pos,] = to_row_vector(input[n,]);
		  pos += 1;
		  }
		}
	  return repmat;
	}

	matrix repeat_rv_to_matrix(row_vector input, int K) {
	  int M = num_elements(input);
	  matrix[K, M] repmat;
	  int pos = 1;

	  for (k in 1:K) {
		repmat[pos,] = input;
		pos += 1;
	  }
	  return repmat;
	}


    // NOTE: The length of the returned real array must match the length of the state input into the function.
    // Source: https://mc-stan.org/docs/2_22/functions-reference/functions-ode-solver.html
    
	real[] ODE_states(real time,     // Time
					  real[] y,      // System state {susceptible,infected,recovered}
					  real[] theta,  // Parameters
					  real[] x_r,    // Real-type data
					  int[] x_i      // Integer-type data
					  ){
	
	  int A       = x_i[1];  // Number of age groups
	  int n_obs   = x_i[2];  // Length of the time series
	  int n_difeq = x_i[3];  // Number of differential equations in the system
	  
	  real dy_dt[A * n_difeq]; // SEIR (ignoring R) then C 
	  real f_inf[A];           // Force of infection
      real init[2 * A];        // Initial values at (S, E1) compartmens

      real age_dist[A]    = x_r[(2 * n_obs + 1):(2 * n_obs + A)];  // Population per age group

      // Estimated parameters:
      real contact[A * A] = theta[1:(A * A)]; // Sampled contact matrix in vectorized format.
	                                          // First A values correspond to number of contact between age class 1 and other classes, etc.
	  real gamma          = theta[A * A + 1]; // Recovery rate
	  real tau            = theta[A * A + 2]; // Incubation rate
																										 
      real pi             = theta[A * A + 3];                 // Number of cases at t0
	  real beta[A]        = theta[(A * A + 4):(A*A + A + 3)]; // Effective contact rate

      // Compartments:
      for (i in 1:A){
        
		init[i]         = age_dist[i]  * (1-pi); // Initial states - Susceptibles
        init[A + i]     = age_dist[i]  *  pi;    // Initial states - Exposed 1
		
	    // Force of infection by age class (SEEIIR model):
		f_inf[i] = sum( to_vector(beta).* ( to_vector( y[(3*A+1):(4*A)] ) + 
		                                    to_vector( y[(4*A+1):(5*A)] ) 
		                                   ) ./ to_vector(age_dist) .* to_vector(contact[(A*(i-1)+1):(i*A)]) ); 
		   
        				
        // S: susceptible
        dy_dt[i] = - f_inf[i] * ( y[i] + init[i] ); 
		
		// E1: incubating (not yet infectious)
        dy_dt[A + i] = f_inf[i] * ( y[i] + init[i] ) - tau  *  ( y[A + i] + init[A + i] );
 
 		// E2: incubating (not yet infectious)
        dy_dt[2 * A + i] = tau * ( ( y[A + i] + init[A + i] ) - y[2 * A + i]  );
 
        // I1: infectious
        dy_dt[3 * A + i] = tau  * y[2 * A + i] - gamma  *  y[3 * A + i];

        // I2: infectious
        dy_dt[4 * A + i] = gamma  *  ( y[3 * A + i] - gamma  *  y[4 * A + i] );

        // C: cumulative number of infections by date of disease onset
        dy_dt[(n_difeq-1) * A + i] = tau * y[2 * A + i];
		
       }// End for

      return dy_dt;
    }// End SEIR function

    // Integration using the trapezoidal rule:
	real[,] integrate_ode_trapezoidal(real[] y_initial, 
									  real initial_time, 
									  real[] times,  // Vector of time indexes				  
									  real[] theta,  // Parameters
					                  real[] x_r,    // Real-type data
					                  int[] x_i      // Integer-type data
									  )
    {
                                        
      real h;
      vector[size(y_initial)] dy_dt_initial_time;
      vector[size(y_initial)] dy_dt_t;
      vector[size(y_initial)] k;

      real y_approx[size(times),size(y_initial)];
	  
      int A       = x_i[1];  // Number of age groups
	  int n_obs   = x_i[2];
      real theta_ODE[A*A + A + 3];

	  real left_t[n_obs]  = x_r[1:n_obs];              // Left and right time bounds for the calculation of the time-dependent incidence rate
      real right_t[n_obs] = x_r[(n_obs+1):(2 * n_obs)];// Left and right time bounds for the calculation of the time-dependent incidence rat
      real beta_N_temp[A*n_obs] = theta[(A*A + 3):(A*n_obs + A*A + 2)];

      // Define the parameter vector that enters the ode_states():		 
      theta_ODE[1:(A * A)] = theta[1:(A * A)];         // Vectorized contact matrix
      theta_ODE[A * A + 1] = theta[A * A + 1];         // gamma = Recovery rate
      theta_ODE[A * A + 2] = theta[A*n_obs + A*A + 4]; // tau   = Incubation rate
      theta_ODE[A * A + 3] = theta[A*n_obs + A*A + 3]; // pi
	  
      for (t in 0:(size(times)-1)) {
        //print("t: ", t);
		//print("beta0: ", theta_ODE[A * A + 4]);
		//print("theta[A*A + t + 2]: ", theta[A*A + t + 2]);
		//print("theta[A*A + n_obs + 2]: ", theta[A * A + n_obs + 2]);
		//print("theta_ODE[A * A + 4]: ", theta_ODE[A * A + 4]);

		if(t == 0){
		
		   for (j in 1:A) theta_ODE[A * A + 3 + j] = theta[A * A + 2]; // beta0

		   h = times[1] - initial_time;
		   dy_dt_initial_time = to_vector(ODE_states(initial_time, y_initial, theta_ODE, x_r, x_i));
		   k = h*dy_dt_initial_time;

		   y_approx[t+1,] = to_array_1d(
		                      to_vector(y_initial) +
							  h*(dy_dt_initial_time +
							     to_vector(ODE_states(times[1],
												      to_array_1d(to_vector(y_initial) + k),
												      theta_ODE, x_r, x_i)))/2);
			
		} else {
			h = (times[t+1] - times[t]);
            
			for (j in 1:A){
				 // Assign the effective contact rate parameter at the last time point
				if (t == (size(times) - 1) ) theta_ODE[A * A + 3 + j] = beta_N_temp[n_obs * (j - 1) +  t + 1];
				else if(t >= left_t[t] && t <= right_t[t]) theta_ODE[A * A + 3 + j] = beta_N_temp[n_obs * (j - 1) +  t];
			}// End for
			
		    dy_dt_t = to_vector(ODE_states(times[t], y_approx[t], theta_ODE, x_r, x_i));

            k = h*dy_dt_t;

			y_approx[t+1,] = to_array_1d(
			                    to_vector(y_approx[t,]) +
								h*(dy_dt_t + to_vector(ODE_states(times[t+1],
															      to_array_1d(to_vector(y_approx[t,]) + k),
															      theta_ODE, x_r, x_i)
															      )) /2 );
		
		}// End if
      }// End for

      return y_approx;

    } // End trapezoidal rule function
}

data {							  
	// Structure:								  
	int<lower = 1> A;                 // Number of age groups
	int<lower = 1> n_obs;             // Length of analysis period

	int<lower = 1> n_pop;             // Population
	
	int<lower = 1, upper = 7> ecr_changes;
	int n_changes;
	int n_remainder;
	
	real age_dist[A];                 // Age distribution of the general population
	vector[A] pop_diag;               // Inverse of population for each age group

	int<lower = 1> n_difeq;           // Number of differential equations (S,I,C)

	vector[A] L_cm[A];                // Lower triangular matrix, stemming from the Cholesky decomposition of the observed Contact matrix
	real<lower = 0> ifr_age[A];       // Infection-fatality rate per age group

	real t0;                          // Initial time point (zero)
	real ts[n_obs];                   // Time bins
	real<lower=0> left_t[n_obs];      // Left time limit
	real<lower=0> right_t[n_obs];     // Right time limit
	vector<lower = 0>[n_obs] I_D;     // Discretized infection to death distribution.          
	row_vector[A] E_deathsByAge_day1; // Age-group deaths at day 1 of the analysis

	// Data to fit:
	//int y_deaths[n_obs,A];            // Mortality data (new daily deaths), per age group
									// NOTE: the function "matrix<lower = 0>[n_obs,A] y_deaths;" led to an error message in the 
									// negative_binomial_2() sampling statement, that the y_deaths object was not of type "int".

	// Fixed parameters:
	real dE;
	real dI;            

	real eta0;                        // Initial transmission rate
	real x_init[A];
	real p_sigmaCM;                

    real<lower = 0, upper = 1> pi;
	real<lower = 0> sigmaBM[A];       // Standard deviation of GBM
	real<lower = 0> phiD;             // Likelihood variance parameter

	// Debugging:
	//int inference; // 0: simulating from priors; 1: fit to data
}

transformed data {
  vector<lower = 0>[n_obs] I_D_rev; // Reversed discretized infection-to-death distribution
  
  int x_i[3];
  real x_r[2 * n_obs + A];
  real<lower = 0> gamma;
  real<lower = 0> tau;
  real<lower = 0> beta0;                     
  
  real init[A * n_difeq] = rep_array(0.0, A * n_difeq); // Initial conditions for the (S,E,I,C) compartments
  
  vector[A] ones_vector_A = rep_vector(1.0, A);
  
  vector[(A * (A + 1)) / 2] L_vector = rep_vector(0, (A * (A + 1)) / 2);

  for(i in 1:n_obs) I_D_rev[i] = I_D[n_obs - i + 1];

  x_i[1] = A;
  x_i[2] = n_obs;
  x_i[3] = n_difeq;
   
  x_r[1:n_obs]                         = left_t;
  x_r[(n_obs+1):(2 * n_obs)]           = right_t;
  x_r[(2 * n_obs + 1):(2 * n_obs + A)] = age_dist;
  
  gamma     = 2.0/dI;
  tau       = 2.0/dE;
  beta0     = exp(eta0);
}

parameters {
  real x_noise[(n_changes - 1)*A];
  vector[(A * (A + 1)) / 2] L_raw;          // Vectorized version of the L matrix. Used to apply a NCP to calculate the sampled contact matrix 
}

transformed parameters{
  matrix[n_changes, A] x_trajectory;
  matrix<lower = 0>[n_obs, A] beta_trajectory;              
  real<lower = 0> beta_N[n_obs*A];
  
  real theta[A*A + A*n_obs + 4];             // Vector of ODE parameters
  real state_solutions[n_obs, A * n_difeq];  // Solution from the ODE solver
  matrix[n_obs, A] comp_C;			         // Store the calculated values for the dummy ODE compartment
  
  matrix<lower = 0>[n_obs, A] E_casesByAge;  // Expected infections per group
  matrix<lower = 0>[n_obs, A] E_deathsByAge; // Expected deaths per age group
  
  matrix[A, A] cm_sym;
  matrix[A, A] cm_sample;
  
  // Transformed parameters for the contact matrix (Non-central parameterisation):
  matrix[A, A] L_raw_mat             = to_triangular(L_raw, A);
  matrix[A, A] L                     = to_triangular(L_vector, A);
  matrix[n_changes - 1, A] x_noise_mat = to_matrix(x_noise, n_changes - 1, A);
   
  for(col in 1:A) for(row in col:A) L[row,col] = L_cm[row,col] + (p_sigmaCM * L_cm[row,col]) *  L_raw_mat[row,col];
  cm_sym    = tcrossprod(L);
  cm_sample = diag_pre_multiply(pop_diag, cm_sym);
  
  // Transformed parameters for the GBM (Non-central parameterisation):                           			
  x_trajectory[1,] = to_row_vector(x_init);
   
  for (t in 2:n_changes) for (j in 1:A) x_trajectory[t,j] = x_trajectory[t-1,j] + sigmaBM[j] * x_noise_mat[t-1,j];
   
  if (ecr_changes == 1) {
	beta_trajectory = exp(x_trajectory);

  } else {
    beta_trajectory = append_row( repeat_matrix(       exp( x_trajectory[1:(n_changes-1),] ), ecr_changes),
                                  repeat_rv_to_matrix( exp( x_trajectory[n_changes] ),        n_remainder)
                                 );
   }

  beta_N = to_vector_colwise(beta_trajectory);

  // Change of format for integrate_ode_euler/ integrate_ode_rk45/ integrate_ode_bdf:
  theta[1:(A * A)]                     = to_vector_rowwise(cm_sample);                      
  theta[A * A + 1]                     = gamma;                     
  theta[A * A + 2]                     = beta0;
  theta[(A*A + 3):(A*n_obs + A*A + 2)] = beta_N; 
  theta[A*n_obs + A*A + 3]             = pi; 
  theta[A*n_obs + A*A + 4]             = tau; 

  // Solution to the ODE system:
  state_solutions = integrate_ode_trapezoidal(init,   // initial states
											  t0,     // initial_time, 
											  ts,     // real times
											  theta,  // parameters
											  x_r,    // real data
											  x_i     // integer data
											  );
  
  // Calculate new daily Expected cases and Expected deaths for each age group:
  for (i in 1:n_obs) {
	  
	if(i == 1) E_deathsByAge[i,] = E_deathsByAge_day1;

    for (j in 1:A){

     // Format ODE results
	 comp_C[i,j] = state_solutions[i,(n_difeq-1) * A + j]  *  n_pop;

     // Alternative option:
     E_casesByAge[i,j] = comp_C[i,j] - (i == 1 ? 0 :  ( comp_C[i,j] > comp_C[i-1,j] ? comp_C[i-1,j] : 0) );

     // Expected deaths by calendar day and age group:
	 if(i != 1) E_deathsByAge[i,j] =  ifr_age[j] * dot_product(head(E_casesByAge[,j],i-1), tail(I_D_rev, i-1));
    }// End for

  }//End for
}

model {
  x_noise ~ std_normal();
  L_raw   ~ std_normal();
}

generated quantities {
	matrix[n_obs, A] y_sim;
	matrix[n_obs, A] Susceptibles; 
	
	for(t in 1:n_obs) {
		for (j in 1:A) {
			Susceptibles[t,j] = ( state_solutions[t,j] + (age_dist[j]  * (1-pi)) )* n_pop ;
			y_sim[t,j] = neg_binomial_2_rng( E_deathsByAge[t,j], E_deathsByAge[t,j]/phiD);
		}// End for
	}// End for
}
