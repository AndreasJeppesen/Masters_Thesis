functions{
  // Number of unique elements in array:
  // int n_unique(array[] int x){
  //   int counter;
  //   array[num_elements(x)] int sorted = sort_asc(x); // sort x, ascending order
  //   for (i in 1:num_elements(x)){
  //     if (i == 1){
  //       counter = counter + 1; // first element is always new
  //     }
  //     else if (sorted[i] != sorted[i-1]){ // if this element is not identical to previous
  //       counter = counter + 1;       // then it is new
  //     }
  //   }
  //   return (counter);
  // }
  // RNG for truncated normal distribution:
  real normal_lb_rng(real mu, real sigma, real lb) {
    real p = normal_cdf(lb | mu, sigma);  // cdf for bounds
    real u = uniform_rng(p, 1);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf for value
  }
}

// The input data
data {
  int<lower=1> ntrials; // number of trials
  
  array[ntrials] int button_pressed; // whether the button was pressed or not
  array[ntrials] int target_present; // whether a target stimulus was present or not
  array[ntrials] int correct; // Correct reaction or not
  array[ntrials] real RT; // reaction times
  int ncorrectRTs;
  array[ntrials] int falsealarm; // whether trial type was false alarm or random
  
  real min_rt;
  
  
  real a_prior_m;
  real a_prior_sd;

  real b_prior_m;
  real b_prior_sd;
  
  real v_prior_m;
  real v_prior_sd;
  
  real t_prior_m;
  real t_prior_sd;
  
  real FA_prior_m;  // parameter FA: fixed effect of "false alarm" trial vs random trials.
  real FA_prior_sd;
  
}

transformed data {
  // Make index vectors for where the correct and incorrect RTs are
  array[ntrials] int correctIndex;
  int position1;
  position1 = 1;
  for (i in 1:ntrials){
    if (button_pressed[i] == 1){
      if (correct[i] == 1){
      correctIndex[position1] = i; // save this trial index, i, at position1 in the correctIndex list
      position1 = position1 + 1; // update position
      }
    }
  } // end trial loop
  
  
  
  
}


parameters {
  // Individual level parameters on real axis
  real log_a;
  real logit_b;
  real log_v;
  real cont_t;
  real log_FA;
}

transformed parameters {
  // Transformed parameters on their actual scales
  real a;
  real b;
  real v;
  real t;
  real FA;
  
  array[ntrials] real bias; // "Bias" for reacting correctly; depends on presence of target stimulus
  array[ntrials] real<lower=0, upper=1> theta; // the accuracy parameter from CDF approximation. Probability of reacting correctly at all.
  array[ntrials] real EEA; // The Efficiency of Evidende Accumulation given by v, d, s, and item_effect
  
  

  a = exp(log_a);
  b = inv_logit(logit_b);
  v = exp(log_v);
  t = inv_logit(cont_t) * (min_rt - 0.0002) + 0.0001;
  FA = exp(log_FA);
  
  // EEA, Bias and Theta
  
    
  // EEA is calculated based on v, the difficulty (by cond.), the false-alarm effect, and the random by-item effect.
  for (i in 1:ntrials){
    EEA[i] = v - FA*falsealarm[i]; // add effect of FA 
  } // end trial loop
  
  
  // The bias, b, is the bias for pressing button.
  // Which way it affects the prob. of a correct reaction depends on whether a target stimulus is present or not.
  // The  probability of hitting upper boundary (theta) is from the RWiener package documentation.
  for (i in 1:ntrials){
    if (target_present[i] == 1){
      bias[i] = b;
    }
    else{
      bias[i] = 1-b; // if target stim not present, you still have bias for pressing and hence lower prob of reacting correctly
    }
    
    if (EEA[i] == 0){
      theta[i] = bias[i]; // if there is no direction to the drift rate, the probablity of reacting correctly is only controlled by the bias
    }

    if (EEA[i] != 0){
      theta[i] = 1 - ( (1 - exp(-2*EEA[i]*a*(1-bias[i]))) / (exp(2*EEA[i]*a*bias[i]) - exp(-2*EEA[i]*a*(1-bias[i]))) );
    }
      
  } // end trial loop
    

}

// The model to be estimated
model {
  // Priors
  target += normal_lpdf(log_v | v_prior_m, v_prior_sd); // Drift rate prior. Uninformed.
  target += normal_lpdf(log_a | a_prior_m, a_prior_sd); // Boundary seperation prior. Maybe needs to have more probability mass at higher values?
  target += normal_lpdf(logit_b | b_prior_m, b_prior_sd); // Bias prior. Uninformed.
  target += normal_lpdf(cont_t | t_prior_m, t_prior_sd); // Non-decision time prior. Most probability mass between 0 and 1 seconds.
  target += normal_lpdf(log_FA | FA_prior_m, FA_prior_sd);
  
  // Likelihood
  target += bernoulli_lpmf(correct | theta); // where they correct in general?
  target += wiener_lpdf(RT[correctIndex[1:ncorrectRTs]] | a, t, b, EEA[correctIndex[1:ncorrectRTs]]); // RT on correct trials 
  // No incorrect RTs in the RVP
  
  
}

generated quantities {
  // Priors
  real log_v_prior;
  real v_prior;
  real log_a_prior;
  real a_prior;
  real logit_b_prior;
  real b_prior;
  real cont_t_prior;
  real t_prior;
  real log_FA_prior;
  real FA_prior;
  
  log_v_prior = normal_rng(v_prior_m, v_prior_sd);
  log_a_prior = normal_rng(a_prior_m, a_prior_sd);
  logit_b_prior = normal_rng(b_prior_m, b_prior_sd);
  cont_t_prior = normal_rng(t_prior_m, t_prior_sd);
  log_FA_prior = normal_rng(FA_prior_m, FA_prior_sd);
  
  // Transformations:
  v_prior = exp(log_v_prior);
  a_prior = exp(log_a_prior);
  b_prior = inv_logit(logit_b_prior);
  t_prior = inv_logit(cont_t_prior) * (min_rt - 0.0002) + 0.0001;
  FA_prior = exp(log_FA_prior);
  
  // Prior & Posterior predictions of Hit Rate
  array[ntrials] real EEA_prior;
  array[ntrials] real Hit_Rate_prior;
  array[ntrials] real Hit_Rate_posterior;
  
  for (i in 1:ntrials){
    EEA_prior[i] = v_prior - FA_prior * falsealarm[i];
    // Hit Rate will correspond to the probability of reacting correctly to a target stimulus, i.e.:
    Hit_Rate_prior[i] = 1 - ( (1 - exp(-2*EEA_prior[i]*a_prior*(1-b_prior))) / (exp(2*EEA_prior[i]*a_prior*b_prior) - exp(-2*EEA_prior[i]*a_prior*(1-b_prior))) );
    Hit_Rate_posterior[i] = 1 - ( (1 - exp(-2*EEA[i]*a*(1-b))) / (exp(2*EEA[i]*a*b) - exp(-2*EEA[i]*a*(1-b))) );
  } // end trial loop
  
}

