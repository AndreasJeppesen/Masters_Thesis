functions {
  // Functions for parallelization of summation
  real partial_sum_bernoulli(array[] int y_slice, int start, int end, array[] real theta) {
    return(bernoulli_lpmf(y_slice | theta[start:end]));
  }
  real partial_sum_wiener(array[] real RT_slice, int start, int end, array[] real caut, real t, real b, array[] real EEA) {
    return(wiener_lpdf(RT_slice | caut[start:end], t, b, EEA[start:end]));
  }
}

// The input data
data {
  int<lower=1> nagents;
  array[nagents] int<lower=1> ntrials; // number of trials pr agent
  int<lower=10> maxntrials; // the maximum number of trials a participant has (in total)
  int<lower=1> nconditions; // number of conditions
  array[maxntrials,nagents] int button_pressed; // whether the button was pressed or not
  array[maxntrials,nagents] int target_present; // whether a target stimulus was present or not
  array[maxntrials,nagents] int falsealarm; // whether trial type was false alarm or random
  array[maxntrials,nagents] int correct; // Correct reaction or not
  array[maxntrials,nagents] real RT; // reaction times
  array[nagents] int ncorrectRTs;  // number of correct and incorrect RTs pr agent
  array[nagents] int nincorrectRTs;
  
  // Dummy variable denoting which condition we are in
  array[maxntrials, nconditions, nagents] int condition_dum;
  
  array[nconditions] real d;
  
  array[nagents] real min_rt;
  
  real a_prior_m;
  real a_prior_sd;

  real b_prior_m;
  real b_prior_sd;
  
  real v_prior_m;
  real v_prior_sd;
  
  real t_prior_m;
  real t_prior_sd;
  
  real FA_prior_m;
  real FA_prior_sd;
  
  real w_prior_m;
  real w_prior_sd;
  
  // real d_prior_m;
  // real d_prior_sd;
  
}

transformed data {
  // Make index vectors for where the correct and incorrect RTs are
  array[maxntrials, nagents] int correctIndex;
  array[maxntrials, nagents] int incorrectIndex;
  int position1;
  int position0;
  for (s in 1:nagents){
    position1 = 1;
    position0 = 1;
    for (i in 1:ntrials[s]){
      if (button_pressed[i,s] == 1){
        if (correct[i,s] == 1){
          correctIndex[position1,s] = i; // save this trial index, i, at position1 in the correctIndex list
          position1 = position1 + 1; // update position
        }
        if (correct[i,s] == 0){
          incorrectIndex[position0, s] = i;
          position0 = position0 + 1;
        }
      }
    } // end trial loop
  } // end agent loop
  
}

parameters {
  // parameters on real axis
  array[nagents] real log_a;
  array[nagents] real logit_b;
  array[nagents] real log_v;
  array[nagents] real cont_t;
  array[nagents] real log_FA;
  array[nagents] real log_w;
  
  //array[nconditions] real log_d;
  
}

transformed parameters {
  // Parameters on their actual scale, i.e. after transformations
  array[nagents] real<upper=10000> a;
  array[nagents] real b;
  array[nagents] real v;
  array[nagents] real t;
  array[nagents] real FA;
  array[nagents] real<upper=10000> w;
  //array[nconditions] real d;

  
  // for (c in 1:nconditions){
  //   d[c] = exp(log_d[c]);
  // }
  
  for (s in 1:nagents){
    a[s] = exp(log_a[s]); // a has lower bound 0
    b[s] = inv_logit(logit_b[s]); // b is between 0 and 1
    v[s] = exp(log_v[s]); // v has lower bound 0
    t[s] = inv_logit(cont_t[s]) * (min_rt[s] - 0.0002) + 0.0001; // We bind t between 0 and minimum RT
    FA[s] = exp(log_FA[s]);
    w[s] = exp(log_w[s]);
  }
  
  
  
  ////// EEA, Bias, Theta, and caution
  
  array[maxntrials, nagents] real bias; // "Bias" for reacting correctly; depends on presence of target stimulus
  array[maxntrials, nagents] real<lower=0, upper=1> theta; // the accuracy parameter from CDF approximation. Probability of reacting correctly at all.
  array[maxntrials, nagents] real EEA; // The Efficiency of Evidende Accumulation given by v, d, s
  array[maxntrials, nagents] real caut; // caution given by the difficulty, d, and the individual param, a.
  array[maxntrials, nagents] real EEA_given_d; // help-variable to create caut
  
  
  for (s in 1:nagents){
    for (i in 1:ntrials[s]){
      
      EEA_given_d[i,s] = v[s]; // initialize EEA_given_d
      for (c in 1:nconditions){
          // EEA_given_d[i,s] = EEA_given_d[i,s] + v[s]*(0.95^d[c])*condition_dum[i,c,s]; // include difficulty-of-condition effect
          EEA_given_d[i,s] = EEA_given_d[i,s] - 0.2 * d[c] * condition_dum[i,c,s];
      }
      if(EEA_given_d[i,s] < 0){ // If difficulty is too high, the smallest EEA can get is 0, not negative.
        EEA_given_d[i,s] = 0;
      }
      
      EEA[i,s] = EEA_given_d[i,s] - FA[s] * falsealarm[i,s]; // Include effect of FA trial if relevant
      
      // Compute caut
      // caut[i,s] = a[s] * (1 + ((v[s]-EEA_given_d[i,s])/v[s]) * w[s] );
      caut[i,s] = a[s] + w[s] * 0.2 * dot_product(d, condition_dum[i, ,s]);
      
      // Trial-wise NDT and b
      // t_pr_trial[i,s] = dot_product(t[s, ], task_dum[i, ,s]);
      // b_pr_trial[i,s] = dot_product(b[s, ], task_dum[i, ,s]);
    
    } // end trial loop
    
    // The bias, b, is the bias for pressing button.
    // Which way it affects the prob. of a correct reaction depends on whether a target stimulus is present or not.
    // The  probability of hitting upper boundary (theta) is from the RWiener package documentation.
    for (i in 1:ntrials[s]){
      if (target_present[i,s] == 1){
        bias[i,s] = b[s];
      }
      else{
        bias[i,s] = 1 - b[s]; // if target stim not present, you still have bias for pressing and hence lower prob of reacting correctly
      }
      
      if (EEA[i,s] <= 0.001 && EEA[i,s] >= -0.001){
        theta[i,s] = bias[i,s]; // if there is no direction to the drift rate, the probablity of reacting correctly is only controlled by the bias
      }
      else if (caut[i,s] <= 0.001){
        theta[i,s] = bias[i,s];
      }
  
      else {
        theta[i,s] = 1 - ( (1 - exp(-2*EEA[i,s]*caut[i,s]*(1-bias[i,s]))) / (exp(2*EEA[i,s]*caut[i,s]*bias[i,s]) - exp(-2*EEA[i,s]*caut[i,s]*(1-bias[i,s]))) );
      }
        
    } // end trial loop
    
  } // end agent loop
  
  
}

model {
  
  // Priors
  target += normal_lpdf(log_v | v_prior_m, v_prior_sd); // Drift rate prior.
  target += normal_lpdf(log_a | a_prior_m, a_prior_sd); // Boundary seperation prior.
  target += normal_lpdf(logit_b | b_prior_m, b_prior_sd); // Bias prior.
  target += normal_lpdf(cont_t | t_prior_m, t_prior_sd); // Non-decision time prior.
  target += normal_lpdf(log_FA | FA_prior_m, FA_prior_sd); // FA effect prior
  target += normal_lpdf(log_w | w_prior_m, w_prior_sd); // w prior
  
  // for (c in 1:nconditions){
  //   target += normal_lpdf(log_d[c] | d_prior_m, d_prior_sd);
  // } 
  
  int grainsize = 1;
  
  for (s in 1:nagents){
    target += reduce_sum(partial_sum_bernoulli, correct[,s],
                            grainsize,
                            theta[,s]); // where they correct in general?
                            
    target += reduce_sum(partial_sum_wiener,
                            RT[correctIndex[1:ncorrectRTs[s], s],s],
                            grainsize,
                            caut[correctIndex[1:ncorrectRTs[s], s],s],
                            t[s], b[s], EEA[correctIndex[1:ncorrectRTs[s], s],s]);
                            
    target += reduce_sum(partial_sum_wiener,
                            RT[incorrectIndex[1:nincorrectRTs[s], s], s],
                            grainsize,
                            caut[incorrectIndex[1:nincorrectRTs[s], s],s],
                            t[s], b[s], -EEA[incorrectIndex[1:nincorrectRTs[s], s],s]);
  }
  
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
  array[nagents] real t_prior;
  real log_FA_prior;
  real FA_prior;
  real log_w_prior;
  real w_prior;
  // array[nconditions] real log_d_prior;
  // array[nconditions] real d_prior;
  
  log_v_prior = normal_rng(v_prior_m, v_prior_sd);
  log_a_prior = normal_rng(a_prior_m, a_prior_sd);
  logit_b_prior = normal_rng(b_prior_m, b_prior_sd);
  cont_t_prior = normal_rng(t_prior_m, t_prior_sd);
  log_FA_prior = normal_rng(FA_prior_m, FA_prior_sd);
  log_w_prior = normal_rng(w_prior_m, w_prior_sd);
  
  // for (c in 1:nconditions){
  //   log_d_prior[c] = normal_rng(d_prior_m, d_prior_sd);
  // }
  
  // Transformations:
  v_prior = exp(log_v_prior);
  a_prior = exp(log_a_prior);
  b_prior = inv_logit(logit_b_prior);
  FA_prior = exp(log_FA_prior);
  w_prior = exp(log_w_prior);
  // for (c in 1:nconditions){
  //   d_prior[c] = exp(log_d_prior[c]);
  // }
  for (s in 1:nagents){
    t_prior[s] = inv_logit(cont_t_prior) * (min_rt[s] - 0.0002) + 0.0001;
  }
  
  // Prior & Posterior predictions of Hit Rate
  array[maxntrials] real EEA_prior;
  array[maxntrials] real EEA_given_d_prior;
  array[maxntrials] real caut_prior;
  array[maxntrials] real Hit_Rate_prior;
  // array[maxntrials] real Hit_Rate_posterior;
  
  for (i in 1:ntrials[1]){ // Use ntrials for agent 1
    EEA_given_d_prior[i] = 0; // initialize EEA_given_d_prior
    for (c in 1:nconditions){
      // EEA_given_d_prior[i] = EEA_given_d_prior[i] + v_prior*(0.95^d_prior[c])*condition_dum[i,c,1]; // include difficulty-of-condition effect, use dummy from agent 1
      EEA_given_d_prior[i] = EEA_given_d_prior[i] - 0.2 * d[c] * condition_dum[i,c,1];
      
    }
    
    EEA_prior[i] = EEA_given_d_prior[i] - FA_prior * falsealarm[i,1]; // Include effect of FA trial if relevant
    
    // Compute caut
    // caut_prior[i] = a_prior * (1 + ((v_prior-EEA_given_d_prior[i])/v_prior) * w_prior);
    caut_prior[i] = a_prior + w_prior * 0.2 * dot_product(d, condition_dum[i, ,1]);
    
    // Hit Rate will correspond to the probability of reacting correctly to a target stimulus, i.e. bias = b:
    Hit_Rate_prior[i] = 1 - ( (1 - exp(-2*EEA_prior[i]*caut_prior[i]*(1-b_prior))) / (exp(2*EEA_prior[i]*caut_prior[i]*b_prior) - exp(-2*EEA_prior[i]*caut_prior[i]*(1-b_prior))) );
    
    // Hit Rate posterior will vary between participants and can be seen on a theta[i] posterior for a trial with a target
    //Hit_Rate_posterior[i] = 1 - ( (1 - exp(-2*EEA[i]*caut[i]*(1-b))) / (exp(2*EEA[i]*caut[i]*b) - exp(-2*EEA[i]*caut[i]*(1-b))) );
  } // end trial loop
  
  
}