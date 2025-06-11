// Task general diffusion model


functions {
  // Functions for parallelization of summation
  real partial_sum_bernoulli(array[] int y_slice, int start, int end, array[] real theta) {
    return(bernoulli_lpmf(y_slice | theta[start:end]));
  }
  real partial_sum_wiener(array[] real RT_slice, int start, int end, array[] real caut, array[] real t, array[] real b, array[] real EEA) {
    return(wiener_lpdf(RT_slice | caut[start:end], t[start:end], b[start:end], EEA[start:end]));
  }
  
  // RNG for truncated normal distribution:
  real normal_lb_rng(real mu, real sigma, real lb) {
    real p = normal_cdf(lb | mu, sigma);  // cdf for bounds
    real u = uniform_rng(p, 1);
    return (sigma * inv_Phi(u)) + mu;  // inverse cdf for value
  }
}


data {
  
  int<lower=1> nagents; // number of agents
  array[nagents] int<lower=1> ntrials; // number of trials pr agent
  int<lower=1> nconditions; // number of conditions (in total)
  int<lower=1> ntasks; // number of different tasks
  int<lower=10> maxntrials; // the maximum number of trials a participant has (in total)
  
  array[maxntrials,nagents] int button_pressed; // whether the button was pressed or not
  array[maxntrials,nagents] int target_present; // whether a target stimulus was present or not
  array[maxntrials,nagents] int falsealarm; // whether trial type was false alarm or random
  array[maxntrials,nagents] int correct; // Correct reaction or not
  array[maxntrials,nagents] real RT; // reaction times
  array[nagents] int ncorrectRTs;  // number of correct and incorrect RTs pr agent
  array[nagents] int nincorrectRTs;
  
  // Dummy variable denoting which condition we are in
  array[maxntrials, nconditions, nagents] int condition_dum;
  
  array[nconditions] real d; // The difficulties of the conditions
  
  // Dummy variable denoting which task we are in
  array[maxntrials, ntasks, nagents] int task_dum;
  
  // Dummy variable denoting if incorrect RT data is present in this trial (in RVP task it is not logged)
  array[maxntrials, nagents] int incorrectRTs_present;
  
  array[nagents, ntasks] real min_rt; // Minimum RT observed for each agent on each task
  
  matrix[ntasks,ntasks] A; // covar matrix for the correlated params
  
  
  // Prior values
  real a_M_prior_m;
  real a_M_prior_sd;
  // real a_SD_prior_m;
  // real a_SD_prior_sd;
  real b_M_prior_m;
  real b_M_prior_sd;
  // real b_SD_prior_m;
  // real b_SD_prior_sd;
  real v_M_prior_m;
  real v_M_prior_sd;
  // real v_SD_prior_m;
  // real v_SD_prior_sd;
  real t_M_prior_m;
  real t_M_prior_sd;
  // real t_SD_prior_m;
  // real t_SD_prior_sd;
  
  // real d_prior_m;
  // real d_prior_sd;
  real w_prior_m;
  real w_prior_sd;
  real FA_prior_m;
  
  real a_SD_input;
  real b_SD_input;
  real v_SD_input;
  real t_SD_input;
  
}

transformed data{
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
          if (incorrectRTs_present[i,s] == 1){
            incorrectIndex[position0, s] = i;
            position0 = position0 + 1;
          }
        }
      }
    } // end trial loop
  } // end agent loop
}

parameters {
  // Individual lvl parameters
  array[nagents] real log_a_M;
  // array[nagents] real<lower=0> log_a_SD;
  array[nagents] real logit_b_M;
  // array[nagents] real<lower=0> logit_b_SD;
  array[nagents] real log_v_M;
  // array[nagents] real<lower=0> log_v_SD;
  array[nagents] real log_t_M;
  // array[nagents] real<lower=0> t_SD;
  
  array[nagents] real log_w; // does not vary between tasks - partially controls how much a varies between tasks
  
  
  // Task- and individual specific params
  array[nagents, ntasks] real log_FA;
  array[nagents, ntasks] real log_a;
  array[nagents, ntasks] real logit_b;
  array[nagents, ntasks] real log_v;
  array[nagents, ntasks] real cont_t;
  
  // Condition specific param
  // array[nconditions] real log_d;
  
}

transformed parameters {
  
  /////// Reparameterizations:
  
  // Individual means on actual scale
  array[nagents] real a_M;
  array[nagents] real b_M;
  array[nagents] real v_M;
  array[nagents] real t_M;
  
  array[nagents] real<lower=0> log_a_SD;
  array[nagents] real<lower=0> logit_b_SD;
  array[nagents] real<lower=0> log_v_SD;
  array[nagents] real<lower=0> t_SD;
  
  
  // Individual values for specific task on actual scale
  array[nagents, ntasks] real a;
  array[nagents, ntasks] real b;
  array[nagents, ntasks] real v;
  array[nagents, ntasks] real t;
  array[nagents, ntasks] real FA;
  array[nagents] real w;
  
  // array[nconditions] real d;
  
  // Non-decision time, individual means. Now bounded below by zero
  for (s in 1:nagents){
    t_M[s] = exp(log_t_M[s]);
    w[s] = exp(log_w[s]);
  }
  
  
  for (s in 1:nagents){
    a_M[s] = exp(log_a_M[s]);
    b_M[s] = inv_logit(logit_b_M[s]);
    v_M[s] = exp(log_v_M[s]);
    
    // Transform input SD's to what they should be on the sampling scales
    log_a_SD[s] = a_SD_input / a_M[s];
    logit_b_SD[s] = b_SD_input;
    log_v_SD[s] = v_SD_input / v_M[s];
    t_SD[s] = t_SD_input;
    
    
    for (task in 1:ntasks){
      a[s,task] = exp(log_a[s,task]);
      b[s,task] = inv_logit(logit_b[s,task]);
      v[s,task] = exp(log_v[s,task]);
      t[s,task] = inv_logit(cont_t[s,task]) * (min_rt[s,task] - 0.0002) + 0.0001; // We bind t on specific task between 0 and minimum RT for that [agent,task]
      FA[s,task] = exp(log_FA[s,task]);
    }
  }
  
  // for (c in 1:nconditions){
  //   d[c] = exp(log_d[c]);
  // }
 
  
  ////// EEA, Bias, Theta, and caution
  array[maxntrials, nagents] real bias; // "Bias" for reacting correctly; depends on presence of target stimulus
  array[maxntrials, nagents] real<lower=0, upper=1> theta; // the accuracy parameter from CDF approximation. Probability of reacting correctly at all.
  array[maxntrials, nagents] real EEA; // The Efficiency of Evidende Accumulation given by v, d, s
  array[maxntrials, nagents] real caut; // caution given by the difficulty, d, and the individual param, a.
  array[maxntrials, nagents] real EEA_given_d; // help-variable to create caut
  array[maxntrials, nagents] real t_pr_trial; // trial-wise NDT
  array[maxntrials, nagents] real b_pr_trial; // trial-wise b
  
  for (s in 1:nagents){
    for (i in 1:ntrials[s]){
      
      EEA_given_d[i,s] = 0; // initialize EEA_given_d
      caut[i,s] = 0; // init caut
      for (c in 1:nconditions){
        for (task in 1:ntasks){
          // EEA_given_d[i,s] = EEA_given_d[i,s] + v[s,task]*task_dum[i,task,s]*(0.95^d[c])*condition_dum[i,c,s]; // include difficulty-of-condition effect
          EEA_given_d[i,s] = EEA_given_d[i,s] + (v[s,task] - 0.2 * d[c]) * task_dum[i,task,s] * condition_dum[i,c,s];
          caut[i,s] = caut[i,s] + (a[s,task] + w[s] * 0.2 * d[c]) * task_dum[i,task,s] * condition_dum[i,c,s];
        }
      }
      if (EEA_given_d[i,s] < 0){
        EEA_given_d[i,s] = 0;
      }
      
      EEA[i,s] = EEA_given_d[i,s]; // init EEA at EEA_given_d
      for (task in 1:ntasks){
        EEA[i,s] = EEA[i,s] - FA[s,task]*task_dum[i,task,s]*falsealarm[i,s]; // add only when task_dum & false_alarm == 1
      }
      
      
      // Trial-wise NDT and b
      t_pr_trial[i,s] = dot_product(t[s, ], task_dum[i, ,s]);
      b_pr_trial[i,s] = dot_product(b[s, ], task_dum[i, ,s]);
    
    } // end trial loop
    
    // The bias, b, is the bias for pressing button.
    // Which way it affects the prob. of a correct reaction depends on whether a target stimulus is present or not.
    // The  probability of hitting upper boundary (theta) is from the RWiener package documentation.
    for (i in 1:ntrials[s]){
      if (target_present[i,s] == 1){
        bias[i,s] = dot_product(b[s, ], task_dum[i, ,s]);
      }
      else{
        bias[i,s] = 1 - dot_product(b[s, ], task_dum[i, ,s]); // if target stim not present, you still have bias for pressing and hence lower prob of reacting correctly
      }
      
      if (EEA[i,s] <= 0.001 && EEA[i,s] >= - 0.001){
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
  ///// Priors
  
  // Individual means and fluctutations (SDs)
  for (s in 1:nagents){
    target += normal_lpdf(log_a_M[s] | a_M_prior_m, a_M_prior_sd);
    // target += normal_lpdf(log_a_SD[s] | a_SD_prior_m, a_SD_prior_sd) - 
    //   normal_lccdf(0 | a_SD_prior_m, a_SD_prior_sd);
      
    target += normal_lpdf(logit_b_M[s] | b_M_prior_m, b_M_prior_sd);
    // target += normal_lpdf(logit_b_SD[s] | b_SD_prior_m, b_SD_prior_sd) - 
    //   normal_lccdf(0 | b_SD_prior_m, b_SD_prior_sd);
      
    target += normal_lpdf(log_v_M[s] | v_M_prior_m, v_M_prior_sd);
    // target += normal_lpdf(log_v_SD[s] | v_SD_prior_m, v_SD_prior_sd) - 
    //   normal_lccdf(0 | v_SD_prior_m, v_SD_prior_sd);
    
    target += normal_lpdf(log_t_M[s] | t_M_prior_m, t_M_prior_sd);
    // target += normal_lpdf(t_SD[s] | t_SD_prior_m, t_SD_prior_sd) - 
    //   normal_lccdf(0 | t_SD_prior_m, t_SD_prior_sd);
    
    target += normal_lpdf(log_w[s] | w_prior_m, w_prior_sd);
  }
  
  
  // Difficulty of conditions
  // for (c in 1:nconditions){
  //   target += normal_lpdf(log_d[c] | d_prior_m, d_prior_sd);
  // }
  
  // Task-specific FA-effect; joint prior to correlate them across tasks within each agent
  for (s in 1:nagents){
    target += multi_normal_cholesky_lpdf(to_vector(log_FA[s, ]) | rep_vector(FA_prior_m, ntasks), A);
  }
  
  
  ///// Draw task-specific parameters for each agent
  for (task in 1:ntasks){
    target += normal_lpdf(log_a[ ,task] | log_a_M, log_a_SD); // Boundary sep, a
    target += normal_lpdf(logit_b[ ,task] | logit_b_M, logit_b_SD); // Bias, b
    target += normal_lpdf(log_v[ ,task] | log_v_M, log_v_SD); // drift rate v
    target += normal_lpdf(t[ ,task] | t_M, t_SD); // NDT, t
  }
  
  ///// Likelihood
  int grainsize = 1;
  
  for (s in 1:nagents){
    // target += bernoulli_lpmf(correct[,s] | theta[,s]); // where they correct in general?
    // target += wiener_lpdf(RT[correctIndex[1:ncorrectRTs[s], s],s] |
    //     caut[correctIndex[1:ncorrectRTs[s], s],s],
    //     t_pr_trial[correctIndex[1:ncorrectRTs[s], s],s],
    //     b_pr_trial[correctIndex[1:ncorrectRTs[s], s],s],
    //     EEA[correctIndex[1:ncorrectRTs[s], s],s]); // RT on correct trials
    // target += wiener_lpdf(RT[incorrectIndex[1:nincorrectRTs[s], s],s] |
    //     caut[incorrectIndex[1:nincorrectRTs[s], s],s],
    //     t_pr_trial[incorrectIndex[1:nincorrectRTs[s], s],s],
    //     b_pr_trial[incorrectIndex[1:nincorrectRTs[s], s],s],
    //     -EEA[incorrectIndex[1:nincorrectRTs[s], s],s]); // RT on incorrect trials
        
    // With reduce sum
    target += reduce_sum(partial_sum_bernoulli, correct[,s], grainsize, theta[,s]);
    target += reduce_sum(partial_sum_wiener,
                            RT[correctIndex[1:ncorrectRTs[s], s],s],
                            grainsize,
                            caut[correctIndex[1:ncorrectRTs[s], s],s],
                            t_pr_trial[correctIndex[1:ncorrectRTs[s], s],s],
                            b_pr_trial[correctIndex[1:ncorrectRTs[s], s],s],
                            EEA[correctIndex[1:ncorrectRTs[s], s],s]
                            );
    target += reduce_sum(partial_sum_wiener,
                            RT[incorrectIndex[1:nincorrectRTs[s], s],s],
                            grainsize,
                            caut[incorrectIndex[1:nincorrectRTs[s], s],s],
                            t_pr_trial[incorrectIndex[1:nincorrectRTs[s], s],s],
                            b_pr_trial[incorrectIndex[1:nincorrectRTs[s], s],s],
                            -EEA[incorrectIndex[1:nincorrectRTs[s], s],s]
                            );
                                
  }
  
}


generated quantities {
  real log_a_M_prior;
  // real<lower=0> log_a_SD_prior;
  real logit_b_M_prior;
  // real<lower=0> logit_b_SD_prior;
  real log_v_M_prior;
  // real<lower=0> log_v_SD_prior;
  real log_t_M_prior;
  // real<lower=0> t_SD_prior;
  // array[nconditions] real log_d_prior; // difficulties of conditions
  vector[ntasks] log_FA_prior; // effect of "false alarm"
  real log_w_prior; // effect of difficulty on caution
  
  // Generate on sampling scales
  log_a_M_prior = normal_rng(a_M_prior_m, a_M_prior_sd);
  // log_a_SD_prior = normal_lb_rng(a_SD_prior_m, a_SD_prior_sd, 0);
  logit_b_M_prior = normal_rng(b_M_prior_m, b_M_prior_sd);
  // logit_b_SD_prior = normal_lb_rng(b_SD_prior_m, b_SD_prior_sd, 0);
  log_v_M_prior = normal_rng(v_M_prior_m, v_M_prior_sd);
  // log_v_SD_prior = normal_lb_rng(v_SD_prior_m, v_SD_prior_sd, 0);
  log_t_M_prior = normal_rng(t_M_prior_m, t_M_prior_sd);
  // t_SD_prior = normal_lb_rng(t_SD_prior_m, t_SD_prior_sd, 0);
  // for (c in 1:nconditions){
  //   log_d_prior[c] = normal_rng(d_prior_m, d_prior_sd);
  // }
  log_FA_prior = multi_normal_cholesky_rng(rep_vector(FA_prior_m, ntasks), A);
  log_w_prior = normal_rng(w_prior_m, w_prior_sd);
  
  // Transform means to actual scales
  real a_M_prior;
  real b_M_prior;
  real v_M_prior;
  real t_M_prior;
  vector[ntasks] FA_prior;
  real w_prior;
  // array[nconditions] real d_prior; // difficulties of conditions
  
  a_M_prior = exp(log_a_M_prior);
  b_M_prior = inv_logit(logit_b_M_prior);
  v_M_prior = exp(log_v_M_prior);
  t_M_prior = exp(log_t_M_prior);
  w_prior = exp(log_w_prior);
  for (task in 1:ntasks){
    FA_prior[task] = exp(log_FA_prior[task]);
  }
  // for (c in 1:nconditions){
  //   d_prior[c] = exp(log_d_prior[c]);
  // }
}
