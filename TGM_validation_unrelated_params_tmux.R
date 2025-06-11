setwd("/work/AndreasHjorthJeppesen#0476/R")



source("Setup_and_init/activate_renv.R")

library(pacman)
pacman::p_load(tidyverse, cmdstanr, posterior, ggridges, MASS)







nagents <- 30
ntrials_pr_task <- c(450, 360, 300) # CPT-IP, Conners', RVP
nconditions_pr_task <- c(3, 1, 1)
ntargets_pr_task <- c(90, 324, 24)
nFAtrials_RVP <- 5
CPT_IP_stimuli <- read.csv("CPT_IP_stimuli.csv")
ntasks <- 3



corr_matrix <- matrix(c(1, 0.9, 0.9,
                        0.9, 1, 0.9,
                        0.9, 0.9, 1), nrow = 3) # Corr matrix for the model to use
variances <- diag(1, nrow = 3, ncol = 3) # In this case, Identity matrix
covar_matrix <- variances %*% corr_matrix %*% variances


d <- array(c(0.5,2,4,
             0.1,NA,NA,
             1,NA,NA), dim = c(ntasks, max(nconditions_pr_task)))
d_vec <- array(d, dim = 9)
d_vec <- d_vec[!is.na(d_vec)]


datetime <- "2025-04-02 105828.565123"

sim_dat <- read.csv(paste("TGM_unr_params_large_SD_sim_dat", paste0(datetime, ".csv")))




file <- file.path("task_general_diffusion_mod_linear_caut_fixed_d_and_SD_partialsum.stan")

TGM_mod <- cmdstan_model(file, cpp_options = list(stan_threads = TRUE))



a_M_prior_m <- 0
a_M_prior_sd <- 1
a_SD_prior_m <- 0
a_SD_prior_sd <- 0.3
b_M_prior_m <- 0
b_M_prior_sd <- 1
b_SD_prior_m <- 0
b_SD_prior_sd <- 0.3
v_M_prior_m <- 0
v_M_prior_sd <- 1
v_SD_prior_m <- 0
v_SD_prior_sd <- 0.3
t_M_prior_m <- 0
t_M_prior_sd <- 1
t_SD_prior_m <- 0
t_SD_prior_sd <- 0.2

d_prior_m <- 0
d_prior_sd <- 1
w_prior_m <- 0
w_prior_sd <- 1
FA_prior_m <- 0




ntrials <- array(NA, dim = nagents)
for (s in 1:nagents){
  ntrials[s] <- length(sim_dat$trial[sim_dat$ID == s]) # number of trials pr participant
}

maxntrials <- max(ntrials) # Maximal number of trials a participant has completed
nconditions <- sum(nconditions_pr_task) # Total number of conditions

# Empty arrays to fill
target_present <- array(0, dim = c(maxntrials, nagents))
button_pressed <- array(0, dim = c(maxntrials, nagents))
correct <- array(0, dim = c(maxntrials, nagents))
RT <- array(0, dim = c(maxntrials, nagents))
falsealarm <- array(0, dim = c(maxntrials, nagents))
incorrectRTs_present <- array(0, dim = c(maxntrials, nagents))
ncorrectRTs <- array(0, dim = nagents)
nincorrectRTs <- array(0, dim = nagents)
condition_dum <- array(0, dim = c(maxntrials, nconditions, nagents))
task_dum <- array(0, dim = c(maxntrials, ntasks, nagents))
min_rt <- array(0, dim = c(nagents, ntasks))



for (s in 1:nagents){
  target_present[1:ntrials[s],s] <- sim_dat$target_present[sim_dat$ID == s]
  button_pressed[1:ntrials[s],s] <- sim_dat$button_pressed[sim_dat$ID == s]
  correct[1:ntrials[s],s] <- sim_dat$correct[sim_dat$ID == s]
  RT[1:ntrials[s],s] <- sim_dat$RT[sim_dat$ID == s]
  falsealarm[1:ntrials[s],s] <- sim_dat$FA_trial[sim_dat$ID == s]
  incorrectRTs_present[1:ntrials[s],s] <- sim_dat$incorrectRTs_present[sim_dat$ID == s]
  ncorrectRTs[s] <- sum(sim_dat$button_pressed[sim_dat$correct == 1 &
                                                 sim_dat$RT > 0 &
                                                 sim_dat$ID == s])
  nincorrectRTs[s] <- sum(sim_dat$button_pressed[sim_dat$correct == 0 &
                                                   sim_dat$RT > 0 &
                                                   sim_dat$ID == s])
  
  for (i in 1:ntrials[s]){
    condition_dum[i, sim_dat$condition[sim_dat$ID == s][i], s] <- 1
    task_dum[i, sim_dat$task_num[sim_dat$ID == s][i], s] <- 1
  }
  for (task in 1:ntasks){
    min_rt[s,task] <- min(sim_dat$RT[sim_dat$RT > 0 & sim_dat$ID == s & sim_dat$task_num == task])
  }
}

A_T <- chol(covar_matrix) # Cholesky decomposition

A <- t(A_T)

iter_warmup = 1000
iter_sampling = 3000
chains = 4




infer_params <- data.frame()

for (s in 17:21){
  print(paste("ITERATION", s))
  
  data <- list(nagents = 1,
               ntrials = ntrials[s],
               maxntrials = maxntrials,
               nconditions = nconditions,
               ntasks = ntasks,
               target_present = array(target_present[,s], dim = c(maxntrials, 1)),
               button_pressed = array(button_pressed[,s], dim = c(maxntrials, 1)),
               correct = array(correct[,s], dim = c(maxntrials, 1)),
               RT = array(RT[,s], dim = c(maxntrials, 1)),
               falsealarm = array(falsealarm[,s], dim = c(maxntrials, 1)),
               incorrectRTs_present = array(incorrectRTs_present[,s], dim = c(maxntrials, 1)),
               ncorrectRTs = ncorrectRTs[s],
               nincorrectRTs = nincorrectRTs[s],
               condition_dum = array(condition_dum[ , ,s], dim = c(maxntrials, nconditions, 1)),
               task_dum = array(task_dum[ , ,s], dim = c(maxntrials, ntasks, 1)),
               min_rt = array(min_rt[s,], dim = c(1, ntasks)),
               A = A,
               d = d_vec,
               
               a_M_prior_m = a_M_prior_m,
               a_M_prior_sd = a_M_prior_sd,
               # a_SD_prior_m = a_SD_prior_m,
               # a_SD_prior_sd = a_SD_prior_sd,
               b_M_prior_m = b_M_prior_m,
               b_M_prior_sd = b_M_prior_sd,
               # b_SD_prior_m = b_SD_prior_m,
               # b_SD_prior_sd = b_SD_prior_sd,
               v_M_prior_m = v_M_prior_m,
               v_M_prior_sd = v_M_prior_sd,
               # v_SD_prior_m = v_SD_prior_m,
               # v_SD_prior_sd = v_SD_prior_sd,
               t_M_prior_m = t_M_prior_m,
               t_M_prior_sd = t_M_prior_sd,
               # t_SD_prior_m = t_SD_prior_m,
               # t_SD_prior_sd = t_SD_prior_sd,
               
               # d_prior_m = d_prior_m,
               # d_prior_sd = d_prior_sd,
               w_prior_m = w_prior_m,
               w_prior_sd = w_prior_sd,
               FA_prior_m = FA_prior_m,
               
               a_SD_input = .06,
               b_SD_input = .2,
               v_SD_input = .06,
               t_SD_input = .006
  )
  
  
  
  #-- Fit model
  
  temp_samples <- TGM_mod$sample(
    data = data,
    seed = 123,
    chains = chains,
    parallel_chains = chains,
    threads_per_chain = 4,
    init = 0,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 10,
    max_treedepth = 10,
    adapt_delta = 0.99
  )
  
  temp_draws_df <- as_draws_df(temp_samples$draws())
  
  temp_trim_df <- temp_draws_df %>%
    dplyr::select(!starts_with("EEA")) %>% 
    dplyr::select(!starts_with("theta")) %>% 
    dplyr::select(!starts_with("bias")) %>% 
    dplyr::select(!starts_with("caut")) %>% 
    dplyr::select(!starts_with("t_pr_trial")) %>% 
    dplyr::select(!starts_with("b_pr_trial"))
  
  temp_trim_df$ID <- s
  
  write.csv(temp_trim_df, sprintf("draws_trim_TGM_unrelated_params_fixed_d_and_large_SD_agent_%s_%s.csv", s, datetime))
  
  infer_params <- rbind(infer_params, temp_trim_df)
}


# write.csv(infer_params, paste("draws_trim_TGM_unrelated_params_fixed_d_and_SD", paste0(datetime, ".csv")))


