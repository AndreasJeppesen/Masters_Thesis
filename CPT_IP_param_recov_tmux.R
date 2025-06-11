setwd("/work/AndreasHjorthJeppesen#0476/R")



source("Setup_and_init/activate_renv.R")


library(pacman)
pacman::p_load(tidyverse, cmdstanr, posterior, ggridges)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}




source("CPT_IP_w_diff_and_linear_caut_wiener_sim_f.R")
ntrials <- 450
stimuli <- read.csv("CPT_IP_stimuli.csv")
nconditions <- 3
nagents = 30


true_params <- read.csv("TGM_fixed_d_and_SD_true_params 2025-03-10 125332.675184.csv")
ntrials_pr_task <- c(450, 360, 300) # CPT-IP, Conners', RVP
nconditions_pr_task <- c(3, 1, 1)
ntargets_pr_task <- c(90, 324, 24)
nFAtrials_RVP <- 5
ntasks <- 3
d <- array(c(0.5,2,4,
             0.1,NA,NA,
             1,NA,NA), dim = c(ntasks, max(nconditions_pr_task)))
CPT_d <- d[,1]

sim_dat <- read.csv("TGM_fixed_d_and_SD_sim_dat 2025-03-10 125332.675184.csv")
sim_dat <- sim_dat %>% filter(task_name == "CPT-IP")


file <- file.path("CPT_IP_diffusion_model_partialsum(2).stan")

wiener_CPT_IP_mod = cmdstan_model(file, cpp_options = list(stan_threads = TRUE))




# Set parameter values
# a <- runif(nagents, 1, 1.7)
# b <- runif(nagents, 0.3, 0.7)
# v <- runif(nagents, 0.8, 1.7)
# t <- runif(nagents, 0.2, 0.25)
# d <- c(0.1, 2, 4)
# FA <- runif(nagents, -1.5, -0.5)
# w <- runif(nagents, 0, 2)

# Simulate data
# sim_dat <- data.frame()
# for (s in 1:nagents){
#   temp_sim <- CPT_IP_w_diff_and_caut_wiener_sim_f(item_type = stimuli$stim_type,
#                                                   v[s], a[s], b[s], t[s], d, FA[s], w[s])
#   temp_sim$ID = s
#   sim_dat <- rbind(sim_dat, temp_sim)
# }
# 
# # Save true param values to csv
# true_params <- data.frame(s = seq(1:nagents), a, b, v, t, w, FA)

datetime <- as.character(Sys.time())
datetime <- gsub(":", "", datetime)

#write.csv(true_params, paste("CPT_IP_true_params_fixed_d", paste0(datetime, ".csv")))



a_prior_m <- 0
a_prior_sd <- 1
b_prior_m <- 0
b_prior_sd <- 1
v_prior_m <- 0
v_prior_sd <- 1
t_prior_m <- 0
t_prior_sd <- 1
d_prior_m <- 0
d_prior_sd <- 1
FA_prior_m <- 0
FA_prior_sd <- 1
w_prior_m <- 0
w_prior_sd <- 1



# Prepare data for Stan
ntrials <- array(NA, dim = nagents)
for (s in 1:nagents){
  ntrials[s] <- length(sim_dat$trial[sim_dat$ID == s]) # number of trials pr participant
}
maxntrials <- max(ntrials) # Maximal number of trials a participant has completed

# Empty arrays to fill
target_present <- array(0, dim = c(maxntrials, nagents))
button_pressed <- array(0, dim = c(maxntrials, nagents))
correct <- array(0, dim = c(maxntrials, nagents))
RT <- array(0, dim = c(maxntrials, nagents))
falsealarm <- array(0, dim = c(maxntrials, nagents))
ncorrectRTs <- array(0, dim = nagents)
nincorrectRTs <- array(0, dim = nagents)
condition_dum <- array(0, dim = c(maxntrials, nconditions, nagents))
min_rt <- array(NA, dim = nagents)

for (s in 1:nagents){
  target_present[1:ntrials[s],s] <- sim_dat$target_present[sim_dat$ID == s]
  button_pressed[1:ntrials[s],s] <- sim_dat$button_pressed[sim_dat$ID == s]
  correct[1:ntrials[s],s] <- sim_dat$correct[sim_dat$ID == s]
  RT[1:ntrials[s],s] <- sim_dat$RT[sim_dat$ID == s]
  falsealarm[1:ntrials[s],s] <- sim_dat$FA_trial[sim_dat$ID == s]
  ncorrectRTs[s] <- sum(sim_dat$button_pressed[sim_dat$correct == 1 &
                                                 sim_dat$RT > 0 &
                                                 sim_dat$ID == s])
  nincorrectRTs[s] <- sum(sim_dat$button_pressed[sim_dat$correct == 0 &
                                                   sim_dat$RT > 0 &
                                                   sim_dat$ID == s])
  min_rt[s] <- min(sim_dat$RT[sim_dat$RT > 0 & sim_dat$ID == s])
  
  for (i in 1:ntrials[s]){
    condition_dum[i, sim_dat$condition[sim_dat$ID == s][i], s] <- 1
  }
}


iter_warmup = 1000
iter_sampling = 3000
chains = 4

infer_params <- data.frame()

for (s in 1:nagents){
  print(paste("ITERATION", s))
  
  data <- list(nagents = 1,
               ntrials = ntrials[s],
               maxntrials = maxntrials,
               nconditions = nconditions,
               target_present = array(target_present[,s], dim = c(ntrials[s], 1)),
               button_pressed = array(button_pressed[,s], dim = c(ntrials[s], 1)),
               correct = array(correct[,s], dim = c(ntrials[s], 1)),
               RT = array(RT[,s], dim = c(ntrials[s], 1)),
               falsealarm = array(falsealarm[,s], dim = c(ntrials[s], 1)),
               ncorrectRTs = ncorrectRTs[s],
               nincorrectRTs = nincorrectRTs[s],
               condition_dum = array(condition_dum[ , ,s], dim = c(ntrials[s], nconditions, 1)),
               min_rt = min_rt[s],
               d = CPT_d,
               a_prior_m = a_prior_m,
               a_prior_sd = a_prior_sd,
               b_prior_m = b_prior_m,
               b_prior_sd = b_prior_sd,
               v_prior_m = v_prior_m,
               v_prior_sd = v_prior_sd,
               t_prior_m = t_prior_m,
               t_prior_sd = t_prior_sd,
               FA_prior_m = FA_prior_m,
               FA_prior_sd = FA_prior_sd,
               # d_prior_m = d_prior_m,
               # d_prior_sd = d_prior_sd,
               w_prior_m = w_prior_m,
               w_prior_sd = w_prior_sd
  )
  
  # Start sampling
  temp_samples <- wiener_CPT_IP_mod$sample(
    data = data,
    seed = 123,
    chains = chains,
    parallel_chains = chains,
    threads_per_chain = 4,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    init = 0,
    refresh = 100,
    max_treedepth = 10,
    adapt_delta = 0.99
  )
  
  temp_draws_df <- as_draws_df(temp_samples$draws())
  
  temp_trim_df <- temp_draws_df %>%
    dplyr::select(!starts_with("EEA")) %>% 
    # dplyr::select(!starts_with("theta")) %>% 
    dplyr::select(!starts_with("bias")) %>% 
    dplyr::select(!starts_with("caut"))
  
  temp_trim_df$ID <- s
  
  write.csv(temp_trim_df, sprintf("draws_trim_CPT_IP_%s_agent_%s_init0.csv", datetime, s))
  # infer_params <- rbind(infer_params, temp_trim_df)
}

# write.csv(infer_params, paste("infer_params_CPT_IP_fixed_d", paste0(datetime, ".csv")))


# dtrim <- draws_df %>%
#   dplyr::select(
#     !starts_with("bias")
#   ) %>% dplyr::select(
#     !starts_with("EEA")
#   ) %>% dplyr::select(
#     !starts_with("caut")
#   ) %>% dplyr::select(
#     !starts_with("theta")
#   )
# 
# write.csv(dtrim, paste("draws_trim_CPT_IP_fixed_d", paste0(datetime, ".csv")))
