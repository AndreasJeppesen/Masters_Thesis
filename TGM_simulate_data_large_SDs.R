setwd("/work/AndreasHjorthJeppesen#0476/R")



source("Setup_and_init/activate_renv.R")



##Load packages

library(pacman)
pacman::p_load(tidyverse, cmdstanr, posterior, ggridges, MASS)



##Load data simulation functions

# For CPT-IP
source("CPT_IP_w_diff_and_linear_caut_wiener_sim_f.R")
# For Conners' CPT
source("Conners_TGM_sim_f.R")
# For RVP
source("RVP_TGM_sim_f.R")



#Simulation environment

nagents <- 30
ntrials_pr_task <- c(450, 360, 300) # CPT-IP, Conners', RVP
nconditions_pr_task <- c(3, 1, 1)
ntargets_pr_task <- c(90, 324, 24)
nFAtrials_RVP <- 5
CPT_IP_stimuli <- read.csv("CPT_IP_stimuli.csv")
ntasks <- 3



#Simulate data

#--Set true param values

# Set true param values
a_m <- runif(nagents, 0.8, 2)
a_sd <- runif(nagents, 0.06, 0.12)
b_m <- runif(nagents, 0.3, 0.7)
b_sd <- runif(nagents,  0.06, 0.12)
v_m <- runif(nagents, 0.8, 2)
v_sd <- runif(nagents, 0.06, 0.12)
t_m <- runif(nagents, 0.2, 0.25)
t_sd <- runif(nagents, 0.010, 0.020)

FA_means <- c(0, 0, 0)
corr_matrix <- matrix(c(1, 0.6, 0.6,
                        0.6, 1, 0.6,
                        0.6, 0.6, 1), nrow = 3)
variances <- diag(1, nrow = 3, ncol = 3) # In this case, Identity matrix
covar_matrix <- variances %*% corr_matrix %*% variances

FA <- mvrnorm(nagents, mu = FA_means, Sigma = corr_matrix) # Generate 3 correlated variables with m=0 and sd=1
FA[,1] <- FA[,1] * 0.2 - 1.2 # change the mean and sd of each variable, so we dont get positive values of FA
FA[,2] <- FA[,2] * 0.2 - 1.5
FA[,3] <- FA[,3] * 0.2 - 0.6

w <- runif(nagents, 0, 2)

d <- array(c(0.5,2,4,
             0.1,NA,NA,
             1,NA,NA), dim = c(ntasks, max(nconditions_pr_task)))
d_vec <- array(d, dim = 9)
d_vec <- d_vec[!is.na(d_vec)]


#-- Generate data

# Task level param values
true_as <- array(NA, dim = c(nagents, ntasks))
true_bs <- array(NA, dim = c(nagents, ntasks))
true_vs <- array(NA, dim = c(nagents, ntasks))
true_ts <- array(NA, dim = c(nagents, ntasks))

sim_dat <- data.frame()

for (s in 1:nagents){
  # draw param values for each of the tasks
  true_as[s,] <- rnorm(ntasks, a_m[s], a_sd[s])
  true_bs[s,] <- rnorm(ntasks, b_m[s], b_sd[s])
  true_vs[s,] <- rnorm(ntasks, v_m[s], v_sd[s])
  true_ts[s,] <- rnorm(ntasks, t_m[s], t_sd[s])
  
  # Generate CPT-IP data
  CPT_IP_sim <- CPT_IP_w_diff_and_caut_wiener_sim_f(item_type = CPT_IP_stimuli$stim_type,
                                                    v = true_vs[s,1], a = true_as[s,1],
                                                    b = true_bs[s,1], tau = true_ts[s,1],
                                                    d = d[,1], FA = FA[s,1], w = w[s])
  
  CPT_IP_sim <- CPT_IP_sim %>% subset(select = -item_type)
  CPT_IP_sim$ID <- s
  CPT_IP_sim$task_name <- "CPT-IP"
  CPT_IP_sim$task_num <- 1
  sim_dat <- rbind(sim_dat, CPT_IP_sim)
  
  # Generate Conners' CPT data
  conners_sim <- Conners_TGM_sim_f(ntrials_pr_task[2], ntargets_pr_task[2],
                                   true_vs[s,2], true_as[s,2],
                                   true_bs[s,2], true_ts[s,2],
                                   d[1,2], FA[s,2], w[s])
  conners_sim$ID <- s
  conners_sim$task_name <- "Conners"
  conners_sim$task_num <- 2
  conners_sim$condition <- 4
  sim_dat <- rbind(sim_dat, conners_sim)
  
  # Generate RVP data
  RVP_sim <- RVP_TGM_sim_f(ntrials_pr_task[3], ntargets_pr_task[3], nFAtrials_RVP,
                           true_vs[s,3], true_as[s,3],
                           true_bs[s,3], true_ts[s,3],
                           d[1,3], FA[s,3], w[s])
  RVP_sim$RT <- ifelse(RVP_sim$correct == 1, RVP_sim$RT, 0) # No incorrect RTs in RVP
  RVP_sim$ID <- s
  RVP_sim$task_name <- "RVP"
  RVP_sim$task_num <- 3
  RVP_sim$condition <- 5
  sim_dat <- rbind(sim_dat, RVP_sim)
  
}
sim_dat$incorrectRTs_present <- ifelse(sim_dat$task_name == "RVP", 0, 1)


# Save true param values to csv
true_params <- data.frame(s = seq(1:nagents), a_m, a_sd, b_m, b_sd, v_m, v_sd, t_m, t_sd, w)

for (task in 1:ntasks){
  tnum <- as.character(task)
  temp <- data.frame(matrix(NA, nrow = nagents, ncol = 5))
  colnames(temp) <- c(paste0("a", tnum),
                      paste0("b", tnum),
                      paste0("v", tnum),
                      paste0("t", tnum),
                      paste0("FA", tnum))
  temp[,1] <- true_as[,task]
  temp[,2] <- true_bs[,task]
  temp[,3] <- true_vs[,task]
  temp[,4] <- true_ts[,task]
  temp[,5] <- FA[,task]
  true_params <- cbind(true_params, temp)
}

datetime <- as.character(Sys.time())
datetime <- gsub(":", "", datetime)

write.csv(true_params, paste("TGM_unr_params_large_SD_true_params", paste0(datetime, ".csv")))
write.csv(sim_dat, paste("TGM_unr_params_large_SD_sim_dat", paste0(datetime, ".csv")))
