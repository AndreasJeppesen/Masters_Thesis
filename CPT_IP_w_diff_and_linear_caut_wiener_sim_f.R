CPT_IP_w_diff_and_caut_wiener_sim_f <- function(item_type, v, a, b, tau, d, FA, w) {
  library(RWiener)
  FA_items <- c(2,3,4,5)
  ntrials <- length(item_type)
  
  # Setup the data frame to fill based on the test environment
  dat <- data.frame("trial" = 1:ntrials,
                    "item_type" = item_type,
                    "target_present" = rep(NA, ntrials),
                    "correct" = rep(NA, ntrials),
                    "button_pressed" = rep(NA, ntrials),
                    "RT" = rep(NA, ntrials),
                    "condition" = c(rep(1, 150), rep(2, 150), rep(3, 150)),
                    "FA_trial" = rep(NA, ntrials)
  )
  dat$target_present <- ifelse(dat$item_type == 6, 1, 0)
  dat$FA_trial <- ifelse(dat$item_type %in% FA_items, 1, 0)
  
  # Simulate data trial by trial
  for (i in 1:ntrials){
    c <- dat$condition[i] # which condition are we in
    EEA_given_d <- v - 0.2 * d[c]
    EEA_given_d <- ifelse(EEA_given_d < 0, 0, EEA_given_d)
    EEA <- EEA_given_d + dat$FA_trial[i] * FA
    # caut <- a * (1 + ((v-EEA_given_d)/v) ) * w
    caut <- a + w * 0.2 * d[c]
    
    # Create the bias depending on target presence here (since we interpret "upper" as correct):
    if (dat$target_present[i] == 1){
      bias <- b
      
      reaction <- rwiener(1, alpha = caut, tau = tau, beta = bias, delta = EEA) # generate one reaction using the Wiener model
    }
    
    else { # else if target not present
      bias <- 1-b  
      reaction <- rwiener(1, alpha = caut, tau = tau, beta = bias, delta = EEA) # generate one reaction using the Wiener model
    }
    
    dat$correct[i] <- ifelse(reaction$resp[1] == "upper", 1, 0)   # Interpret "upper" reactions as correct
    
    dat$button_pressed[i] <- ifelse(dat$target_present[i] == 1 & dat$correct[i] == 1, 1, # target present & correct => button press
                                    ifelse(dat$target_present[i] == 0 & dat$correct[i] == 0, 1, # target absent & incorrect => button p
                                           0)) # Else, no button press
    
    dat$RT[i] <- ifelse(dat$button_pressed[i] == 1, reaction$q[1], 0) # if button pressed, save RT
    
  } # end trial loop
  
  return(dat)
}