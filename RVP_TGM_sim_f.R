RVP_TGM_sim_f <- function(ntrials, ntargets, nFAtrials, v, a, b, tau, d, FA, w) {
  library(RWiener)
  # Create test environment
  dat <- data.frame("trial" = 1:ntrials,
                    "target_present" = sample(c(rep(1, ntargets), rep(0, ntrials-ntargets))), # distribute targets randomly
                    "FA_trial" = rep(NA, ntrials),
                    "correct" = rep(NA, ntrials),
                    "button_pressed" = rep(NA, ntrials),
                    "RT" = rep(NA, ntrials)
  )
  
  # Simulate data trial by trial
  for (i in 1:ntrials){
    
    # Assign FA_trial to some non-target trials
    dat$FA_trial[i] <- ifelse(dat$target[i] != 1, rbinom(1, 1, prob = (nFAtrials/(ntrials-ntargets))), 0)
    
    # EEA_given_d <- v*0.95^d
    EEA_given_d <- v - 0.2 *d
    EEA_given_d <- ifelse(EEA_given_d < 0, 0, EEA_given_d)
    EEA <- EEA_given_d + FA * dat$FA_trial[i]
    # caut <- a * (1 + ((v-EEA_given_d)/v) ) * w
    caut <- a + w * 0.2 * d
    
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
    
    dat$RT[i] <- ifelse(dat$button_pressed[i] == 1 & dat$correct[i] == 1, reaction$q[1], 0) # if button pressed & correct, save RT
    
  } # trial loop end
  
  return(dat)
}