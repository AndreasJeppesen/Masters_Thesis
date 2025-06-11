Conners_TGM_sim_f <- function(ntrials, ntargets, v, a, b, tau, d, FA, w) {
  library(RWiener)
  
  n_conditions <- length(d) # number of conditions
  output_data <- data.frame() # empty df to fill
  
  for (c in 1:n_conditions){
    
    # Create test environment (for this condition)
    dat <- data.frame("trial" = 1:ntrials[c],
                      "target_present" = sample(c(rep(1, ntargets[c]), rep(0, ntrials[c]-ntargets[c]))), # distribute targets randomly
                      "correct" = rep(NA, ntrials[c]),
                      "button_pressed" = rep(NA, ntrials[c]),
                      "RT" = rep(NA, ntrials[c]),
                      "condition" = c,
                      "FA_trial" = rep(0, ntrials[c])
    )
    dat$FA_trial <- ifelse(dat$target_present == 0, 1, 0) # All non-target trials are considered FA-trials in the TGM setup
    # Simulate data trial by trial
    for (i in 1:ntrials[c]){
      
      EEA_given_d <- v - 0.2 *d[c]
      EEA_given_d <- ifelse(EEA_given_d < 0, 0, EEA_given_d)
      EEA <- EEA_given_d + FA * dat$FA_trial[i]
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
    output_data <- rbind(output_data, dat)
  } # end condition loop
  
  
  
  return(output_data)
}