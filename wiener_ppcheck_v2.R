wiener_ppcheck <- function(RT_data,
                           ntrials_pr_condition,
                           a_samples,
                           b_samples,
                           v_samples,
                           t_samples,
                           w_samples,
                           nconditions = 1,
                           d = 0,
                           niter = 100,
                           prior_only = FALSE,
                           xlim = c(0,4),
                           ylim = c(0,4)
                           ){
  library(RWiener)
  nsamples <- length(a_samples)
  subset_idx <- floor(runif(niter, min = 1, max = nsamples+1)) # Draw 100 random iterations to visualize
  EEA <- array(NA, dim = c(nconditions, niter))
  caut <- array(NA, dim = c(nconditions, niter))
  
  pred_RT <- array(NA, dim = c(sum(ntrials_pr_condition), niter))
  
  for (c in 1:nconditions){
    EEA[c, ] <- v_samples[subset_idx] - 0.2 * d[c]
    EEA[c, ] <- ifelse(EEA[c, ] < 0, 0, EEA[c, ])
    caut[c, ] <- a_samples[subset_idx] + 0.2 * w_samples[subset_idx] * d[c]
  }
  
  
  for (i in 1:niter){
    print(paste("Simulating iteration", i))
    pred_RT_temp <- c()
    
    for (c in 1:nconditions){
      # Simulate data
      N <- ntrials_pr_condition[c]
      reactions <- rwiener(N,
                           alpha = caut[c, i],
                           tau = t_samples[subset_idx[i]],
                           beta = b_samples[subset_idx[i]],
                           delta = EEA[c, i])
      # Save correct RTs
      corr_RTs <- reactions$q[reactions$resp == "upper"]
      
      # append to pred_RT_temp
      pred_RT_temp <- c(pred_RT_temp, corr_RTs)
    } # end condition loop
    
    # save in pred_RT
    pred_RT[1:length(pred_RT_temp), i] <- pred_RT_temp
  } # end iter loop
  
  
  # Visualize RT data and simulated RTs
  
  plot_title <- ifelse(prior_only, "Prior Predictive Check", "Posterior Predictive Check")
  
  plot(density(pred_RT[,1][!is.na(pred_RT[,1])], adjust = 0.7), col = "lightblue", xlim = xlim, ylim = ylim,
       main = plot_title, xlab = "RT (s)", ylab = "Density")
  for (i in 2:niter){
    lines(density(pred_RT[,i][!is.na(pred_RT[,i])], adjust = 0.7), col = "lightblue")
  }
  lines(density(RT_data, adjust = 0.7), col = "black", lwd = 2)
  legend("topright", legend = c("y", "pred"), col = c("black", "lightblue"), lty = 1)
}