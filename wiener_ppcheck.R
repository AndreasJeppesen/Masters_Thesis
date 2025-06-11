wiener_ppcheck <- function(RT_data, draws_df, prior_only = FALSE, nconditions = 1, n_draws, xlim = c(0,4), ylim = c(0,4)){
  
  x <- seq(from = 0, to = 4, length.out = 101) # x help line to calculate densities on
  
  if (nconditions == 1){
  
    densities <- array(NA, dim = c(101, n_draws))
    
    
    if (prior_only){ # Prior Predictive check
      for (i in 1:n_draws){
        # sample from the prior distribution
        row_num <- floor(runif(1, min = 1, max = length(draws_df$v + 1)))
        # Calulate densities:
        densities[, i] <- dwiener(x, draws_df$a_prior[row_num], draws_df$t_prior[row_num],
                                  draws_df$b_prior[row_num], draws_df$v_prior[row_num],
                                  rep("upper", length(x)))
      }
    } else {  # Posterior Predictive check
      for (i in 1:n_draws){
        # sample from the posterior distribution
        row_num <- floor(runif(1, min = 1, max = length(draws_df$v + 1)))
        # Calulate densities:
        densities[, i] <- dwiener(x, draws_df$a[row_num], draws_df$t[row_num],
                                  draws_df$b[row_num], draws_df$v[row_num],
                                  rep("upper", length(x)))
      }
    }
    
    plot_title <- ifelse(prior_only, "Prior Predictive Check", "Posterior Predictive Check")
    
    plot(x = x, y = densities[,1], type = "l", col = "lightblue", xlim = xlim, ylim = ylim,
         main = plot_title, xlab = "RT (s)", ylab = "Density")
    for (i in 2:n_draws){
      lines(x = x, y = densities[,i], type = "l", col = "lightblue")
    }
    lines(density(RT_data, adjust = 0.7), col = "black", lwd = 2)
    legend("topright", legend = c("y", "pred"), col = c("black", "lightblue"), lty = 1)
    
    
    
  } else {  # if multiple conditions:
    
      densities <- array(NA, dim = c(101, n_draws))
      
      n_draws_pr_difficulty <- floor(n_draws / nconditions)
      
      
      if (prior_only){ # Prior Predictive check
        # Get the (unknown number of) difficulty samples in a seperate df:
        d_samples <- select(draws_df, starts_with("d_prior["))
        
        for (c in 1:nconditions){
          for (i in (n_draws_pr_difficulty * (c-1) + 1):(n_draws_pr_difficulty * c)){ # calculate n_draws_pr_difficulty densities at a time
            # sample from the prior distribution
            row_num <- floor(runif(1, min = 1, max = length(draws_df$v + 1)))
            # Calulate densities:
            densities[, i] <- dwiener(x, draws_df$a_prior[row_num], draws_df$t_prior[row_num],
                                      draws_df$b_prior[row_num],
                                      draws_df$v_prior[row_num] * 0.95 ^ as.numeric(d_samples[row_num, c]),
                                      rep("upper", length(x)))
          }
        } # end condition loop
        
        
      } else {  # Posterior Predictive check
        
        d_samples <- select(draws_df, starts_with("d["))
        
        for (c in 1:nconditions){
          for (i in (n_draws_pr_difficulty * (c-1) + 1):(n_draws_pr_difficulty * c)){ # calculate n_draws_pr_difficulty densities at a time
            # sample from the posterior distribution
            row_num <- floor(runif(1, min = 1, max = length(draws_df$v + 1)))
            # Calulate densities:
            densities[, i] <- dwiener(x, draws_df$a[row_num], draws_df$t[row_num],
                                      draws_df$b[row_num],
                                      draws_df$v[row_num] * 0.95 ^ as.numeric(d_samples[row_num, c]),
                                      rep("upper", length(x)))
          }
        } # end condition loop
      }
      
      plot_title <- ifelse(prior_only, "Prior Predictive Check", "Posterior Predictive Check")
      
      plot(x = x, y = densities[,1], type = "l", col = "lightblue", xlim = xlim, ylim = ylim,
           main = plot_title, xlab = "RT (s)", ylab = "Density")
      for (i in 2:n_draws){
        lines(x = x, y = densities[,i], type = "l", col = "lightblue")
      }
      lines(density(RT_data, adjust = 0.7), col = "black", lwd = 2)
      legend("topright", legend = c("y", "pred"), col = c("black", "lightblue"), lty = 1)
    
  }
  
  
}