drift_process_plot <- function(a, b, v, niter, xlim = NULL){
  drift_sim <- data.frame()
  
  for (i in 1:niter){
    drift_sim_temp = data.frame("time" = 0, "total" = a * b, "iteration" = i) # initiate new temp data
    
    total = a*b # reset total and time
    time = 0
    
    while (total > 0 & total < a){ # Start drift process
      
      drift = rnorm(1, v/1000, 0.1 / 5) # Sample the drift
      
      total = total + drift # add drift to total
      
      time = time + 1 / 1000
      
      drift_sim_temp = rbind(drift_sim_temp, data.frame("time" = time, "total" = total, "iteration" = i))
    }
    drift_sim <- rbind(drift_sim, drift_sim_temp) # Add iteration to 
  }
  
  if (is.null(xlim)){
    drift_sim %>% ggplot(aes(time, total, color = as.factor(iteration))) +
      geom_line() +
      ylim(-0.01, a+0.1) +
      geom_segment(aes(x = 0, xend = (a-(a*b))/v, y = a * b, yend = a), color = "black", size = 1) +
      geom_segment(aes(x = 0, xend = max(time)+0.1, y = a, yend = a), color = "red", linetype = "dashed", size = 1) +
      geom_segment(aes(x = 0, xend = max(time)+0.1, y = 0, yend = 0), color = "red", linetype = "dashed", size = 1) +
      labs(title = "Drift Processes", x = "Time (s)", y = "Total", color = "Drift number") +
      annotate(geom = "text", x = -0.02-(0.07*max(time)), y = a + 0.01, label = "a", color = "black", size = 5) +
      annotate(geom = "text", x = -0.02-(0.07*max(time)), y = (a * b) + 0.01, label = "z", color = "black", size = 5) +
      annotate(geom = "text", x = (a-(a*b))/v, y = a+(0.05*a), label = paste("v =", v), color = "black", size = 4) +
      theme_minimal()
  } else {
    drift_sim %>% ggplot(aes(time, total, color = as.factor(iteration))) +
      geom_line() +
      xlim(xlim[1], xlim[2]) +
      ylim(-0.01, a+0.1) +
      geom_segment(aes(x = 0, xend = (a-(a*b))/v, y = a * b, yend = a), color = "black", size = 1) +
      geom_segment(aes(x = 0, xend = xlim[2], y = a, yend = a), color = "red", linetype = "dashed", size = 1) +
      geom_segment(aes(x = 0, xend = xlim[2], y = 0, yend = 0), color = "red", linetype = "dashed", size = 1) +
      labs(title = "Drift Processes", x = "Time (s)", y = "Total", color = "Drift number") +
      annotate(geom = "text", x = -0.02-(0.07*xlim[2]), y = a + 0.01, label = "a", color = "black", size = 5) +
      annotate(geom = "text", x = -0.02-(0.07*xlim[2]), y = (a * b) + 0.01, label = "z", color = "black", size = 5) +
      annotate(geom = "text", x = (a-(a*b))/v, y = a+(0.05*a), label = paste("v =", v), color = "black", size = 4) +
      theme_minimal()
  }
  
}