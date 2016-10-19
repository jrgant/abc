# ABC rejection sampler of Pritchard et al. 1999
# as outlined in Toni et al. 2009
# 
# First example: 
# Deterministic Lotka-Volterra

# Pull in the defined data, and project functions
source('toni09_lv_preamble.R')

# Set parameters, threshold value and number of data generation steps
epsilon <- sse(sample_data[,c("x", "y")], true_soln[,c("x", "y")])
N <- 1000000 # be careful, 1 mill will take ages

# Set up containers for results, counters, etc
accepted_thetas <- list()
counter <- 1

# Start the stopwatch
start_time <- Sys.time()

for (i in 1:N){
    
    # Sample from pi(theta), prior for parameter vector (a, b)
    theta_star <- runif(2, -10, 10)
    names(theta_star) <- c("a", "b")
    
    # Simulate a dataset using these parameters (at the specified times)
    try(out <- lsoda(func = LV, 
            y = init, 
            parms = theta_star, 
            atol = 1e-3,
            times = times), 
            silent = TRUE)
    dfs <- as.data.frame(out)
    
    # Only proceed if lsoda didn't throw and error
    if(NROW(dfs) == length(times)){
        
        # Find the distance (sse) from the simulated points to the data
        x_star <- dfs[sample_rows,c("x", "y")]
        distance <- sse(sample_data[,c("x", "y")], x_star)
        
        if(distance < epsilon){
            cat("a ")
            cat(i, " ")
            accepted_thetas[[counter]] <- theta_star
            counter <- counter + 1
        }
    }
}

# Print the final time step
end_time <- Sys.time()
(total_time <- end_time - start_time)
# Time difference of 2.41 hours (for 1 million data generation steps)

# Turn into one dataset
output <- do.call(rbind, accepted_thetas)

# Plot a hist of a or b (see fig 1b of Toni et al 2009)
#hist(output[,1], xlim = c(0.7, 1.4))
#hist(output[,2], xlim = c(0.7, 1.4))

#write.csv(output, file.path(main_dir, "0_2_1mill_runs.csv"))

# Acceptance rate of (Toni et al. (2009) reported 7E-5): 
# (length(accepted_thetas)/N)
# 0.005343
