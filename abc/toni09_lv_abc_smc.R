# ABC SMC from Toni et al. (2009)
# 
# It's not clear what to take as the proposal distribution ... ?
# With some clever optimisation one could gen. the 1st popn outside the loop

# Boolean for whether figures should be output/saved
DO_PLOTTING <- TRUE

set.seed(19)

# Pull in the defined data, project functions, required packages
source('toni09_lv_preamble.R')

# Initialize epsilon values (first run through just saves points from the prior)
epsilons <- c(Inf, 30.0, 16.0, 6.0, 5.0, 4.3)
N.populations <- length(epsilons) # denoted T in Toni et al. (2009)

# Number of particles
N.particles <- 1000

# Bivariate kernel used in the weighting calculation
require(mvtnorm) # could also use MASS::mvrnorm for random draws
dK <- function(theta1, theta2){
    return(dmvnorm(theta1, mean = theta2))
}

rK <- function(mean, sigma){
    return(rmvnorm(1, mean = mean, sigma = sigma))
}


# Set population indicator (t) to 1
t <- 1

# List of previous population of parameter values and weights
population <- list()
weights <- list()
data.gen.steps <- c()

# Start the timer to record how long this takes
start_time <- Sys.time()

while(t <= N.populations){
    
    # Set particle indicator (i) to 1 (since R starts indexing from 1)
    i <- 1
    population[[t]] <- list()
    weights[[t]] <- rep(NA, N.particles) # pre-allocate weights vector
    
    # Record how many data generation steps there are
    data.gen.counter <- 0
    while(i <= N.particles){
        
        # Reset the distance
        distance <- Inf
        
        while(distance >= epsilons[t]){
            # Increment the counter for the number of data generation steps
            data.gen.counter <- data.gen.counter + 1
            
            cat('.')
            if(t == 1){
                theta_star_star <- runif(2, min = -10, max = 10)
                names(theta_star_star) <- c("a", "b")
            }else{
                # Sample a new theta_star from the previous popn
                # using previous weights
                theta_star <- sample(population[[t-1]], 1, 
                        prob = weights[[t-1]])
                theta_star <- unlist(theta_star)
                
                # Perturb theta* to get theta**
                theta_star_star <- rK(mean = theta_star, sigma = 0.25*diag(2))
                names(theta_star_star) <- c("a", "b")
            }
            # Need to wrap the above in a while look to make 
            # sure the parameters are appropriate: 
            # i.e. make sure pi(theta_star_star) > 0 
            
            # Data generation step
            try(out <- lsoda(func = LV,
                    y = init,
                    parms = theta_star_star,
                    atol = 1e-3,
                    times = c(0, sample_data$time)),
                silent = TRUE)
            dfs <- as.data.frame(out)
            
            # Only proceed if lsoda didn't throw and error
            if(NROW(dfs) == (NROW(sample_data) + 1)){
                
                # Find the distance (sse) from the 
                # simulated points to the data
                x_star <- dfs[-1, c("x", "y")]
                x0 <- sample_data[, c("x", "y")]
                
                distance <- sse(x0, x_star)
                
            }
        }
        
        # Set the current particle to theta_star_star
        population[[t]][[i]] <- theta_star_star
        
        # Calculate the weight for this particle
        if (t == 1){
            weights[[t]][i] <- 1
        }else{
            pi_theta <- 0.0025 # 1/(20*20)
            
            thetas <- population[[t-1]]
            w <- weights[[t-1]]
            Ks <- sapply(thetas, function(x) dK(x, theta_star_star))
            denominator <- sum(w * Ks)
            
            weights[[t]][i] <- pi_theta * denominator
        }
        cat("*", t, "/", i)
        i <- i + 1
    } # while(i <= N.particles)
    
    # Record the number of data generation steps taken
    data.gen.steps <- c(data.gen.steps, data.gen.counter)
    
    # Normalize the weights; increment t
    weights[[t]] <- weights[[t]]/sum(weights[[t]], na.rm = TRUE)
    t <- t + 1
    
}# while(t <= N.populations)

end_time <- Sys.time()
(total_time <- end_time - start_time)
# Time difference of 14.27762 mins (35.75432 mins using mv norm dist)

# Cumulative number of data generation steps
cumsum(data.gen.steps)

# Make a long dataset of the final population of particles
full <- list()
for (t in 1:N.populations){
    final.pop <- do.call(rbind, population[[t]])
    long.final.pop <- melt(final.pop)
    long.final.pop$population <- t
    full[[t]] <- long.final.pop
}
long.all.pop <- do.call(rbind, full)
names(long.all.pop) <- c("particleID", "parameter", "estimate", "population")


# Plot results?  
if(DO_PLOTTING){
    par(mfrow = c(length(epsilons), 2), mar = c(1,1,1,1))
    axlimits <- list(c(-10, 10), c(-2, 3), c(0, 2), c(0.8, 1.4), c(0.8, 1.4), c(0.8, 1.4))
    bxlimits <- list(c(-10, 10), c(-10, 10), c(0, 3.5), c(0.6, 1.5), c(0.6, 1.5), c(0.6, 1.4))
    
    bins <- c(10, 10, 10, 8, 5, 4)
    for(i in 1:N.populations){
        breaks <- bins[i]
        for (p in c("a", "b")){
            
            if(p == "a"){
                xlims <- axlimits[[i]]
            }else{
                xlims <- bxlimits[[i]]
            }
            
            c1 <- with(long.all.pop, parameter == p)
            c2 <- with(long.all.pop, population == i)
                hist(long.all.pop[c1&c2,"estimate"], bty = 'n', 
                    ylim = c(0, 400), xlim = xlims, 
                    main = "", breaks = breaks)
        }
    }
    
    
    output_plot <- ggplot(long.all.pop, 
            aes(estimate, fill = factor(parameter))) + 
            geom_histogram(bins = 20) +
            facet_grid(population ~ parameter, scales = "free")
    #output_plot
    pdf('../output/figure_1b.pdf')
    print(output_plot)
    dev.off()
}
