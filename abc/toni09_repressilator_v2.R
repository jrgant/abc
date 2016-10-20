# Repressilator (Elowitz and Leibler, 2000)
# as outlined in Toni et al. 2009
# 
# Questions ... 
# 
# "We assume that the initial conditions are known."
# How many data points were drawn in the first instance as the test data? 
# What were the epsilon values used?  
# How many data generation steps were performed (guess this depends upon the epsilons)?  
# What perturbation kernel was used?  Does this matter that much?  
# How to determine the number of particles to use?  

require(reshape2) # for reshape2::melt()
require(deSolve) # for deSolve::lsoda()
require(ggplot2) # only if plotting 

# Define square of sum of squared errors
sse <- function(x, xd){
    return( sqrt(sum((x - xd)^2)))
}

set.seed(19)


# Logical for whether to output plotting
DO_PLOTTING <- FALSE

# Differential equations defining the system
Repressilator <- function(Time, State, Pars){
    
    with(as.list(c(State, Pars)), {
    
    dm1 <- -m1 + alpha/(1+p3^n) + alpha0
    dp1 <- -beta * (p1 - m1)
    dm2 <- -m2 + alpha/(1+p1^n) + alpha0
    dp2 <- -beta * (p2 - m2)
    dm3 <- -m3 + alpha/(1+p2^n) + alpha0
    dp3 <- -beta * (p3 - m3)
    
    return(list(c(dm1, dp1, dm2, dp2, dm3, dp3)))
    }
    )
}

# Parameters
pars <- c(  alpha0 = 1, 
            n = 2, 
            beta = 5,
            alpha = 1000)

# Initial conditions
init <- c(  m1 = 0, 
            p1 = 2, 
            m2 = 0,
            p2 = 1,
            m3 = 0,
            p3 = 3)

# Time steps
times <- seq(0, 20, length.out = 5000)

# Solve, coerce into data.frame
	out <- lsoda(func = Repressilator, y = init, parms = pars, times = times)
rep.df <- as.data.frame(out)

# Reshape into long format
long.rep.df <- melt(rep.df, id.vars = c("time"))

if(DO_PLOTTING){
    ggplot(long.rep.df, aes(x = time, y = value, 
        group = factor(variable), colour = factor(variable))) + 
    geom_line()
}

# Generate some data (choosing 19 data points)
sample_n <- 19 
sampled_rows <- sample(NROW(rep.df), sample_n, replace = TRUE)
sample.rep.df <- rep.df[sampled_rows,]
x0 <- sample.rep.df[,c("m1", "m2", "m3")]

# Add noise to create the simulated data
noise <- matrix(rnorm(3*sample_n, mean = 0, sd = 5), sample_n, 3)
sample.rep.df[,c("m1", "m2", "m3")] <- sample.rep.df[,c("m1", "m2", "m3")]+noise
sample.rep.df <- sample.rep.df[,c("time", "m1", "m2", "m3")]
x <- sample.rep.df[,c("m1", "m2", "m3")]
epsilon_0 <- sse(x0, x)


# "Infer" the "pseudo-posterior" distribution from Figure 4 (a)
# This should give an indicator for values of epsilons
epsilons_post <- c()
 for(i in 1:10^3) {
 	theta <- c(runif(1, min = 0.7, max = 1.3), runif(1, min = 1.8, max = 2.3), runif(1, min = 3.5, max = 8), runif(1, min = 800, max = 1400))
 	names(theta) <- c("alpha0", "n", "beta", "alpha")
	out <- lsoda(func = Repressilator, y = init, parms = theta, times = times)
            	dfs <- as.data.frame(out)
            
                x <- dfs[sampled_rows,c("m1", "m2", "m3")]
                x0 <- sample.rep.df[,c("m1", "m2", "m3")]
                
                dist <- sse(x0, x)
                epsilons_post <- c(epsilons_post, dist)	
}

if(DO_PLOTTING){ hist(epsilons_post) }

# Initialize epsilon values based on the above quantiles and sse between data and actual solution
epsilons <- c(quantile(epsilons_post, c(0.75, 0.5, 0.3, 0.2, 0.1, 0.05, 0.025, 0.01), names=FALSE), epsilon_0 + 2)
print("Epsilons: ")
print(epsilons)

N.populations <- length(epsilons) # denoted T in Toni et al. (2009)

# Number of particles
N.particles <- 1000

# perturbation ranges
sigma <- c(2, 2, 5, 100)

# Perturbation kernel, adds some noise based on the above
Kt <- function(x, sigma){
    return(x + c(runif(1, min = -sigma[[1]], max = sigma[[1]]), runif(1, min = -sigma[[2]], max = sigma[[2]]), runif(1, min = -sigma[[3]], max = sigma[[3]]), runif(1, min = -sigma[[4]], max = sigma[[4]])))
}

# Set population indicator (t) to 1
t <- 1

# List of previous population of parameter values and weights
population <- list()
data.gen.steps <- c()

# Start the timer to record how long this takes
start_time <- Sys.time()

while(t <= N.populations){
    
    # Set particle indicator (i) to 1 (since R starts indexing from 1)
    i <- 1
    population[[t]] <- list()
    
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
                theta_star_star <- c(runif(1, min = -2, max = 10), runif(1, min = 0, max = 10), runif(1, min = 0, max = 20), runif(1, min = 500, max = 2500))            
                names(theta_star_star) <- c("alpha0", "n", "beta", "alpha")
            }else{
                # Sample a new theta_star from the previous popn
                theta_star <- sample(population[[t-1]], 1)
                theta_star <- unlist(theta_star)
                
                # Perturb theta* to get theta**
                theta_star_star <- Kt(theta_star, sigma)
                names(theta_star_star) <- c("alpha0", "n", "beta", "alpha")
            }

            
            if((theta_star_star[[1]]>=-2) && (theta_star_star[[1]]<=10) && (theta_star_star[[2]]>=0) && (theta_star_star[[2]]<=10) && (theta_star_star[[3]]>=-5) && (theta_star_star[[3]]<=20) && (theta_star_star[[4]]>=500) && (theta_star_star[[4]]<=2500)) {
            	# Data generation step
            	out <- lsoda(func = Repressilator, y = init, parms = theta_star_star, times = times)
            	dfs <- as.data.frame(out)
            
                 # Only proceed if lsoda didn't throw and error or NaN
            if((NROW(dfs) == length(times)) && (is.na(dfs[length(times),"m1"])==FALSE) ){
                # Find the distance (sse) from the 
                # simulated points to the data
  			   x_star <- dfs[sampled_rows,c("m1", "m2", "m3")]
                x0 <- sample.rep.df[,c("m1", "m2", "m3")]
                

                
                distance <- sse(x0, x_star)
            }

           
        	} 
        }
        
        # Set the current particle to theta_star_star
        population[[t]][[i]] <- theta_star_star
        
         cat("*", t, "/", i)
        i <- i + 1
    } # while(i <= N.particles)
    
    # Record the number of data generation steps taken
    data.gen.steps <- c(data.gen.steps, data.gen.counter)
    
     t <- t + 1
    
}# while(t <= N.populations)

end_time <- Sys.time()
(total_time <- end_time - start_time)




