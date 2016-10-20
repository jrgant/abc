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
times <- seq(0, 10, length.out = 1000)

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

# Generate some data (choosing 8 data points)
sample_n <- 8 
sampled_rows <- sample(NROW(rep.df), 8, replace = TRUE)
sample.rep.df <- rep.df[sampled_rows,]

# Add noise to create the simulated data
noise <- matrix(rnorm(3*sample_n, mean = 0, sd = 5), sample_n, 3)
sample.rep.df[,c("m1", "m2", "m3")] <- sample.rep.df[,c("m1", "m2", "m3")]+noise
sample.rep.df <- sample.rep.df[,c("time", "m1", "m2", "m3")]

# Define the prior distributions
pi_alpha0 <- runif(1, min = -2, max = 10)
pi_n <- runif(1, min = 0, max = 10)
pi_beta <- runif(1, min = -5, max = 20)
pi_alpha <- runif(1, 500, 2500)

epsilons <- c(9000, 4000, 2000, 1000)
