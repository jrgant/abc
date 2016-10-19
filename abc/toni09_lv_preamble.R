# Script to reproduce the output from Toni et al., 2009
# First analysis using Lotka-Volterra equations; this script generates the 
# data for the analysis and defines the functions to be used.  

require(reshape2) # for reshape2::melt()
require(deSolve) # for deSolve::lsoda()
require(ggplot2) # only if plotting 

# Turn on plotting?  
DO_PLOTTING <- FALSE

# Define sum of squared errors
sse <- function(x, xd){
    return( sum((x - xd)^2))
    #return( sum((x[,1] - xd[,1])^2 + (x[,2] - xd[,2])^2) )
}

# Differential equations defining the system
LV <- function(Time, State, Pars){
    
    with(as.list(c(State, Pars)), {
    
    dx <- a*x - x*y
    dy <- b*x*y - y
    
    return(list(c(dx, dy)))
    }
    )
}

# Parameters
pars <- c(a = 1.0, b = 1.0)

# Initial conditions (from eyeing off the graph)
init <- c(x = 1.0, y = 0.5)

# Time steps to simulate
all_times <- seq(0, 15, length.out = 151)

# Solve and coerce the output into a data.frame
lv_soln <- lsoda(func = LV, y = init, parms = pars, times = all_times)
lv.df <- as.data.frame(lv_soln)

# Read in the sample data
sample_data <- read.csv("../data/toni09_lv_data.csv")

# Evaluate the true solution at the x positions of interest
sample_times <- c(0, sample_data$time)
true_soln <- lsoda(func = LV, y = init, parms = pars, times = sample_times)
true_soln <- as.data.frame(true_soln[-1,])

#true_soln <- lv.df[lv.df$time == sample_data$time,]

# Calculate SSE for the sample data
sse(sample_data[,c("x", "y")], true_soln[,c("x", "y")])
#>> [1] 4.313869


# Plot the solution with sample points overlaid
if (DO_PLOTTING){
    png("../output/figure_1.png", width = 480, height = 480, res = 90)
    with(lv.df, plot(time, x, 
        col = 'grey', type = 'l', bty = 'n',
        xlim = c(0, 15), ylim = c(0, 4), xlab = "", ylab = ""))
    with(lv.df, lines(time, y, col = 'black', type = 'l', lty = 2))
    
    with(sample_data, points(time, x, pch = 21, col = 'grey', cex = 0.7))
    with(sample_data, points(time, y, pch = 24, col = 'black', cex = 0.7))
    dev.off()
}
