#how can tau be continuous when our data time points are discrete?
#tau prior is currently discrete tau_star<-ceiling(runif(1,-0.5,5))

#S0 is being perturbed into being continuous...



source('toni09_model_section_preamble.R')
toni09_tdc_data <- read.csv("toni09_tdc_data.csv")

require(deSolve)
require(mvtnorm)

# Boolean for whether figures should be output/saved
DO_PLOTTING <- TRUE

# Initialize epsilon values 
#epsilons<-c(Inf,100, 90, 80, 73, 70, 60, 50, 40, 30, 25, 20, 16, 15, 14, 13.8)

epsilons<-c(Inf,100, 90,80,73)

N.populations <- length(epsilons) # denoted T in Toni et al. (2009)

# Number of particles
N.particles <- 1000

# List of previous population of parameter values and weights
#popultation and weights need to be lists of lists
#first index is t, second index is m
population <- rep(list(list()),N.populations )
weights <- rep(list(list()),N.populations )
m_vec<-rep(list())
data.gen.steps <- c()

t<-1

# Start the timer to record how long this takes
start_time <- Sys.time()


while(t <= N.populations){
  # Set particle indicator (i) to 1 (since R starts indexing from 1)
  i <- 1
  
   population[[t]][[1]] <- list();population[[t]][[2]] <- list()
   population[[t]][[3]] <- list();population[[t]][[4]] <- list()
  # 
  # # pre-allocate weights vector
   weights[[t]][[1]] <-list(); weights[[t]][[2]] <- list(); 
   weights[[t]][[3]] <- list(); weights[[t]][[4]] <- list(); 
  
  m_vec[[t]]<-list()
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
        #Sample which model you will use
         m_star<-ceiling(runif(1,0,4))
        
        #Find the associated parameter vector with m_star
         theta_star_star<-theta_star_fn(m_star)

      }else{
        #Sample which model you will use
        #If a model dies out, then we just sample from the others
        repeat{
          m_star<-ceiling(runif(1,0,4))
          if(length(population[[t-1]][[m_star]])!=0) break
        } 
        
        # Sample a new theta_star from the previous popn
        # using previous weights
        
        #sample from the previous population
    
        length_pop<-length(population[[t-1]][[m_star]])
        
        weights_pop<-weights[[t-1]][[m_star]]
        
        #this samples the index
        which_theta_star <- sample(1:length_pop, 1, 
                             prob = weights_pop )
        
        #extract correct sample
        theta_star <- population[[t-1]][[m_star]][[which_theta_star]]

        # Perturb theta* to get theta**
        theta_star_star <-rK_fn(m_star,theta_star)
        
      }
      # Need to wrap the above in a while look to make 
      # sure the parameters are appropriate: 
      # i.e. make sure pi(theta_star_star) > 0 
      
      # Data generation step
      # Select which model based on m_star
      #obtain the correct model for m_star
      current_model<-get(paste('model.',m_star,sep=''))
      
      #make a named parameter vector and state vector
      par_vec<-make_par_vec(t,m_star,theta_star_star)
      y_vec<-make_y_vec(m_star,theta_star_star)
    
      #the model is ran by calling a function
      #because model 2 is run with dede and a different time vector
      dfs<-run_model(m_star,current_model,y_vec,par_vec)
    
      # Only proceed if lsoda didn't throw and error
      if(NROW(dfs) == (NROW(toni09_tdc_data$I))){
        
        # Find the distance (sse) from the 
        # data (only I and R)
        
        x_star <- dfs[,c('I','R')]
        x0 <- toni09_tdc_data[,c('I','R')]

        distance <- sse(x0, x_star)
        
        # to catch model 2 chucking out NAs
        if(is.na(distance)) distance<-Inf

      }
    }
    # Set the current particle to theta_star_star
    #appending to the list
    population[[t]][[m_star]][[length(population[[t]][[m_star]])+1]]<-theta_star_star
    m_vec[[t]][[length(m_vec[[t]])+1]]<-m_star
    
    # Calculate the weight for this particle
    if (t == 1){
      weights[[t]][[m_star]] <- c(unlist(weights[[t]][[m_star]]),1)
    }else{
      pi_theta <- pi_theta_fn(m_star,theta_star_star)
      w <- weights[[t-1]][[m_star]]
      thetas <- population[[t-1]][[m_star]]
      Ks <- sapply(thetas, function(x) dK_fn(m_star,x, theta_star_star))
      denominator <- sum(w * Ks)
      

      weights[[t]][[m_star]] <- c(unlist(weights[[t]][[m_star]]),pi_theta / denominator)
    }
    cat("*", t, "/", i)
    i <- i + 1
  } # while(i <= N.particles)
   
  # Record the number of data generation steps taken
  data.gen.steps <- c(data.gen.steps, data.gen.counter)

    # Normalize the weights for all m; increment t
  for(m in 1:4){
    if(length(weights[[t]][[m]])==0)  weights[[t]][[m]] <- 1
    else weights[[t]][[m]] <- weights[[t]][[m]]/sum(weights[[t]][[m]], na.rm = TRUE)
  }
  if(sum(sapply(x<-1:4,function(x) length(population[[t]][[x]])))!=N.particles) break
  
  t <- t + 1
  
  }


if(DO_PLOTTING){
  
  #Figure 7
  no_m<-4
  par(mfrow=c(3,5))
  for(i in 1:N.populations) barplot(get_count(no_m,i,population),ylim=c(0,N.particles))


}

#Figure 8 (when all epislons have been run)
#model 3
final_m3<-matrix(unlist(population[[N.populations]][[3]]),ncol=4,byrow=T)
par(mfrow=c(2,2))
hist(final_m3[,1],main='S0')
hist(final_m3[,2],main='gamma')
hist(final_m3[,3],main='v')
hist(final_m3[,4],main='delta')


