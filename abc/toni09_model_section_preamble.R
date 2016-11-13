
  model.1 <- function( t, x, parameters ) { 
    ##Inputs:
    #t time 
    #x state variables
    # parameters parameter vector
    S <- x[1]
    I <- x[2]
    R <- x[3]
    #The with statement means we can use the parameter names as in parameters
    with( as.list(parameters), { 
      dS <-alpha-gamma*S*I-d*S
      dI <- gamma*S*I-v*I-d*I
      dR <-v*I-d*R
      res <- c(dS,dI,dR)
      list(res)
    }
    )
  }
  
   # out <- lsoda(func =model.1,
   #                  y = c(S=84,I=1,R=0),
   #                  parms = c(alpha=0,gamma=0.005,d=0,v=0.25),
   #                  atol = 1e-3,
   #                  times = seq(1, 21,by=1))
  
  model.2 <- function( t, x, parameters ) {
    ##Inputs:
    #t time
    #x state variables
    # parameters parameter vector
    S <- x[1]
    I <- x[2]
    R <- x[3]
    #The with statement means we can use the parameter names as in parameters
    with( as.list(parameters), {
      if(t < tau){
        I.lag<-0
      }else{
        I.lag<-lagvalue(t-tau,2)
      }
      dS <-alpha-gamma*S*I.lag-d*S
      dI <- gamma*S*I.lag-v*I-d*I
      dR <-v*I-d*R
      res <- c(dS,dI,dR)
      list(res)
    }
    )
  }

    
      # out <-dede(y = c(S=38,I=1,R=0),times = seq(-1, 21,by=1), func = model.2,
      #                parms = c(alpha=0,gamma=0.025,d=0,v=0.2,tau=1))
      #                 
    
  #this is model 3
  model.3 <- function( t, x, parameters ) { 
    ##Inputs:
    #t time 
    #x state variables
    # parameters parameter vector
    S <- x[1]
    L <- x[2]
    I <- x[3]
    R <- x[4]
    #The with statement means we can use the parameter names as in parameters
    with( as.list(parameters), { 
      dS <-alpha-gamma*S*I-d*S
      dL<- gamma*S*I-delta*L-d*L
      dI <- delta*L-v*I-d*I
      dR <-v*I-d*R
      res <- c(dS,dL,dI,dR)
      list(res)
    }
    )
  }
  
  # out <- lsoda(func =model.3,
  #                  y = c(S=84,L=0,I=1,R=0),
  #                  parms = c(alpha=0,gamma=0.005,delta=4,d=0,v=0.25),
  #                  atol = 1e-3,
  #                  times = seq(1, 21,by=1))

  model.4 <- function( t, x, parameters ) { 
    ##Inputs:
    #t time 
    #x state variables
    # parameters parameter vector
    S <- x[1]
    I <- x[2]
    R <- x[3]
    #The with statement means we can use the parameter names as in parameters
    with( as.list(parameters), { 
      dS <-alpha-gamma*S*I-d*S+e*R
      dI <- gamma*S*I-v*I-d*I
      dR <-v*I-(d+e)*R
      res <- c(dS,dI,dR)
      list(res)
    }
    )
  }
  
   # out <- lsoda(func =model.4,
   #                  y = c(S=84,I=1,R=0),
   #                  parms = c(alpha=0,gamma=0.005,d=0,v=0.25,e=0.01),
   #                  atol = 1e-3,
   #                  times = seq(1, 21,by=1))
  
  theta_star_fn<-function(m_star){
    #Function to sample from prior vector associated with m_star
    #need a catch here in case m_star is out of bounds
    
    #Input
    #m_star model number 1, 2, 3 or 4
    
    #Sample parameters in every model
    S0_star<-ceiling(runif(1,36,100))
    gamma_star<-runif(1,0,3)
    v_star<-runif(1,0,3)

    #Additional parameters
    theta_star<-c(S0=S0_star,gamma=gamma_star,v=v_star)
    
    if(m_star==2){
      tau_star<-ceiling(runif(1,-0.5,5))
      theta_star<-c(theta_star,tau=tau_star)
    }else if(m_star==3){
      delta_star<-runif(1,-0.5,5)
      theta_star<-c(theta_star,delta=delta_star)  
    }else if(m_star==4){
      e_star<-runif(1,-0.5,5)
      theta_star<-c(theta_star,e=e_star)  
    }
    return(theta_star)
  }
  
  
  rK_fn <- function(m_star,theta_star){
    #Function to perturb values 
    mean<-theta_star

    #all additional paramters have sigma=1
    #so we only have two possible sigma vectors
    
    if(m_star==1){
      sigma<-c(3,0.3,0.3)
    }else if(m_star!=1){
      sigma<-c(3,0.3,0.3,1)
    }
    
    sigma_diag<-sigma^2*diag(length(sigma))
    
    perturbed<-rmvnorm(1, mean = mean, sigma = sigma_diag)
    perturbed[1]<-ceiling(perturbed[1])
    
    return(perturbed)
  }
  
  make_par_vec<-function(t,m_star,theta_star){
    #Function to make a named parameter vector
    if(t==1) par_vec<-c(theta_star[-1],alpha=0,d=0)
    else{
      par_vec<-c(gamma=theta_star[2],v=theta_star[3],alpha=0,d=0)
      
      if(m_star==2){
        par_vec<-c(gamma=theta_star[2],v=theta_star[3],tau=theta_star[4],alpha=0,d=0)
        
      }else if(m_star==3){
        par_vec<-c(gamma=theta_star[2],v=theta_star[3],delta=theta_star[4],alpha=0,d=0)
        
      }else if(m_star==4){
        par_vec<-c(gamma=theta_star[2],v=theta_star[3],e=theta_star[4],alpha=0,d=0)
        
      }
    }
    return(par_vec)
  }
  
 
  make_y_vec<-function(m_star,theta_star){
    #Function to make the y vector a model 3 has an additional state
    
    y_vec<-c(S=theta_star_star[1],I=1,R=0)
    
    if(m_star==3) y_vec<-c(S=theta_star_star[1],L=0,I=1,R=0)
    
    return(y_vec)
  }  
  
  run_model<-function(m_star,current_model,y_vec,par_vec){
    #Function to run model with dede or lsoda
    if(m_star==2){
      out<-tryCatch(dede(y = y_vec,times = seq(-par_vec[3], 21,by=1), func = current_model,
                                  parms =  par_vec),
          silent = TRUE,error=function(err) return(matrix(NA,ncol=3,nrow=5)))
      
      upper<-par_vec[3]+1
      dfs <- as.data.frame(out[-seq(1:upper),])
    }else{
    out<-tryCatch(lsoda(func = current_model,
                     y = y_vec,
                     parms = par_vec,
                     times = seq(1, 21,by=1)),
        silent = TRUE,error=function(err) return(matrix(NA,ncol=3,nrow=5)))
      dfs <- as.data.frame(out)
    }
   
  }
  
  dK_fn<-function(m_star,theta1,theta2){
    #Function to calculate density of values
    sig<-c(3,0.3,0.3)
     
    #all additional parameters have sigma=1
   if(m_star!=1) sig<-c(sig,1)
      
    sigma_diag<-sig^2*diag(length(sig))
    
    return(dmvnorm(theta1, mean = theta2, sigma = sigma_diag))
    
  }
  
  #dK_fn(4,population[[1]][[4]][[3]],population[[2]][[4]][[1]])
  
  pi_theta_fn<-function(m_star,theta_star_star){
    #Function to calculate the prior weight of current parameters
    #Multiplied
    
    #need to do by name to avoid incorrect calculations really...
    
    prior_S0<-dunif(theta_star_star[1],36,100)
    
    prior_gamma<-dunif(theta_star_star[2],0,3)
 
    prior_v<-dunif(theta_star_star[3],0,3)
    
    
    prior_theta<-prior_S0*prior_gamma*prior_v
    
    if(m_star==2){
      prior_tau<-dunif(theta_star_star[4],-0.5,5)
      prior_theta<-prior_theta*prior_tau
    }else if(m_star==3){
      prior_delta<-dunif(theta_star_star[4],-0.5,5)
      prior_theta<-prior_theta*prior_delta 
    }else if(m_star==4){
      prior_e<-dunif(theta_star_star[4],-0.5,5)
      prior_theta<-prior_theta*prior_e
    }
    
    return(prior_theta)
  }
  
  
  sse <- function(x, xd){
    # Define sum of squared errors
    
    return( sqrt(sum((x[,1] - xd[,1])^2+(x[,2] - xd[,2])^2)))
  }
  
  
  get_count<-function(no_m,t,population){
    #To get the count of m for each t
    #used to make bar plot
    
    #make a vector from 1 to the number of models
    m_val<-1:no_m
    #how many particles were accepted for that model
    count<-sapply(m_val,function(m_val) length(population[[t]][[m_val]]))
    return(count)
  }
