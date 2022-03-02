model{
  
  for(i in 1:N){
    
    y_fp[i] ~ dnorm(mu[i],tau_obs)

    mu[i] <- beta[1] + beta[2]*ctd[i] 
    
    ctd[i] <- ifelse(index_ctd_obs[i]=0,y_ctd[i],last_ctd_obs)
    last_ctd_obs <- ifelse(index_ctd_obs[i]=0,y_ctd[i],last_ctd_obs)
    
  }
  
  #### Priors
  tau_obs ~ dgamma(0.001,0.001)
  beta ~ dmnorm(as.vector(c(0,0,0,0)),solve(diag(1E-03,4)))
  
}

