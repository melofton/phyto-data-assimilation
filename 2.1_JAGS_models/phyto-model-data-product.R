model{
  
  for(i in 1:N){
    
    (a[1]*y_exo_chla + a[2]*y_ctd + a[3]*y_fp + a[4]*y_count) / 4 ~ dnorm(mu[i],tau_obs)
    
  }
  
  #### Priors
  tau_obs ~ dgamma(0.001,0.001)
  
  a ~ dmnorm(as.vector(c(0,0,0,0)),solve(diag(1E-03,4)))
  
}