model{

  for(i in 1:N){
    
    (a_g[1]*y_exo_chla + a_g[2]*y_ctd + a_g[3]*y_fp_g + a_g[4]*y_count_g) / 4 ~ dnorm(mu_g[i],tau_obs_g)
    
    (a_c[1]*y_exo_chla + a_c[2]*y_exo_phy + a_c[3]*y_ctd + a_c[4]*y_fp_c + a_c[5]*y_count_c) / 4 ~ dnorm(mu_c[i],tau_obs_c)
    
    (a_d[1]*y_exo_chla + a_d[2]*y_ctd + a_d[3]*y_fp_d + a_d[4]*y_count_d) / 4 ~ dnorm(mu_d[i],tau_obs_d)
    
    (a_f[1]*y_exo_chla + a_f[2]*y_ctd + a_f[3]*y_fp_f + a_f[4]*y_count_f) / 4 ~ dnorm(mu_f[i],tau_obs_f)

  }
  
  #### Priors
  tau_obs_g ~ dgamma(0.001,0.001)
  tau_obs_c ~ dgamma(0.001,0.001)
  tau_obs_d ~ dgamma(0.001,0.001)
  tau_obs_f ~ dgamma(0.001,0.001)
  
  a_g ~ dmnorm(as.vector(c(0,0,0,0)),solve(diag(1E-03,4)))
  a_c ~ dmnorm(as.vector(c(0,0,0,0,0)),solve(diag(1E-03,5)))
  a_d ~ dmnorm(as.vector(c(0,0,0,0)),solve(diag(1E-03,4)))
  a_f ~ dmnorm(as.vector(c(0,0,0,0)),solve(diag(1E-03,4)))

}