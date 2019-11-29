library(rootSolve) 
#******************Evolutionary invasion analysis**************************

#***********************COMBINED MODEL - ACUTE AND PERSISTENT INFECTION*************************************

#************* mutant virus fitness function without delay in budding or fixed time to apoptosis************
mutfit <- function(mu_mv=0.1                  # mutant virus clearance rate
                         ,mu_mc = 1/120       # mutant virus cell death rate
                         ,beta_m = 10^-9      # mutant probability of infecting susceptible cells
                         ,gamma_m = 0      # mutant virus yeild at apoptosis
                         ,lambda_m = 100
                         ,alpha_m =0      # mutant virus apoptosis rate
                         ,mu_v=0.1            # resident virus clearance rate  
                         ,mu_c = 1/120        # resident infected cell death rate 
                         ,beta = 10^-9        # resident probability of infecting susceptible cells
                         ,alpha = 1/24        # resident virus apoptosis rate
                         ,gamma = 2400 
                         ,lambda = 100){      # resident virus yeild at apoptosis 
  fit1 <- (-(alpha_m+mu_mc+mu_mv) + sqrt( (alpha_m+mu_mc+mu_mv)^2 - 4*(alpha_m*mu_mv + mu_mv*mu_mc - (beta_m*(gamma_m*alpha_m+lambda_m)*mu_c*(mu_v + alpha) )/((gamma*alpha + lambda)*beta)) ) ) / 2
  fit2 <- (-(alpha_m+mu_mc+mu_mv) - sqrt( (alpha_m+mu_mc+mu_mv)^2 - 4*(alpha_m*mu_mv + mu_mv*mu_mc - (beta_m*(gamma_m*alpha_m+lambda_m)*mu_c*(mu_v + alpha))/((gamma*alpha + lambda)*beta) ) ) ) / 2
  return(fit1)
}
#******************************************************************************************

#****************With fixed time to apoptosis and budding delay*************************************

mutfit.delay <- function(mu_mv=0.1            # mutant virus clearance rate
                         ,mu_mc = 1/120       # mutant virus cell death rate
                         ,beta_m = 10^-9      # mutant probability of infecting susceptible cells
                         ,gamma_m = 2400     # mutant virus yeild at apoptosis
                         ,lambda_m = 100      # mutant virus budding rate
                         ,tau_a_m = 24      # mutant virus fixed time to apoptosis
                         ,tau_b_m = 2       # delay in budding
                         ,mu_v=0.1            # resident virus clearance rate  
                         ,mu_c = 1/120        # resident infected cell death rate 
                         ,beta = 10^-9        # resident probability of infecting susceptible cells
                         ,tau_a = 24        # resident virus fixed time to apoptosis
                         ,gamma = 20
                         ,lambda = 1        
                         ,tau_b=2
                         ,r=0.03){            # susceptible cell growth rate  
  sigma_m <- exp(-mu_mc*tau_a_m)
  sigmab_m <- exp(-mu_mc*tau_b_m)
  sigma  <- exp(-mu_c*tau_a)
  sigmab <- exp(-mu_c*tau_b)
  
  Shat <-   mu_c*( (mu_c / (lambda*beta*(1-sigma)*sigmab)) + (1/(gamma*beta*sigma) )      )

  a <- (1  + gamma_m*beta_m*Shat*sigma_m*tau_a_m + beta_m*Shat*sigma_m*tau_a_m*lambda_m*sigmab_m*tau_b_m)
  b <- (mu_mc + mu_mv - gamma_m*beta_m*Shat*sigma_m + gamma_m*beta_m*Shat*sigma_m*tau_a_m*mu_mc  + lambda_m*sigmab_m*tau_b_m*beta_m*Shat - lambda_m*sigmab_m*beta_m*Shat*sigma_m*tau_a_m  - lambda_m*sigmab_m*tau_b_m*beta_m*Shat*sigma_m)
  c <- (mu_mc*mu_mv - gamma_m*beta_m*Shat*sigma_m*mu_mc - lambda_m*sigmab_m*beta_m*Shat + lambda_m*sigmab_m*beta_m*Shat*sigma_m)
 fit1 <- (- b + sqrt(b^2 - 4*a*c))/ (2*a)
 
  
 return(fit1)

}
#********************************vectorize************************************
mutfit.v <- Vectorize(mutfit)                     # enable multiple values to be input for a parameter
mutfit.delay.v <- Vectorize(mutfit.delay)         # enable multiple values to be input for a parameter
#**************************************************************************



# re-write mutant virus fitness function for input to unitroot.all function
  
#*************Without delay in budding and fixed time to apoptosis*****************
fun.nodelay <- function(x
                ,mu_mv=0.1           # mutant virus clearance rate
                ,mu_mc = 1/120       # mutant virus cell death rate
                ,beta_m = 10^-6      # mutant probability of infecting susceptible cells
                ,gamma_m = 0      # mutant virus yeild at apoptosis
                ,lambda_m = 100
                ,alpha_m = 0      # mutant virus apoptosis rate
                ,mu_v=0.1            # resident virus clearance rate  
                ,mu_c = 1/120        # resident infected cell death rate 
                ,beta = 10^-6        # resident probability of infecting susceptible cells
                ,alpha = 1/24        # resident virus apoptosis rate
                ,gamma = 2400
                ,lambda = 100
){
  
  if(gamma=="NA"){
    
    return( ( -(0+mu_mc+mu_mv) + sqrt( (0+mu_mc+mu_mv)^2 - 4*(0*mu_mv + mu_mv*mu_mc - ( (beta_m*(0*0+lambda_m)*mu_c*(mu_v+alpha))/( (x*alpha+0)*beta) ) ) ) )  / 2 )
  }
}

#*************With delay in budding and fixed time to apoptosis*****************
fun.delay <- function(x
                ,mu_mv=0.1            # mutant virus clearance rate
                ,mu_mc = 1/120       # mutant virus cell death rate
                ,beta_m = 10^-6      # mutant probability of infecting susceptible cells
                ,gamma_m = 2400      # mutant virus yeild at apoptosis
                ,lambda_m = 1
                ,tau_a_m = 24       # mutant virus apoptosis 
                ,tau_b_m = 24 
                ,mu_v=0.1            # resident virus clearance rate  
                ,mu_c = 1/120        # resident infected cell death rate 
                ,beta = 10^-6        # resident probability of infecting susceptible cells
                ,gamma = 2400
                ,lambda = 1
                ,tau_a = 24        # resident virus apoptosis rate
                ,tau_b = 24
                ,r=0.03
){
  
  if(lambda=="NA"){
    
   # gamma <- 1
    sigma_m <- exp(-mu_mc*tau_a_m)
    sigmab_m <- exp(-mu_mc*tau_b_m)
    sigma  <- exp(-mu_c*tau_a)
    sigmab <- exp(-mu_c*tau_b)
    
    Shat <-   mu_c*( (mu_c / (x*beta*(1-sigma)*sigmab)) + (1/(gamma*beta*sigma) )      )
    
    a <- (1  + gamma_m*beta_m*Shat*sigma_m*tau_a_m + beta_m*Shat*sigma_m*tau_a_m*lambda_m*sigmab_m*tau_b_m)
    b <- (mu_mc + mu_mv - gamma_m*beta_m*Shat*sigma_m + gamma_m*beta_m*Shat*sigma_m*tau_a_m*mu_mc  + lambda_m*sigmab_m*tau_b_m*beta_m*Shat - lambda_m*sigmab_m*beta_m*Shat*sigma_m*tau_a_m  - lambda_m*sigmab_m*tau_b_m*beta_m*Shat*sigma_m)
    c <- (mu_mc*mu_mv - gamma_m*beta_m*Shat*sigma_m*mu_mc - lambda_m*sigmab_m*beta_m*Shat + lambda_m*sigmab_m*beta_m*Shat*sigma_m)
    fit1 <- (- b + sqrt(b^2 - 4*a*c))/ (2*a)
    
    
  }
  return(fit1)
}




mutfit.delay(mu_mv=0.1            # mutant virus clearance rate
                         ,mu_mc = 1/120       # mutant virus cell death rate
                         ,beta_m = 10^-9      # mutant probability of infecting susceptible cells
                         ,gamma_m = 200    # mutant virus yeild at apoptosis
                         ,lambda_m = 10      # mutant virus budding rate
                         ,tau_a_m = 24      # mutant virus fixed time to apoptosis
                         ,tau_b_m = 24       # delay in budding
                         ,mu_v=0.1            # resident virus clearance rate  
                         ,mu_c = 1/120        # resident infected cell death rate 
                         ,beta = 10^-9        # resident probability of infecting susceptible cells
                         ,gamma =200
                         ,lambda = 10 
                         ,tau_a = 24        # resident virus fixed time to apoptosis
                         ,tau_b=24
                         ,r=0.03)
