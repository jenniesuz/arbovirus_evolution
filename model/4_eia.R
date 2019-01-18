library(rootSolve) 
#******************Evolutionary invasion analysis**************************

#***********************COMBINED MODEL ACUTE AND PERSISTENT INFECTION*************************************

#************* mutant virus fitness function without delay in budding or fixed time to apoptosis************
mutfit <- function(mu_mv=0.1                  # mutant virus clearance rate
                         ,mu_mc = 1/120       # mutant virus cell death rate
                         ,beta_m = 10^-9      # mutant probability of infecting susceptible cells
                         ,gamma_m = 2400      # mutant virus yeild at apoptosis
                         ,lambda_m = 100
                         ,alpha_m = 1/24      # mutant virus apoptosis rate
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
                         ,beta_m = 10^-9     # mutant probability of infecting susceptible cells
                         ,gamma_m = 2400      # mutant virus yeild at apoptosis
                         ,lambda_m = 100
                         ,f.alpha_m = 24      # fixed mutant virus time to apoptosis
                         ,b_delay_m = 2           # delay in budding
                         ,mu_v=0.1            # resident virus clearance rate  
                         ,mu_c = 1/120        # resident infected cell death rate 
                         ,beta = 10^-9        # resident probability of infecting susceptible cells
                         ,f.alpha = 24          # resident virus fixed time to apoptosis
                         ,gamma = 2400 
                         ,lambda = 100
                         ,b_delay=2){    
  rho_m <- exp(-mu_mc*f.alpha_m)
  rhob_m <- exp(-mu_mc*b_delay_m)
  rho  <- exp(-mu_c*f.alpha)
  rhob <- exp(-mu_c*b_delay)
  
  #Shat <- (mu_v / ( beta*( (lambda - lambda*exp(-mu_c*f.alpha))/mu_c + gamma*exp(-mu_c*f.alpha) ))  ) 
  Shat <- ( mu_c*mu_v*(0.03-mu_c) ) / ( beta*gamma*rho*mu_c*(0.03+mu_c) + ( beta*lambda*rhob*(0.03*rho-mu_c*rho-0.03-mu_c) ) )
  #part.1 <- (2*(1+beta_m*Shat*f.alpha_m*rho_m*(gamma_m + b_delay_m*lambda_m)) )
  
  #part.2 <- (- mu_mc - mu_mv - beta_m*Shat*lambda_m*b_delay_m + beta_m*Shat*rho_m*(gamma_m - f.alpha_m*mu_mc*gamma_m + f.alpha_m*lambda_m + b_delay_m*lambda_m))
  
  #part.3 <-  (-4*(1+beta_m*Shat*f.alpha_m*rho_m*(gamma_m+b_delay_m*lambda_m)) *(-beta_m*Shat*lambda_m + mu_mc*mu_mv - beta_m*Shat*rho_m*(mu_mc*gamma_m-lambda_m)))
  
  #part.4 <-  (mu_mc+mu_mv+beta_m*Shat*b_delay_m*lambda_m - beta_m*Shat*rho_m*(gamma_m - f.alpha_m*mu_mc*gamma_m + f.alpha_m*lambda_m + b_delay_m*lambda_m)) 
  
  #fit1 <- (part.2+sqrt(part.4^2 + part.3) )/part.1
  
  a <- (1 + beta_m*Shat*f.alpha_m*rho_m*(gamma_m+b_delay_m*rhob_m*lambda_m))
  b <- (mu_mc + mu_mv + beta_m*Shat*(-gamma_m*rho_m + gamma_m*rho_m*f.alpha_m*mu_mc + lambda_m*b_delay_m*rhob_m - rho_m*rhob_m*b_delay_m - lambda_m*f.alpha_m*rhob_m*rho_m))
  c <- (mu_mc*mu_mv + beta_m*Shat*(- gamma_m*rho_m*mu_mc - lambda_m*rhob_m + rho_m*lambda*rhob_m))
  fit1 <- (- b + sqrt(b^2 - 4*a*c))/2*a
 #return(fit1)

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
                ,gamma_m = 2400      # mutant virus yeild at apoptosis
                ,lambda_m = 100
                ,alpha_m = 1/24      # mutant virus apoptosis rate
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
                ,gamma_m = 0      # mutant virus yeild at apoptosis
                ,lambda_m = 100
                ,f.alpha_m = 2000      # mutant virus apoptosis rate
                ,b_delay_m = 12  
                ,mu_v=0.1            # resident virus clearance rate  
                ,mu_c = 1/120        # resident infected cell death rate 
                ,beta = 10^-6        # resident probability of infecting susceptible cells
                ,f.alpha = 24        # resident virus apoptosis rate
                ,gamma = 2400
                ,lambda = 0
                ,b_delay = 2
){
  
  if(gamma=="NA"){
    
    gamma_m <- 0
    f.alpha_m <- 1000
    rho_m <- exp(-mu_mc*f.alpha_m)
    
    Shat <- mu_v / ( beta*( (lambda - lambda*exp(-mu_c*f.alpha))/mu_c + x*exp(-mu_c*f.alpha) ))   
    
    part.1 <- 1 / 2*(1+beta_m*Shat*f.alpha_m*rho_m*(gamma_m + b_delay_m*lambda_m) ) 
    
    part.2 <- - mu_mc - mu_mv - beta_m*Shat*lambda_m*b_delay_m + beta_m*Shat*rho_m*(gamma_m - f.alpha_m*mu_mc*gamma_m + f.alpha_m*lambda_m + b_delay_m*lambda_m)
    
    part.3 <-  4*(1+beta_m*Shat*f.alpha_m*rho_m*(gamma_m+b_delay_m*lambda_m)) *(-beta_m*Shat*lambda_m + mu_mc*mu_mv - beta_m*Shat*rho_m*(mu_mc*gamma_m-lambda_m))
    
    part.4 <-  mu_mc+mu_mv+beta_m*Shat*b_delay_m*lambda_m - beta_m*Shat*rho_m*(gamma_m - f.alpha_m*mu_mc*gamma_m + f.alpha_m*lambda_m + b_delay_m*lambda_m) 
    fit1 <- part.1*(part.2+sqrt(part.4^2 - part.3) )
    return(fit1)
    
    }
}


