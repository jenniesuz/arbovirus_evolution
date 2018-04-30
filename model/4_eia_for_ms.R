library(rootSolve) 
library("ggplot2")
#******************Evolutionary invasion analysis**************************

#***********************COMBINED MODEL ACUTE AND PERSISTENT INFECTION*************************************

#************* mutant virus fitness function without delay in budding or fixed time to apoptosis************
mutfit <- function(mu_mv=0.1            # mutant virus clearance rate
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
                         ,lambda = 100){     # resident virus yeild at apoptosis 
  fit1 <- (-(alpha_m+mu_mc+mu_mv) + sqrt( (alpha_m+mu_mc+mu_mv)^2 - 4*(alpha_m*mu_mv + mu_mv*mu_mc - (beta_m*(gamma_m*alpha_m+lambda_m)*mu_c*(mu_v + alpha) )/((gamma*alpha + lambda)*beta)) ) ) / 2
  fit2 <- (-(alpha_m+mu_mc+mu_mv) - sqrt( (alpha_m+mu_mc+mu_mv)^2 - 4*(alpha_m*mu_mv + mu_mv*mu_mc - (beta_m*(gamma_m*alpha_m+lambda_m)*mu_c*(mu_v + alpha))/((gamma*alpha + lambda)*beta) ) ) ) / 2
  return(fit1)
}
#******************************************************************************************

mutfit.v <- Vectorize(mutfit)         # enable multiple values to be input for a parameter

lam <- seq(10,500,0.1) 
lam.mut.a24 <- mutfit.v(lambda_m=lam,alpha_m=0,gamma_m=0,lambda=0,gamma=2400,alpha=1/24)
mut.fitness.dat <- cbind.data.frame(lam,lam.mut.a24)

#*****************Plot without delay/ fixed time***************************
#pdf(file="fig_ms_2.pdf",width=5,height=4)
ggplot(mut.fitness.dat, aes(x=lam,y=lam.mut.a24)) +
  geom_line(aes(x=lam,y=lam.mut.a24)) +
  labs(y="Mutant virus fitness", x=expression(paste("Mutant virus budding rate (",italic(lambda)["m"],")")),")")

#*****************************************************************

mutfit.delay <- function(mu_mv=0.1            # mutant virus clearance rate
                   ,mu_mc = 1/120       # mutant virus cell death rate
                   ,beta_m = 10^-6      # mutant probability of infecting susceptible cells
                   ,gamma_m = 2400      # mutant virus yeild at apoptosis
                   ,lambda_m = 100
                   ,f.alpha_m = 24      # fixed mutant virus time to apoptosis
                   ,b_delay_m = 2           # delay in budding
                   ,mu_v=0.1            # resident virus clearance rate  
                   ,mu_c = 1/120        # resident infected cell death rate 
                   ,beta = 10^-6        # resident probability of infecting susceptible cells
                   ,f.alpha = 24          # resident virus fixed time to apoptosis
                   ,gamma = 2400 
                   ,lambda = 100
                   ,b_delay=2){    
  rho_m <- exp(-mu_mc*f.alpha_m)
  
  Shat <- mu_v / ( beta*( (lambda - lambda*exp(-mu_c*f.alpha))/mu_c + gamma*exp(-mu_c*f.alpha) ))   
  
  part.1 <- (1 / 2*(1+beta_m*Shat*f.alpha_m*rho_m*(gamma_m + b_delay_m*lambda_m) ) )
  
  part.2 <- ( -mu_mc - mu_mv - beta_m*Shat*lambda_m*b_delay_m + beta_m*Shat*rho_m*(gamma_m - f.alpha_m*mu_mc*gamma_m + f.alpha_m*lambda_m + b_delay_m*lambda_m) )
  
  part.3 <-  -4*(1+beta_m*Shat*f.alpha_m*rho_m*(gamma_m+b_delay_m*lambda_m)) *(-beta_m*Shat*lambda_m + mu_mc*mu_mv - beta_m*Shat*rho_m*(mu_mc*gamma_m-lambda_m))
  
  part.4 <- ( mu_mc+mu_mv+beta_m*Shat*b_delay_m*lambda_m - beta_m*Shat*rho_m*(gamma_m - f.alpha_m*mu_mc*gamma_m + f.alpha_m*lambda_m + b_delay_m*lambda_m) )
  
  fit1 <- part.1*(part.2+sqrt(part.3 + part.4^2) )
  return(fit1)
}

#************calculate for without delay/ fixed time*********************
mutfit.delay.v <- Vectorize(mutfit.delay)         # enable multiple values to be input for a parameter

lam <- seq(10,500,0.1) 
lam.mut.a24 <- mutfit.delay.v(lambda_m=lam,f.alpha_m=24,gamma_m=0,b_delay_m=12,lambda=0,gamma=2400,f.alpha=24)
mut.fitness.dat$delay <- lam.mut.a24



#*****************Plot without delay/ fixed time***************************
#pdf(file="fig_ms_2.pdf",width=5,height=4)
ggplot(mut.fitness.dat, aes(x=lam,y=delay)) +
  geom_line(aes(x=lam,y=lam.mut.a24)) 
  labs(y="Mutant virus fitness", x=expression(paste("Mutant virus budding rate (",italic(lambda)["m"],")")),")")
#dev.off()
#*********************************************************************************************





# re-write mutant virus fitness function for input to unitroot.all function
  
#*************Without delay in budding and fixed time to apoptosis*****************
fun <- function(x
                ,mu_mv=0.1            # mutant virus clearance rate
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


lambda_m <-  seq(1,100,0.1)# range of values for mutant virus probabillity of susceptible cell infection
fit.bound1 <- sapply(lambda_m,function(y){
  uniroot.all(fun,c(0,10^3*24),lambda_m=y,gamma="NA")
})

bounds <- cbind.data.frame(lambda_m,"fit.bound1"=fit.bound1/24)

bounds[,2]/bounds[,1] # 6 *

#pdf(file="fig_ms_2.pdf",width=5,height=4)
ggplot(bounds, aes(x=fit.bound1,y=lambda_m)) +
  geom_line(aes(x=fit.bound1,y=lambda_m)) +
  labs(y=expression(paste("Mutant virus budding rate ( ",lambda["m"],")" )), x=expression(paste("Resident virus yeild at apopotosis * apopotosis rate (",gamma, alpha,")")))
#dev.off()

#*************With delay in budding and fixed time to apoptosis*****************
fun <- function(x
                ,mu_mv=0.1            # mutant virus clearance rate
                ,mu_mc = 1/120       # mutant virus cell death rate
                ,beta_m = 10^-6      # mutant probability of infecting susceptible cells
                ,gamma_m = 2400      # mutant virus yeild at apoptosis
                ,lambda_m = 100
                ,f.alpha_m = 24      # mutant virus apoptosis rate
                ,b_delay_m = 2  
                ,mu_v=0.1            # resident virus clearance rate  
                ,mu_c = 1/120        # resident infected cell death rate 
                ,beta = 10^-6        # resident probability of infecting susceptible cells
                ,f.alpha = 24        # resident virus apoptosis rate
                ,gamma = 2400
                ,lambda = 100
                ,b_delay = 2
){
  
  if(gamma=="NA"){
    
    gamma_m <- 0
    f.alpha_m <- 0
    rho_m <- exp(-mu_mc*f.alpha_m)
    
    Shat <- mu_v / ( beta*( (lambda - lambda*exp(-mu_c*f.alpha))/mu_c + gamma*exp(-mu_c*f.alpha) ))   
    
    part.1 <- (1 / 2*(1+beta_m*Shat*f.alpha_m*rho_m*(gamma_m + b_delay_m*lambda_m) ) )
    
    part.2 <- ( -mu_mc - mu_mv - beta_m*Shat*lambda_m*b_delay_m + beta_m*Shat*rho_m*(gamma_m - f.alpha_m*mu_mc*gamma_m + f.alpha_m*lambda_m + b_delay_m*lambda_m) )
    
    part.3 <-  -4*(1+beta_m*Shat*f.alpha_m*rho_m*(gamma_m+b_delay_m*lambda_m)) *(-beta_m*Shat*lambda_m + mu_mc*mu_mv - beta_m*Shat*rho_m*(mu_mc*gamma_m-lambda_m))
    
    part.4 <- ( mu_mc+mu_mv+beta_m*Shat*b_delay_m*lambda_m - beta_m*Shat*rho_m*(gamma_m - f.alpha_m*mu_mc*gamma_m + f.alpha_m*lambda_m + b_delay_m*lambda_m) )
    
    return(part.1*(part.2+sqrt(part.3 + part.4^2) ))
    
    }
}


lambda_m <-  seq(1,100,0.1)# range of values for mutant virus probabillity of susceptible cell infection
fit.bound2 <- sapply(lambda_m,function(y){
  uniroot.all(fun,c(1,10^2),lambda_m=y,gamma="NA")
})

bounds <- cbind.data.frame(lambda_m,"fit.bound2"=fit.bound2)

bounds[,2]/bounds[,1] # 6 *

#pdf(file="fig_ms_2.pdf",width=5,height=4)
ggplot(bounds, aes(x=fit.bound2,y=lambda_m)) +
  geom_line(aes(x=fit.bound2,y=lambda_m)) +
  labs(y=expression(paste("Mutant virus budding rate ( ",lambda["m"],")" )), x=expression(paste("Resident virus yeild at apopotosis * apopotosis rate (",gamma, alpha,")")))
#dev.off()


