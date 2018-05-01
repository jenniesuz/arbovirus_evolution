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
#****************With fixed time to apoptosis and budding delay*************************************

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
  
  part.1 <- 1 / 2*(1+beta_m*Shat*f.alpha_m*rho_m*(gamma_m + b_delay_m*lambda_m) ) 
  
  part.2 <- - mu_mc - mu_mv - beta_m*Shat*lambda_m*b_delay_m + beta_m*Shat*rho_m*(gamma_m - f.alpha_m*mu_mc*gamma_m + f.alpha_m*lambda_m + b_delay_m*lambda_m)
  
  part.3 <-  4*(1+beta_m*Shat*f.alpha_m*rho_m*(gamma_m+b_delay_m*lambda_m)) *(-beta_m*Shat*lambda_m + mu_mc*mu_mv - beta_m*Shat*rho_m*(mu_mc*gamma_m-lambda_m))
  
  part.4 <-  mu_mc+mu_mv+beta_m*Shat*b_delay_m*lambda_m - beta_m*Shat*rho_m*(gamma_m - f.alpha_m*mu_mc*gamma_m + f.alpha_m*lambda_m + b_delay_m*lambda_m) 
  
  fit1 <- part.1*(part.2+sqrt(part.4^2 - part.3) )
  return(fit1)
}
#********************************vectorize************************************

mutfit.v <- Vectorize(mutfit)         # enable multiple values to be input for a parameter
mutfit.delay.v <- Vectorize(mutfit.delay)         # enable multiple values to be input for a parameter
#**************************************************************************


#***************acute resident virus, range of budding rates for mutant persistent virus*****
lam <- seq(1,100,0.1) 

lam.mut <- mutfit.v(lambda_m=lam,alpha_m=0,gamma_m=0,lambda=0,gamma=2400,alpha=1/24)
lam.mut <- cbind.data.frame(lam,lam.mut)

lam.delay <- mutfit.delay.v(lambda_m=lam
                              ,f.alpha_m=1000
                              ,gamma_m=0
                              ,b_delay_m=2
                              ,lambda=0
                              ,gamma=2400
                              ,f.alpha=24)
lam.mut$delay <- lam.delay
#*************acute resident virus, range of budding delays for mutant persistent virus****
delay <- seq(1,500,0.1)
delay.mut <- mutfit.delay.v(lambda_m=50
                            ,f.alpha_m=1000
                            ,gamma_m=0
                            ,b_delay_m=delay
                            ,lambda=0
                            ,gamma=2400
                            ,f.alpha=24)
delay.mut <- cbind.data.frame(delay,delay.mut)
#*****************Plot both***************************
#pdf(file="fig_ms_2.pdf",width=5,height=4)
ggplot(lam.mut, aes(x=lam,y=delay)) +
  geom_line(aes(x=lam,y=delay,color="Delay")) +
  geom_line(aes(x=lam,y=lam.mut,color="No delay"),linetype=2) +
  scale_color_manual(name="",values=c("No delay"="black","Delay"="black")) +
  guides(color=guide_legend(override.aes=list(linetype=c(1,2)))) +
  labs(y="Mutant virus fitness", x=expression(paste("Mutant virus budding rate (",italic(lambda)["m"],")")),")")
#dev.off()
#*********************************************************************************************
ggplot(delay.mut, aes(x=delay,y=delay.mut)) +
  geom_line(aes(x=delay,y=delay.mut)) +
  labs(y="Mutant virus fitness"
       , x=expression(paste("Mutant virus budding delay (",italic(tau),"')")))
#*******************************************************




# re-write mutant virus fitness function for input to unitroot.all function
  
#*************Without delay in budding and fixed time to apoptosis*****************
fun1 <- function(x
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

#*************With delay in budding and fixed time to apoptosis*****************
fun2 <- function(x
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

#******************************************************************************************

lambda_m <-  seq(1,100,0.1)       # range of values for mutant virus probabillity of susceptible cell infection
fit.bound1 <- sapply(lambda_m,function(y){
  uniroot.all(fun1,c(0,10^3*24),lambda_m=y,gamma="NA")
})

fit.bound2 <- sapply(lambda_m,function(y){
  return(uniroot.all(fun2,c(1,10^3*24),lambda_m=y,gamma="NA"))
})

lam.roots <- cbind.data.frame(lambda_m,"fit.bound1"=fit.bound1/24,"fit.bound2"=as.numeric(fit.bound2)/24)

pdf(file="fig_ms_2.pdf",width=5,height=4)
ggplot(lam.roots, aes(x=fit.bound1,y=lambda_m)) +
  geom_line(aes(x=fit.bound1,y=lambda_m,color="Model 1"),linetype=1) +
  geom_line(aes(x=fit.bound2,y=lambda_m,color="Model 2"),linetype=2) +
  scale_color_manual(name="",values=c("Model 1"="black","Model 2"="black")) +
  guides(color=guide_legend(override.aes=list(linetype=c(1,2)))) +
  labs(y=expression(paste("Mutant virus budding rate (",italic(lambda)["m"],")"))
                          , x=expression(paste("Resident virus yeild at apoptosis/time to apoptosis (",italic(gamma),"/",italic(tau),")")))
dev.off()


