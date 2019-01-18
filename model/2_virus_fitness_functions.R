

#**********************************Virus fitness functions***********************************

#****************Model with budding and apoptosis, no delays*******************
fit.nodelay <- function(mu_v=0.1                # virus clearance rate
                        ,mu_i = 1/120       # infected cell death rate
                        ,beta = 10^-6       # probability of single virion infecting susceptible cell
                        ,lambda=10          # release rate of virus from infected cells
                        ,alpha=1/24         # apoptosis rate
                        ,gamma=50           # virus yield at apoptosis
                        ,S = 10^6){         # susceptible cells
  fit1 <- (- (mu_i+mu_v+alpha) + sqrt( (mu_i+mu_v+alpha)^2 - 4*(mu_v*(alpha+mu_i) - beta*S*(gamma*alpha+lambda) ) ))  / 2
  fit2 <- (- (mu_i+mu_v+alpha) - sqrt( (mu_i+mu_v+alpha)^2 - 4*(mu_v*(alpha+mu_i) - beta*S*(gamma*alpha+lambda) ) )) / 2
  return(fit1)
}
#*****allow to take multiple values for a parameter******
fit.nodelay.v <- Vectorize(fit.nodelay)
#***********************************************************************************

#****************Model with budding and apoptosis, with delays**********************
fit.delay <- function(mu_v=0.1              # virus clearance rate
                        ,mu_i = 1/120       # infected cell death rate
                        ,beta = 10^-6       # probability of single virion infecting susceptible cell
                        ,lambda=10          # release rate of virus from infected cells
                        ,gamma=50           # virus yield at apoptosis
                        ,tau.a=24           # time to apoptosis  
                        ,tau.b=24           # time to start of budding
                        ,S = 10^6){         # susceptible cells
  
  sigma <- exp(-mu_i*tau.a)
  
  part.1 <- ( 2*(1 + beta*S*tau.a*sigma*(gamma+tau.b*lambda)))
  part.2 <- (- mu_i - mu_v - beta*S*lambda*tau.b + beta*S*sigma*(gamma - tau.a*mu_i*gamma + tau.a*lambda + tau.b*lambda) )
  part.3 <- (-4*(1 + beta*S*sigma*tau.a*(gamma + tau.b*lambda))*(-beta*S*lambda + mu_i*mu_v - beta*S*sigma*(mu_i*gamma - lambda)))
  part.4 <- (mu_i + mu_v + beta*S*tau.b*lambda - beta*S*sigma*(gamma - tau.a*mu_i*gamma + tau.a*lambda + tau.b*lambda))^2
  
  fit1 <- (part.2 + sqrt(part.3 + part.4))/part.1
  fit2 <- (part.2 - sqrt(part.3 + part.4))/part.1
  return(fit1)
}
#*****allow to take multiple values for a parameter******
fit.delay.v <- Vectorize(fit.delay)

#*********************************************************************************************************
