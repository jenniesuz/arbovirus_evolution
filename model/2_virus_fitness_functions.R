

#**********************************Virus fitness functions***********************************

#****************Model with budding and apoptosis, no delays*******************
fit.nodelay <- function(mu_v=0.1                # virus clearance rate
                        ,mu_c = 1/120       # infected cell death rate
                        ,beta = 10^-6       # probability of single virion infecting susceptible cell
                        ,lambda=10          # release rate of virus from infected cells
                        ,alpha=1/24         # apoptosis rate
                        ,gamma=50           # virus yield at apoptosis
                        ,S = 10^6){         # susceptible cells
  fit1 <- (- (mu_c+mu_v+alpha) + sqrt( (mu_c+mu_v+alpha)^2 - 4*(mu_v*(alpha+mu_c) - beta*S*(gamma*alpha+lambda) ) ))  / 2
  fit2 <- (- (mu_c+mu_v+alpha) - sqrt( (mu_c+mu_v+alpha)^2 - 4*(mu_v*(alpha+mu_c) - beta*S*(gamma*alpha+lambda) ) )) / 2
  return(fit1)
}
#*****allow to take multiple values for a parameter******
fit.nodelay.v <- Vectorize(fit.nodelay)
#***********************************************************************************


#****************Model with budding and apoptosis, with delays**********************
fit.delay <- function(mu_v=0.1           # virus clearance rate
                        ,mu_c = 1/120       # infected cell death rate
                        ,beta = 10^-6       # probability of single virion infecting susceptible cell
                        ,lambda=10          # release rate of virus from infected cells
                        ,gamma=50           # virus yield at apoptosis
                        ,tau.a=24           # time to apoptosis  
                        ,tau.b=24           # time to start of budding
                        ,S = 10^6){         # susceptible cells
  
     sigma  <- exp(-mu_c*tau.a)
     sigmab <- exp(-mu_c*tau.b)

  a <- (1  + gamma*beta*S*sigma*tau.a + beta*S*sigma*tau.a*lambda*sigmab*tau.b)
   b <- (mu_c + mu_v - gamma*beta*S*sigma + gamma*beta*S*sigma*tau.a*mu_c  + lambda*sigmab*tau.b*beta*S - lambda*sigmab*beta*S*sigma*tau.a  - lambda*sigmab*tau.b*beta*S*sigma)
   c <- (mu_c*mu_v - gamma*beta*S*sigma*mu_c - lambda*sigmab*beta*S + lambda*sigmab*beta*S*sigma)
  fit1 <- (- b + sqrt(b^2 - 4*a*c))/ (2*a)
  
  # without accounting for infected cell death for budding delay
 #   a <- (1 + beta*S*tau.a*sigma*gamma + beta*S*tau.a*sigma*tau.b*lambda)
#   b <- (mu_c + mu_v + beta*S*lambda*tau.b - beta*S*sigma*gamma + beta*S*sigma*gamma*tau.a*mu_c - beta*S*sigma*lambda*tau.a - beta*S*sigma*lambda*tau.b)
  #  c <- (mu_c*mu_v - beta*S*lambda - beta*S*sigma*gamma*mu_c +  beta*S*sigma*lambda)
  # fit1 <- (- b + sqrt(b^2 - 4*a*c))/ (2*a)

  return(fit1)
}
#*****allow to take multiple values for a parameter******
fit.delay.v <- Vectorize(fit.delay)

#*********************************************************************************************************


# 
# #****************Acute virus with delay**********************
# fit.a.delay <- function(mu_v=0.1           # virus clearance rate
#                         ,mu_c = 1/120       # infected cell death rate
#                         ,beta = 10^-6       # probability of single virion infecting susceptible cell
#                         ,gamma=50           # virus yield at apoptosis
#                         ,tau.a=24           # time to apoptosis  
#                         ,S = 10^6){         # susceptible cells
#   
#   sigma  <- exp(-mu_c*tau.a)
#   
#   a <- (1 + beta*S*tau.a*sigma*gamma)
#   b <- (mu_c + mu_v + gamma*beta*S*sigma*(tau.a*mu_c - 1))
#   c <- (mu_c*mu_v - beta*S*sigma*gamma*mu_c)
#   fit1 <- (- b + sqrt(b^2 - 4*a*c)) / (2*a)
#   
#   return(fit1)
# }
# #*****allow to take multiple values for a parameter******
# fit.a.delay.v <- Vectorize(fit.a.delay)
# 
# #*********************************************************************************************************
