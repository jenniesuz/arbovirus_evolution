

#***********************************Force of selection*********************************
library(Deriv)

#*******************************Without delays******************************************
# derivative w.r.t gamma
Deriv(  expression( (-(mu_i + mu_v + alpha)  + sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda)) )  )/2) ,"gamma")

# derivative w.r.t beta
Deriv(  expression(  (- (mu_i + mu_v + alpha)  + sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda)))   )/2) ,"beta")

# derivative w.r.t alpha
Deriv(  expression(  (-(mu_i + mu_v + alpha)  + sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda)))   )/2) ,"alpha")

# derivative w.r.t. lambda
Deriv(  expression(  (-(mu_i + mu_v + alpha)  + sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda)))   )/2) ,"lambda")

#*******Force of selection with respect to gamma********************
do.dg <- function(gamma= 50
                  ,lambda=0
                  ,alpha = 1/24
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1){
  return(  ( alpha*beta*S / sqrt( (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda ))  ) )
  )
}
do.dg.v <- Vectorize(do.dg)
#*******Force of selection with respect to beta********************
do.db <- function(gamma= 50
                  ,lambda=0
                  ,alpha = 1/24
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1 ){ 
  return( ( (gamma*alpha+lambda)*S / sqrt( (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda ))  ) )
  )
}
do.db.v <- Vectorize(do.db)
#********Force of selection with respect to alpha*************
do.da <- function(gamma= 50
                  ,lambda=0
                  ,alpha = 1/24
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1 ){ 
  return( 
    0.5*(0.5*((2*(alpha+mu_i+mu_v) - 4*(mu_v - beta*gamma*S) )/sqrt( (alpha+mu_i+mu_v)^2 - 4*(mu_v*(alpha + mu_i) - beta*S*(gamma*alpha+lambda)) )  )-1)
  )
}

do.da.v <- Vectorize(do.da)
#*******Force of selection with respect to lambda********************
do.dl <- function(gamma= 0
                  ,lambda=100
                  ,alpha = 0
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1){
  return(  ( beta*S / sqrt( (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda ))  ) )
  )
}
do.dl.v <- Vectorize(do.dg)
#*************************************************************************************************


#*******************************With delays******************************************
# derivative w.r.t gamma
Deriv(expression(  (- (mu_c + mu_v - gamma*beta*S*exp(-mu_c*tau.a) + gamma*beta*S*exp(-mu_c*tau.a)*tau.a*mu_c  + lambda*exp(-mu_c*tau.b)*tau.b*beta*S - lambda*exp(-mu_c*tau.b)*beta*S*exp(-mu_c*tau.a)*tau.a  - lambda*exp(-mu_c*tau.b)*tau.b*beta*S*exp(-mu_c*tau.a)) + sqrt((mu_c + mu_v - gamma*beta*S*exp(-mu_c*tau.a) + gamma*beta*S*exp(-mu_c*tau.a)*tau.a*mu_c  + lambda*exp(-mu_c*tau.b)*tau.b*beta*S - lambda*exp(-mu_c*tau.b)*beta*S*exp(-mu_c*tau.a)*tau.a  - lambda*exp(-mu_c*tau.b)*tau.b*beta*S*exp(-mu_c*tau.a))^2 - 4*(1  + gamma*beta*S*exp(-mu_c*tau.a)*tau.a + beta*S*exp(-mu_c*tau.a)*tau.a*lambda*exp(-mu_c*tau.b)*tau.b)*(mu_c*mu_v - gamma*beta*S*exp(-mu_c*tau.a)*mu_c - lambda*exp(-mu_c*tau.b)*beta*S + lambda*exp(-mu_c*tau.b)*beta*S*exp(-mu_c*tau.a))))/ (2*(1  + gamma*beta*S*exp(-mu_c*tau.a)*tau.a + beta*S*exp(-mu_c*tau.a)*tau.a*lambda*exp(-mu_c*tau.b)*tau.b))
),"gamma")

# derivative w.r.t lambda
Deriv(expression( ),"lambda")

# derivative w.r.t beta
Deriv(expression( ),"beta")

# derivative w.r.t tau.a
Deriv(expression( ),"tau.a")

# derivative w.r.t tau.b
Deriv(expression( ),"tau.b")
# 
# 
# #*******Force of selection with respect to gamma********************
# do.dg <- function(gamma= 50
#                   ,lambda=0
#                   ,beta = 10^-6
#                   ,tau.a=24
#                   ,tau.b=24
#                   ,S = 10^6
#                   ,mu_i = 1/120
#                   ,mu_v = 0.1){
#   return(  ( alpha*beta*S / sqrt( (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda ))  ) )
#   )
# }
# do.dg.v <- Vectorize(do.dg)
# #***************************************************************
