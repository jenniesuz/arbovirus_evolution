

#***********************************force of selection*********************************
library(Deriv)

#***********************PERSISTENT INFECTION WITHOUT LATENT PERIOD************************************
#**********************************
# derivative w.r.t lambda
Deriv(  expression(  - (( (mu_i + mu_v)  - sqrt(  (mu_i + mu_v)^2 - 4*(mu_i*mu_v - lambda*beta*S) )  )/2)) ,"lambda")

# derivative w.r.t beta
Deriv(  expression(  - (( (mu_i + mu_v)  - sqrt(  (mu_i + mu_v)^2 - 4*(mu_i*mu_v - lambda*beta*S) )  )/2)) ,"beta")
#**************************************************************

#*******force of selection with respect to lambda********************
do.dl <- function(lambda = 10^2
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1){
  return(  beta*S/ sqrt( (mu_i + mu_v)^2 - 4*(mu_i*mu_v - beta*S*lambda)  ) 
    )
}
do.dl.v <- Vectorize(do.dl)
#*********force of selection with respect to beta*****************
do.db <- function(lambda = 10^2
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1 ){ 
  return( lambda*S / sqrt( (mu_i + mu_v)^2 - 4*(mu_i*mu_v - beta*S*lambda)  ) 
  )
}
do.db.v <- Vectorize(do.db)
#****************************************************
lam <- seq(1,1000,0.1) #seq(0,10^2.5,2)
bet <- seq(10^-10,10^-8,10^-10) #seq(10^-10,10^-5,10^-8)
t <- seq(0.5,10,0.1)

#********Plots************
pdf(file="fig2.pdf",width=6,height=4)
par(mfcol=c(1,2),mar=c(4,4,1,1),cex=0.6)
foi.lam <- do.dl.v(lambda=lam)

plot(lam,foi.lam
     ,bty="n"
     ,type="l"
     ,xlab=expression(paste("Virus budding rate (",lambda,")"))
     ,ylab=expression(paste("Force of selection"))
     ,main="a"
)

foi.bet <- do.db.v(beta=bet)

plot(bet,foi.bet
     ,bty="n"
     ,type="l"
     ,xlab=expression(paste("Virus infection rate (",beta,")"))
     ,ylab=" "
     ,main="b"
)
dev.off()
#**************************************************



#***********************PERSISTENT INFECTION WITH LATENT PERIOD************************************

#**********************************
# derivative w.r.t lambda
Deriv(  expression(  - (( (mu_i + mu_v - beta*S*lambda*tau)  - sqrt(  (mu_i + mu_v - beta*S*lambda*tau)^2 - 4*(mu_i*mu_v - beta*S*lambda) )  )/2)) ,"lambda")

# derivative w.r.t beta
Deriv(  expression(  - (( (mu_i + mu_v)  - sqrt(  (mu_i + mu_v - beta*S*lambda*tau)^2 - 4*(mu_i*mu_v - beta*S*lambda*tau - lambda*beta*S) )  )/2)) ,"beta")

# derivative w.r.t tau
Deriv(  expression(  - (( (mu_i + mu_v - beta*S*lambda*tau)  - sqrt(  (mu_i + mu_v - beta*S*lambda*tau)^2 - 4*(mu_i*mu_v - lambda*beta*S) )  )/2)) ,"tau")


#***********************************

#*******force of selection functions********************
do.dl <- function(lambda = 100
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1){
  return(  ( beta*S / sqrt( (mu_i + mu_v)^2 - 4*(mu_i*mu_v - beta*S*lambda)  ) )
  )
}

do.dl.v <- Vectorize(do.dl)

do.db <- function(lambda = 100
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1 ){ 
  return( ( lambda*S / sqrt( (mu_i + mu_v)^2 - 4*(mu_i*mu_v - beta*S*lambda)  ) )
  )
}

do.db.v <- Vectorize(do.db)


#****************************************************







#***********************ACUTE INFECTION************************************

#**********************************
# derivative w.r.t gamma
Deriv(  expression(  - (( (mu_i + mu_v + alpha)  - sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*gamma*alpha) )  )/2)) ,"gamma")

# derivative w.r.t beta
Deriv(  expression(  - (( (mu_i + mu_v + alpha)  - sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*gamma*alpha) )  )/2)) ,"beta")

# derivative w.r.t alpha
Deriv(  expression(  - (( (mu_i + mu_v + alpha)  - sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*gamma*alpha) )  )/2)) ,"alpha")

#***********************************


#*******Force of selection with respect to gamma********************
do.dg <- function(gamma= 100
                  ,alpha = 1/24
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1){
  return(  ( alpha*beta*S / sqrt( (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*gamma*alpha)  ) )
  )
}
do.dg.v <- Vectorize(do.dg)
#*******Force of selection with respect to beta********************
do.db <- function(gamma= 100
                  ,alpha = 1/24
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1 ){ 
  return( ( gamma*alpha*S / sqrt( (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*gamma*alpha)  ) )
  )
}
do.db.v <- Vectorize(do.db)
#********Force of selection with respect to alpha*************
do.da <- function(gamma= 100
                  ,alpha = 1/24
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1 ){ 
  return( 
    -0.5*(1-0.5*(2*(alpha+mu_i+mu_v) - 4*(mu_v - beta*gamma*S) )/sqrt( (alpha+mu_i+mu_v)^2 - 4*(mu_v*(alpha + mu_i) - beta*gamma*alpha*S) )  )
  )
}

do.da.v <- Vectorize(do.da)
#****************************************************

#*********************RANGES FOR PARAMETER VALUES***************************
gam <- seq(10,10^3,10) #seq(0,10^2.5,2)
bet <- seq(10^-10,10^-8,10^-10) #seq(10^-10,10^-5,10^-8)
alph <- seq(1/120,1/3,1/200)
t <- seq(0.5,10,0.1)
#***************************************************************************


#******************Plot****************************************
foi.alph <- do.da.v(alpha=alph)
foi.gam <- do.dg.v(gamma=gam)
foi.bet <- do.db.v(beta=bet)

pdf(file="fig5.pdf",width=4,height=4)
par(mfcol=c(2,2),mar=c(4,4,1,1),cex=0.6)

plot(alph,foi.alph
     ,bty="n"
     ,type="l"
     ,xlab=expression(paste("Apoptosis rate (",alpha,")"))
     ,ylab="Force of selection"
     ,main="a")

plot(gam,foi.gam
     ,bty="n"
     ,type="l"
     ,xlab=expression(paste("Virus yeild (",gamma,")"))
     ,ylab="Force of selection"
     ,main="b"
     )

plot(bet,foi.bet
     ,bty="n"
     ,type="l"
     ,xlab=expression(paste("Virus infection rate (",beta,")"))
     ,ylab="Force of selection"
     ,main="c")

dev.off()
#***************************************************************************