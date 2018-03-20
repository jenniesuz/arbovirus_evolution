

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
bet <- seq(10^-10,10^-5,10^-8)#seq(10^-10,10^-5,10^-10)  #
t <- seq(0.5,10,0.1)
S <- 10^6
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

plot(bet*S,foi.bet
     ,bty="n"
     ,type="l"
     ,xlab=expression(paste(beta,"S"))
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
do.dg <- function(gamma= 2400
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
do.db <- function(gamma= 2400
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
do.da <- function(gamma= 2400
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
gam <- lam*24 #seq(0,10^2.5,2)
bet <- seq(10^-10,10^-5,10^-10) #seq(10^-10,10^-5,10^-8)
alph <- seq(1/72,1/3,1/200)

gamforalph <-  1/alph*100  

t <- seq(0.5,10,0.1)
#***************************************************************************


#******************Plot****************************************
foi.alph <- do.da.v(alpha=alph)
foi.alphwithgam <- do.da.v(alpha=alph,gamma=gamforalph)

foi.gam <- do.dg.v(gamma=gam)
foi.bet <- do.db.v(beta=bet)

pdf(file="fig5.pdf",width=4,height=4)
par(mfcol=c(2,2),mar=c(4,4,1,1),cex=0.6)
plot(gam,foi.gam
     ,bty="n"
     ,type="l"
     ,xlab=expression(paste("Virus yeild (",gamma,")"))
     ,ylab="Force of selection"
     ,main="a"
)

plot(bet*S,foi.bet
     ,bty="n"
     ,type="l"
     ,xlab=expression(paste(beta,italic("S")))
     ,ylab="Force of selection"
     #,log="y"
     ,main="b")

plot(1/alph
     ,foi.alphwithgam
     ,bty="n"
     ,type="l"
     ,xlab=expression(paste("Inverse of the apoptosis rate (1/",alpha,")"))
     ,ylab="Force of selection"
     ,main="c")

plot(1/alph
     ,foi.alph
     ,bty="n"
     ,type="l"
     ,xlab=expression(paste("Inverse of the apoptosis rate (1/",alpha,")"))
     ,ylab="Force of selection"
     ,main="d")




dev.off()
#***************************************************************************






#*******************************FULL MODEL******************************************

#**********************************
# derivative w.r.t gamma
Deriv(  expression(  - (( (mu_i + mu_v + alpha)  - sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda)) )  )/2)) ,"gamma")

# derivative w.r.t beta
Deriv(  expression(  - (( (mu_i + mu_v + alpha)  - sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda)))   )/2)) ,"beta")

# derivative w.r.t alpha
Deriv(  expression(  - (( (mu_i + mu_v + alpha)  - sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda)))   )/2)) ,"alpha")

# derivative w.r.t. lambda
Deriv(  expression(  - (( (mu_i + mu_v + alpha)  - sqrt(  (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda)))   )/2)) ,"lambda")

#***********************************


#*******Force of selection with respect to gamma********************
do.dg <- function(gamma= 2400
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
do.db <- function(gamma= 2400
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
do.da <- function(gamma= 2400
                  ,lambda=0
                  ,alpha = 1/24
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1 ){ 
  return( 
    -0.5*(1-0.5*(2*(alpha+mu_i+mu_v) - 4*(mu_v - beta*gamma*S) )/sqrt( (alpha+mu_i+mu_v)^2 - 4*(mu_v*(alpha + mu_i) - beta*S*(gamma*alpha+lambda)) )  )
  )
}

do.da.v <- Vectorize(do.da)
#****************************************************
#*******Force of selection with respect to lambda********************
do.dl <- function(gamma= 0
                  ,lambda=10^2
                  ,alpha = 0
                  ,beta = 10^-6
                  ,S = 10^6
                  ,mu_i = 1/120
                  ,mu_v = 0.1){
  return(  ( beta*S / sqrt( (mu_i + mu_v + alpha)^2 - 4*(mu_i*mu_v + mu_v*alpha - beta*S*(gamma*alpha+lambda ))  ) )
  )
}
do.dl.v <- Vectorize(do.dg)


#*********************RANGES FOR PARAMETER VALUES***************************
lam <- seq(1,1000,0.1)

gam <- lam*24  #seq(0,10^2.5,2)
bet <- seq(10^-10,10^-5,10^-10) #seq(10^-10,10^-5,10^-8)
alph <- seq(1/72,1/3,1/200)
S<-10^6

gamforalph <-  1/alph*100  

t <- seq(0.5,10,0.1)
#***************************************************************************


#******************Plot****************************************
foi.alph <- do.da.v(alpha=alph)
foi.alphwithgam <- do.da.v(alpha=alph,gamma=gamforalph)

foi.gam <- do.dg.v(gamma=gam)
foi.betnolam <- do.db.v(beta=bet)
foi.betnogam <- do.db.v(beta=bet,lambda=10^2,gamma=0,alpha=0)

foi.lam <- do.dl.v(lambda=lam)


pdf(file="fig8.pdf",width=5,height=3)
par(mfcol=c(1,3),mar=c(4,4,1,1),cex=0.5)
plot(gam
     ,foi.gam
     ,bty="n"
     ,type="l"
     ,col="red"
     #,xlab=expression(paste("Virus yeild (",gamma,")"))
     ,xlab="Variable of interest"
     ,ylab="Force of selection"
     ,main="a"
)
legend("topleft",legend=c("No apoptosis","No budding"),col=c("red","blue"),bty="n",lty=c(1,2))
par(new=T)
plot(lam
     ,foi.lam
     ,bty="n"
     ,type="l"
     ,col="blue"
     ,lty=2
     ,yaxt="n"
     ,xaxt="n"
     ,xlab=" "
     #,xlab=expression(paste("Virus yeild (",lamb,")"))
     ,ylab=" "
     ,main="a"
)

plot(bet*S,foi.betnolam
     ,bty="n"
     ,type="l"
     ,col="red"
     ,xlab=expression(paste(beta,italic("S")))
     ,ylab=" "
     #,log="y"
     ,main="b")
par(new=T)
plot(bet*S,foi.betnogam
     ,bty="n"
     ,type="l"
     ,lty=2
     ,col="blue"
     ,xlab=expression(paste(beta,italic("S")))
     ,ylab=" "
     ,yaxt="n"
     ,xaxt="n"
     #,log="y"
     ,main="b")
plot(1/alph
     ,foi.alphwithgam
     ,bty="n"
     ,col="blue"
     ,lty=2
     ,type="l"
     ,xlab=expression(paste("Inverse of the apoptosis rate (1/",alpha,")"))
     ,ylab=" "
     ,main="c")
par(new=T)
plot(1/alph
     ,foi.alph
     ,xaxt="n"
     ,yaxt="n"
     ,col="blue"
     ,lty=1
     ,bty="n"
     ,type="l"
     ,xlab=" "
     ,ylab=" "
     #,xlab=expression(paste("Inverse of the apoptosis rate (1/",alpha,")"))
    # ,ylab="Force of selection"
     ,main=" ")
legend("topleft",legend=c("Constant virus yield","Variable virus yield"),col=c("blue","blue"),bty="n",lty=c(1,2))


dev.off()
#***************************************************************************





























