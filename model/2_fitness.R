library(ggplot2)

# This script contains the virus fitness functions for persistent and acute infections #


#********************************PLOT FUNCTION**************************************
plot.func <- function(var,vals,name,plotno,ylabel=expression(paste("Virus fitness (",omega,")"))){
  plot(var
       ,vals
       ,bty="n"
       ,type="l"
       ,xlab=name
       ,ylab=ylabel
       ,main=plotno
       ,cex.main=1.2
  )
}
#**************************************************************************************



#****************PERSISTENT INFECTION WITH NO DELAY/ LATENT PERIOD*******************
#******************************fitness function*************************************
persist.fit <- function(mu_v=0.1            # virus clearance rate
                        ,mu_i = 1/120       # infected cell death rate
                        ,beta = 10^-6       # probability of single virion infecting susceptible cell
                        ,lambda=10^2        # release rate of virus from infected cells
                        ,S = 10^6){         # susceptible cells
  
  fit1 <- ((-(mu_i+mu_v) + sqrt( (mu_i+mu_v)^2 - 4*(mu_v*mu_i - beta*S*lambda) ) )  / 2)
  fit2 <- ((-(mu_i+mu_v) - sqrt( (mu_i+mu_v)^2 - 4*(mu_v*mu_i - beta*S*lambda) ) ) / 2)
  return(fit1)
}
#*****allow to take multiple values for a parameter******
persist.fit.v <- Vectorize(persist.fit)
#*********************************************************************************

#*************PERSISTENT INFECTION WITH LATENT PERIOD******************************
#***************************fitness function*******************
p.latent.fit <- function(mu_v=0.1
                         ,mu_c = 1/120 #1/120   
                         ,beta = 10^-6       # probability of infecting susceptible cells
                         ,lambda=10^2        # release rate of virus
                         ,tau=3
                         ,S = 10^6){
  
  fit1 <- ((-(mu_c+mu_v - beta*S*lambda*tau) + sqrt( (mu_c+mu_v - beta*S*lambda*tau)^2 - 4*(mu_v*mu_c - beta*S*lambda) ) )  / 2)
  fit2 <- ((-(mu_c+mu_v - beta*S*lambda*tau) - sqrt( (mu_c+mu_v - beta*S*lambda*tau)^2 - 4*(mu_v*mu_c - beta*S*lambda) ) ) / 2)
  return(fit1)
}
#********************************************
p.latent.fit.v <- Vectorize(p.latent.fit)
#*********************************************************************************

#*********************RANGES FOR PARAMETER VALUES***************************
lam <- seq(1,1000,0.1) #seq(0,10^2.5,2)
bet <- seq(10^-10,10^-5,10^-10) #seq(10^-10,10^-5,10^-8)
t <- seq(0.5,10,0.1)
S <- 10^6
#**************************************************************************

#********PLOT NO LATENT for range of beta and lambda values***********
vals.lam <- persist.fit.v(lambda=lam)
vals.bet.p <- persist.fit.v(beta=bet)

pdf(file="fig1.pdf",width=6,height=4)
par(mfcol=c(1,2),mar=c(4,4,1,1),cex=0.5)
plot.func(var=lam
          ,vals=vals.lam
          ,name=expression(paste("Virus budding rate (",lambda,")"))
          ,plotno="a"
)
plot.func(var=bet*S
          ,vals=vals.bet.p
          ,name=expression(paste("Product of virus infection rate and number of susceptible cells (", beta,italic("S"),")"))
          ,plotno="b"
          ,ylabel="")

dev.off()
#****************************************************************


#**********************PLOT latent period*************************************
vals.lam <- p.latent.fit.v(lambda=lam,tau=3)
vals.bet <- p.latent.fit.v(beta=bet,tau=2)
vals.tau <- p.latent.fit.v(tau=t,mu_c=1/500)

pdf(file="figx.pdf",width=6,height=3)
par(mfcol=c(2,2),mar=c(4,4,1,1),cex=0.6)
plot.func(var=lam
          ,vals=vals.lam
          ,name=expression(paste("Virus budding rate (",lambda,")"))
          ,plotno="a")
plot.func(var=bet
          ,vals=vals.bet
          ,name=expression(paste("Virus infection rate (",beta,")"))
          ,plotno="b")
plot.func(var=t
          ,vals=vals.tau
          ,name=expression(paste("Delay (",tau,")"))
          ,plotno="c")
dev.off()
#**************************************************************************





#********************ACUTE INFECTION NO LATENT PERIOD****************
#**********************fitness function**************************
acute.fit <- function(mu_v =0.1
                      ,mu_c = 1/120   
                      ,beta = 10^-6       # probability of infecting susceptible cells
                      ,gamma = 2400       # virus yeild at apoptosis - set so equivalent to lambda
                      ,alpha = 1/24       # apoptosis rate
                      ,K = 10^6){
  
  fit1 <- (-(alpha+mu_c+mu_v) + sqrt( (alpha+mu_c+mu_v)^2 - 4*(mu_v*alpha + mu_v*mu_c - beta*K*gamma*alpha) ) ) / 2
  fit2 <- (-(alpha+mu_c+mu_v) - sqrt( (alpha+mu_c+mu_v)^2 - 4*(mu_v*alpha + mu_v*mu_c - beta*K*gamma*alpha) ) ) / 2
  return(fit1)
}
#*****************************************************************
acute.fit.v <- Vectorize(acute.fit)
#**************************************************************

#*******************ACUTE INFECTION WITH LATENT PERIOD******************
#***************************fitness function*******************
acute.fit.latent <- function(mu_v=0.1
                      ,mu_c = 1/120   
                      ,beta = 10^-6       
                      ,gamma= 2400        
                      ,alpha=1/24
                      ,tau=3
                      ,S = 10^6){
  
  fit1 <- ((-(mu_c + alpha + mu_v - tau*alpha*mu_v - beta*S*gamma*alpha) + sqrt( (mu_c + alpha + mu_v - tau*alpha*mu_v - beta*S*gamma*alpha)^2 - 4*( (1 - tau*alpha) * (mu_v*mu_c + mu_v*alpha - beta*S*gamma*alpha)  ) ) )  / (2 - 2*tau*alpha) )
  fit2 <- ((-(mu_c + alpha + mu_v - tau*alpha*mu_v - beta*S*gamma*alpha) - sqrt( (mu_c + alpha + mu_v - tau*alpha*mu_v - beta*S*gamma*alpha)^2 - 4*((1 - tau*alpha) * (mu_v*mu_c + mu_v*alpha - beta*S*gamma*alpha)) ) ) / (2 - 2*tau*alpha) )
  return(fit1)
}
#********************************************
acute.fit.latent.v <- Vectorize(acute.fit.latent)
#*********************************************************************************

#*********************RANGES FOR PARAMETER VALUES***************************
gam <- lam*24 #seq(0,10^2.5,2)

bet <- seq(10^-10,10^-5,10^-10) #seq(10^-10,10^-5,10^-8)

alph <- seq(1/72,1/2,1/200)
gamforalph <- 1/alph*100        # keep as equivalent virus production per unit time as persistent infection

t <- seq(0.5,10,0.1)
#**************************************************************************

#********PLOT NO LATENT for range of values***********
vals.gam <- acute.fit.v(gamma=gam)        # alpha really influences initial increase for given gamma
vals.bet <- acute.fit.v(beta=bet)
vals.alp.gam <- acute.fit.v(alpha=alph,gamma=gamforalph)
vals.alp <- acute.fit.v(alpha=alph)

pdf(file="fig4.pdf",width=4,height=4)
par(mfcol=c(2,2),mar=c(4,4,1,1),cex=0.5)
plot.func(var=gam
          ,vals=vals.gam
          ,name=expression(paste("Virus yeild at apoptosis (",gamma,")"))
          ,plotno="a")

plot.func(var=bet*S
          ,vals=vals.bet
          ,name=expression(paste(beta,italic("S")))
          ,plotno="b")

plot.func(var=1/alph
          ,vals=vals.alp.gam
          ,name=expression(paste("Inverse of virus apoptosis rate (1/",alpha,")"))
          ,plotno="c")

plot.func(var=1/alph
          ,vals=vals.alp
          ,name=expression(paste("Inverse of virus apoptosis rate (1/",alpha,")"))
          ,plotno="d")

dev.off()
#****************************************************************

#********PLOT LATENT for range of beta and lambda values***********
vals.gam <- acute.fit.latent.v(gamma=gam)        # alpha really influences initial increase for given gamma
vals.bet <- acute.fit.latent.v(beta=bet)
vals.alp <- acute.fit.latent.v(alpha=alph)
vals.tau <- acute.fit.latent.v(tau=t,mu_c=1/500)

pdf(file="figz.pdf",width=6,height=3)
par(mfcol=c(2,2),mar=c(4,4,1,1),cex=0.6)
plot.func(var=gam
          ,vals=vals.gam
          ,name=expression(paste("Virus yeild at apoptosis (",gamma,")"))
          ,plotno="a")

plot.func(var=bet
          ,vals=vals.bet
          ,name=expression(paste("Virus infection rate (",beta,")"))
          ,plotno="b")

plot.func(var=bet
          ,vals=vals.bet
          ,name=expression(paste("Virus apoptosis rate (",alpha,")"))
          ,plotno="c")

plot.func(var=t
          ,vals=vals.tau
          ,name=expression(paste("Virus apoptosis rate (",alpha,")"))
          ,plotno="d")


dev.off()

#***********************************************************************************************************************************











#****************MODEL PERMITTING BOTH TYPES OF INFECTION WITH NO DELAY/ LATENT PERIOD*******************
#******************************fitness function*************************************
mod.fit <- function(mu_v=0.1            # virus clearance rate
                        ,mu_i = 1/120       # infected cell death rate
                        ,beta = 10^-6       # probability of single virion infecting susceptible cell
                        ,lambda=10^2        # release rate of virus from infected cells
                        ,alpha=1/24
                        ,gamma=2400
                        ,S = 10^6){         # susceptible cells
  
  fit1 <- (-(mu_i+mu_v + alpha) + sqrt( (mu_i+mu_v+alpha)^2 - 4*(mu_v*(alpha+mu_i) - beta*S*(gamma*alpha+lambda) ) ))  / 2
  fit2 <- (-(mu_i+mu_v + alpha) - sqrt( (mu_i+mu_v+alpha)^2 - 4*(mu_v*(alpha+mu_i) - beta*S*(gamma*alpha+lambda) ) )) / 2
  return(fit1)
}
#*****allow to take multiple values for a parameter******
mod.fit.v <- Vectorize(mod.fit)
#*********************************************************************************
lam <- seq(1,1000,0.1) #seq(0,10^2.5,2)
bet <- seq(10^-10,10^-5,10^-10) #seq(10^-10,10^-5,10^-8)
gam <- lam*24  #seq(0,10^2.5,2)
alph <- seq(1/72,1/2,1/200)

gamforalph <- 1/alph*10^2


S = 10^6

#********PLOT NO LATENT for range of values***********

#*****first assume zero apoptosis rate and zero yeild
vals.bet.noalph <- mod.fit.v(beta=bet,alpha=0,gamma=0)      # infection rate
vals.lam.noalph <- mod.fit.v(lambda=lam,alpha=0,gamma=0)    # budding rate


#******no budding rate 
vals.bet.nolam <- mod.fit.v(beta=bet,lambda=0) 

vals.gam.nolam <- mod.fit.v(gamma=gam,lambda=0)     # yeild at apoptosis

vals.gam.nolam.alp <- mod.fit.v(gamma=gam,lambda=0)

vals.alp.nolam.gam <- mod.fit.v(alpha=alph,gamma=gamforalph,lambda=0)    # apoptosis rate
vals.alp.nolam <-mod.fit.v(alpha=alph,lambda=0)
#********************************************************************


pdf(file="fig7.pdf",width=5,height=3)

par(mfcol=c(1,3),mar=c(4,4,1,1),cex=0.5)
plot(bet*S
    ,vals.bet.noalph
    ,xlab=expression(paste(beta,italic("S")))
    ,ylab=expression(paste("Virus fitness (",omega,")"))
    ,main="a"
    ,type="l"
    ,col="red"
    ,bty="n"
    ,ylim=c(0,40))
legend("topleft",legend=c("No apoptosis","No budding"),col=c("red","blue"),bty="n",lty=c(1,2))
par(new=T)
plot(bet*S
          ,vals.bet.nolam
          ,col="blue"
          ,ylab=""
          ,xlab=""
          ,type="l"
          ,lty=2
          ,bty="n"
          ,yaxt="n"
          ,xaxt="n"
          ,ylim=c(0,40)
     )


plot(lam
     ,vals.lam.noalph
     ,xlab=expression(paste(lambda," or ",gamma, alpha))
     ,ylab=" "
     ,main="b"
     ,type="l"
     ,col="red"
     ,bty="n"
      ,ylim=c(0,40)
   )

par(new=T)
plot(gam
     ,vals.gam.nolam
     ,xlab=""
     ,ylab=" "
     ,type="l"
     ,col="blue"
     ,yaxt="n"
     ,xaxt="n"
     ,bty="n"
     ,lty=2
     ,ylim=c(0,40)
     )

plot(1/alph
     ,vals.alp.nolam
     ,xlab=expression(paste("1/",alpha))
     ,ylab=" "
     ,main="c"
     ,type="l"
     ,col="blue"
     ,bty="n"
     ,ylim=c(0,40)
     )
par(new=T)
plot(1/alph
     ,vals.alp.nolam.gam
     ,xlab=" "
     ,ylab=" "
     ,type="l"
     ,col="blue"
     ,yaxt="n"
     ,xaxt="n"
     ,bty="n"
     ,lty=2
     ,ylim=c(0,40)
)
legend("topleft",legend=c("Constant virus yield","Variable virus yield"),col=c("blue","blue"),bty="n",lty=c(1,2))

dev.off()
#***************************************************************


