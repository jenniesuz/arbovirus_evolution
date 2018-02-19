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
  
  fit1 <- -(((mu_i+mu_v) + sqrt( (mu_i+mu_v)^2 - 4*(mu_v*mu_i - beta*S*lambda) ) )  / 2)
  fit2 <- -(((mu_i+mu_v) - sqrt( (mu_i+mu_v)^2 - 4*(mu_v*mu_i - beta*S*lambda) ) ) / 2)
  return(fit2)
}
#*****allow to take multiple values for a parameter******
persist.fit.v <- Vectorize(persist.fit)
#*********************************************************************************

#*************PERSISTENT INFECTION WITH LATENT PERIOD******************************
#***************************fitness function*******************
p.latent.fit <- function(mu_v=0.1
                         ,mu_c = 1/240 #1/120   
                         ,beta = 10^-6       # probability of infecting susceptible cells
                         ,lambda=10^2        # release rate of virus
                         ,tau=3
                         ,S = 10^6){
  
  fit1 <- -(((mu_c+mu_v - beta*S*lambda*tau) + sqrt( (mu_c+mu_v - beta*S*lambda*tau)^2 - 4*(mu_v*mu_c - beta*S*lambda) ) )  / 2)
  fit2 <- -(((mu_c+mu_v - beta*S*lambda*tau) - sqrt( (mu_c+mu_v - beta*S*lambda*tau)^2 - 4*(mu_v*mu_c - beta*S*lambda) ) ) / 2)
  return(fit2)
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
vals.bet <- persist.fit.v(beta=bet)

pdf(file="fig1.pdf",width=6,height=4)
par(mfcol=c(1,2),mar=c(4,4,1,1),cex=0.7)
plot.func(var=lam
          ,vals=vals.lam
          ,name=expression(paste("Virus budding rate (",lambda,")"))
          ,plotno="a"
)
plot.func(var=bet*S
          ,vals=vals.bet
          ,name=expression(paste(beta,italic("S")))
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
                      ,gamma = 10^3       # virus yeild at apoptosis
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
                      ,gamma=10^2         
                      ,alpha=1/24
                      ,tau=3
                      ,S = 10^6){
  
  fit1 <- -(((mu_c + alpha + mu_v - tau*alpha*mu_v - beta*S*gamma*alpha) + sqrt( (mu_c + alpha + mu_v - tau*alpha*mu_v - beta*S*gamma*alpha)^2 - 4*( (1 - tau*alpha) * (mu_v*mu_c + mu_v*alpha - beta*S*gamma*alpha)  ) ) )  / (2 - 2*tau*alpha) )
  fit2 <- -(((mu_c + alpha + mu_v - tau*alpha*mu_v - beta*S*gamma*alpha) - sqrt( (mu_c + alpha + mu_v - tau*alpha*mu_v - beta*S*gamma*alpha)^2 - 4*((1 - tau*alpha) * (mu_v*mu_c + mu_v*alpha - beta*S*gamma*alpha)) ) ) / (2 - 2*tau*alpha) )
  return(fit2)
}
#********************************************
acute.fit.latent.v <- Vectorize(acute.fit.latent)
#*********************************************************************************

#*********************RANGES FOR PARAMETER VALUES***************************
gam <- seq(0,10^4,100) #seq(0,10^2.5,2)
bet <- seq(10^-10,10^-5,10^-10) #seq(10^-10,10^-5,10^-8)
alph <- seq(1/72,1/2,1/200)
t <- seq(0.5,10,0.1)
#**************************************************************************

#********PLOT NO LATENT for range of values***********
vals.gam <- acute.fit.v(gamma=gam)        # alpha really influences initial increase for given gamma
vals.bet <- acute.fit.v(beta=bet)
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
          ,vals=vals.alp
          ,name=expression(paste("Inverse of virus apoptosis rate (1/",alpha,")"))
          ,plotno="c")

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

