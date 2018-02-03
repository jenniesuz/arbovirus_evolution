

#***********FITNESS FUNCTIONS**********

#****************PERSISTENT INFECTION WITHOUT DELAY/ LATENT PERIOD*******************
# Persistent infection, no delay in infected cells producing virus
#*******************fitness function***************
persist.fit <- function(mu_v=0.1            # virus clearance rate
                        ,mu_i = 1/120       # infected cell death rate
                        ,beta = 10^-6       # probability of single virion infecting susceptile cell
                        ,lambda=10^2        # release rate of virus from infected cells
                        ,S = 10^6){         # steady-state level of susceptible cells
  
  fit1 <- -(((mu_i+mu_v) + sqrt( (mu_i+mu_v)^2 - 4*(mu_v*mu_i - beta*S*lambda) ) )  / 2)
  fit2 <- -(((mu_i+mu_v) - sqrt( (mu_i+mu_v)^2 - 4*(mu_v*mu_i - beta*S*lambda) ) ) / 2)
  return(fit2)
}
#*****allow to take multiple values for a parameter******
persist.fit.v <- Vectorize(persist.fit)
#*********************************************************************************



#*************PERSISTENT INFECTION WITH LATENT PERIOD******************************
#*********fitness function*******************
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

#************************

#*********************RANGES FOR PARAMETER VALUES***************************
lam <- seq(1,1000,0.1)#seq(0,10^2.5,2)
bet <- seq(10^-10,10^-5,10^-10)#seq(10^-10,10^-5,10^-8)
t <- seq(0.5,10,0.1)
#**************************************************************************

par(mfcol=c(1,2))
#****************RUN FUNCS AND PLOT FITNESS AS A FUNCTION OF PARAMETERS*************
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

#***********************no latent period***************************
vals.lam <- persist.fit.v(lambda=lam)
vals.bet <- persist.fit.v(beta=bet)

pdf(file="fig1.pdf",width=6,height=4)
par(mfcol=c(1,2),mar=c(4,4,1,1),cex=0.6)
plot.func(var=lam
          ,vals=vals.lam
          ,name=expression(paste("Virus budding rate (",lambda,")"))
          ,plotno="a"
          )
plot.func(var=bet
          ,vals=vals.bet
          ,name=expression(paste("Virus infection rate (",beta,")"))
          ,plotno="b"
          ,ylabel="")

dev.off()
#****************************************************************


#**********************latent period*************************************
vals.lam <- p.latent.fit.v(lambda=lam,tau=3)
vals.bet <- p.latent.fit.v(beta=bet,tau=2)
vals.tau <- p.latent.fit.v(tau=t,mu_c=1/500)

par(mfcol=c(2,2))

plot.func(var=lam,vals=vals.lam,name="Virus budding rate")
plot.func(var=bet,vals=vals.bet,name="Probability of infection")
plot.func(var=t,vals=vals.tau,name="Delay")
