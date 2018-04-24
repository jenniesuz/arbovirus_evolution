library(ggplot2)

# This script contains the virus fitness functions for persistent and acute infections #


#****************MODEL PERMITTING BOTH BUDDING AND APOPTOSIS WITH NO DELAY/ LATENT PERIOD*******************
#******************************fitness function*************************************
mod.fit <- function(mu_v=0.1            # virus clearance rate
                        ,mu_i = 1/120       # infected cell death rate
                        ,beta = 10^-6       # probability of single virion infecting susceptible cell
                        ,lambda=10        # release rate of virus from infected cells
                        ,alpha=1/24
                        ,production=50
                        ,S = 10^6){         # susceptible cells
  if(alpha>0){gamma <- 1/alpha*production - 1/alpha*lambda}else{gamma <- 0}
  fit1 <- (- (mu_i+mu_v + alpha) + sqrt( (mu_i+mu_v+alpha)^2 - 4*(mu_v*(alpha+mu_i) - beta*S*(gamma*alpha+lambda) ) ))  / 2
  fit2 <- (- (mu_i+mu_v + alpha) - sqrt( (mu_i+mu_v+alpha)^2 - 4*(mu_v*(alpha+mu_i) - beta*S*(gamma*alpha+lambda) ) )) / 2
  return(fit1)
}
#*****allow to take multiple values for a parameter******
mod.fit.v <- Vectorize(mod.fit)
#*********************************************************************************
lam <- seq(1,50,0.1) #seq(0,10^2.5,2)
bet <- seq(10^-10,10^-5,10^-10) #seq(10^-10,10^-5,10^-8)
alph <- seq(1/40,1/2,1/200)
S = 10^6

#********PLOT NO LATENT for range of values***********
vals.lam <- mod.fit.v(lambda=lam,alpha=0.01)
vals.alph <- mod.fit.v(alpha=alph)

results <- cbind.data.frame(alph,vals.alph)

ggplot(results, aes(x=1/alph,y=vals.alph)) +
  geom_line(aes(x=1/alph,y=vals.alph)) +
  labs( x=expression(paste(1/alpha)),y="Fitness") 

















pdf(file="fig7.pdf",width=4,height=4)

par(mfcol=c(2,2),mar=c(4,4,1,1),cex=0.5)
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
     ,xlab=expression(paste(lambda))
     ,ylab=" "
     ,main="b"
     ,type="l"
     ,col="red"
     ,bty="n"
      ,ylim=c(0,40)
   )


plot(gam
     ,vals.gam.nolam
     ,xlab=expression(paste(gamma))
     ,ylab=expression(paste("Virus fitness (",omega,")"))
     ,type="l"
     ,col="blue"
     ,main="c"
     ,bty="n"
     ,lty=2
     ,ylim=c(0,40)
     )

plot(1/alph
     ,vals.alp.nolam
     ,xlab=expression(paste("1/",alpha))
     ,ylab=" "
     ,main="d"
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


