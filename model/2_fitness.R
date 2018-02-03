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
#**************************************************************************

#********PLOT NO LATENT plot for range of beta and lambda values***********
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














#******************lambda dependent on mu_c**************************************
#*********fitness function*******************
persist.fit <- function(mu_v=0.1
                        ,mu_c = 1/120   
                        ,beta = 10^-6       # probability of infecting susceptible cells
                        #,lambda=10^2        # release rate of virus
                        ,K = 10^6){
  lambda <- (10*exp(30*mu_c))
  fit1 <- -(((mu_c+mu_v) + sqrt( (mu_c+mu_v)^2 - 4*(mu_v*mu_c - beta*K*lambda) ) )  / 2)
  fit2 <- -(((mu_c+mu_v) - sqrt( (mu_c+mu_v)^2 - 4*(mu_v*mu_c - beta*K*lambda) ) ) / 2)
  return(fit2)
}
#********************************************
persist.fit.v <- Vectorize(persist.fit)

#********plot for rangseq of beta and lambda values***********
m<-seq(0,1/5,1/200)
test<-persist.fit.v(mu_c=m)

plot((10*exp(30*m))
     ,test,type="l"
     ,bty="n"
     ,xlab=expression(paste(lambda," - virus budding rate"))
     ,ylab="Virus fitness"
     ,main="Persistent infection"
     ,cex.main=0.8
)

#***********************************************************************************








# Acute infection, no delay in virus production, density-dependent susceptible cell
# replication


acute.fit <- function(mu_v =0.1
                      ,mu_c = 1/120   
                      ,beta = 10^-6       # probability of infecting susceptible cells
                      ,gamma = 10^3       # virus yeild at apoptosis
                      ,alpha = 1/24        # apoptosis rate
                      ,K = 10^6){
  
  fit1 <- (-(alpha+mu_c+mu_v) + sqrt( (alpha+mu_c+mu_v)^2 - 4*(mu_v*alpha + mu_v*mu_c - beta*K*gamma*alpha) ) ) / 2
  fit2 <- (-(alpha+mu_c+mu_v) - sqrt( (alpha+mu_c+mu_v)^2 - 4*(mu_v*alpha + mu_v*mu_c - beta*K*gamma*alpha) ) ) / 2
  return(fit1)
}

acute.fit.v <- Vectorize(acute.fit)

test <- acute.fit.v(gamma=seq(0,10^4,100))        # alpha really influences initial increase for given gamma
plot(seq(0,10^4,100)
     ,test,type="l"
     ,bty="n"
     ,xlab=expression(paste(gamma," - virus yeild at apoptosis"))
     ,ylab="Virus fitness"
     ,main="Acute infection"
     ,cex.main=0.8
)


test <- acute.fit.v(beta=seq(10^-10,10^-5,10^-8))
plot(seq(10^-10,10^-5,10^-8)*10^6
     ,test,bty="n"
     ,type="l"
     ,xlab=expression(paste(beta,italic(S),", for constant ",italic(S)))
     ,ylab="Virus fitness"
     ,main="Acute infection"
     ,cex.main=0.8
)



test <- acute.fit.v(alpha=seq(1/120,1/3,1/200))        # alpha really influences initial increase for given gamma
plot(seq(1/120,1/3,1/200)
     ,test,bty="n",type="l"
     ,xlab=expression(alpha)
     ,ylab="Virus fitness"
     ,main="Acute infection"
     ,cex.main=0.8)



dev.off()

#**********************************************************************************************************************************

#*******************ACUTE INFECTION WITH DELAY/ LATENT PERIOD******************

#*********fitness function*******************
acute.fit <- function(mu_v=0.1
                        ,mu_c = 1/120   
                        ,beta = 10^-6       # probability of infecting susceptible cells
                        ,gamma=10^2        # release rate of virus
                        ,alpha=1/24
                        ,tau=3
                        ,S = 10^6){
  
  fit1 <- -(((mu_c + alpha + mu_v - tau*alpha*mu_v + beta*S*gamma*alpha) + sqrt( (mu_c + alpha + mu_v - tau*alpha*mu_v + beta*S*gamma*alpha)^2 - 4*( (1 - tau*alpha) * (mu_v*mu_c + mu_v*alpha - beta*S*gamma*alpha)  ) ) )  / (2 - 2*tau*alpha) )
  fit2 <- -(((mu_c + alpha + mu_v - tau*alpha*mu_v + beta*S*gamma*alpha) - sqrt( (mu_c + alpha + mu_v - tau*alpha*mu_v + beta*S*gamma*alpha)^2 - 4*((1 - tau*alpha) * (mu_v*mu_c + mu_v*alpha - beta*S*gamma*alpha)) ) ) / (2 - 2*tau*alpha) )
  return(fit2)
}
#********************************************
acute.fit.v <- Vectorize(acute.fit)
#*********************************************************************************
par(mfcol=c(2,2),mar=c(4,4,4,1))

gamma<-seq(0,10^4,100)
tau <- 3
tau <- seq(1,10,0.1)

#********plot for range of beta and lambda values***********
var <-seq(0,10^4,100)
test <- acute.fit.v(gamma=var)

var <- seq(0.5,10,0.5)
test <- acute.fit.v(tau=var)

var <- seq(1/120,1/3,1/200)
test <- acute.fit.v(alpha=var)

var <- seq(10^-10,10^-5,10^-8)
test <- acute.fit.v(beta=var)

#tiff("fitness_func_plots.tiff", height =4, width = 4, units = 'in', compression="lzw", res=400)

plot(var
     ,test,type="l"
     ,bty="n"
     #,xlab=expression(paste(lambda," - release rate of virus / hr"))
     ,ylab="Virus fitness"
     ,main="Acute infection"
     ,cex.main=0.8
)
#**************************


#***********************************************************************************************************************************

















#*************************************************
#******************gamma dependent on alpha**************************************
#*********fitness function*******************
acute.fit <- function(mu_v =0.1
                      ,mu_c = 1/120   
                      ,beta = 10^-6       # probability of infecting susceptible cells
                      #,gamma = 10^3       # virus yeild at apoptosis
                      ,alpha = 1/24        # apoptosis rate
                      ,K = 10^6){
  gamma <- (1000*exp(-20*alpha))
  fit1 <- (-(alpha+mu_c+mu_v) + sqrt( (alpha+mu_c+mu_v)^2 - 4*(mu_v*alpha + mu_v*mu_c - beta*K*gamma*alpha) ) ) / 2
  fit2 <- (-(alpha+mu_c+mu_v) - sqrt( (alpha+mu_c+mu_v)^2 - 4*(mu_v*alpha + mu_v*mu_c - beta*K*gamma*alpha) ) ) / 2
  return(fit1)
}

acute.fit.v <- Vectorize(acute.fit)


alpha <- seq(1/120,1/3,1/200)
gamma <- (1000*exp(-20*alpha))
plot(alpha,gamma)

#********plot for rangseq of beta and lambda values***********

test<-acute.fit.v(alpha=alpha)

plot((1000*exp(-20*alpha))
     ,test,type="l"
     ,bty="n"
     ,xlab=expression(paste(gamma," - virus yeild"))
     ,ylab="Virus fitness"
     ,main="Acute infection"
     ,cex.main=0.8
)
#*******************************************************************


a <- seq(1/120,1,1/25)
g <- seq(0,10^4,100)

acute.param.values <- expand.grid("mu_v" = 0.1
                            ,"mu_c" = 1/120
                            ,"beta" = 10^-10  
                            ,"gamma" = g
                            ,"alpha" = a
                            ,"S" = 10^6
                            
)

#**** acute infection****
test <- sapply(1:length(acute.param.values$mu_v),function(x){
  temp <- acute.param.values[x,]
  res <- acute.fit(temp[1],temp[2],temp[3],temp[4],temp[5],temp[6])
  return(as.vector(res))
})
test <-unlist(test)

acute.res <- cbind.data.frame(acute.param.values,res=test)
acute.res$alpha <- round(acute.res$alpha,2)
acute.res$gamma <- round(acute.res$gamma,2)


acute.res$quart <- cut(acute.res$alpha,quantile(acute.res$alpha))


#***********************************************************
par(mfcol=c(1,1))
ggplot(acute.res,aes(x=gamma,y=res,color=factor(quart))) +
geom_point(shape=16,size=3) +
theme(legend.position=c(0.8,0.85)
,legend.background=element_rect(colour="lightgrey"),
panel.background = element_rect(fill="lightgrey")) +
labs(x=expression(gamma),y="Fitness") 
scale_colour_manual(values=c("#fdcc8a", "#fc8d59", "#e34a33", "#b30000"),
name="Quartiles of alpha",
labels=c("0.01-0.09","0.09-0.17","0.17-0.25","0.25=0.33"))
#**********************************************************


persist.param.values <- expand.grid("mu_v" = 0.1
                                    ,"mu_c" = 1/120
                                    ,"beta" = seq(10^-20,10^-2,10^-4)
                                    ,"lambda" = c(1,10,100,1000)
                                    ,"S" = 10^6
                                    
)
