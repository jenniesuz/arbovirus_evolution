
#******************Evolutionary invasion analysis**************************


#**************************PERSISTENT INFECTION***********************************
#*************persistent mutant virus fitness function************
persist.mutfit <- function(mu_mv=0.1            # mutant virus clearance rate
                           ,mu_mc = 1/120       # mutant virus cell death rate
                           ,beta_m = 10^-6      # mutant probability of infecting susceptible cells
                           ,lambda_m=10^2       # release rate of mutant virus
                           ,mu_v=0.1            # resident virus clearance rate  
                           ,mu_c = 1/120        # resident infected cell death rate 
                           ,beta = 10^-6        # resident probability of infecting susceptible cells
                           ,lambda = 10^2 ){    # release rate of resident virus
  
  fit1 <- - (((mu_mc+mu_mv) + sqrt( (mu_mc+mu_mv)^2 - 4*(mu_mv*mu_mc - (beta_m*lambda_m*mu_c*mu_v)/(lambda*beta)) ) ) / 2)
  fit2 <- - (((mu_mc+mu_mv) - sqrt( (mu_mc+mu_mv)^2 - 4*(mu_mv*mu_mc - (beta_m*lambda_m*mu_c*mu_v)/(lambda*beta) ) ) ) / 2)
  return(fit2)
}

persist.mutfit.v <- Vectorize(persist.mutfit)         # enable multiple values to be input for a parameter


#******************************************************************

#***************test function******************************
persist.mutfit(lambda_m=10^2)

test <- persist.mutfit.v(lambda_m=seq(0,10^3,10))

plot(seq(0,10^3,10)
     ,test,type="l"
     ,bty="n"
     ,xlab=expression(paste(lambda," - release rate of virus"))
     ,ylab="Virus fitness"
     ,main="Persistent infection"
     ,cex.main=0.8
)
abline(h=0)

test <- persist.mutfit.v(beta_m=seq(10^-10,10^-5,10^-8))
plot(seq(10^-10,10^-5,10^-8)*10^6
     ,test,bty="n"
     ,type="l"
     ,xlab=expression(paste(beta,italic(S),", for constant ",italic(S)))
     ,ylab="Virus fitness"
     ,main="Persistent infection"
     ,cex.main=0.8
)
abline(h=0)
#**********************************************************************************

#****************conditions favouring mutant virus evolution***********************************
library(rootSolve) 

# look at where mutant virus fitness function is zero across a range of parameter values relative to
# those of the resident virus

# re-write mutant virus fitness function for input to unitroot.all function
fun <- function(x
                ,mu_mv=0.1
                ,mu_mc = 1/120   
                ,beta_m = 10^-6      
                ,lambda_m = 10^2       
                ,mu_v=0.1
                ,mu_c = 1/120   
                ,beta = 10^-6      
                ,lambda = 10^2){
  if(lambda_m=="NA"){
  return(-(((mu_mc+mu_mv) - sqrt( (mu_mc+mu_mv)^2 - 4*(mu_mv*mu_mc - (beta_m*x*mu_c*mu_v)/(lambda*beta)) ) ) / 2))
  }else{
    if(beta_m=="NA"){
      return(-(((mu_mc+mu_mv) - sqrt( (mu_mc+mu_mv)^2 - 4*(mu_mv*mu_mc - (x*lambda_m*mu_c*mu_v)/(lambda*beta)) ) ) / 2))
   }
  }

}

#************fitness boundary for lambda and beta***************

beta_m <- seq(10^-10,10^-6,10^-10) # range of values for mutant virus probabillity of susceptible cell infection

fit.bound <- sapply(beta_m,function(y){
  uniroot.all(fun,c(0,10^3),beta_m=y,lambda_m="NA")
})

pdf(file="fig3.pdf",width=6,height=4)
par(mfcol=c(1,1),mar=c(4,4,1,1),cex=0.6)
plot(beta_m
     ,fit.bound
     ,type="l"
     ,ylab=expression(lambda)
     ,xlab=expression(beta)
     ,bty="n")
dev.off()
#*********************************************************


#************fitness boundary for lambda and mu_c***************

mu_mc <- seq(1/240,1/10,1/600)  # range of values for mutant virus probabillity of susceptible cell infection

test <- sapply(mu_mc,function(y){
  uniroot.all(fun,c(0,10^3),mu_mc=y,lambda_m="NA")
})

plot(mu_mc
     ,test
     ,type="l"
     ,ylab=expression(lambda)
     ,xlab=expression(paste(mu[I])))
#**************************************************************
#***********************************************************************************


#***********************ACUTE INFECTION*************************************
#*************acute mutant virus fitness function************
acute.mutfit <- function(mu_mv=0.1            # mutant virus clearance rate
                           ,mu_mc = 1/120       # mutant virus cell death rate
                           ,beta_m = 10^-6      # mutant probability of infecting susceptible cells
                           ,gamma_m = 2400      # mutant virus yeild at apoptosis
                           ,alpha_m = 1/24      # mutant virus apoptosis rate
                           ,mu_v=0.1            # resident virus clearance rate  
                           ,mu_c = 1/120        # resident infected cell death rate 
                           ,beta = 10^-6        # resident probability of infecting susceptible cells
                           ,alpha = 1/24        # resident virus apoptosis rate
                           ,gamma = 2400 ){     # resident virus yeild at apoptosis 
  
  fit1 <- - (((alpha_m+mu_mc+mu_mv) + sqrt( (alpha_m+mu_mc+mu_mv)^2 - 4*(alpha_m*mu_mv + mu_mv*mu_mc - (beta_m*gamma_m*alpha_m*mu_c*mu_v)/(gamma*alpha*beta)) ) ) / 2)
  fit2 <- - (((alpha_m+mu_mc+mu_mv) - sqrt( (alpha_m+mu_mc+mu_mv)^2 - 4*(alpha_m*mu_mv + mu_mv*mu_mc - (beta_m*gamma_m*alpha_m*mu_c*mu_v)/(gamma*alpha*beta) ) ) ) / 2)
  return(fit2)
}

acute.mutfit.v <- Vectorize(acute.mutfit)         # enable multiple values to be input for a parameter


# re-write mutant virus fitness function for input to unitroot.all function
fun <- function(x
                ,mu_mv=0.1            # mutant virus clearance rate
                ,mu_mc = 1/120       # mutant virus cell death rate
                ,beta_m = 10^-6      # mutant probability of infecting susceptible cells
                ,gamma_m = 2400      # mutant virus yeild at apoptosis
                ,alpha_m = 1/24      # mutant virus apoptosis rate
                ,mu_v=0.1            # resident virus clearance rate  
                ,mu_c = 1/120        # resident infected cell death rate 
                ,beta = 10^-6        # resident probability of infecting susceptible cells
                ,alpha = 1/24        # resident virus apoptosis rate
                ,gamma = 2400
                ){
  
  if(gamma_m=="NA"){
    return(-(((alpha_m+mu_mc+mu_mv) - sqrt( (alpha_m+mu_mc+mu_mv)^2 - 4*(alpha_m*mu_mv + mu_mv*mu_mc - (beta_m*x*alpha_m*mu_c*mu_v)/(gamma*alpha*beta) ) ) ) / 2))
  }else{
    if(beta_m=="NA"){
      return(-(((alpha_m+mu_mc+mu_mv) - sqrt( (alpha_m+mu_mc+mu_mv)^2 - 4*(alpha_m*mu_mv + mu_mv*mu_mc - (x*gamma_m*alpha_m*mu_c*mu_v)/(gamma*alpha*beta) ) ) ) / 2))
      }else{
        if(alpha_m=="NA"){
          return( -(((x+mu_mc+mu_mv) - sqrt( (x+mu_mc+mu_mv)^2 - 4*(x*mu_mv + mu_mv*mu_mc - (beta_m*gamma_m*x*mu_c*mu_v)/(gamma*alpha*beta) ) ) ) / 2))
      
          }
  }
  }
  
}



#************fitness boundary for gamma and beta***************

beta_m <- seq(10^-10,10^-6,10^-10) # range of values for mutant virus probabillity of susceptible cell infection
alpha_m <- seq(1/120,1/2,1/200)
fit.bound1 <- sapply(beta_m,function(y){
  uniroot.all(fun,c(0,10^3*24),beta_m=y,gamma_m="NA")
})
fit.bound2 <- sapply(alpha_m,function(y){
  uniroot.all(fun,c(10^-10,10^-2),alpha_m=y,beta_m="NA")
})

pdf(file="fig6.pdf",width=6,height=4)
par(mfcol=c(1,2),mar=c(4,4,1,1),cex=0.6)
plot(beta_m
     ,fit.bound1
     ,type="l"
     ,ylab=expression(gamma)
     ,xlab=expression(beta)
     ,main="a"
     ,bty="n")

plot(alpha_m
     ,fit.bound2
     ,type="l"
     ,ylab=expression(beta)
     ,xlab=expression(alpha)
     ,main="b"
     ,bty="n")

dev.off()
#*********************************************************
