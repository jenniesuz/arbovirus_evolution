
require(deSolve) 

# model in hours
# amount of free virus here is not influenced by viral cell entry which it should be

#*********************PARAMETERS*************************************
params <- function(     
  
  c_death =  1/120          # background cell death rate
  
  ,double =  24             # 24-48h population doubling time # if population doubles in 24-48hs doubling time = ln(2)/ln(1+r)
  
  ,inf = 10^-10              # prob virion infection / hr 
  
  ,apoptosis = 1/24         # apoptosis rate
  
  ,v_death = 0.1            # free virus clearance rate
  
  ,budding = 100
  
  ,K=10^6
  
)
return(as.list(environment()))
#***************************************************************



#**************************************************************
initial <- c(S = 10^6      # number of susceptible cells at start
             
             ,I_m = 0        # number of infected cells at start
             
             ,I_r = 0
             
             ,V_r = (10^6)*0.1      # number of virions at start MOI of 1:1
             
             ,V_m = (10^6)*0.1 
)

times <- seq(0,72,1)
#**************************************************************



#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  r <- exp(log(2)/double) - 1             # calculate cell growth rate given doubling time
  
  #r <- 1/120
  deriv <- rep(NA,5)
  
  #deriv[1] <- (r-c_death)*S*(1-((S)/K)) - inf*S*V    # include density-dependence in susceptible
  
  deriv[1] <- c_death*(S+I) - inf*S*V - c_death*S
  
  deriv[2] <- inf*S*V - c_death*I - apoptosis*I # infected cells
  
  deriv[3] <- max.virions*apoptosis*I - v_death*V               # free virus
  
  deriv[4] <- max.virions*apoptosis*I - v_death*V               # free virus
  
  deriv[5] <- max.virions*apoptosis*I - v_death*V               # free virus
  
  return(list(deriv))
})
#*************************************************************



#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#****************************************************************


test <- simPop(parms=params(apoptosis=1))


#*************PLOT**************************************

par(mfcol=c(2,3),cex=0.8)

plot(test$time
     #,log10(test$I+1)
     ,test$I
     ,bty="n"
     ,type="l"
     ,col="red"
     ,xlab="Time (hours) post infection"
     ,ylab=expression(paste('Log'[10], " number of cells"))
     ,main=paste("cell death = ",round(params()$c_death,3)
                 ,"double =",params()$double 
                 ,"inf rate =",params()$inf
                 ,"\n"
                 ,"eclipse =",params()$eclipse
                 ,"virus death =",params()$v_death
                 ,"virus release =", params()$release)
     
     ,xlim=c(0,max(times))
     # ,ylim=c(1,log10(max(test$S)))
)
abline(h=0)
#par(new=T)
plot(test$time
     #,log10(test$S+1)
     ,test$S
     ,bty="n"
     ,type="l"
     ,col="blue"
     #,xaxt="n"
     #,yaxt="n"
     ,xlab=" "
     ,ylab=" "
     ,xlim=c(0,max(times))
     #,ylim=c(1,log10(max(test$S)))
)



plot(test$time
     ,log10(test$V)
     ,bty="n"
     ,type="l"
     ,col="black"
     # ,log="y"
     ,ylim=c(2,10)
     ,main=params()$g
     ,xlim=c(0,max(times))
     ,xlab="Time (hours) post infection"
     ,ylab=expression(paste('Log'[10], " virus"))
)

#**************************************************

