
require(deSolve) 

# model in hours
# amount of free virus here is not influenced by viral cell entry which it should be

#*********************PARAMETERS*************************************
params <- function(     
  
  c_death =  1/120          # background cell death rate
  
  ,double =  24             # 24-48h population doubling time # if population doubles in 24-48hs doubling time = ln(2)/ln(1+r)
  
  ,inf = 10^-8             # prob virion infection / hr 
  
  ,apoptosis = 1/24         # apoptosis rate
  
  ,v_death = 0.1            # free virus clearance rate
  
  ,budding = 10            # budding rate - also used to calculate virus yield at apoptosis
  
  ,K=10^6
  
)
return(as.list(environment()))
#***************************************************************



#**************************************************************
initial <- c(S = 10^6      # number of susceptible cells at start
             
             ,I_m = 0      # number of cells infected with mutant virus at start
             
             ,I_r = 0      # number of cells infected with 'resident' virus at start
             
             ,V_m = (10^6)*0.1      # number of  mutant virions at start MOI of 1:1
             
             ,V_r = (10^6)*0.1      # number of resident virions at start
)

times <- seq(0,72,1)
#**************************************************************



#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  r <- exp(log(2)/double) - 1             # calculate cell growth rate given doubling time
  
  #r <- 1/120
  deriv <- rep(NA,5)
  
 # deriv[1] <- (r-c_death)*S*(1-((S)/K)) - - inf*S*V_m - inf*S*V_r       # include density-dependence in susceptible
  
  deriv[1] <- r*S - inf*S*V_m - inf*S*V_r - c_death*S      # S
  
  deriv[2] <- inf*S*V_m - c_death*I_m                      # cells infected with mutant budding virus
  
  deriv[3] <- inf*S*V_r - c_death*I_r - apoptosis*I_r      # cells infected with resident acute virus
  
  deriv[4] <- budding*I_m - v_death*V_m                    # free m virus
  
  deriv[5] <- (budding*(1/apoptosis))*apoptosis*I_r - v_death*V_r    # free r virus
  
  return(list(deriv))
})
#*************************************************************



#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#****************************************************************


test <- simPop(parms=params(apoptosis=1/24))


#*************PLOT**************************************
pdf(file="fig11.pdf",width=5,height=4)

par(mfcol=c(1,2),cex=0.5)

plot(test$time
     #,log10(test$I_m+1)
     ,test$I_m
     ,bty="n"
     ,type="l"
     ,col="red"
     ,ylim=c(0,max(test$I_m))
    # ,ylim=c(0,250)
     ,xlab="Time (hours) post infection"
     ,ylab=expression(paste( " Number of infected cells"))
     ,xlim=c(0,max(times))
     # ,ylim=c(1,log10(max(test$S)))
)
legend("topright",legend=c("Persistent virus","Acute virus"),bty="n",col=c("red","blue"),lty=1)
par(new=T)
plot(test$time
     #,log10(test$I_r+1)
     ,test$I_r
     ,bty="n"
     ,type="l"
     ,col="blue"
     ,ylim=c(0,max(test$I_m))
     #,xaxt="n"
     #,yaxt="n"
     ,xlab=" "
     ,ylab=" "
     ,xlim=c(0,max(times))
     #,ylim=c(1,log10(max(test$S)))
)



plot(test$time
     ,log10(test$V_m)
     ,bty="n"
     ,type="l"
     ,col="red"
     # ,log="y"
     ,ylim=c(2,10)
     ,main=params()$g
     ,xlim=c(0,max(times))
     ,xlab="Time (hours) post infection"
     ,ylab=expression(paste('Log'[10], " Virus"))
)
par(new=T)
plot(test$time
     #,log10(test$S+1)
     ,log10(test$V_r)
     ,bty="n"
     ,type="l"
     ,col="blue"
     ,ylim=c(2,10)
     #,xaxt="n"
     #,yaxt="n"
     ,xlab=" "
     ,ylab=" "
     ,xlim=c(0,max(times))
     #,ylim=c(1,log10(max(test$S)))
)
dev.off()
#**************************************************

