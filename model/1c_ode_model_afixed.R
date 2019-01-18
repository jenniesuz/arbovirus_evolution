
# Model of between cell virus infection fixed time to apoptosis

library(PBSddesolve) 
# model in hours
#*********************PARAMETERS*************************************
params <- function(     
  
  c_death =  1/120          # background cell death rate, assuming expected cell survival of 5 days
  
  ,double =  24             # 24-48h population doubling time # if population doubles in 24-48hs doubling time = ln(2)/ln(1+r)
  
  ,inf = 10^-8              # prob virion infection / hr - 
  
  ,v_death = 0.1            # free virus clearance rate
  
  ,budding = 10             # budding rate - also used to calculate virus yield at apoptosis
  
  ,apoptosis = 24            # time to apoptosis for acute infection
  
  #,K=10^6                  # number of susceptible cells at 'carrying capacity' not currently used
)
return(as.list(environment()))
#***************************************************************

#*****************INITIAL CONDITIONS*****************************
initial <- c(S = 10^6      # number of susceptible cells at start
             
             ,I_p = 0      # number of cells infected with persistent virus at start
             
             ,I_a = 0      # number of cells infected with acute virus at start
             
             ,V_p = (10^6)*0.1      # number of  persistent virions at start MOI of 1:1
             
             ,V_a = (10^6)*0.1      # number of acute virions at start - assuming as in competition experiment both added in equal amounts
)

times <- seq(0,72,1)               # times to solve at
#**************************************************************



#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  
  tlag <- tt - apoptosis                    # previous time-point
  
  if (tlag <= 0){
    
    lag <- c(0,0,0,0,0)               # if before time 1 then use initial conditions
    
  }else {
    
    lag <- pastvalue(tlag)  
    
  }
  
  
  r <- exp(log(2)/double) - 1             # calculate cell growth rate given doubling time
  
  deriv <- rep(NA,5)
  
 # deriv[1] <- (r-c_death)*S*(1-((S)/K)) - - inf*S*V_m - inf*S*V_r       # include density-dependence in susceptible
  
  deriv[1] <- r*S - inf*S*V_p - inf*S*V_a - c_death*S      # S
  
  deriv[2] <- inf*S*V_p - c_death*I_p                      # cells infected with persistent virus
  
  deriv[3] <- inf*S*V_a - c_death*I_a - inf*lag[1]*lag[5]*exp(c_death*apoptosis)     # cells infected with acute virus
  
  deriv[4] <- budding*I_p - v_death*V_p  # free persistent virus
  
  deriv[5] <- (budding*(apoptosis))*inf*lag[1]*lag[5]*exp(c_death*apoptosis) - v_death*V_a    # free acute virus
  
  return(list(deriv))
})
#*************************************************************

#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(dde(init, tseq, modFunction, parms=parms,hbsize=50000))
  return(simDat)
}
#***************************************************************
acute.delay <- simPop(parms=params(apoptosis=40))






#*************PLOTS**************************************
par(mfcol=c(1,2),cex=0.5)

plot(acute.delay$time
     ,acute.delay$I_p
     ,bty="n"
     ,type="l"
     ,col="red"
     ,ylim=c(1,max(acute.delay$I_p))
     ,xlab="Time (hours) post infection"
     ,ylab=expression(paste( " Number of infected cells"))
     ,xlim=c(0,max(times))
     ,log="y"
)
legend("topright",legend=c("Persistent virus","Acute virus"),bty="n",col=c("red","blue"),lty=1)
par(new=T)
plot(acute.delay$time
     ,acute.delay$I_a
     ,bty="n"
     ,type="l"
     ,col="blue"
     ,ylim=c(1,max(acute.delay$I_p))
     ,xaxt="n"
     ,yaxt="n"
     ,xlab=" "
     ,ylab=" "
     ,xlim=c(0,max(times))
     ,log="y"
)



plot(acute.delay$time
     ,log10(acute.delay$V_p)
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
plot(acute.delay$time
     #,log10(acute.delay$S+1)
     ,log10(acute.delay$V_a)
     ,bty="n"
     ,type="l"
     ,col="blue"
     ,ylim=c(2,10)
     #,xaxt="n"
     #,yaxt="n"
     ,xlab=" "
     ,ylab=" "
     ,xlim=c(0,max(times))
     #,ylim=c(1,log10(max(acute.delay$S)))
)
#**************************************************

