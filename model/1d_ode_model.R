
# Model of between cell virus infection with persistent delay
library("ggplot2")
library("gridExtra")
library(PBSddesolve) 
# model in hours
#*********************PARAMETERS*************************************
params <- function(     
  
  c_death =  1/120          # background cell death rate, assuming expected cell survival of 5 days
  
  ,double =  24             # 24-48h population doubling time # if population doubles in 24-48hs doubling time = ln(2)/ln(1+r)
  
  ,inf = 10^-8              # prob virion infection / hr - 
  
  ,v_death = 0.1            # free virus clearance rate
  
  ,apoptosis = 24         # apoptosis rate - yeild at apoptosis calculated as 1/apoptosis * budding
  
 #  ,budding_p = 10           # budding rate - also used to calculate virus yield at apoptosis
  
 #  ,budding_a = 10           #
  
  ,delay_p = 4              # delay in budding rate for persistent virus
  
  ,delay_a = 4

 # ,K=10^6
)
return(as.list(environment()))
#***************************************************************

#*****************INITIAL CONDITIONS*****************************
initial <- c(S = 10^6      # number of susceptible cells at start
             
             ,I_p = 0      # number of cells infected with persistent virus at start
             
             ,I_a = 0      # number of cells infected with acute virus at start
             
             ,V_p = (10^6)*0.01      # number of  persistent virions at start MOI of 1:1
             
             ,V_a = (10^6)*0.01      # number of acute virions at start - assuming as in competition experiment both added in equal amounts
)

times <- seq(0,72,1)               # times to solve at
#**************************************************************



#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  budding_p <- exp(delay_p/2)
  
  rep <- exp(delay_a/2)
  budding_a <- (rep/5)
  yield <- (apoptosis)*(rep/5)*4
  
  tlag_a <- tt - delay_a                    # previous time-point
  if (tlag_a <= 0){
    lag_a <- c(0,0,0,0,0)               # if before time 1 then use initial conditions
  }else {
    lag_a <- pastvalue(tlag_a)  
  }
  
  tlag_p <- tt - delay_p                    # previous time-point
  if (tlag_p <= 0){
    lag_p <- c(0,0,0,0,0)               # if before time 1 then use initial conditions
  }else {
    lag_p <- pastvalue(tlag_p)  
  }
  
  r <- exp(log(2)/double) - 1             # calculate cell growth rate given doubling time
  
  deriv <- rep(NA,5)
  
#  deriv[1] <- (r-c_death)*S*(1-((S)/K)) - inf*S*V_p - inf*S*V_a      # include density-dependence in susceptible
  
  deriv[1] <- r*S - inf*S*V_p - inf*S*V_a - c_death*S      # S
  
  deriv[2] <- inf*S*V_p - c_death*I_p                      # cells infected with persistent virus
  
  deriv[3] <- inf*S*V_a - c_death*I_a - inf*lag_a[1]*lag_a[5]*exp(-c_death*apoptosis)    # cells infected with acute virus
  
  deriv[4] <- budding_p*(lag_p[2] - lag_p[2]*c_death*delay_p) - v_death*V_p  # free persistent virus
  
  deriv[5] <-  budding_a*(lag_a[3] - lag_a[3]*c_death*delay_a) + yield*inf*lag_a[1]*lag_a[5]*exp(-c_death*apoptosis) - v_death*V_a    # free acute virus
  
  return(list(deriv))
})
#*************************************************************

#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(dde(init, tseq, modFunction, parms=parms,hbsize=50000))
  return(simDat)
}
#****************************************************************

# apoptosis can't equal 0
persist.delay <- simPop(parms=params(delay_a=8,delay_p=4,apoptosis=1/12))

par(mfcol=c(1,2))
p1 <- ggplot(persist.delay, aes(x=time,y=log10(V_p))) +
  geom_line(aes(x=time,y=log10(V_a),color="Acute")) +
  geom_line(aes(x=time,y=log10(V_p),color="Persistent"),linetype=2) +
  scale_color_manual(name="",values=c("Acute"="black","Persistent"="black")) +
  guides(color=guide_legend(override.aes=list(linetype=c(1,2)))) +
  labs( y=expression(paste("Number of virions (",log[10],")")),x="Time (hours)") +
  theme(legend.position="none")

p2 <- ggplot(persist.delay, aes(x=time,y=I_p)) +
  geom_line(aes(x=time,y=I_a,color="Acute")) +
  geom_line(aes(x=time,y=I_p,color="Persistent"),linetype=2) +
  scale_color_manual(name="",values=c("Acute"="black","Persistent"="black")) +
  guides(color=guide_legend(override.aes=list(linetype=c(1,2)))) +
  labs( y="Number of infected cells",x="Time (hours)") +
  theme(legend.position="none")

grid.arrange(p1, p2, nrow = 1)
#*************PLOTS**************************************


par(mfcol=c(1,2),cex=0.5)

plot(persist.delay$time
     ,persist.delay$I_p
     ,bty="n"
     ,type="l"
     ,col="red"
     ,ylim=c(1,max(persist.delay$I_p))
     ,xlab="Time (hours) post infection"
     ,ylab=expression(paste( " Number of infected cells"))
     ,xlim=c(0,max(times))
     ,log="y"
)
legend("topright",legend=c("Persistent virus","Acute virus"),bty="n",col=c("red","blue"),lty=1)
par(new=T)
plot(persist.delay$time
     ,persist.delay$I_a
     ,bty="n"
     ,type="l"
     ,col="blue"
     ,ylim=c(1,max(persist.delay$I_p))
     ,xaxt="n"
     ,yaxt="n"
     ,xlab=" "
     ,ylab=" "
     ,xlim=c(0,max(times))
     ,log="y"
)



plot(persist.delay$time
     ,log10(persist.delay$V_p)
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
plot(persist.delay$time
     #,log10(persist.delay$S+1)
     ,log10(persist.delay$V_a)
     ,bty="n"
     ,type="l"
     ,col="blue"
     ,ylim=c(2,10)
     #,xaxt="n"
     #,yaxt="n"
     ,xlab=" "
     ,ylab=" "
     ,xlim=c(0,max(times))
     #,ylim=c(1,log10(max(persist.delay$S)))
)
#**************************************************

