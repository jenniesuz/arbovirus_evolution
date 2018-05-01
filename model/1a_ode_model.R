
# Model of between cell virus infection without delays

require(deSolve) 
# model in hours
#*********************PARAMETERS*************************************
params <- function(     
  
  c_death =  1/120          # background cell death rate, assuming expected cell survival of 5 days
  
  ,double =  24             # 24-48h population doubling time # if population doubles in 24-48hs doubling time = ln(2)/ln(1+r)
  
  ,inf = 10^-8              # prob virion infection / hr - 
  
  ,apoptosis = 1/24         # apoptosis rate - yeild at apoptosis calculated as 1/apoptosis * budding
  
  ,v_death = 0.1            # free virus clearance rate
  
  ,budding = 10             # budding rate - also used to calculate virus yield at apoptosis
  
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
  
  r <- exp(log(2)/double) - 1             # calculate cell growth rate given doubling time
  
  deriv <- rep(NA,5)
  
 # deriv[1] <- (r-c_death)*S*(1-((S)/K)) - - inf*S*V_m - inf*S*V_r       # include density-dependence in susceptible
  
  deriv[1] <- r*S - inf*S*V_p - inf*S*V_a - c_death*S      # S
  
  deriv[2] <- inf*S*V_p - c_death*I_p                      # cells infected with persistent virus
  
  deriv[3] <- inf*S*V_a - c_death*I_a - apoptosis*I_a      # cells infected with acute virus
  
  deriv[4] <- budding*I_p - v_death*V_p                    # free persistent virus
  
  deriv[5] <- (budding*(1/apoptosis))*apoptosis*I_a - v_death*V_a    # free acute virus
  
  return(list(deriv))
})
#*************************************************************



#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#****************************************************************


no.delays <- simPop(parms=params(apoptosis=1/24))



library("ggplot2")
#pdf(file="fig_ms_1.pdf",width=5,height=4)
p1<-ggplot(no.delays, aes(x=time,y=log10(V_p))) +
 geom_line(aes(x=time,y=log10(V_a))) +
 geom_line(aes(x=time,y=log10(V_p)),linetype=2) +
# scale_color_manual(name="",values=c("Acute"="black","Persistent"="black")) +
# guides(color=guide_legend(override.aes=list(linetype=c(1,2)))) +
 labs( y=expression(paste("Number of virions (",log[10],")")),x="Time (hours)") 
#dev.off()

# 
# #*************PLOTS**************************************
# pdf(file="fig11.pdf",width=5,height=4)
# 
# par(mfcol=c(1,2),cex=0.5)
# 
# plot(no.delays$time
#      ,no.delays$I_p
#      ,bty="n"
#      ,type="l"
#      ,col="red"
#      ,ylim=c(1,max(no.delays$I_p))
#      ,xlab="Time (hours) post infection"
#      ,ylab=expression(paste( " Number of infected cells"))
#      ,xlim=c(0,max(times))
#      ,log="y"
# )
# legend("topright",legend=c("Persistent virus","Acute virus"),bty="n",col=c("red","blue"),lty=1)
# par(new=T)
# plot(no.delays$time
#      ,no.delays$I_a
#      ,bty="n"
#      ,type="l"
#      ,col="blue"
#      ,ylim=c(1,max(no.delays$I_p))
#      ,xaxt="n"
#      ,yaxt="n"
#      ,xlab=" "
#      ,ylab=" "
#      ,xlim=c(0,max(times))
#      ,log="y"
# )
# 
# 
# 
# plot(no.delays$time
#      ,log10(no.delays$V_p)
#      ,bty="n"
#      ,type="l"
#      ,col="red"
#      # ,log="y"
#      ,ylim=c(2,10)
#      ,main=params()$g
#      ,xlim=c(0,max(times))
#      ,xlab="Time (hours) post infection"
#      ,ylab=expression(paste('Log'[10], " Virus"))
# )
# par(new=T)
# plot(no.delays$time
#      #,log10(no.delays$S+1)
#      ,log10(no.delays$V_a)
#      ,bty="n"
#      ,type="l"
#      ,col="blue"
#      ,ylim=c(2,10)
#      #,xaxt="n"
#      #,yaxt="n"
#      ,xlab=" "
#      ,ylab=" "
#      ,xlim=c(0,max(times))
#      #,ylim=c(1,log10(max(no.delays$S)))
# )
# dev.off()
# #**************************************************
# 
