source("1a_ode_model.R")
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
  
  ,delay_p = 12              # delay in budding rate for persistent virus
  
  ,delay_a = 12

 # ,K=10^6
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
  budding_p <- 10 #exp(delay_p/2)
  
  rep <- 10 #exp(delay_a/2)
  budding_a <- 0 #(rep/5)
  yield <- rep*apoptosis #(apoptosis)*(rep/5)*4
  
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
  
  deriv[5] <-  yield*inf*lag_a[1]*lag_a[5]*exp(-c_death*apoptosis) - v_death*V_a #+budding_a*(lag_a[3] - lag_a[3]*c_death*delay_a)    # free acute virus
  
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

combined.mods <- cbind.data.frame(time=c(no.delays$time,persist.delay$time)
                                  ,p=c(no.delays$V_p,persist.delay$V_p)
                                  ,a=c(no.delays$V_a,persist.delay$V_a)
                                  ,mod=c(rep("Model 1",length(no.delays$time)),rep("Model 2",length(no.delays$time)))
                                 )

pdf(file="fig_ms_1.pdf",width=5,height=4)
ggplot(combined.mods, aes(x=time,y=p)) +
  geom_line(aes(x=time,y=log10(p))) +
  geom_line(aes(x=time,y=log10(a)),linetype=2) +
  facet_wrap( ~ mod,labeller = label_wrap_gen(width=30,multi_line=FALSE)) +
  labs( y=expression(paste("Number of virions (",log[10],")")),x="Time (hours)") +
  theme(legend.position="none")
dev.off()


#p2 <- ggplot(persist.delay, aes(x=time,y=log10(V_p))) +
#  geom_line(aes(x=time,y=log10(V_a),color="Acute")) +
#  geom_line(aes(x=time,y=log10(V_p),color="Persistent"),linetype=2) +
#  scale_color_manual(name="",values=c("Acute"="black","Persistent"="black")) +
#  guides(color=guide_legend(override.aes=list(linetype=c(1,2)))) +
#  labs( y=expression(paste("Number of virions (",log[10],")")),x="Time (hours)") +
#  theme(legend.position="none")

#p2 <- ggplot(persist.delay, aes(x=time,y=I_p)) +
#  geom_line(aes(x=time,y=I_a,color="Acute")) +
#  geom_line(aes(x=time,y=I_p,color="Persistent"),linetype=2) +
#  scale_color_manual(name="",values=c("Acute"="black","Persistent"="black")) +
#  guides(color=guide_legend(override.aes=list(linetype=c(1,2)))) +
#  labs( y="Number of infected cells",x="Time (hours)") +
#  theme(legend.position="none")

#grid.arrange(p1, p2, nrow = 1)
#*************PLOTS**************************************
#pdf(file="fig_ms_1.pdf",width=5,height=4)
#grid.arrange(p1, p2, nrow = 1)

#dev.off()

