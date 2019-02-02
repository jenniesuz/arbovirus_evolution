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
  
  ,apoptosis = 24           # apoptosis 
  
  ,apoptosis_m = 24
  
  ,budding = 10           # budding rate - also used to calculate virus yield at apoptosis
  
  ,budding_m = 10           #
  
  ,delay = 12              # delay in budding rate for persistent virus
  
  ,delay_m = 12
  
  ,yield = 10*24
  
  ,yield_m = 10*24

)
return(as.list(environment()))
#***************************************************************

#*****************INITIAL CONDITIONS*****************************
initial <- c(S = 10^6      # number of susceptible cells at start
             
             ,I = 0      # number of cells infected with virus at start
             
             ,I_m = 0      # number of cells infected with mutant virus at start
             
             ,V =  (10^6)*0.1      # number of  persistent virions at start MOI of 1:1
             
             ,V_m = (10^6)*0.1      # number of acute virions at start - assuming as in competition experiment both added in equal amounts
)

times <- seq(0,240,1)               # times to solve at
#**************************************************************

#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
 tlag <- tt - delay                    # previous time-point
  if (tlag <= 0){
    lag <- c(0,0,0,0,0)                # if before time 1 then use initial conditions
  }else {
    lag <- pastvalue(tlag)  
  }
  
  tlag_m <- tt - delay_m                   # previous time-point
  if (tlag_m <= 0){
    lag_m <- c(0,0,0,0,0)                # if before time 1 then use initial conditions
  }else {
    lag_m <- pastvalue(tlag_m)  
  }
  
  tlag1 <- tt - apoptosis                    # previous time-point
  if (tlag1 <= 0){
    lag1 <- c(0,0,0,0,0)                # if before time 1 then use initial conditions
  }else {
    lag1 <- pastvalue(tlag1)  
  }
  
  tlag1_m <- tt - apoptosis_m                   # previous time-point
  if (tlag1_m <= 0){
    lag1_m <- c(0,0,0,0,0)                # if before time 1 then use initial conditions
  }else {
    lag1_m <- pastvalue(tlag1_m)  
  }
  
  r <- exp(log(2)/double) - 1             # calculate cell growth rate given doubling time
  
  deriv <- rep(NA,5)
  
  deriv[1] <- r*S - inf*S*V - inf*S*V_m - c_death*S      # S
  
  deriv[2] <- inf*S*V - c_death*I - inf*lag1[1]*lag1[4]*exp(-c_death*apoptosis)                     
  
  deriv[3] <- inf*S*V_m - c_death*I_m - inf*lag1_m[1]*lag1_m[5]*exp(-c_death*apoptosis_m)  
  
  deriv[4] <- budding*(lag[2] - lag[2]*c_death*delay) + yield*inf*lag1[1]*lag1[4]*exp(-c_death*apoptosis) - v_death*V 
  
  deriv[5] <-  budding*(lag_m[3] - lag_m[3]*c_death*delay_m) + yield_m*inf*lag1_m[1]*lag1_m[5]*exp(-c_death*apoptosis_m) - v_death*V_m 
  
  return(list(deriv))
})
#*************************************************************

#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(dde(init, tseq, modFunction, parms=parms,hbsize=50000))
  return(simDat)
}
#****************************************************************


delays <- simPop(tseq= seq(0,120,100)
                 ,parms=params(c_death=1/120
                                 ,apoptosis=2
                                 ,yield=10
                                 ,budding=0
                                 ,budding_m=0
                                 ,delay=2
                                 ,delay_m=60
                                 ,apoptosis_m=60
                                 ,yield_m=10))

library("ggplot2")
#pdf(file="fig_ms_1.pdf",width=5,height=4)
p1<-ggplot(delays, aes(x=time,y=log10(V_m))) +
  geom_line(aes(x=time,y=log10(V))) +
  geom_line(aes(x=time,y=log10(V_m)),linetype=2) +
  scale_color_manual(name="",values=c("resident"="black","mutant"="black")) +
  labs(y=expression(paste("Number of virions (",log[10],")")),x="Time (hours)") +
  theme_set(theme_bw())  +
  theme
p1





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

