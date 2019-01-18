# Between cell virus transmission
# models the presence of two types of virus - acute and persistent
# assuming delay in budding virus and
# fixed time to apoptosis

# model in hours
library(PBSddesolve) 
source("1a_ode_model_nodelays.R") # for initial conditions and times to solve at

#*********************PARAMETERS*************************************
params.b <- function(     
  
  c_death =  1/120          # background cell death rate, assuming expected cell survival of 5 days
  
  ,double =  24             # 24-48h population doubling time 
  
  ,inf = 10^-8              # probability at least one virion infects a cell assumed same for both virus types
  
  ,budding_p = 10           # budding rate 
  
  ,delay_p = 12             # delay in budding for persistent virus
  
  ,delay_a = 24             # time to apoptosis for acute virus
  
  ,yield_a = 10*24          # yield at apoptosis for acute virus

  ,v_death = 0.1            # free virus clearance rate assumed same for both virus types
  
)
return(as.list(environment()))
#***************************************************************

#****************MODEL*****************************************
mod.b <- function(tt,yy,parms) with(c(parms,as.list(yy)), {

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
  
  deriv[4] <- budding_p*lag_p[2]*exp(-c_death*delay_p) - v_death*V_p  # free persistent virus
  
  deriv[5] <-  yield_a*inf*lag_a[1]*lag_a[5]*exp(-c_death*apoptosis) - v_death*V_a #+budding_a*(lag_a[3] - lag_a[3]*c_death*delay_a)    # free acute virus
  
  return(list(deriv))
})
#*************************************************************

#**************SIMULATE***************************************
simPop.b <- function(init=initial, tseq = times, modFunction=mod.b, parms = params.b()) {
  simDat <- as.data.frame(dde(init, tseq, modFunction, parms=parms,hbsize=50000))
  return(simDat)
}
#****************************************************************



# 
# #***********************************************************************
# pd1 <- simPop(parms=params(delay_a=24
#                            ,delay_p=24
#                            ,budding_p=10
#                            ,yield_a=10*10
#                            ,apoptosis=24))
# 
# pd2 <- simPop(parms=params(delay_a=24
#                            ,delay_p=24
#                            ,budding_p=30
#                            ,yield_a=10*30
#                            ,apoptosis=24))
# 
# 
# pd3 <- simPop(parms=params(delay_a=24
#                            ,delay_p=24
#                            ,budding_p=60
#                            ,yield_a=10*60
#                            ,apoptosis=24))
# 
# pd4 <- simPop(parms=params(delay_a=24
#                            ,delay_p=24
#                            ,budding_p=90
#                            ,yield_a=10*90
#                            ,apoptosis=24))
# 
# pd5 <- simPop(parms=params(delay_a=24
#                            ,delay_p=24
#                            ,budding_p=120
#                            ,yield_a=10*120
#                            ,apoptosis=24))
# 
# combined.mods <- rbind.data.frame(pd1,pd2,pd3,pd4,pd5)
# combined.mods$nam <- c(rep("pd1",length(pd1[,1]))
#                        ,rep("pd2",length(pd2[,1]))
#                        ,rep("pd3",length(pd3[,1]))
#                        ,rep("pd4",length(pd4[,1]))
#                        ,rep("pd5",length(pd4[,1])))
# 
# #pdf(file="fig_ms_1.pdf",width=5,height=4)
# ggplot(combined.mods, aes(x=time,y=V_p)) +
#   geom_line(aes(x=time,y=log10(V_p),col=nam)) +
#   geom_line(aes(x=time,y=log10(V_a),col=nam),linetype=2) +
#   scale_color_manual(values=as.character(cols)) +
#   theme_set(theme_bw())  +
#   theme
# #dev.off()





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

