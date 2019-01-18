
# Model of between cell virus infection without delays
library(ggplot2)
library(gridExtra)
library(deSolve)
theme <-   theme(panel.border = element_blank()
                 ,axis.line = element_line(color = 'black')
                 ,text=element_text(size=6)
                 ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
                 ,axis.text=element_text(size=5)
                 ,legend.key.size = unit(0.8,"line")
                 ,legend.background = element_blank()
                 ,legend.text=element_text(size=5)
                 ,legend.position =c(0.83,0.9)
                 ,legend.title = element_blank()
)
 
# model in hours
#*********************PARAMETERS*************************************
params <- function(     
  
  c_death =  1/120          # background cell death rate, assuming expected cell survival of 5 days
  
  ,double =  24             # 24-48h population doubling time # if population doubles in 24-48hs doubling time = ln(2)/ln(1+r)
  
  ,inf = 10^-8              # prob virion infection / hr - 
  
  ,v_death = 0.1            # free virus clearance rate
  
  ,apoptosis = 1/24         # apoptosis rate - yeild at apoptosis calculated as 1/apoptosis * budding
  
  ,budding = 10             # budding rate 
  
  ,yield= 10*24
  
  ,apoptosis_m = 1/24         # apoptosis rate - yeild at apoptosis calculated as 1/apoptosis * budding
  
  ,budding_m = 10             # budding rate 
  
  ,yield_m = 10*24
)
return(as.list(environment()))
#***************************************************************

#*****************INITIAL CONDITIONS*****************************
initial <- c(S = 10^6      # number of susceptible cells at start
             
             ,I = 0      # number of cells infected with persistent virus at start
             
             ,I_m = 0      # number of cells infected with acute virus at start
             
             ,V = (10^6)*0.1      # number of  persistent virions at start MOI of 1:1
             
             ,V_m = (10^6)*0.1      # number of acute virions at start - assuming as in competition experiment both added in equal amounts
)

times <- seq(0,72,1)               # times to solve at
#**************************************************************



#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  r <- exp(log(2)/double) - 1             # calculate cell growth rate given doubling time
  
  deriv <- rep(NA,5)
  
  deriv[1] <- r*S - inf*S*V_m - inf*S*V - c_death*S      # S
  
  deriv[2] <- inf*S*V - c_death*I - apoptosis*I           # cells infected with virus
  
  deriv[3] <- inf*S*V_m - c_death*I_m - apoptosis_m*I_m      # cells infected with mutant virus
  
  deriv[4] <- budding*I + yield*apoptosis*I - v_death*V                    # free virus
  
  deriv[5] <- budding_m*I_m + yield_m*apoptosis_m*I_m - v_death*V_m    # free mutant virus
  
  return(list(deriv))
})
#*************************************************************



#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#****************************************************************


no.delays <- simPop(parms=params(c_death=1/120
                                 ,apoptosis=1/24
                                 ,yield=10*24
                                 ,budding=0
                                 ,budding_m=0
                                 ,apoptosis_m=1/240
                                 ,yield_m=10*240))

library("ggplot2")
#pdf(file="fig_ms_1.pdf",width=5,height=4)
p1<-ggplot(no.delays, aes(x=time,y=log10(V_m))) +
 geom_line(aes(x=time,y=log10(V))) +
 geom_line(aes(x=time,y=log10(V_m)),linetype=2) +
 scale_color_manual(name="",values=c("resident"="black","mutant"="black")) +
 labs(y=expression(paste("Number of virions (",log[10],")")),x="Time (hours)") +
 theme_set(theme_bw())  +
 theme
p1
#dev.off()
