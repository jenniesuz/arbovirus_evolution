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
    lag_a <- c(0,0,0,0,0)                   # if before time 1 then use initial conditions
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
  
  deriv[1] <- r*S - inf*S*V_p - inf*S*V_a - c_death*S      # Susceptible cells
  
  deriv[2] <- inf*S*V_p - c_death*I_p                      # cells infected with persistent virus
  
  deriv[3] <- inf*S*V_a - c_death*I_a - inf*lag_a[1]*lag_a[5]*exp(-c_death*delay_a)    # cells infected with acute virus
  
  deriv[4] <- budding_p*lag_p[2]*exp(-c_death*delay_p) - v_death*V_p  # free persistent virus
  
  deriv[5] <-  yield_a*inf*lag_a[1]*lag_a[5]*exp(-c_death*delay_a) - v_death*V_a     # free acute virus
  
  return(list(deriv))
})
#*************************************************************

#**************SIMULATE***************************************
simPop.b <- function(init=initial, tseq = times, modFunction=mod.b, parms = params.b()) {
  simDat <- as.data.frame(dde(init, tseq, modFunction, parms=parms,hbsize=50000))
  return(simDat)
}
#****************************************************************


