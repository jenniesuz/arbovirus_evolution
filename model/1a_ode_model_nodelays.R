# Between cell virus transmission
# models the presence of two types of virus - acute and persistent
# assuming no delay in budding virus and
# constant rate of apoptosis

# model in hours
#*********************PARAMETERS*************************************
params.nd <- function(     
  
  c_death =  1/120          # background cell death rate, assuming expected cell survival of 5 days
  
  ,double =  24             # 24-48h cell population doubling time 
  
  ,inf = 10^-8              # probability at least one virion infects a cell
  
  ,apoptosis = 1/24         # apoptosis rate for acute virus
  
  ,yield= 10*24             # yeild at apoptosis for acute virus
  
  ,budding = 10             # budding rate for persistent virus 
  
  ,v_death = 0.1            # free virus clearance rate - assumed equal for both virus types
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
mod.nd <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  
  r <- exp(log(2)/double) - 1             # calculate cell growth rate given doubling time
  
  deriv <- rep(NA,5)
  
  deriv[1] <- r*S - inf*S*V_p - inf*S*V_a - c_death*S      # susceptible cells
  
  deriv[2] <- inf*S*V_p - c_death*I_p                      # cells infected with persistent virus
  
  deriv[3] <- inf*S*V_a - c_death*I_a - apoptosis*I_a      # cells infected with acute virus
  
  deriv[4] <- budding*I_p - v_death*V_p                    # free persistent virus
  
  deriv[5] <- yield*apoptosis*I_a - v_death*V_a            # free acute virus
  
  return(list(deriv))
})
#*************************************************************

#**************SIMULATE***************************************
simPop.nd <- function(init=initial, tseq = times, modFunction=mod.nd, parms = params.nd()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#****************************************************************

