
# Model of persistent viral dynamics

require(deSolve)          

# model in hours


#*********************PARAMETERS*************************************
params <- function(     
  
   c_death =  1/120         # background cell death rate - assumed same for susceptible and infected cells
  
  ,double =  24             # 24-48h population doubling time # if population doubles in 24-48hs doubling time = ln(2)/ln(1+r)
  
  ,inf = 10^-8             # virus infection rate

  ,eclipse = 12             # duration of eclipse phase before infected cells start producing virus
  
  ,v_death = 0.1            # free virus clearance rate
  
  ,release = 10^2            # virus release rate from infected cells
  
  ,K=10^6

  )
return(as.list(environment()))
#***************************************************************



#*****************INITIAL CONDITIONS****************************
initial <- c(S = 10^6   # number of susceptible cells at start
             
             ,I = 0     # number of infected cells at start
             
             ,V = 10^6)  # number of free virions at start

times <- seq(0,72,1)
#**************************************************************



#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(parms,as.list(yy)), {
  

    
    r <- exp(log(2)/double) - 1             # calculate cell growth rate given doubling time
    
    
  
    deriv <- rep(NA,3)
    
    #deriv[1] <- r*S - inf*S*V  - c_death*S                                 # susceptible cells
    deriv[1] <- (r-c_death)*S*(1-((S)/K)) - inf*S*V
    
    deriv[2] <- inf*S*V - c_death*I              # infected cells
    
    deriv[3] <- release*I*(1/eclipse)  - v_death*V                        # free virus
    
  
  return(list(deriv))
})
#*************************************************************



#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#****************************************************************



#************PLOT FUN************************************

plot.func <- function(dat=test){

plot(dat$time
     #,log10(dat$I+1)
     ,dat$I
     ,bty="n"
     ,type="l"
     ,col="red"
     ,xlab="Time (hours) post infection"
     #,ylab=expression(paste('Log'[10], " number of cells"))
     ,ylab="Number of cells"
     ,main=paste("cell death = ",round(params()$c_death,3)
                 ,"double =",params()$double 
                 ,"inf rate =",params()$inf
                 ,"\n"
                 ,"eclipse =",params()$eclipse
                 ,"virus death =",params()$v_death
                 ,"virus release =", params()$release)
     ,xlim=c(0,max(times))
     ,ylim=c(min(dat$S),max(dat$I))
)

par(new=T)
plot(dat$time
     #,log10(dat$S+1)
     ,dat$S
     ,bty="n"
     ,type="l"
     ,col="blue"
     ,xaxt="n"
     ,yaxt="n"
     ,xlab=" "
     ,ylab=" "
     ,xlim=c(0,max(times))
     ,ylim=c(min(dat$S),max(dat$S))
     #,ylim=c(1,log10(max(dat$S)))
)

legend("bottomright",legend=c("Infected","Susceptible"),bty="n",col=c("red","blue"),lty=1)
# 
# plot(dat$time
#      ,(dat$I)/dat$S
#      ,bty="n"
#      ,type="l"
#      ,col="green"
#      ,xlab="Time (hours)"
#      ,ylab="% cells infected"
#      ,xlim=c(0,max(times))
#      ,ylim=c(0,max(dat$I/dat$S))
# )
# # 

plot(dat$time
     ,log10(dat$V)
     ,bty="n"
     ,type="l"
     ,col="black"
     # ,log="y"
     ,ylim=c(min(log10(dat$V)),max(log10(dat$V)))
     ,xlim=c(0,max(times))
     ,xlab="Time (hours) post infection"
     ,ylab=expression(paste('Log'[10], " virus"))
)
}

#**************************************************



#***********Run and plot*********************
par(mfcol=c(2,2),cex=0.6)

times<-seq(0,72,1)
test <- simPop()
plot.func(test)

test <- simPop(parms=params(inf=10^-9))
plot.func(test)
  
  
test <- simPop(parms=params(inf=10^-6,v_death=0.6,c_death=1/72))
plot.func(test)
  #*********************************************



#**************Run model for range of values of each variable*******************

subsParms <- function(inf.rate, release.rate, fixed.params = params()){
  model.output <- simPop(parms = params(inf=inf.rate,release=release.rate) )
  return(log10(max(model.output$V)))
}


simPop_VEC <- Vectorize(subsParms, list("inf.rate","release.rate"))


inf.seq <- seq(10^-14, 10^-2, l = 10)
inf.seq


release.seq <- seq(0.1, 1, l = 20)
release.seq


mat <- outer(inf.seq, release.seq, simPop_VEC) 


#**********************************************************

res <- lapply(release.seq,function(x){
  test <- simPop(parms=params(release = x))
  
  
  plot(test$time
       ,(test$I)/test$S*100
       ,bty="n"
       ,type="l"
       ,col="green"
       ,xlab="Time (hours)"
       ,ylab="% cells infected"
       ,xlim=c(0,max(times))
       ,ylim=c(0,100)
  )
  
  
  plot(test$time
       ,log10(test$V)
       ,bty="n"
       ,type="l"
       ,col="black"
       # ,log="y"
       ,ylim=c(min(log10(test$V)),max(log10(test$V)))
       ,xlim=c(0,max(times))
       ,xlab="Time (hours) post infection"
       ,ylab="Log10 virus"
  )
})




#**************************************************************************

