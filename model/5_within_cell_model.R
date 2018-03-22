
require(deSolve) 


#*********************PARAMETERS*************************************
params <- function(     
  
   k1 = 1/10         # time taken mRNA produced 
  
  ,k2 = 1/10        # time taken mRNA used in production of mature virus 
  
  ,p = 60
  
  ,k2form = "constant"

)
return(as.list(environment()))
#***************************************************************

#**************************************************************
initial <- c( M = 1     # + rna
              
             ,V = 0     # virions
)

times <- seq(0,100,1)
#**************************************************************
par(mfcol=c(2,3))
curve(params()$k2*x,from=0,to=100,bty="n")
curve(params()$k2*exp(x/20),from=0,to=100,bty="n")

#****************MODEL*****************************************
mod <- function(tt,yy,parms) with(c(tt,parms,as.list(yy)), {

  if(k2form=="constant"){
      k2 <- k2
    }else{
      if(k2form=="exp"){
        k2 <- k2*exp(tt/90)
    }else{
      if(k2form=="linear"){
        k2 <- k2*tt
      }
    }
  }

  
  deriv <- rep(NA,2)
  
  deriv[1] <- k1*M - k2*M/p
  
  deriv[2] <- k2*M/p

  return(list(deriv))
})
#*************************************************************



#**************SIMULATE***************************************
simPop <- function(init=initial, tseq = times, modFunction=mod, parms = params()) {
  simDat <- as.data.frame(ode(init, tseq, modFunction, parms=parms))
  return(simDat)
}
#****************************************************************
pdf(file="fig12.pdf",width=5,height=4)

plot.func <- function(x=test$time,y=test$M,z=test$V,k2=1/params()$k2,name="Linear"){
  
   plot(x
       ,y
       ,bty="n"
       ,type="l"
       ,col="red"
       ,xlab="time"
       ,main=name
      # ,ylim=c(0,20000)
       # ,log="y"
       ,ylab="Within-cell mRNA")
  
  plot(x
       ,z
       ,type="l"
       ,col="red"
     #  ,ylim=c(0,3000)
       ,xlab="time"
       ,ylab="Within-cell mature virus"
       ,main=name
        #,log="y"
       ,bty="n")
  
}

par(mfcol=c(2,3),cex=0.5)

test <- simPop(parms=params(k2form="exp"))
plot.func(name="Exponential")

test <- simPop(parms=params(k2form="linear"))
plot.func(name="Linear")

test <- simPop(parms=params(k2form="constant"))
plot.func(name="Constant")
dev.off()
