library(ggplot2)
library(gridExtra)
source("4_eia.R")   # script containing evolutionary invasion analysis functions

#****************colours for plots***********8
cols <- c("#e41a1c"
          ,"#377eb8"
          ,"#4daf4a"
          ,"#984ea3"
          ,"#ff7f00"
          ,"#ffff33"
          ,"#a65628" 
          ,"#f781bf"
          ,"#999999"
          ,"#cab2d6"
)
#***************plot commands******************
theme <-   theme(panel.border = element_blank()
                 ,axis.line = element_line(color = 'black')
                 ,text=element_text(size=8)
                 ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
                 ,axis.text=element_text(size=8)
                 ,legend.key.size = unit(0.8,"line")
                 ,legend.background = element_blank()
                 ,legend.text=element_text(size=7)
                 ,legend.position =c(0.2,0.9)
                 ,legend.title = element_blank()
)
#*****************************************************

#*********growth of persisent virus in present of acute as a function of mutant budding rate**8
lam <- seq(1,6000,0.1)  # range of budding rates for mutant virus

lam.mut <- mutfit.v(lambda_m=lam    # assuming no delays
                    ,beta_m=10^-9
                    ,beta=10^-9
                    ,alpha_m=0
                    ,gamma_m=0
                    ,lambda=0
                    ,gamma=6000
                    ,alpha=1/24)
lam.mut <- cbind.data.frame(lam,lam.mut)

lam.delay <- mutfit.delay.v(lambda_m=lam   # assuming budding delay and fixed time to apoptosis
                            ,lambda=0
                            ,beta_m=10^-9
                            ,beta=10^-9
                            ,tau_a_m=1000
                            ,tau_a=24
                            ,gamma_m=0
                            ,gamma=6000 
                            ,tau_b_m=24
                            ,tau_b=24
                            
                            )

lam.mut <- rbind.data.frame(lam.mut,cbind.data.frame(lam=lam,lam.mut=lam.delay))
lam.mut$nam <- c(rep("No delay",length(lam)),rep("Delay",length(lam)))

#pdf(file="fig_9.pdf",width=5,height=4)
ggplot(lam.mut, aes(x=lam,y=lam.mut)) +
   geom_line(aes(x=lam,y=lam.mut,color=nam)) +
  scale_color_manual(values=as.character(cols)) +
  labs(y="Mutant virus fitness", x=expression(paste("Mutant virus budding rate (",italic(lambda)["m"],")")),")") +
  theme_set(theme_bw())  +
  theme
#dev.off()
#*******************************************************************************

#*******************************************************************************

lambda_m <-  seq(1,400,0.05)       # range of values for mutant virus 
fit.bound1 <- sapply(lambda_m,function(y){
 return( uniroot.all(fun.nodelay,c(1,10000),lambda_m=y,gamma="NA",mu_mc=1/24,mu_c=1/24))
})

fit.bound2 <- sapply(lambda_m,function(y){
  return( uniroot.all(fun.nodelay,c(1,10000),lambda_m=y,gamma="NA",mu_mc=1/120,mu_c=1/120))
})

fit.bound3 <- sapply(lambda_m,function(y){
  return(uniroot.all(fun.delay,c(1,10000),lambda_m=y,gamma="NA"))
})

lam.roots <- cbind.data.frame(lambda_m,"fit.bound1"=as.numeric(fit.bound1),"fit.bound2"=as.numeric(fit.bound2))

pdf(file="fig_10.pdf",width=5,height=4)
ggplot(lam.roots, aes(x=fit.bound1,y=lambda_m)) +
  geom_ribbon(aes(x=fit.bound1/24,ymin=0,ymax=lambda_m),fill="grey") +
  #geom_line(aes(x=fit.bound2/24,y=lambda_m),color="red") +
  labs(y=expression(paste("Mutant virus budding rate (",italic(lambda)["m"],")"))
       , x=expression(paste("Resident virus yield at apoptosis * apoptosis rate  (",italic(gamma),italic(alpha),")"))) +
  theme_set(theme_bw())  +
  theme
dev.off()
#**********************************************************************************
