library(ggplot2)
library(gridExtra)
source("4_eia.R")
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

theme <-   theme(panel.border = element_blank()
                 ,axis.line = element_line(color = 'black')
                 ,text=element_text(size=6)
                 ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
                 ,axis.text=element_text(size=5)
                 ,legend.key.size = unit(0.8,"line")
                 ,legend.background = element_blank()
                 ,legend.text=element_text(size=5)
                 ,legend.position =c(0.4,0.9)
                 ,legend.title = element_blank()
)


lam <- seq(1,500,0.1)

lam.mut <- mutfit.v(lambda_m=lam,beta_m=10^-9,beta=10^-9,alpha_m=0,gamma_m=0,lambda=0,gamma=6000,alpha=1/24)
lam.mut <- cbind.data.frame(lam,lam.mut)

lam.delay <- mutfit.delay.v(lambda_m=lam
                            ,lambda=0
                            ,beta_m=10^-9
                            ,beta=10^-9
                            ,f.alpha_m=500
                            ,f.alpha=24
                            ,gamma_m=0
                            ,gamma=6000 
                            ,b_delay_m=24
                            ,b_delay=24
                            
                            )

lam.mut <- rbind.data.frame(lam.mut,cbind.data.frame(lam=lam,lam.mut=lam.delay))
lam.mut$nam <- c(rep("no delay",length(lam)),rep("delay",length(lam)))

ggplot(lam.mut, aes(x=lam,y=lam.mut)) +
   geom_line(aes(x=lam,y=lam.mut,color=nam)) +
  scale_color_manual(values=as.character(cols)) +
  labs(y="Mutant virus fitness", x=expression(paste("Mutant virus budding rate (",italic(lambda)["m"],")")),")") +
  theme_set(theme_bw())  +
  theme


#*******************************************************************************



#******************************************************************************************

lambda_m <-  seq(1,100,0.1)       # range of values for mutant virus 
fit.bound1 <- sapply(lambda_m,function(y){
  uniroot.all(fun1,c(0,10^3*24),lambda_m=y,gamma="NA")
})

fit.bound2 <- sapply(lambda_m,function(y){
  return(uniroot.all(fun2,c(1,10^3*24),lambda_m=y,gamma="NA"))
})

lam.roots <- cbind.data.frame(lambda_m,"fit.bound1"=fit.bound1/24,"fit.bound2"=as.numeric(fit.bound2)/24)

pdf(file="fig_3.pdf",width=5,height=4)
ggplot(lam.roots, aes(x=fit.bound1,y=lambda_m)) +
  geom_line(aes(x=fit.bound1,y=lambda_m,color="Model 1"),linetype=1) +
  geom_line(aes(x=fit.bound2,y=lambda_m,color="Model 2"),linetype=2) +
  scale_color_manual(name="",values=c("Model 1"="black","Model 2"="black")) +
  guides(color=guide_legend(override.aes=list(linetype=c(1,2)))) +
  labs(y=expression(paste("Mutant virus budding rate (",italic(lambda)["m"],")"))
       , x=expression(paste("Resident virus yeild at apoptosis/time to apoptosis (",italic(gamma),"/",italic(tau),")"))) +
  theme_set(theme_bw())  +
  theme
dev.off()

