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

#*********growth of acute  virus in present of persistent as a function of yield********
gam <- seq(1,2000,0.1)  # range of budding rates for mutant virus

gam.mut <- mutfit.v(lambda_m=1   # assuming no delays
                    ,beta_m=10^-9
                    ,beta=10^-9
                    ,gamma_m=gam
                    ,lambda=50
                    ,gamma=0
                    ,alpha_m=1/24
                    ,alpha=0)
gam.mut <- cbind.data.frame(gam,gam.mut)

gam.delay <- mutfit.delay.v(mu_mv=0.1            # mutant virus clearance rate
                            ,mu_v=0.1            # resident virus clearance rate  
                            ,mu_mc = 1/120       # mutant virus cell death rate
                            ,mu_c = 1/120        # resident cell death rate   
                            ,beta_m = 10^-9      # mutant probability of infecting susceptible cells
                            ,beta = 10^-9        # resident probability of infecting susceptible cells
                            
                            ,gamma = 1*24        # resident virus yeild at apoptosis
                            ,gamma_m = gam       # mutant virus yeild at apoptosis
                            ,lambda = 50         # resident virus budding rate
                            ,lambda_m = 0        # mutant virus budding rate
                            ,tau_a = 24*5          # resident virus fixed time to apoptosis
                            ,tau_a_m = 24        # mutant virus fixed time to apoptosis
                            ,tau_b= 12           # resident delay in budding
                            ,tau_b_m = 12        # mutant delay in budding
                          )

gam.mut <- rbind.data.frame(gam.mut,cbind.data.frame(gam=gam,gam.mut=gam.delay))
gam.mut$nam <- c(rep("No delay",length(gam)),rep("Delay",length(gam)))

#pdf(file="fig_9.pdf",width=5,height=4)
ggplot(gam.mut, aes(x=gam,y=gam.mut)) +
   geom_line(aes(x=gam,y=gam.mut,color=nam)) +
  scale_color_manual(values=as.character(cols)) +
  labs(y="Mutant virus fitness", x=expression(paste("Mutant virus yield (",italic(gamma)["m"],")")),")") +
  theme_set(theme_bw())  +
  theme
#dev.off()


mutfit.delay.v(mu_mv=0.1            # mutant virus clearance rate
               ,mu_v=0.1            # resident virus clearance rate  
               ,mu_mc = 1/120       # mutant virus cell death rate
               ,mu_c = 1/120        # resident cell death rate   
               ,beta_m = 10^-9      # mutant probability of infecting susceptible cells
               ,beta = 10^-9        # resident probability of infecting susceptible cells
               
               ,gamma = 1*24        # resident virus yeild at apoptosis
               ,gamma_m = 1*24       # mutant virus yeild at apoptosis
               ,lambda = 50         # resident virus budding rate
               ,lambda_m = 50        # mutant virus budding rate
               ,tau_a = 24          # resident virus fixed time to apoptosis
               ,tau_a_m = 24        # mutant virus fixed time to apoptosis
               ,tau_b= 12           # resident delay in budding
               ,tau_b_m = 12        # mutant delay in budding
)

#*******************************************************************************

#*******************************************************************************
# 
 gamma_m <-  seq(1,10000,1)       # range of values for mutant virus 
# fit.bound1 <- sapply(gamma_m,function(y){
#  return( uniroot.all(fun.nodelay
#                      ,c(1,10000)
#                      ,gamma_m=y
#                      ,gamma="NA"
#                      ,mu_mc=1/24
#                      ,mu_c=1/24))
# })
# 
# fit.bound2 <- sapply(lambda_m,function(y){
#   return( uniroot.all(fun.nodelay
#                       ,c(1,10000)
#                       ,lambda_m=y
#                       ,gamma="NA"
#                       ,mu_mc=1/120
#                       ,mu_c=1/120))
# })
# 
#gam.roots <- cbind.data.frame(lambda_m,"fit.bound1"=as.numeric(fit.bound1),"fit.bound2"=as.numeric(fit.bound2))

# #pdf(file="fig_10.pdf",width=5,height=4)
# ggplot(lam.roots, aes(x=fit.bound1,y=lambda_m)) +
#   geom_ribbon(aes(x=fit.bound1/24,ymin=0,ymax=lambda_m),fill="grey") +
#   #geom_line(aes(x=fit.bound2/24,y=lambda_m),color="red") +
#   labs(y=expression(paste("Mutant virus budding rate (",italic(lambda)["m"],")"))
#        , x=expression(paste("Resident virus yield at apoptosis * apoptosis rate  (",italic(gamma),italic(alpha),")"))) +
#   theme_set(theme_bw())  +
#   theme
# #dev.off()
# 
fit.bound3 <- lapply(gamma_m,function(y){
  return(uniroot.all(fun.delay
                     ,c(1,10000000)
                     ,mu_v=0.1  
                     ,mu_mv=0.1            # mutant virus clearance rate
                     ,mu_mc = 1/120       # mutant virus cell death rate
                     ,mu_c = 1/120        # resident infected cell death rate 
                     ,beta = 10^-6        # resident probability of infecting susceptible cells
                     ,beta_m = 10^-6      # mutant probability of infecting susceptible cells
                     ,gamma = 1           # resident virus yield at apoptosis
                     ,gamma_m=y
                     ,lambda="NA"         # resident budding
                     ,lambda_m = 1        # mutant budding
                     ,tau_a = 24          # resident virus apoptosis 
                     ,tau_a_m = 24        # mutant virus apoptosis 
                     ,tau_b_m = 12  
                     ,tau_b = 12
 ))
})


uniroot.all(fun.delay
            ,c(1,10000000)
            ,mu_v=0.1  
            ,mu_mv=0.1            # mutant virus clearance rate
            ,mu_mc = 1/120       # mutant virus cell death rate
            ,mu_c = 1/120        # resident infected cell death rate 
            ,beta = 10^-6        # resident probability of infecting susceptible cells
            ,beta_m = 10^-6      # mutant probability of infecting susceptible cells
            ,gamma = 200           # resident virus yield at apoptosis
            ,gamma_m=200
            ,lambda="NA"         # resident budding
            ,lambda_m = 20        # mutant budding
            ,tau_a = 24          # resident virus apoptosis 
            ,tau_a_m = 24        # mutant virus apoptosis 
            ,tau_b_m = 12  
            ,tau_b = 12
)
#**********************************************************************************
ggplot(lam.roots, aes(x=fit.bound3,y=lambda_m)) +
  geom_ribbon(aes(x=fit.bound3,ymin=0,ymax=lambda_m),fill="grey") +
  labs(y=expression(paste("Mutant virus budding rate (",italic(lambda)["m"],")"))
       , x=expression(paste("Resident virus yield at apoptosis * apoptosis rate  (",italic(gamma),italic(alpha),")"))) +
  theme_set(theme_bw())  +
  theme