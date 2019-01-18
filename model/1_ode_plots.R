library(ggplot2)
library(gridExtra)
library(deSolve)

source("1a_ode_model_nodelays.R")
source("1d_ode_model_both.R")

theme <-   theme(panel.border = element_blank()
                 ,axis.line = element_line(color = 'black')
                 ,text=element_text(size=9)
                 ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
                 ,axis.text=element_text(size=8)
                 ,legend.key.size = unit(0.8,"line")
                 ,legend.background = element_blank()
                 ,legend.text=element_text(size=8)
                 ,legend.position =c(0.83,0.5)
                 ,legend.title = element_blank()
)


#************colours for plots***************************
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
#***********************************************


#****************run models**********************
no.delays <- simPop.nd(parms=params.nd(apoptosis=1/24
                                     ,yield=6000
                                     ,budding=200
                                     ,inf=10^-9))


delays <- simPop.b(parms=params.b(delay_a=24
                                     ,delay_p=24
                                     ,budding_p=100
                                     ,yield_a=6000
                                     ,apoptosis=24
                                     ,v_death=0.1
                                     ,c_death=1/120
                                     ,inf=10^-9))
#*************************************************


#*****plot**********
combined.mods <- cbind.data.frame(time=c(no.delays$time,no.delays$time,delays$time,delays$time)
                                  ,results=c(no.delays$V_p,no.delays$V_a,delays$V_p,delays$V_a)
                                  ,virustype=c(rep("Persistent",length(no.delays$time))
                                                   ,rep("Acute",length(no.delays$time))
                                                   ,rep("Persistent",length(no.delays$time))
                                                   ,rep("Acute",length(no.delays$time)))
                                  ,mod=c(rep("No delays",length(no.delays$time)*2),rep("Budding delay & fixed time to apoptosis",length(no.delays$time)*2))
)

#pdf(file="fig_ms_1.pdf",width=5,height=4)
ggplot(combined.mods, aes(x=time,y=results)) +
  geom_line(aes(x=time,y=log10(results),color=virustype)) +
  scale_color_manual(values=as.character(cols)) +
  facet_wrap( ~ mod,labeller = label_wrap_gen(width=30,multi_line=FALSE)) +
  labs( y=expression(paste("Number of virions (",log[10],")")),x="Time (hours)") +
  theme_set(theme_bw())  +
  theme
#dev.off()

