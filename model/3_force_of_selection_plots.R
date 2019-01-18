library(ggplot2)
library(gridExtra)
source("3_force_of_selection_functions.R")

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

#*******************Force of selection*********************************
#*********************RANGES FOR PARAMETER VALUES***************************
lam <- seq(0,50,0.1) #seq(0,10^2.5,2)
bet <- seq(10^-8,10^-5,10^-8) #seq(10^-10,10^-5,10^-8)
alph <- seq(0,1/2,1/200)
gam <- seq(0,1200,10)
gamforalph <-  1/alph*100 
alphforgam <- 1/gam*100


#***************************************************************************
#******************Plot****************************************

vals.lam <- do.dl.v(lambda=lam)
vals.gam <- do.dg.v(gamma=gam)
vals.alph <- do.da.v(alpha=alph)
vals.bet <- do.db.v(beta=bet)
vals.alphwithgam <- do.da.v(alpha=alph,gamma=gamforalph)
vals.gamwithalph <- do.dg.v(gamma=gam,alpha=alphforgam)


pdf(file="fig2.pdf",width=4,height=4)
l <- ggplot(cbind.data.frame(lam,vals.lam), aes(x=lam,y=vals.lam)) +
  geom_line() +
  labs( x=expression(paste(lambda)),y="Force of selection", title="a") +
  theme_set(theme_bw())  +
  theme
g <- ggplot(cbind.data.frame(gam,vals.gam), aes(x=gam,y=vals.gam)) +
  geom_line() +
  labs( x=expression(paste(gamma)),y=" ", title="b") +
  theme_set(theme_bw())  +
  theme
a <- ggplot(cbind.data.frame(alph,vals.alph), aes(x=alph,y=vals.alph)) +
  geom_line() +
  labs( x=expression(paste(alpha)),y="Force of selection", title="c") +
  theme_set(theme_bw())  +
  theme

ag <- ggplot(cbind.data.frame(alph,vals.alphwithgam), aes(x=alph,y=vals.alphwithgam)) +
  geom_line() +
  labs( x=expression(paste(alpha)),y=" ", title="d") +
  theme_set(theme_bw())  +
  theme

ga <- ggplot(cbind.data.frame(gam,vals.gamwithalph), aes(x=gam,y=vals.gamwithalph)) +
  geom_line() +
  labs( x=expression(paste(gamma)),y="Force of selection", title="e") +
  theme_set(theme_bw())  +
  theme

b <- ggplot(cbind.data.frame(bet,vals.bet), aes(x=bet,y=vals.bet)) +
  geom_line() +
  labs( x=expression(paste(beta)),y=" ", title="f") +
  theme_set(theme_bw())  +
  theme

grid.arrange(l, g, a, ag, ga, b, nrow = 3)


dev.off()
#***************************************************************************


