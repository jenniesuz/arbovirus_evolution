library(ggplot2)
library(gridExtra)
source("2_virus_fitness_functions.R")

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


#****************Effect of parameter values on fitness function witout delays**********************
lam <- seq(1,500,0.1) 
bet <- seq(0,10^-5,10^-10)
alph <- seq(1/72,1/2,1/200)
gam <- seq(0,1200,10)
S = 10^6

vals.lam <- fit.nodelay.v(lambda=lam,gamma=0,alpha=0)
vals.alph <- fit.nodelay.v(alpha=alph, lambda=0)
vals.gam <- fit.nodelay.v(gamma=gam,lambda=0,alpha=1/24)
vals.bet <- fit.nodelay.v(beta=bet,gamma=0,alpha=0)

#****************fitness function plots for model without delays****************
pdf(file="fig1.pdf",width=4,height=4)
l <- ggplot(cbind.data.frame(lam,vals.lam), aes(x=lam,y=vals.lam)) +
  geom_line() +
 # ylim(0,7) +
  labs( x=expression(paste(lambda)),y="Fitness", title="a") +
  theme_set(theme_bw())  +
  theme
g <- ggplot(cbind.data.frame(gam,vals.gam), aes(x=gam,y=vals.gam)) +
  geom_line() +
  ylim(0,7) +
  labs( x=expression(paste(gamma)),y=" ", title="b") +
  theme_set(theme_bw())  +
  theme
a <- ggplot(cbind.data.frame(alph,vals.alph), aes(x=alph,y=vals.alph)) +
  geom_line() +
  labs( x=expression(paste(alpha)),y="Fitness", title="c") +
  theme_set(theme_bw())  +
  theme
b <- ggplot(cbind.data.frame(bet,vals.bet), aes(x=bet,y=vals.bet)) +
  geom_line() +
  labs( x=expression(paste(beta)),y=" ", title="d") +
  theme_set(theme_bw())  +
  theme

grid.arrange(l, g, a, b, nrow = 2)

dev.off()
#***********************************************************************

lam <- seq(1,50,0.1) 
bet <- seq(0,10^-5,10^-10)
alph <- seq(1/72,1/2,1/200)
gam <- seq(0,50,0.1)
t.a <- 1/alph
t.b <- 1/alph
S = 10^6


vals.lam <- fit.delay.v(lambda=lam,gamma=0,tau.a=72,tau.b=1)
vals.tau.a <- fit.delay.v(tau.a=t.a,gamma=100, lambda=0,tau.b=1)
vals.tau.b <- fit.delay.v(tau.b=t.b, lambda=10,gamma=0,tau.a=72)
vals.gam <- fit.delay.v(gamma=gam,lambda=0,tau.b=72,tau.a=1)
vals.bet <- fit.delay.v(beta=bet)

pdf(file="fig2.pdf",width=4,height=4)
l <- ggplot(cbind.data.frame(lam,vals.lam), aes(x=lam,y=vals.lam)) +
  geom_line() +
  labs( x=expression(paste(lambda)),y="Fitness", title="a") +
  theme_set(theme_bw())  +
  theme
g <- ggplot(cbind.data.frame(gam,vals.gam), aes(x=gam,y=vals.gam)) +
  geom_line() +
  labs( x=expression(paste(gamma)),y=" ", title="b") +
  theme_set(theme_bw())  +
  theme
ta <- ggplot(cbind.data.frame(t.a,vals.tau.a), aes(x=t.a,y=vals.tau.a)) +
  geom_line() +
  labs( x=expression(paste(tau," delay in apoptosis")),y="Fitness", title="c") +
  theme_set(theme_bw())  +
  theme
tb <- ggplot(cbind.data.frame(t.b,vals.tau.b), aes(x=t.b,y=vals.tau.b)) +
  geom_line() +
  labs( x=expression(paste(tau,"' delay in budding")),y=" ", title="d") +
  theme_set(theme_bw())  +
  theme
b <- ggplot(cbind.data.frame(bet,vals.bet), aes(x=bet,y=vals.bet)) +
  geom_line() +
  labs( x=expression(paste(beta)),y="Fitness", title="e") +
  theme_set(theme_bw())  +
  theme

grid.arrange(l, g, ta,tb, b, nrow = 3)

dev.off()
#*****************************************************************************************
