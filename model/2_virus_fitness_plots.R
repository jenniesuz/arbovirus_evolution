library(ggplot2)
library(gridExtra)
source("2_virus_fitness_functions.R")

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
#*********************************************
theme <-   theme(panel.border = element_blank()
                 ,axis.line = element_line(color = 'black')
                 ,text=element_text(size=5)
                 ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
                 ,axis.text=element_text(size=5)
                 ,legend.key.size = unit(0.3,"line")
                 ,legend.background = element_blank()
                 ,legend.text=element_text(size=3)
                 ,legend.position =c(0.2,0.7)
                ,legend.title = element_blank()
)


#****************Effect of parameter values on fitness function witout delays**********************
lam <- seq(1,500,0.1) 
bet <- seq(0,10^-5,10^-8)
alph <- seq(1/72,1/2,1/200)
gam <- seq(1,500*24,100)
muc <- seq(1/500,1/2,1/200)
S = 10^6


vals.lam <- fit.nodelay.v(lambda=lam,gamma=0,alpha=0)
vals.alph <- fit.nodelay.v(alpha=alph, lambda=0)
vals.gam <- fit.nodelay.v(gamma=gam,lambda=0,alpha=1/24)
vals.bet <- fit.nodelay.v(beta=bet,gamma=0,alpha=0)
vals.mucp <- fit.nodelay.v(mu_c=muc,gamma=0,alpha=0)
vals.muca <- fit.nodelay.v(mu_c=muc,lambda=0)

#****************fitness function plots for model without delays****************
pdf(file="fig_2.pdf",width=4,height=4)
l <- ggplot(cbind.data.frame(lam,vals.lam), aes(x=lam,y=vals.lam)) +
  geom_line() +
  ylim(0,20) +
  labs( x=expression(paste(lambda)),y=expression(paste("Virus fitness (",italic(omega),")")), title="a)") +
  theme_set(theme_bw())  +
  theme
g <- ggplot(cbind.data.frame(gam,vals.gam), aes(x=gam,y=vals.gam)) +
  geom_line() +
  ylim(0,20) +
  labs( x=expression(paste(gamma)),y=" ", title="b)") +
  theme_set(theme_bw())  +
  theme
a <- ggplot(cbind.data.frame(alph,vals.alph), aes(x=alph,y=vals.alph)) +
  geom_line() +
  ylim(0,20) +
  labs( x=expression(paste(alpha)),y=expression(paste("Virus fitness (",italic(omega),")")), title="c)") +
  theme_set(theme_bw())  +
  theme
b <- ggplot(cbind.data.frame(bet,vals.bet), aes(x=bet,y=vals.bet)) +
  geom_line() +
  ylim(0,20) +
  labs( x=expression(paste(beta)),y=" ", title="d)") +
  theme_set(theme_bw())  +
  theme
mp <- ggplot(cbind.data.frame(muc,vals.mucp), aes(x=muc,y=vals.mucp)) +
  geom_line() +
  ylim(0,20) +
  labs( x=expression(paste(mu),"I"),y=" ", title="e)") +
  theme_set(theme_bw())  +
  theme
ma <- ggplot(cbind.data.frame(muc,vals.muca), aes(x=muc,y=vals.muca)) +
  geom_line() +
  ylim(0,20) +
  labs( x=expression(paste(mu),"I"),y=" ", title="f)") +
  theme_set(theme_bw())  +
  theme

grid.arrange(l, g, a, b,mp,ma, nrow = 3)

dev.off()

#*********************************************************************************
#*********************************************************************************
vfs.plot <- function(dat,title,xlab=expression(italic(lambda))){
  pt <- ggplot(dat, aes(x=p,y=val)) +
    geom_line(aes(x=p,y=val,color=as.factor(Param))) +
    scale_color_manual(values=as.character(cols)) +
    labs( y=expression(paste("Virus fitness (",italic(omega),")")),x=xlab) +
    ggtitle(title) +
    theme_set(theme_bw())  +
    theme
  return(pt)
}
#*******Effect of other parameters on fitness as a function of lambda**********************************
#********effect of infection probability****************
lambet <- expand.grid(p=lam,Param=c(1e-06, 2e-06, 3e-06, 4e-06, 5e-06, 6e-06,7e-06, 8e-06, 9e-06,1e-05)*S)
lambet.vals <- fit.nodelay.v(lambda=lambet[1],beta=lambet[2]/S,gamma=0,alpha=0)
lambet$val <- lambet.vals
#********effect of cell death rate******
lammuc <- expand.grid(p=lam,Param=1/c(1/24,1/48,1/72,1/96,1/120,1/144))
lammuc.vals <- fit.nodelay.v(lambda=lammuc[1],mu_c=1/lammuc[2],gamma=0,alpha=0)
lammuc$val <- lammuc.vals
#******effect of virus clearance rate****
lammuv <- expand.grid(p=lam,Param=seq(0,0.8,0.1))
lammuv.vals <- fit.nodelay.v(lambda=lammuv[1],mu_v=lammuv[2],gamma=0,alpha=0)
lammuv$val <- lammuv.vals

#***plots******
p1 <- vfs.plot(dat=lambet,title="a) Effect of probability of infection \n (shown as beta*S)")
p2 <- vfs.plot(dat=lammuc,title="b) Effect of cell life span")
p3 <- vfs.plot(dat=lammuv,title="c) Effect of virus death rate")

pdf(file="fig_3.pdf",width=4,height=4)
grid.arrange(p1,p2,p3,nrow=2, ncol=2)
dev.off()

#*******Effect of other parameters on fitness as a function of gamma**********************************
#********effect of infection probability****************
gambet <- expand.grid(p=gam,Param=c(1e-06, 2e-06, 3e-06, 4e-06, 5e-06, 6e-06,7e-06, 8e-06, 9e-06,1e-05)*S)
gambet.vals <- fit.nodelay.v(gamma=gambet[1],beta=gambet[2]/S,lambda=0)
gambet$val <- gambet.vals
#********effect of cell death rate******
gammuc <- expand.grid(p=gam,Param=1/c(1/24,1/48,1/72,1/96,1/120,1/144))
gammuc.vals <- fit.nodelay.v(gamma=gammuc[1],mu_c=1/gammuc[2],lambda=0)
gammuc$val <- gammuc.vals
#******effect of virus clearance rate****
gammuv <- expand.grid(p=gam,Param=seq(0,0.8,0.1))
gammuv.vals <- fit.nodelay.v(gamma=gammuv[1],mu_v=gammuv[2],lambda=0)
gammuv$val <- gammuv.vals
#******effect of apoptosis rate****
gamal <- expand.grid(p=gam,Param=1/c(1/6,1/12,1/18,1/24,1/30,1/36,1/42,1/48))
gamal.vals <- fit.nodelay.v(gamma=gamal[1],alpha=1/gamal[2],lambda=0)
gamal$val <- gamal.vals


#***plots******
p4 <- vfs.plot(dat=gambet,title="a) Effect of probability of infection \n (shown as beta*S)",xlab=expression(italic(gamma)))
p5 <- vfs.plot(dat=gammuc,title="b) Effect of cell life span (hrs)",xlab=expression(italic(gamma)))
p6 <- vfs.plot(dat=gammuv,title="c) Effect of virus death rate",xlab=expression(italic(gamma)))
p7 <- vfs.plot(dat=gamal,title="d) Effect of time to apoptosis",xlab=expression(italic(gamma)))


pdf(file="fig_4.pdf",width=4,height=4)
grid.arrange(p4,p5,p6,p7,nrow=2, ncol=2)
dev.off()

#*************************************************************************
#**********************************************************************




#*****************Effect of parameter values on fitness function with delays******************************
lam <- seq(10,200,0.1) 
bet <- seq(10^-6,10^-5,10^-8)
alph <- seq(1/72,1/2,1/300)
gam <- seq(100,50*24,1)
t.a <- 1/alph
t.b <- 1/alph
S = 10^6


vals.lam <- fit.delay.v(lambda=lam,gamma=0,tau.a=20000,tau.b=24)
vals.tau.a <- fit.delay.v(tau.a=t.a,gamma=100, lambda=0,tau.b=120,mu_c=1/120)
vals.tau.b <- fit.delay.v(tau.b=t.b, lambda=100,gamma=0,tau.a=2000)
vals.gam <- fit.delay.v(gamma=gam,lambda=0,tau.a=24,tau.b=120)
vals.bet <- fit.delay.v(beta=bet,lambda=0,gamma=50,tau.a=24,tau.b=120)

pdf(file="fig_5.pdf",width=4,height=4)
l <- ggplot(cbind.data.frame(lam,vals.lam), aes(x=lam,y=vals.lam)) +
  geom_line() +
  labs( x=expression(paste(lambda)),y=expression(paste("Virus fitness (",italic(omega),")")), title="a)") +
  theme_set(theme_bw())  +
  theme
g <- ggplot(cbind.data.frame(gam,vals.gam), aes(x=gam,y=vals.gam)) +
  geom_line() +
  labs( x=expression(paste(gamma)),y=" ", title="b)") +
  theme_set(theme_bw())  +
  theme
ta <- ggplot(cbind.data.frame(t.a,vals.tau.a), aes(x=t.a,y=vals.tau.a)) +
  geom_line() +
  labs( x=expression(paste(tau," delay in apoptosis")),y=expression(paste("Virus fitness (",italic(omega),")")), title="c)") +
  theme_set(theme_bw())  +
  theme
tb <- ggplot(cbind.data.frame(t.b,vals.tau.b), aes(x=t.b,y=vals.tau.b)) +
  geom_line() +
  labs( x=expression(paste(tau,"' delay in budding")),y=" ", title="d)") +
  theme_set(theme_bw())  +
  theme
b <- ggplot(cbind.data.frame(bet,vals.bet), aes(x=bet,y=vals.bet)) +
  geom_line() +
  labs( x=expression(paste(beta)),y=expression(paste("Virus fitness (",italic(omega),")")), title="e)") +
  theme_set(theme_bw())  +
  theme

grid.arrange(l, g, ta,tb, b, nrow = 3)

dev.off()
#*****************************************************************************************
#*******Effect of other parameters on fitness as a function of tau**********************************
#********effect of infection probability****************
taubet <- expand.grid(p=t.a,Param=c(1e-06, 2e-06, 3e-06, 4e-06, 5e-06, 6e-06,7e-06, 8e-06, 9e-06,1e-05)*S)
taubet.vals <- fit.delay.v(tau.a=taubet[1],beta=taubet[2]/S,gamma=50,tau.b=500,lambda=0)
taubet$val <- taubet.vals
#********effect of cell death rate******
taumuc <- expand.grid(p=t.a,Param=1/c(1/24,1/48,1/72,1/96,1/120,1/144))
taumuc.vals <- fit.delay.v(tau.a=taumuc[1],mu_c=1/taumuc[2],gamma=50,tau.b=500,lambda=0)
taumuc$val <- taumuc.vals
#******effect of virus clearance rate****
taumuv <- expand.grid(p=t.a,Param=seq(0,0.8,0.1))
taumuv.vals <- fit.delay.v(tau.a=taumuv[1],mu_v=taumuv[2],gamma=50,tau.b=500,lambda=0)
taumuv$val <- taumuv.vals
#******effect of yield***************
taugam <- expand.grid(p=t.a,Param=c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000))
taugam.vals <- fit.delay.v(tau.a=taugam[1],gamma=taugam[2],tau.b=500,lambda=0)
taugam$val <- taugam.vals

#***plots******
p8 <- vfs.plot(dat=taubet,title="a) Effect of probability of infection \n (shown as beta*S)",xlab=expression(tau))
p9 <- vfs.plot(dat=taumuc,title="b) Effect of cell life span (hrs)",xlab=expression(tau))
p10 <- vfs.plot(dat=taumuv,title="c) Effect of virus death rate",xlab=expression(tau))
p11 <- vfs.plot(dat=taugam,title="d) Effect of yield",xlab=expression(tau))


pdf(file="fig_6.pdf",width=4,height=4)
grid.arrange(p8,p9,p10,p11,nrow=2, ncol=2)
dev.off()





# #*****************Effect of parameter values on acute fitness function with delay******************************
# bet <- seq(10^-6,10^-5,10^-8)
# alph <- seq(1/72,1/2,1/200)
# gam <- seq(0,5000,10)
# t.a <- seq(2,120,10)
# S = 10^6
# 
# vals.tau.a <- fit.a.delay.v(tau.a=t.a,gamma=100,mu_c=0)
# vals.gam <- fit.a.delay.v(gamma=gam,tau.a=2,beta=10^-9)
# vals.bet <- fit.a.delay.v(beta=bet,gamma=50,tau.a=24)
# 
# #pdf(file="fig2.pdf",width=4,height=4)
# g <- ggplot(cbind.data.frame(gam,vals.gam), aes(x=gam,y=vals.gam)) +
#   geom_line() +
#   labs( x=expression(paste(gamma)),y=" ", title="b") +
#   theme_set(theme_bw())  +
#   theme
# ta <- ggplot(cbind.data.frame(t.a,vals.tau.a), aes(x=t.a,y=vals.tau.a)) +
#   geom_line() +
#   labs( x=expression(paste(tau," delay in apoptosis")),y="Fitness", title="c") +
#   theme_set(theme_bw())  +
#   theme
# b <- ggplot(cbind.data.frame(bet,vals.bet), aes(x=bet,y=vals.bet)) +
#   geom_line() +
#   labs( x=expression(paste(beta)),y="Fitness", title="e") +
#   theme_set(theme_bw())  +
#   theme
# 
# grid.arrange(g, ta,b, ncol = 3)
# 
# #dev.off()
# #*****************************************************************************************
# #*************************************************************************************
