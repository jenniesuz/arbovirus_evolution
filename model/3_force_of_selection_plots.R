library(ggplot2)
library(gridExtra)
source("3_force_of_selection_functions.R")


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

theme <-   theme(panel.border = element_blank()
                 ,axis.line = element_line(color = 'black')
                 ,text=element_text(size=6)
                 ,plot.margin=unit(c(0.2,0.3,0.1,0.1), "cm")
                 ,axis.text=element_text(size=5)
                 ,legend.key.size = unit(0.4,"line")
                 ,legend.background = element_blank()
                 ,legend.text=element_text(size=5)
                 ,legend.position =c(0.8,0.8)
                 ,legend.title = element_blank()
)

#*******************Force of selection*********************************
#*********************RANGES FOR PARAMETER VALUES***************************
lam <- seq(1,50,0.1) #seq(0,10^2.5,2)
bet <- seq(10^-8,10^-5,10^-8) #seq(10^-10,10^-5,10^-8)
alph <- seq(0,1/2,1/200)
gam <- seq(0,1200,10)
gamforalph <-  1/alph*100 
alphforgam <- 1/gam*100
S<-10^6


#***************************************************************************
#******************Plot****************************************

vals.lam <- do.dl.v(lambda=lam)
vals.gam <- do.dg.v(gamma=gam)
vals.alph <- do.da.v(alpha=alph)
vals.bet <- do.db.v(beta=bet)
vals.alphwithgam <- do.da.v(alpha=alph,gamma=gamforalph)
vals.gamwithalph <- do.dg.v(gamma=gam,alpha=alphforgam)


pdf(file="fig_7.pdf",width=4,height=4)
l <- ggplot(cbind.data.frame(lam,vals.lam), aes(x=lam,y=vals.lam)) +
  geom_line() +
  labs( x=expression(paste(lambda)),y=expression(paste(italic("d"),italic(omega),"/",italic("d"),italic(lambda))), title="a)") +
  theme_set(theme_bw())  +
  theme
g <- ggplot(cbind.data.frame(gam,vals.gam), aes(x=gam,y=vals.gam)) +
  geom_line() +
  labs( x=expression(paste(gamma)),y=expression(paste(italic("d"),italic(omega),"/",italic("d"),italic(gamma))), title="b)") +
  theme_set(theme_bw())  +
  theme

a <- ggplot(cbind.data.frame(alph,vals.alph), aes(x=alph,y=vals.alph)) +
  geom_line() +
  labs( x=expression(paste(alpha)),y=expression(paste(italic("d"),italic(omega),"/",italic("d"),italic(alpha))), title="c)") +
  theme_set(theme_bw())  +
  theme

ag <- ggplot(cbind.data.frame(alph,vals.alphwithgam), aes(x=alph,y=vals.alphwithgam)) +
  geom_line() +
  labs( x=expression(paste(alpha)),y=expression(paste(italic("d"),italic(omega),"/",italic("d"),italic(alpha))), title="d)") +
  theme_set(theme_bw())  +
  theme

ga <- ggplot(cbind.data.frame(gam,vals.gamwithalph), aes(x=gam,y=vals.gamwithalph)) +
  geom_line() +
  labs( x=expression(paste(gamma)),y=expression(paste(italic("d"),italic(omega),"/",italic("d"),italic(gamma))), title="e)") +
  theme_set(theme_bw())  +
  theme

b <- ggplot(cbind.data.frame(bet,vals.bet), aes(x=bet,y=vals.bet)) +
  geom_line() +
  labs( x=expression(paste(beta)),y=expression(paste(italic("d"),italic(omega),"/",italic("d"),italic(beta))), title="f)") +
  theme_set(theme_bw())  +
  theme

grid.arrange(l, g, a, b, nrow = 2,ncol=2)


dev.off()
#***************************************************************************

#*********************************************************************************
#*********************************************************************************
vfs.plot <- function(dat,title,xlab=expression(italic(lambda)),ylab=expression(paste(italic("d"),italic(omega),"/",italic("d"),italic(gamma)))){
  pt <- ggplot(dat, aes(x=p,y=val)) +
    geom_line(aes(x=p,y=val,color=as.factor(Param))) +
    scale_color_manual(values=as.character(cols)) +
    labs( y=ylab,x=xlab) +
    ggtitle(title) +
    theme_set(theme_bw())  +
    theme
  return(pt)
}
#*******Effect of other parameters on fitness as a function of lambda**********************************
gam <- lam
#********effect of infection probability****************
lambet <- expand.grid(p=lam,Param=c(1e-06, 2e-06, 3e-06, 4e-06, 5e-06, 6e-06,7e-06, 8e-06, 9e-06,1e-05)*S)
lambet.vals <- do.dl.v(lambda=lambet[1],beta=lambet[2]/S)
lambet$val <- lambet.vals
#*********************************************8
gambet <- expand.grid(p=gam,Param=c(1e-06, 2e-06, 3e-06, 4e-06, 5e-06, 6e-06,7e-06, 8e-06, 9e-06,1e-05)*S)
gambet.vals <- do.dg.v(gamma=gambet[1],beta=gambet[2]/S)
gambet$val <- gambet.vals

#******effect of apoptosis rate****
gamal <- expand.grid(p=gam,Param=1/c(1/6,1/12,1/18,1/24,1/30,1/36,1/42,1/48))
gamal.vals <- do.dg.v(gamma=gamal[1],alpha=1/gamal[2])
gamal$val <- gamal.vals


#***plots******
p1 <- vfs.plot(dat=lambet,title="a) Effect of probability of infection \n (shown as beta*S)",ylab=expression(paste(italic("d"),italic(omega),"/",italic("d"),italic(lambda))))
p2 <- vfs.plot(dat=gambet,title="b) Effect of probability of infection \n (shown as beta*S)",xlab=expression(italic(gamma)))
p3 <- vfs.plot(dat=gamal,title="c) Effect of virus apoptosis (shown as 1/alpha)",xlab=expression(italic(gamma)))

pdf(file="fig_8.pdf",width=4,height=4)
grid.arrange(p1,p2,p3,nrow=2, ncol=2)
dev.off()

