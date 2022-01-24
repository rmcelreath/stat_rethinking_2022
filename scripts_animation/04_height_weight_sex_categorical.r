# new height, weight, sex categorical variable example

dens(d2$height[d2$male==1],lwd=3,col=4,xlab="height (cm)")
dens(d2$height[d2$male==0],lwd=3,col=2,add=TRUE)

dens(d2$weight[d2$male==1],lwd=3,col=4,xlab="weight (kg)")
dens(d2$weight[d2$male==0],lwd=3,col=2,add=TRUE)

# W ~ S

data(Howell1)
d <- Howell1
d <- d[ d$age>=18 , ]
dat <- list(
    W = d$weight,
    S = d$male + 1 ) # S=1 female, S=2 male

m_SW <- quap(
    alist(
        W ~ dnorm(mu,sigma),
        mu <- a[S],
        a[S] ~ dnorm(60,10),
        sigma ~ dunif(0,10)
    ), data=dat )

# posterior means
post <- extract.samples(m_SW)
dens( post$a[,1] , xlim=c(39,50) , lwd=3 , col=2 , xlab="posterior mean weight (kg)" )
dens( post$a[,2] , lwd=3 , col=4 , add=TRUE )

# posterior W distributions
W1 <- rnorm( 1000 , post$a[,1] , post$sigma )
W2 <- rnorm( 1000 , post$a[,2] , post$sigma )
dens( W1 , xlim=c(20,70) , ylim=c(0,0.085) , lwd=3 , col=2 , xlab="posterior predicted weight (kg)" )
dens( W2 , lwd=3 , col=4 , add=TRUE )

# causal contrast (in means)
mu_contrast <- post$a[,2] - post$a[,1]
dens( mu_contrast , xlim=c(3,10) , lwd=3 , col=1 , xlab="posterior mean weight contrast (kg)" )

# W contrast
W_contrast <- W2 - W1
dens( W_contrast , xlim=c(-25,35) , lwd=3 , col=1 , xlab="posterior weight contrast (kg)" )

Wdens <- density(W_contrast,adj=0.5)
polygon(c(Wdens$x[Wdens$x>0], max(Wdens$x), 0), c(Wdens$y[Wdens$x>0], 0, 0), col = 4, border = NA )
polygon(c(Wdens$x[Wdens$x<0], 0, min(Wdens$x)), c(Wdens$y[Wdens$x<0], 0, 0), col = 2, border = NA )

# proportion above zero
sum( W_contrast > 0 ) / 1000
# proportion below zero
sum( W_contrast < 0 ) / 1000

# W ~ S + H

data(Howell1)
d <- Howell1
d <- d[ d$age>=18 , ]
dat <- list(
    W = d$weight,
    H = d$height,
    Hbar = mean(d$height),
    S = d$male + 1 ) # S=1 female, S=2 male

m_SHW <- quap(
    alist(
        W ~ dnorm(mu,sigma),
        mu <- a[S] + b[S]*(H-Hbar),
        a[S] ~ dnorm(60,10),
        b[S] ~ dlnorm(0,1),
        sigma ~ dunif(0,10)
    ), data=dat )
    
# full system as SCM

data(Howell1)
d <- Howell1
d <- d[ d$age>=18 , ]
dat <- list(
    W = d$weight,
    H = d$height,
    Hbar = mean(d$height),
    S = d$male + 1 ) # S=1 female, S=2 male

m_SHW_full <- quap(
    alist(

        # weight
        W ~ dnorm(mu,sigma),
        mu <- a[S] + b[S]*(H-Hbar),
        a[S] ~ dnorm(60,10),
        b[S] ~ dlnorm(0,1),
        sigma ~ dunif(0,10),

        # height
        H ~ dnorm(nu,tau),
        nu <- h[S],
        h[S] ~ dnorm(160,10),
        tau ~ dunif(0,10)

    ), data=dat )

# compute total causal effect of S on W
post <- extract.samples(m_SHW_full)
Hbar <- dat$Hbar
n <- 1e4

with( post , {
# simulate W for S=1
H_S1 <- rnorm(n, h[,1] , tau )
W_S1 <- rnorm(n, a[,1] + b[,1]*(H_S1-Hbar) , sigma)
# simulate W for S=2
H_S2 <- rnorm(n, h[,2] , tau)
W_S2 <- rnorm(n, a[,2] + b[,2]*(H_S2-Hbar) , sigma)
# compute contrast
W_do_S <<- W_S2 - W_S1
})

# automated way
HWsim <- sim(m_SHW_full,
             data=list(S=c(1,2)),
             vars=c("H","W"))
W_do_S_auto <- HWsim$W[,2] - HWsim$W[,1]
