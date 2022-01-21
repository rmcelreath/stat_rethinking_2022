# week 2 solutions

# 1

library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]
Hbar <- mean(d2$height)
dat <- list(W=d2$weight,H=d2$height,Hbar=Hbar)

m1 <- quap(
    alist(
        W ~ dnorm( mu , sigma ) ,
        mu <- a + b*( H - Hbar ) ,
        a ~ dnorm( 60 , 10 ) ,
        b ~ dlnorm( 0 , 1 ) ,
        sigma ~ dunif( 0 , 10 )
    ) , data=dat )

dat2 <- list( H=c(140,160,175,100) , Hbar=Hbar )
h_sim <- sim( m1 , data=dat2 )
Ew <- apply(h_sim,2,mean)
h_ci <- apply(h_sim,2,PI,prob=0.89)

datr <- cbind( H=c(140,160,175) , Ew , L89=h_ci[1,] , U89=h_ci[2,] )
round(datr,1)

# 2

library(rethinking)
data(Howell1)
d <- Howell1
d <- d[ d$age < 13 , ]

# sim from priors
n <- 10
a <- rnorm(n,5,1)
b <- rlnorm(n,0,1)
# blank(bty="n")
plot( NULL , xlim=range(d$age) , ylim=range(d$weight) )
for ( i in 1:n ) abline( a[i] , b[i] , lwd=3 , col=2 )

m2 <- quap(
    alist(
        W ~ dnorm( mu , sigma ),
        mu <- a + b*A,
        a ~ dnorm(5,1),
        b ~ dlnorm(0,1),
        sigma ~ dexp(1)
    ), data=list(W=d$weight,A=d$age) )

precis(m2)

# blank(bty="n")
plot( d$age , d$weight , lwd=3, col=2 )
post <- extract.samples(m2)
for ( i in 1:10 ) abline( post$a[i] , post$b[i] , lwd=3 , col=1 )

# 3

library(rethinking)
data(Howell1)
d <- Howell1
d <- d[ d$age < 13 , ]

dat <- list(W=d$weight,A=d$age,S=d$male+1)

m3 <- quap(
    alist(
        W ~ dnorm( mu , sigma ),
        mu <- a[S] + b[S]*A,
        a[S] ~ dnorm(5,1),
        b[S] ~ dlnorm(0,1),
        sigma ~ dexp(1)
    ), data=dat )

# blank(bty="n")
plot( d$age , d$weight , lwd=3, col=ifelse(d$male==1,4,2) , xlab="age (years)" , ylab="weight (kg)" )
Aseq <- 0:12

# girls
muF <- link(m3,data=list(A=Aseq,S=rep(1,13)))
shade( apply(muF,2,PI,0.99) , Aseq , col=col.alpha(2,0.5) )
lines( Aseq , apply(muF,2,mean) , lwd=3 , col=2 )

# boys
muM <- link(m3,data=list(A=Aseq,S=rep(2,13)))
shade( apply(muM,2,PI,0.99) , Aseq , col=col.alpha(4,0.5) )
lines( Aseq , apply(muM,2,mean) , lwd=3 , col=4 )


# contrast at each age
Aseq <- 0:12
mu1 <- sim(m3,data=list(A=Aseq,S=rep(1,13)))
mu2 <- sim(m3,data=list(A=Aseq,S=rep(2,13)))
mu_contrast <- mu1
for ( i in 1:13 ) mu_contrast[,i] <- mu2[,i] - mu1[,i]
plot( NULL , xlim=c(0,13) , ylim=c(-15,15) , xlab="age" , ylab="weight difference (boys-girls)" )

for ( p in c(0.5,0.67,0.89,0.99) )
shade( apply(mu_contrast,2,PI,prob=p) , Aseq )

abline(h=0,lty=2,lwd=2)

for ( i in 1:13 ) points( mu_contrast[1:1000,i] , col=ifelse(mu_contrast[1:1000,i]>0,4,2) , lwd=3 )

# 4 

data(Oxboys)
d <- Oxboys

d$delta <- NA
for ( i in 1:nrow(d) ) {
    if ( d$Occasion[i] > 1 ) d$delta[i] <- d$height[i] - d$height[i-1]
}
d <- d[ !is.na(d$delta) , ]

# simulation from priors
n <- 1e3
alpha <- rnorm(n,0,0.1)
sigma <- rexp(n,3)
delta_sim <- rlnorm(n,alpha,sigma)
dens(delta_sim)

# the model
m4 <- quap(
    alist(
        delta ~ dlnorm( alpha , sigma ),
        alpha ~ dnorm( 0 , 0.1 ),
        sigma ~ dexp( 3 )
    ), data=list(delta=d$delta) )

# compute posterior sum of 8 increments
post <- extract.samples(m4)

dsim <- rlnorm(1e3,post$alpha,post$sigma)
dens(dsim)

inc_sum <- sapply( 1:1000 , function(s) sum(rlnorm(8,post$alpha[s],post$sigma[s])) )
dens(inc_sum)

