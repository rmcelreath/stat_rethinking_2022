# WEEK 4 PROBLEMS

################################
# 1

## R code 6.21
library(rethinking)
d <- sim_happiness( seed=1977 , N_years=1000 )
precis(d)

## R code 6.22
d2 <- d[ d$age>17 , ] # only adults
d2$A <- ( d2$age - 18 ) / ( 65 - 18 )

## R code 6.23
d2$mid <- d2$married + 1
m6.9 <- quap(
    alist(
        happiness ~ dnorm( mu , sigma ),
        mu <- a[mid] + bA*A,
        a[mid] ~ dnorm( 0 , 1 ),
        bA ~ dnorm( 0 , 2 ),
        sigma ~ dexp(1)
    ) , data=d2 )

## R code 6.24
m6.10 <- quap(
    alist(
        happiness ~ dnorm( mu , sigma ),
        mu <- a + bA*A,
        a ~ dnorm( 0 , 1 ),
        bA ~ dnorm( 0 , 2 ),
        sigma ~ dexp(1)
    ) , data=d2 )


compare( m6.9 , m6.10 , func=PSIS )
compare( m6.9 , m6.10 , func=WAIC )

m1 <- quap(
    alist(
        happiness ~ dnorm( mu , sigma ),
        mu <- a[mid],
        a[mid] ~ dnorm( 0 , 1 ),
        sigma ~ dexp(1)
    ) , data=d2 )

################################
# 2 

library(rethinking)
data(foxes)
d <- foxes
d$W <- standardize(d$weight)
d$A <- standardize(d$area)
d$F <- standardize(d$avgfood)
d$G <- standardize(d$groupsize)

tau <- 0.5

m1 <- quap(
    alist(
        W ~ dnorm( mu , sigma ),
        mu <- a + bF*F + bG*G + bA*A,
        a ~ dnorm(0,0.2),
        c(bF,bG,bA) ~ dnorm(0,tau),
        sigma ~ dexp(1)
    ), data=d )

m2 <- quap(
    alist(
        W ~ dnorm( mu , sigma ),
        mu <- a + bF*F + bG*G,
        a ~ dnorm(0,0.2),
        c(bF,bG) ~ dnorm(0,tau),
        sigma ~ dexp(1)
    ), data=d )

m3 <- quap(
    alist(
        W ~ dnorm( mu , sigma ),
        mu <- a + bG*G + bA*A,
        a ~ dnorm(0,0.2),
        c(bG,bA) ~ dnorm(0,tau),
        sigma ~ dexp(1)
    ), data=d )

m4 <- quap(
    alist(
        W ~ dnorm( mu , sigma ),
        mu <- a + bF*F,
        a ~ dnorm(0,0.2),
        bF ~ dnorm(0,tau),
        sigma ~ dexp(1)
    ), data=d )

m5 <- quap(
    alist(
        W ~ dnorm( mu , sigma ),
        mu <- a + bA*A,
        a ~ dnorm(0,0.2),
        bA ~ dnorm(0,tau),
        sigma ~ dexp(1)
), data=d )

compare( m1 , m2 , m3 , m4 , m5 , func=PSIS )

m6 <- quap(
    alist(
        W ~ dnorm( mu , sigma ),
        mu <- a + bA*A + bF*F,
        a ~ dnorm(0,0.2),
        c(bA,bF) ~ dnorm(0,0.5),
        sigma ~ dexp(1)
    ), data=d )

#######
# 3

data(cherry_blossoms)
d <- cherry_blossoms

d$D <- standardize(d$doy)
d$T <- standardize(d$temp)

dd <- d[ complete.cases(d$D,d$T) , ]

m3a <- quap(
    alist(
        D ~ dnorm(mu,sigma),
        mu <- a,
        a ~ dnorm(0,10),
        sigma ~ dexp(1)
    ) , data=list(D=dd$D,T=dd$T) )

m3b <- quap(
    alist(
        D ~ dnorm(mu,sigma),
        mu <- a + b*T,
        a ~ dnorm(0,10),
        b ~ dnorm(0,10),
        sigma ~ dexp(1)
    ) , data=list(D=dd$D,T=dd$T) )

m3c <- quap(
    alist(
        D ~ dnorm(mu,sigma),
        mu <- a + b1*T + b2*T^2,
        a ~ dnorm(0,10),
        c(b1,b2) ~ dnorm(0,10),
        sigma ~ dexp(1)
    ) , data=list(D=dd$D,T=dd$T) )

compare( m3a , m3b , m3c , func=PSIS )

# predict

Tval <- (9 - mean(d$temp,na.rm=TRUE))/sd(d$temp,na.rm=TRUE)
D_sim <- sim( m3b , data=list(T=Tval) )
# put back on natural scale
doy_sim <- D_sim*sd(d$doy,na.rm=TRUE) + mean(d$doy,na.rm=TRUE)
dens( doy_sim , lwd=4 , col=2 , xlab="day in year 1st bloom")
abline(v=89,lty=1)
dens( d$doy , add=TRUE , lwd=3 )

#######
# 4

data(Dinosaurs)
d <- Dinosaurs

dd <- d[ d$sp_id==1 , ]

dat <- list(
    A = dd$age,
    M = dd$mass/max(dd$mass)
)

m4a <- ulam(
    alist(
        M ~ normal(mu,sigma),
        mu <- a + b*A,
        a ~ normal(0,1),
        b ~ normal(0,1),
        sigma ~ exponential(1)
    ) , data=dat , chains=4 , log_lik=TRUE )

plot( dat$A , dat$M , xlab="age" , ylab="mass (normalized)" , lwd=3 , xlim=c(0,16) )
Aseq <- seq(from=0,to=16,len=50)
mu <- link(m4a,data=list(A=Aseq))
lines( Aseq , apply(mu,2,mean) , lwd=3 )
shade( apply(mu,2,PI) , Aseq )

m4b <- ulam(
    alist(
        M ~ normal(mu,sigma),
        mu <- k*(1-exp(-b*A)),
        b ~ exponential(1),
        k ~ normal(1,0.5),
        sigma ~ exponential(1)
    ) , data=dat , chains=4 , log_lik=TRUE )

plot( dat$A , dat$M , xlab="age" , ylab="mass (normalized)" , lwd=3 , xlim=c(0,16) )
mu <- link(m4b,data=list(A=Aseq))
lines( Aseq , apply(mu,2,mean) , lwd=3 )
shade( apply(mu,2,PI) , Aseq )

m4c <- ulam(
    alist(
        M ~ normal(mu,sigma),
        mu <- k*(1-exp(-b*A))^a,
        a ~ exponential(0.1),
        b ~ exponential(1),
        k ~ normal(1,0.5),
        sigma ~ exponential(1)
    ) , data=dat , chains=4 , log_lik=TRUE )

plot( dat$A , dat$M , xlab="age" , ylab="mass (normalized)" , lwd=3 , xlim=c(0,16) )
mu <- link(m4c,data=list(A=Aseq))
lines( Aseq , apply(mu,2,mean) , lwd=3 )
shade( apply(mu,2,PI) , Aseq )

compare( m4a , m4b , m4c , func=PSIS )

# try to fit to all species

# scale max size within each species
d$Ms <- sapply( 1:nrow(d) , function(i) d$mass[i]/max(d$mass[d$sp_id==d$sp_id[i]]) )
dat_all <- list( A=d$age , M=d$Ms , S=d$sp_id )

m4x <- ulam(
    alist(
        M ~ normal(mu,sigma),
        mu <- k[S]*(1-exp(-b[S]*A))^a[S],
        a[S] ~ exponential(0.1),
        b[S] ~ exponential(1),
        k[S] ~ normal(1,0.5),
        sigma ~ exponential(1)
    ) , data=dat_all , chains=4 , log_lik=TRUE )

plot( NULL , xlim=c(0,max(d$age)) , ylim=c(0,1.5) , xlab="dinosaur age" , ylab="mass (normalized)" )
post <- extract.samples(m4x)
Aseq <- seq(from=0,to=16,len=50)
for ( i in 1:max(d$sp_id) ) {
    j <- which(dat_all$S==i)
    points( dat_all$A[j] , dat_all$M[j] , col=i , lwd=3 )
    mu <- link( m4x , data=list(S=rep(i,50),A=Aseq) )
    lines( Aseq , apply(mu,2,mean) , lwd=3 , col=i )
}

# with intervals
plot( NULL , xlim=c(0,max(d$age)) , ylim=c(0,2) , xlab="dinosaur age" , ylab="mass (normalized)" )
post <- extract.samples(m4x)
Aseq <- seq(from=0,to=16,len=50)
for ( i in 1:max(d$sp_id) ) {
    j <- which(dat_all$S==i)
    points( dat_all$A[j] , dat_all$M[j] , col=i , lwd=3 )
    mu <- link( m4x , data=list(S=rep(i,50),A=Aseq) )
    lines( Aseq , apply(mu,2,mean) , lwd=3 , col=i )
    shade( apply(mu,2,PI) , Aseq , col=col.alpha(i,0.3) )
}

# line model
m4z <- ulam(
    alist(
        M ~ normal(mu,sigma),
        mu <- a[S] + b[S]*A,
        a[S] ~ normal(0,1),
        b[S] ~ normal(0,1),
        sigma ~ exponential(1)
    ) , data=dat_all , chains=4 , log_lik=TRUE )

compare( m4z , m4x , func=PSIS )

plot( NULL , xlim=c(0,max(d$age)) , ylim=c(0,1) , xlab="age" , ylab="mass (normalized)" )
post <- extract.samples(m4z)
Aseq <- seq(from=0,to=16,len=50)
for ( i in 1:max(d$sp_id) ) {
    j <- which(dat_all$S==i)
    points( dat_all$A[j] , dat_all$M[j] , col=i , lwd=3 )
    mu <- link( m4z , data=list(S=rep(i,50),A=Aseq) )
    lines( Aseq , apply(mu,2,mean) , lwd=2 , col=i )
}

