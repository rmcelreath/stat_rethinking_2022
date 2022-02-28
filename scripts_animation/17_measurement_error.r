# lecture 17
# measurement error

library(rethinking)
library(animation)
library(ellipse)

# simple parent-child income example
# recall bias on parental income

N <- 500
P <- rnorm(N)
C <- rnorm(N,0.75*P)
Pstar <- rnorm(N, 0.8*P + 0.2*C )

precis(lm(C~P))
precis(lm(C~Pstar))

mCP <- ulam(
    alist(
        C ~ normal(mu,sigma),
        mu <- a + b*P,
        a ~ normal(0,1),
        b ~ normal(0,1),
        sigma ~ exponential(1)
    ) , data=list(P=Pstar,C=C) , chains=4 )

precis(mCP)

post <- extract.samples(mCP)
dens( post$b , lwd=4 , col=2 , xlab="effect of P on C" , xlim=c(0,0.8) )
abline(v=0.75,lty=3)

###########
# divorce example

## R code 15.2
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

# points
plot( d$Divorce ~ d$Marriage , ylim=c(4,15) , xlim=c(12,35) ,
    xlab="Marriage rate" , ylab="Divorce rate" , lwd=3 , col=2 )

# standard errors
for ( i in 1:nrow(d) ) {
    ci <- d$Divorce[i] + c(-1,1)*d$Divorce.SE[i]
    x <- d$Marriage[i]
    lines( c(x,x) , ci , lwd=4 , col=col.alpha(2,0.4) )
    ci <- d$Marriage[i] + c(-1,1)*d$Marriage.SE[i]
    y <- d$Divorce[i]
    lines( ci , c(y,y) , lwd=4 , col=col.alpha(2,0.4) )
}

plot( d$Divorce ~ log(d$Population) , ylim=c(4,15) ,
    xlab="Population (log)" , ylab="Divorce rate" , lwd=3 , col=2 )
for ( i in 1:nrow(d) ) {
    ci <- d$Divorce[i] + c(-1,1)*d$Divorce.SE[i]
    x <- log(d$Population)[i]
    lines( c(x,x) , ci , lwd=6 , col=col.alpha(2,0.5) )
}

## R code 15.3
dlist <- list(
    D_obs = standardize( d$Divorce ),
    D_sd = d$Divorce.SE / sd( d$Divorce ),
    M = standardize( d$Marriage ),
    A = standardize( d$MedianAgeMarriage ),
    N = nrow(d)
)

m15.1 <- ulam(
    alist(
        # model for D* (observed)
        D_obs ~ dnorm( D_true , D_sd ),

        # model for D (unobserved)
        vector[N]:D_true ~ dnorm( mu , sigma ),
        mu <- a + bA*A + bM*M,
        a ~ dnorm(0,0.2),
        bA ~ dnorm(0,0.5),
        bM ~ dnorm(0,0.5),
        sigma ~ dexp(1)
    ) , data=dlist , chains=4 , cores=4 )

precis( m15.1 , depth=2 )

# show the regression and shrinkage

plot( dlist$D_obs ~ dlist$A ,
    xlab="Age of marriage" , ylab="Divorce rate" , lwd=3 , col=grau(0.8) )

post <- extract.samples(m15.1)

Aseq <- seq(from=-3,to=3.5,len=100)
mu <- link(m15.1,data=list(M=rep(0,100),A=Aseq))
shade(apply(mu,2,PI),Aseq,col=col.alpha(2,0.3))

Dest <- apply( post$D_true , 2, mean )
points( dlist$A , Dest , lwd=3 , col=2 )
for ( i in 1:nrow(d) ) {
    x <- dlist$A[i]
    lines( c(x,x) , c(dlist$D_obs[i],Dest[i]) , lwd=3 , col=col.alpha(2,0.5) )
}

mDAM <- ulam(
    alist(
        D_obs ~ dnorm( mu , sigma ),
        mu <- a + bA*A + bM*M,
        a ~ dnorm(0,0.2),
        bA ~ dnorm(0,0.5),
        bM ~ dnorm(0,0.5),
        sigma ~ dexp(1)
    ) , data=dlist , chains=4 , cores=4 )

Aseq <- seq(from=-3,to=3.5,len=100)
mu2 <- link(mDAM,data=list(M=rep(0,100),A=Aseq))
# shade(apply(mu2,2,PI),Aseq,col=col.alpha(1,0.3))
pi2 <- apply(mu2,2,PI)
lines( Aseq , pi2[1,] , lty=3 , col=grau(0.8) )
lines( Aseq , pi2[2,] , lty=3 , col=grau(0.8) )

# animate shrinkage

library(animation)
ani.record(reset=TRUE)
nf <- 50
col1 <- t(col2rgb(1)/255)
col2 <- t(col2rgb(2)/255)
for ( f in 0:nf ) {

    Dvals <- (1-f/nf)*dlist$D_obs + f/nf*Dest

    the_col <- rgb( (1-f/nf)*col1 + f/nf*col2 )
    plot( Dvals ~ dlist$A , xlab="Age of marriage" , ylab="Divorce rate" , lwd=3 , col=col.alpha(the_col,0.8) , ylim=range(dlist$D_obs) )

    shade(apply(mu,2,PI),Aseq,col=col.alpha(2, 0.3*f/nf ))
    shade(apply(mu2,2,PI),Aseq,col=col.alpha(1, 0.3*(1-f/nf) ))

    ani.record()
}

oopts = ani.options(interval = 0.01)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 8 -loop 0 frame*.png divorce_shrink.gif
# convert -delay 4 gp_demo.gif gp_demo.gif
# convert gp_demo.gif -coalesce -fuzz 2% +dither -layers Optimize +map gp_demo_z.gif


## R code 15.5
dlist2 <- list(
    D_obs = standardize( d$Divorce ),
    D_sd = d$Divorce.SE / sd( d$Divorce ),
    M_obs = standardize( d$Marriage ),
    M_sd = d$Marriage.SE / sd( d$Marriage ),
    A = standardize( d$MedianAgeMarriage ),
    N = nrow(d)
)

m15.2 <- ulam(
    alist(
        # D* model (observed)
        D_obs ~ dnorm( D_true , D_sd ),

        # D model (unobserved)
        vector[N]:D_true ~ dnorm( mu , sigma ),
        mu <- a + bA*A + bM*M_true[i],
        a ~ dnorm(0,0.2),
        bA ~ dnorm(0,0.5),
        bM ~ dnorm(0,0.5),
        sigma ~ dexp( 1 ),

        # M* model (observed)
        M_obs ~ dnorm( M_true , M_sd ),

        # M model (unobserved)
        vector[N]:M_true ~ dnorm( nu , tau ),
        nu <- aM + bAM*A,
        aM ~ dnorm(0,0.2),
        bAM ~ dnorm(0,0.5),
        tau ~ dexp( 1 )

    ) , data=dlist2 , chains=4 , cores=4 )

## R code 15.6
post <- extract.samples( m15.2 )
D_true <- apply( post$D_true , 2 , mean )
M_true <- apply( post$M_true , 2 , mean )
plot( dlist2$M_obs , dlist2$D_obs , pch=1 , lwd=3 , col=grau(0.8) ,
    xlab="marriage rate (std)" , ylab="divorce rate (std)" )
points( M_true , D_true , lwd=3 , col=2 )
for ( i in 1:nrow(d) )
    lines( c( dlist2$M_obs[i] , M_true[i] ) , c( dlist2$D_obs[i] , D_true[i] ) , lwd=2 , col=col.alpha(2,0.5) )

post_old <- extract.samples(mDAM)
dens( post_old$bA , lwd=4 , col=grau() , xlab="bA (effect of A)" )
dens( post$bA , lwd=4 , col=2 , add=TRUE )
abline(v=0,lty=3)

dens( post_old$bM , lwd=4 , col=grau() , xlab="bM (effect of M)" , xlim=c(-0.6,1) )
dens( post$bM, lwd=4 , col=2 , add=TRUE )
abline(v=0,lty=3)

########################################################################
# misclassification
# father analysis

# you don't have this data! - cannot be shared for privacy reasons
d <- read.csv( "Himba.csv" )
dat <- as.list(d)
dat$X <- dat$y

# without false positive rate
mXn <- ulam(
    alist(
        X ~ bernoulli(p),
        logit(p) <- a + z[mom_id]*sigma + x[dyad_id]*tau,
        a ~ normal(0,1.5),
        z[mom_id] ~ normal(0,1),
        sigma ~ normal(0,1),
        x[dyad_id] ~ normal(0,1),
        tau ~ normal(0,1)
    ) , data=dat , chains=4 , cores=4 , iter=4000 ,
    constraints=list(sigma="lower=0",tau="lower=0") )

# with false positive rate
dat$f <- 0.05
mX <- ulam(
    alist(
        X|X==1 ~ custom( log( p + (1-p)*f ) ),
        X|X==0 ~ custom( log( (1-p)*(1-f) ) ),
        logit(p) <- a + z[mom_id]*sigma + x[dyad_id]*tau,
        a ~ normal(0,1.5),
        z[mom_id] ~ normal(0,1),
        sigma ~ normal(0,1),
        x[dyad_id] ~ normal(0,1),
        tau ~ normal(0,1)
    ) , data=dat , chains=4 , cores=4 , iter=4000 ,
    constraints=list(sigma="lower=0",tau="lower=0") )

precis(mX)
precis(mXn)

post <- extract.samples( mX )
postn <- extract.samples( mXn )

dens( inv_logit(postn$a) , xlim=c(0,1) , lwd=4 , col=grau() , xlab="probability extra-pair paternity" )
dens( inv_logit(post$a) , add=TRUE , lwd=4 , col=2 )

# numerically kosher version
mX2 <- ulam(
    alist(
        #X|X==1 ~ custom( log( p + (1-p)*f ) ),
        X|X==1 ~ custom( log_sum_exp( log(p) , log1m(p)+log(f) ) ),
        #X|X==0 ~ custom( log( (1-p)*(1-f) ) ),
        X|X==0 ~ custom( log1m(p) + log1m(f) ),
        logit(p) <- a + z[mom_id]*sigma + x[dyad_id]*tau,
        a ~ normal(0,1.5),
        z[mom_id] ~ normal(0,1),
        sigma ~ normal(0,1),
        x[dyad_id] ~ normal(0,1),
        tau ~ normal(0,1)
    ) , data=dat , chains=4 , cores=4 , iter=4000 ,
    constraints=list(sigma="lower=0",tau="lower=0") )

precis(mX)
precis(mX2)

# double kosher
mX3 <- ulam(
    alist(
        #X|X==1 ~ custom( log( p + (1-p)*f ) ),
        X|X==1 ~ custom( log_sum_exp( log_p , log1m_exp(log_p)+log(f) ) ),

        #X|X==0 ~ custom( log( (1-p)*(1-f) ) ),
        X|X==0 ~ custom( log1m_exp(log_p) + log1m(f) ),
        
        log_p <- log_inv_logit( a + z[mom_id]*sigma + x[dyad_id]*tau ),
        a ~ normal(0,1.5),
        z[mom_id] ~ normal(0,1),
        sigma ~ normal(0,1),
        x[dyad_id] ~ normal(0,1),
        tau ~ normal(0,1)
    ) , data=dat , chains=4 , cores=4 , iter=4000 ,
    constraints=list(sigma="lower=0",tau="lower=0") )

