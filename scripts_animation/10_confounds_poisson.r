# more modeling events - sensitivity analysis, latent factors, poisson GLMs

library(rethinking)
library(animation)
library(ellipse)

# continuing from UCBadmit example
# what happens when there is a confound?

set.seed(17)
N <- 2000 # number of applicants
# even gender distribution
G <- sample( 1:2 , size=N , replace=TRUE )
# sample ability, high (1) to average (0)
u <- rbern(N,0.1)
# gender 1 tends to apply to department 1, 2 to 2
# and G=1 with greater ability tend to apply to 2 as well
D <- rbern( N , ifelse( G==1 , u*0.5 , 0.8 ) ) + 1
# matrix of acceptance rates [dept,gender]
accept_rate_u0 <- matrix( c(0.1,0.1,0.1,0.3) , nrow=2 )
accept_rate_u1 <- matrix( c(0.2,0.3,0.2,0.5) , nrow=2 )
# simulate acceptance
p <- sapply( 1:N , function(i) 
    ifelse( u[i]==0 , accept_rate_u0[D[i],G[i]] , accept_rate_u1[D[i],G[i]] ) )
A <- rbern( N , p )

table(G,D)
table(G,A)

dat_sim <- list( A=A , D=D , G=G )

# total effect gender
m1 <- ulam(
    alist(
        A ~ bernoulli(p),
        logit(p) <- a[G],
        a[G] ~ normal(0,1)
    ), data=dat_sim , chains=4 , cores=4 )

post1 <- extract.samples(m1)
post1$fm_contrast <- post1$a[,1] - post1$a[,2]
precis(post1)

# direct effects - now confounded!
m2 <- ulam(
    alist(
        A ~ bernoulli(p),
        logit(p) <- a[G,D],
        matrix[G,D]:a ~ normal(0,1)
    ), data=dat_sim , chains=4 , cores=4 )

precis(m2,3)

post2 <- extract.samples(m2)
post2$fm_contrast_D1 <- post2$a[,1,1] - post2$a[,2,1]
post2$fm_contrast_D2 <- post2$a[,1,2] - post2$a[,2,2]
precis(post2)

dens( post2$fm_contrast_D1 , lwd=4 , col=4 , xlab="F-M contrast in each department" )
dens( post2$fm_contrast_D2 , lwd=4 , col=2 , add=TRUE )
abline(v=0,lty=3)

dens( post2$a[,1,1] , lwd=4 , col=2 , xlim=c(-3,1) )
dens( post2$a[,2,1] , lwd=4 , col=4 , add=TRUE )
dens( post2$fm_contrast_D1 , lwd=4 , add=TRUE )

dens( post2$a[,1,2] , lwd=4 , col=2 , add=TRUE , lty=4 )
dens( post2$a[,2,2] , lwd=4 , col=4 , add=TRUE , lty=4)
dens( post2$fm_contrast_D2 , lwd=4 , add=TRUE , lty=4)

# if we could measure u
dat_sim$u <- u
m3 <- ulam(
    alist(
        A ~ bernoulli(p),
        logit(p) <- a[G,D] + buA*u,
        matrix[G,D]:a ~ normal(0,1),
        buA ~ normal(0,1)
    ), data=dat_sim , constraints=list(buA="lower=0") ,
    chains=4 , cores=4 )

post3 <- extract.samples(m3)
post3$fm_contrast_D1 <- post3$a[,1,1] - post3$a[,2,1]
post3$fm_contrast_D2 <- post3$a[,1,2] - post3$a[,2,2]

precis(post3)

dens( post2$fm_contrast_D1 , lwd=1 , col=4 , xlab="F-M contrast in each department" , xlim=c(-2,1) )
dens( post2$fm_contrast_D2 , lwd=1 , col=2 , add=TRUE )
abline(v=0,lty=3)
dens( post3$fm_contrast_D1 , lwd=4 , col=4 , add=TRUE )
dens( post3$fm_contrast_D2 , lwd=4 , col=2 , add=TRUE )

# sensitivity

dat_sim$D2 <- ifelse( D==2 , 1 , 0 )
dat_sim$b <- c(1,1)
dat_sim$g <- c(1,0)
dat_sim$N <- length(dat_sim$D2)

m3s <- ulam(
    alist( 
        # A model
        A ~ bernoulli(p),
        logit(p) <- a[G,D] + b[G]*u[i],
        matrix[G,D]:a ~ normal(0,1),

        # D model
        D2 ~ bernoulli(q),
        logit(q) <- delta[G] + g[G]*u[i],
        delta[G] ~ normal(0,1),

        # declare unobserved u
        vector[N]:u ~ normal(0,1)
    ), data=dat_sim , chains=4 , cores=4 )

precis(m3s,3,pars=c("a","delta"))

post3s <- extract.samples(m3s)
post3s$fm_contrast_D1 <- post3s$a[,1,1] - post3s$a[,2,1]
post3s$fm_contrast_D2 <- post3s$a[,1,2] - post3s$a[,2,2]

dens( post2$fm_contrast_D1 , lwd=1 , col=4 , xlab="F-M contrast in each department" , xlim=c(-2,1) )
dens( post2$fm_contrast_D2 , lwd=1 , col=2 , add=TRUE )
abline(v=0,lty=3)
dens( post3s$fm_contrast_D1 , lwd=4 , col=4 , add=TRUE )
dens( post3s$fm_contrast_D2 , lwd=4 , col=2 , add=TRUE )

plot( jitter(u) , apply(post3s$u,2,mean) , col=ifelse(G==1,2,4) , lwd=3 )

##################
# now sensitivity on real data
# must convert to logistic regression

data(UCBadmit)
d <- UCBadmit

dat <- list( 
    A = d$admit,
    N = d$applications,
    G = ifelse(d$applicant.gender=="female",1,2),
    D = as.integer(d$dept)
)

# total effect gender
mG <- ulam(
    alist(
        A ~ binomial(N,p),
        logit(p) <- a[G],
        a[G] ~ normal(0,1)
    ), data=dat , chains=4 , cores=4 )

# direct effects
mGD <- ulam(
    alist(
        A ~ binomial(N,p),
        logit(p) <- a[G,D],
        matrix[G,D]:a ~ normal(0,1)
    ), data=dat , chains=4 , cores=4 )

# this code converts UCBadmit to long format (logistic format)
library(tidyverse)
dat_long1 <- uncount( d , admit )
dat_long1$Y <- 1
dat_long1$reject <- NULL
dat_long0 <- uncount( d , reject )
dat_long0$Y <- 0
dat_long0$admit <- NULL
dat_long01 <- rbind( dat_long1 , dat_long0 )
dat_long01$applications <- NULL
dat_long01$admit <- dat_long01$Y
dat_long01$Y <- NULL

# write out long table and include in future version of rethinking package
# write.table(dat_long01,file="UCBadmit_long.csv",sep=";",row.names=FALSE)
# data(UCBadmit_long)

datl <- list(
    A = dat_long01$admit,
    G = ifelse(dat_long01$applicant.gender=="female",1,2),
    D = as.integer(dat_long01$dept)
)

datl$D1 <- ifelse(datl$D==1,1,0)
datl$N <- length(datl$D)
datl$b <- c(1,1)
datl$g <- c(1,0)

mGDu <- ulam(
    alist( 
        # A model
        A ~ bernoulli(p),
        logit(p) <- a[G,D] + b[G]*u[i],
        matrix[G,D]:a ~ normal(0,1),

        # D model
        D1 ~ bernoulli(q),
        logit(q) <- delta[G] + g[G]*u[i],
        delta[G] ~ normal(0,1),

        # declare unobserved u
        vector[N]:u ~ normal(0,1)
    ), data=datl , chains=4 , cores=4 )

precis(mGDu,3,pars=c("delta"))
precis(mGD,3)

postGD <- extract.samples(mGD)
postGD$fm_contrast_D1 <- postGD$a[,1,1] - postGD$a[,2,1]
postGD$fm_contrast_D2 <- postGD$a[,1,2] - postGD$a[,2,2]

postGDu <- extract.samples(mGDu)
postGDu$fm_contrast_D1 <- postGDu$a[,1,1] - postGDu$a[,2,1]
postGDu$fm_contrast_D2 <- postGDu$a[,1,2] - postGDu$a[,2,2]

dens( postGD$fm_contrast_D1 , lwd=2 , col=2 , xlab="F-M contrast in department A" , xlim=c(-1,1.75) )
#dens( postGD$fm_contrast_D2 , lwd=2 , col=4 , add=TRUE )
abline(v=0,lty=3)
dens( postGDu$fm_contrast_D1 , lwd=5 , col=2 , add=TRUE )
#dens( postGDu$fm_contrast_D2 , lwd=5 , col=4 , add=TRUE )

# extract depat A u estimates
postuA <- postGDu$u[,datl$D==1]
u_est <- apply(postuA,2,mean)
u_PI <- apply(postuA,2,PI,0.5)
# sort by posterior mean
Gi <- datl$G[datl$D==1]
o <- order( u_est )

plot( u_est[o] , lwd=1 , pch=16 , col=ifelse(Gi[o]==1,2,4) , xlim=c(1,sum(datl$D==1)) , ylim=c(-1.2,1.6) , xlab="application to department A" , ylab="posterior u" )
for ( i in 1:length(u_est) ) lines( c(i,i) , u_PI[,o[i]] , col=ifelse(Gi[o[i]]==1,2,4) )

# raw Stan version - just as an example

mc <- "
data{
    int A[4526];
    int G[4526];
    int D[4526];
}
parameters{
    matrix[2,6] a;
    vector[6] theta[2];
    vector[4526] u;
    real<lower=0> b;
}
model{
    vector[4526] p;
    vector[6] q;
    u ~ normal(0,1);
    b ~ normal(0,0.5);
    for ( i in 1:2 ) theta[i] ~ normal( 0 , 1 );
    for ( i in 1:4526 ) {
        q = theta[G[i]];
        if ( G[i]==1 ) q[1] = q[1] + u[i];
        D[i] ~ categorical_logit( q );
    }
    to_vector( a ) ~ normal( 0 , 1 );
    for ( i in 1:4526 ) {
        p[i] = a[G[i], D[i]] + u[i];
        p[i] = inv_logit(p[i]);
    }
    A ~ bernoulli( p );
}
"

m5 <- cstan( model_code=mc , data=datl , chains=4 , cores=4 )

precis(mGD,3)
precis(m5,3,pars=c("a","theta","b"))



#############
# proxies

T1 <- rnorm(N,u,0.1)
T2 <- rnorm(N,u,0.5)
T3 <- rnorm(N,u,0.25)

dat_sim$T1 <- T1
dat_sim$T2 <- T2
dat_sim$T3 <- T3
dat_sim$N <- N
dat_sim$u <- NULL

m4 <- ulam(
    alist(
        # A model
        A ~ bernoulli(p),
        logit(p) <- a[G,D] + b*u[i],
        matrix[G,D]:a ~ normal(0,1),
        b ~ normal(0,1),
        # u and T model
        vector[N]:u ~ normal(0,1),
        T1 ~ normal(u,tau[1]),
        T2 ~ normal(u,tau[2]),
        T3 ~ normal(u,tau[3]),
        vector[3]:tau ~ exponential(1)
    ), data=dat_sim , chains=4 , cores=4 , constraints=list(b="lower=0") )

precis(m4,3,pars=c("a","b","tau"))

post4 <- extract.samples(m4)
post4$fm_contrast_D1 <- post4$a[,1,1] - post4$a[,2,1]
post4$fm_contrast_D2 <- post4$a[,1,2] - post4$a[,2,2]

dens( post2$fm_contrast_D2 , lwd=1 , col=2 , xlab="F-M contrast in each department" , xlim=c(-2,1) )
#dens( post2$fm_contrast_D2 , lwd=1 , col=2 , add=TRUE )
abline(v=0,lty=3)
dens( post4$fm_contrast_D2 , lwd=4 , col=2 , add=TRUE )
#dens( post4$fm_contrast_D2 , lwd=4 , col=2 , add=TRUE )

post4 <- extract.samples(m4)
plot( jitter(u) , apply(post4$u,2,mean) , xlab="u (true)" , ylab="posterior mean u" , lwd=2 , col=ifelse(G==1,col.alpha(2,0.6),col.alpha(4,0.6)) , xaxt="n" )
axis(1,at=0:1,labels=0:1)

###########
# poisson regression example

# priors

curve( dlnorm(x,0,10) , from=0 , to=100 , lwd=4 , col=2 , xlab="mean number of tools" , ylab="Density" , ylim=c(0,0.05) )
curve( dlnorm(x,3,0.5) , add=TRUE , lwd=4 , col=4 )

N <- 10
a <- rnorm( N , 0 , 1 )
b <- rnorm( N , 0 , 10 )
plot( NULL , xlim=c(-2,2) , ylim=c(0,100) , xlab="x value" , ylab="expected count" )
for ( i in 1:N ) curve( exp( a[i] + b[i]*x ) , add=TRUE , col=2 , lwd=4 )

N <- 10
a <- rnorm( N , 3 , 0.5 )
b <- rnorm( N , 0 , 0.2 )
plot( NULL , xlim=c(-2,2) , ylim=c(0,100) , xlab="x value" , ylab="expected count" )
for ( i in 1:N ) curve( exp( a[i] + b[i]*x ) , add=TRUE , col=2 , lwd=4 )

# model

library(rethinking)
data(Kline)
d <- Kline
d$P <- scale( log(d$population) )
d$contact_id <- ifelse( d$contact=="high" , 2 , 1 )

dat <- list(
    T = d$total_tools ,
    P = d$P ,
    C = d$contact_id )

# intercept only
m11.9 <- ulam(
    alist(
        T ~ dpois( lambda ),
        log(lambda) <- a,
        a ~ dnorm( 3 , 0.5 )
    ), data=dat , chains=4 , log_lik=TRUE )

# interaction model
m11.10 <- ulam(
    alist(
        T ~ dpois( lambda ),
        log(lambda) <- a[C] + b[C]*P,
        a[C] ~ dnorm( 3 , 0.5 ),
        b[C] ~ dnorm( 0 , 0.2 )
    ), data=dat , chains=4 , log_lik=TRUE )

compare( m11.9 , m11.10 , func=PSIS )

k <- PSIS( m11.10 , pointwise=TRUE )$k
plot( dat$P , dat$T , xlab="log population (std)" , ylab="total tools" ,
    col=ifelse( dat$C==1 , 4 , 2 ) , lwd=4+4*normalize(k) ,
    ylim=c(0,75) , cex=1+normalize(k) )

# set up the horizontal axis values to compute predictions at
P_seq <- seq( from=-1.4 , to=3 , len=100 )

# predictions for C=1 (low contact)
lambda <- link( m11.10 , data=data.frame( P=P_seq , C=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( P_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , P_seq , xpd=TRUE , col=col.alpha(4,0.3) )

# predictions for C=2 (high contact)
lambda <- link( m11.10 , data=data.frame( P=P_seq , C=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( P_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , P_seq , xpd=TRUE , col=col.alpha(2,0.3))

identify( dat$P , dat$T , d$culture )

# natural scale now

plot( d$population , d$total_tools , xlab="population" , ylab="total tools" ,
    col=ifelse( dat$C==1 , 4 , 2 ) , lwd=4+4*normalize(k) ,
    ylim=c(0,75) , cex=1+normalize(k) )
P_seq <- seq( from=-5 , to=3 , length.out=100 )
# 1.53 is sd of log(population)
# 9 is mean of log(population)
pop_seq <- exp( P_seq*1.53 + 9 )
lambda <- link( m11.10 , data=data.frame( P=P_seq , C=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE , col=col.alpha(4,0.3))

lambda <- link( m11.10 , data=data.frame( P=P_seq , C=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE , col=col.alpha(2,0.3) )

identify( d$population , d$total_tools , d$culture )

# innovation/loss model

dat2 <- list( T=d$total_tools, P=d$population, C=d$contact_id )
m11.11 <- ulam(
    alist(
        T ~ dpois( lambda ),
        lambda <- exp(a[C])*P^b[C]/g,
        a[C] ~ dnorm(1,1),
        b[C] ~ dexp(1),
        g ~ dexp(1)
    ), data=dat2 , chains=4 , cores=4 , log_lik=TRUE )

precis(m11.11,2)

plot( d$population , d$total_tools , xlab="population" , ylab="total tools" ,
    col=ifelse( dat$C==1 , 4 , 2 ) , lwd=4+4*normalize(k) ,
    ylim=c(0,75) , cex=1+normalize(k) )
P_seq <- seq( from=-5 , to=3 , length.out=100 )

pop_seq <- exp( P_seq*1.53 + 9 )
lambda <- link( m11.11 , data=data.frame( P=pop_seq , C=1 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=2 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE , col=col.alpha(4,0.3))

lambda <- link( m11.11 , data=data.frame( P=pop_seq , C=2 ) )
lmu <- apply( lambda , 2 , mean )
lci <- apply( lambda , 2 , PI )
lines( pop_seq , lmu , lty=1 , lwd=1.5 )
shade( lci , pop_seq , xpd=TRUE , col=col.alpha(2,0.3) )

identify( d$population , d$total_tools , d$culture )

