# week 6 - multilevel models 1

# 1
# prior predictive varying effects

n <- 1e4
sigma <- rexp(n,1) # show also with 1/10
#sigma <- rnorm(n,1,0.1)
abar <- rnorm(n,0,1)
aT <- rnorm(n,abar,sigma)
dens(inv_logit(aT),xlim=c(0,1),adj=0.1,lwd=4,col=2,xlab="probability survival")

mtext(concat("sigma~exponential(0.1)"))

# 2

library(rethinking)
data(reedfrogs)
d <- reedfrogs

dat <- list(
    S = d$surv,
    D = d$density,
    T = 1:nrow(d),
    P = ifelse( d$pred=="no" , 0L , 1L ),
    G = ifelse( d$size=="small" , 1L , 2L )
)

dat$P <- ifelse( d$pred=="no" , 1L , 2L )
m2 <- ulam(
    alist(
        S ~ binomial( D , p ),
        logit(p) <- a[T] + b[P,G],
        a[T] ~ normal( 0 , sigma ),
        matrix[P,G]:b ~ normal( 0 , 1 ),
        #a_bar ~ normal( 0 , 1.5 ),
        sigma ~ exponential( 1 )
    ), data=dat , chains=4 , cores=4 , log_lik=TRUE )

precis(m2,3,pars=c("b","sigma"))

# 3
# D as effect

dat$Do <- standardize(log(d$density))
m3 <- ulam(
    alist(
        S ~ binomial( D , p ),
        logit(p) <- a[T] + b[P,G] + bD[P]*Do,
        a[T] ~ normal( 0 , sigma ),
        matrix[P,G]:b ~ normal( 0 , 1 ),
        bD[P] ~ normal(0,0.5),
        sigma ~ exponential( 1 )
    ), data=dat , chains=4 , cores=4 , log_lik=TRUE )

precis(m3,3,pars=c("b","bD","sigma"))

compare( m2 , m3 , func=PSIS )

# 4



# setup status quo groups
nreps <- 20
D <- rep( c(10,35,10,35) , times=nreps )
Do <- rep( ( D - mean(d$density) ) / sd(d$density) , times=nreps )
G <- rep( c(1,1,2,2) , times=nreps )
P <- rep( rep(2,4) , times=nreps )

# simulate tanks
post <- extract.samples(m3)
aT <- replicate( length(D) , rnorm(2000,0,post$sigma) )

# simulate survival status quo
p1 <- sapply( 1:length(D) , 
    function(i) 
        inv_logit( 
            aT[,i] + 
            post$b[,P[i],G[i]] + 
            post$bD[,P[i]]*Do[i] ) )

# simulate survival intervention
P <- rep( rep(1,4) , times=nreps )
p2 <- sapply( 1:length(D) , 
    function(i) 
        inv_logit( 
            aT[,i] + 
            post$b[,P[i],G[i]] + 
            post$bD[,P[i]]*Do[i] ) )

dens( p1 , lwd=4 , col=4 , xlab="probability of survival" , xlim=c(0,1) )
mtext("status quo")

dens( p2 , lwd=4 , col=2 , xlab="probability of survival" , xlim=c(0,1) )
mtext("intervention")

dens( p2 - p1 , lwd=4 , col=6 , xlab="change in probability of survival" , xlim=c(-0.2,1) )
abline(v=0,lty=3,lwd=3)
mtext("effect of intervention")

# shuffle and take difference
s1 <- sample( 1:80 , size=80 )
s2 <- sample( 1:80 , size=80 )
dens( p2[,s2] - p1[,s1] , lwd=4 , col=6 , xlab="change in probability of survival" )
abline(v=0,lty=3,lwd=3)
mtext("effect ignoring order")
