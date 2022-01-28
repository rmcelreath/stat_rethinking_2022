# WEEK 3 SOLUTIONS

library(rethinking)
data(foxes)
d <- foxes
d$W <- standardize(d$weight)
d$A <- standardize(d$area)
d$F <- standardize(d$avgfood)
d$G <- standardize(d$groupsize)

# 1 

m1 <- quap(
    alist(
        F ~ dnorm( mu , sigma ),
        mu <- a + bA*A,
        a ~ dnorm(0,0.2),
        bA ~ dnorm(0,0.5),
        sigma ~ dexp(1)
    ), data=d )

# 2

# total effect

m2 <- quap(
    alist(
        W ~ dnorm( mu , sigma ),
        mu <- a + bF*F,
        a ~ dnorm(0,0.2),
        bF ~ dnorm(0,0.5),
        sigma ~ dexp(1)
    ), data=d )

# direct effect of F

m2b <- quap(
    alist(
        W ~ dnorm( mu , sigma ),
        mu <- a + bF*F + bG*G,
        a ~ dnorm(0,0.2),
        c(bF,bG) ~ dnorm(0,0.5),
        sigma ~ dexp(1)
    ), data=d )

# direct effect of F on G

m2c <- quap(
    alist(
        G ~ dnorm( mu , sigma ),
        mu <- a + bF*F,
        a ~ dnorm(0,0.2),
        bF ~ dnorm(0,0.5),
        sigma ~ dexp(1)
    ), data=d )

# 3


library(dagitty)
t2f_dag <- dagitty( "dag {
    X -> Y
    S -> X
    S -> Y
    A -> Y
    A -> X
    A -> S
    S <- U -> Y
}")
adjustmentSets( t2f_dag , exposure="X" , outcome="Y" )


# 4
# Writeasyntheticdatasimulationforthecausal model shown in Problem 3. Be sure to include the unobserved confound in the simulation. Choose any functional relationships that you like—you don’t have to get the epidemiology correct. You just need to honor the causal structure. Then design a regression model to estimate the influence of X on Y and use it on your synthetic data. How large of a sample do you need to reliably estimate P(Y|do(X))? Define “reliably” as you like, but justify your definition.

f <- function(N=100,bX=0) {
    U <- rnorm(N)
    A <- rnorm(N)
    S <- rnorm(N,A+U)
    X <- rnorm(N,A+S)
    Y <- rnorm(N,A+S+bX*X+U)
    return(list(
        A=standardize(A),
        S=standardize(S),
        X=standardize(X),
        Y=standardize(Y)))
}

flist <- alist(
        Y ~ dnorm(mu,exp(log_sigma)),
        mu <- a + bX*X + bS*S + bA*A,
        a ~ dnorm(0,0.2),
        c(bX,bS,bA) ~ dnorm(0,0.5),
        log_sigma ~ dnorm(0,1)
    )

sim_dat <- f(N=10,bX=0)
m4 <- quap( flist , data=sim_dat )
precis(m4)

# function to repeat analysis at specific N
g <- function(N_reps=1e3,N=100,bX=0) {
    r <- mcreplicate( N_reps , mc.cores=8 , {
        sim_dat <- f(N=N,bX=bX)
        m <- quap( flist , data=sim_dat )
        # return width of 89% interval and mean distance from true value
        iw <- as.numeric( precis(m)[2,4] - precis(m)[2,3] )
        md <- as.numeric( abs(precis(m)[2,1] - 0) )
        return(c(md,iw))
    } )
    return(r)
}

x <- g(N=10)

plot( x[2,] , lwd=2 , col=4 , ylim=c(0,max(x)) , xlab="simulation" , ylab="value" )
points( 1:ncol(x) , x[1,] , lwd=2 , col=2 )

abline( h=mean(x[1,]) , lwd=6 , col="white" )
abline( h=mean(x[1,]) , lwd=3 , col=2 )

abline( h=mean(x[2,]) , lwd=6 , col="white" )
abline( h=mean(x[2,]) , lwd=3 , col=4 )

# now plot mean bias and precision across values of N
N_seq <- c(10,20,50,100,500,1000)

y <- sapply( N_seq , function(n) {
    x <- g(N=n)
    return( c( apply(x,1,mean) , apply(x,1,sd) ) ) 
} )

plot( N_seq , y[2,] , lwd=3 , col=4 , type="b" , ylim=c(0,1) , xlab="N" , ylab="value" )
points( N_seq , y[1,] , lwd=3 , col=2 , type="b" )
for ( i in 1:length(N_seq) ) {
    lines( c(N_seq[i],N_seq[i]) , c(y[1,i]+y[3,i],y[1,i]-y[3,i]) , lwd=3 , col=2 )
    lines( c(N_seq[i],N_seq[i]) , c(y[2,i]+y[4,i],y[2,i]-y[4,i]) , lwd=3 , col=4 )
}
