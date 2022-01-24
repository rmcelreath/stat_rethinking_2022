# sim confounding by post-treatment variable

f <- function(n=100,bXZ=1,bZY=1) {
    X <- rnorm(n)
    u <- rnorm(n)
    Z <- rnorm(n, bXZ*X + u)
    Y <- rnorm(n, bZY*Z + u )
    bX <- coef( lm(Y ~ X) )['X']
    bXZ <- coef( lm(Y ~ X + Z) )['X']
    return( c(bX,bXZ) )
}

sim <- mcreplicate( 1e4 , f(bZY=0) , mc.cores=8 )

dens( sim[1,] , lwd=3 , xlab="posterior mean" , xlim=c(-1,0.8) , ylim=c(0,2.7)  )
dens( sim[2,] , lwd=3 , col=2 , add=TRUE )

# case control bias

f <- function(n=100,bXY=1,bYZ=1) {
    X <- rnorm(n)
    Y <- rnorm(n, bXY*X )
    Z <- rnorm(n, bYZ*Y )
    bX <- coef( lm(Y ~ X) )['X']
    bXZ <- coef( lm(Y ~ X + Z) )['X']
    return( c(bX,bXZ) )
}

sim <- mcreplicate( 1e4 , f() , mc.cores=8 )

dens( sim[1,] , lwd=3 , xlab="posterior mean" , xlim=c(0,1.5) , ylim=c(0,5)  )
dens( sim[2,] , lwd=3 , col=2 , add=TRUE )

# precision parasite

f <- function(n=100,bZX=1,bXY=1) {
    Z <- rnorm(n)
    X <- rnorm(n, bZX*Z )
    Y <- rnorm(n, bXY*X )
    bX <- coef( lm(Y ~ X) )['X']
    bXZ <- coef( lm(Y ~ X + Z) )['X']
    return( c(bX,bXZ) )
}

sim <- mcreplicate( 1e4 , f(n=50) , mc.cores=8 )

dens( sim[1,] , lwd=3 , xlab="posterior mean" , xlim=c(0.5,1.5) , ylim=c(0,4.5)  )
dens( sim[2,] , lwd=3 , col=2 , add=TRUE )

# bias amplifier

f <- function(n=100,bZX=1,bXY=1) {
    Z <- rnorm(n)
    u <- rnorm(n)
    X <- rnorm(n, bZX*Z + u )
    Y <- rnorm(n, bXY*X + u )
    bX <- coef( lm(Y ~ X) )['X']
    bXZ <- coef( lm(Y ~ X + Z) )['X']
    return( c(bX,bXZ) )
}

sim <- mcreplicate( 1e4 , f(bXY=0) , mc.cores=8 )

dens( sim[1,] , lwd=3 , xlab="posterior mean" , xlim=c(-0.5,1) , ylim=c(0,5.5)  )
dens( sim[2,] , lwd=3 , col=2 , add=TRUE )

abline_w <- function(...,col=1,lwd=1,dlwd=2) {
    abline(...,col="white",lwd=lwd+dlwd)
    abline(...,col=col,lwd=lwd)
}

n <- 1000
Z <- rbern(n)
u <- rnorm(n)
X <- rnorm(n, 7*Z + u )
Y <- rnorm(n, 0*X + u )

cols <- c( col.alpha(2,0.5) , col.alpha(4,0.5) )
plot( X , Y  , col=cols[Z+1] , lwd=2 )

abline_w( lm(Y~X) , lwd=3 )

abline_w( lm(Y[Z==1]~X[Z==1]) , lwd=3 , col=4 )

abline_w( lm(Y[Z==0]~X[Z==0]) , lwd=3 , col=2 )

