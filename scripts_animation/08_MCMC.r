# new MCMC illustrations

# king markov

num_weeks <- 1e5
positions <- rep(0,num_weeks)
current <- 10
for ( i in 1:num_weeks ) {
  ## record current position
    positions[i] <- current
  ## flip coin to generate proposal
    proposal <- current + sample( c(-1,1) , size=1 )
  ## now make sure he loops around the archipelago
    if ( proposal < 1 ) proposal <- 10
    if ( proposal > 10 ) proposal <- 1
  ## move?
    prob_move <- proposal/current
    current <- ifelse( runif(1) < prob_move , proposal , current )
}


# animate

library(rethinking)
library(animation)

my_circle <- function(x=0,y=0,r=1,angle=0,...) {
    a <- seq(angle, angle + 2 * pi, length = 360)
    #lines( r*cos(a)+x , r*sin(a)+y , ... )
    polygon( r*cos(a)+x , r*sin(a)+y , ... )
}
# plot(NULL,xlim=c(-1,1),ylim=c(-1,1)); my_circle(0,0,0.5,lty=2)
p2c <- function(a) c(cos(a) , sin(a))

blank(h=2,w=4,bty="n")
par(mfrow=c(1,2))

aseq <- seq(from=0,to=2*pi,len=11)
drawbg <- function() {
    plot(NULL,xlim=c(-1.2,1.3),ylim=c(-1.2,1.2),xaxt="n",yaxt="n",xlab="",ylab="")
    for ( i in 1:10 ) {
        xy <- p2c( aseq[i] )
        r <- sqrt((i/40)/pi)
        my_circle( xy[1] , xy[2] , r=r , lwd=3 , col="white" )
        text( xy[1] , xy[2] , i )
    }
}

drawhist <- function(n=1) {
    ymax <- max(10,max( table(positions[1:n]) ))
    plot( NULL , xlim=c(1,10) , ylim=c(0,ymax+2) , xlab="Island" , ylab="Visits" , xaxt="n" )
    axis( 1 , at=1:10 , labels=1:10 )
    for ( i in 1:10 ) {
    if ( i==positions[n] )
        if ( i==positions[n] )
            lines( c(i,i) , c(0,sum(positions[1:n]==i)) , lwd=16 , col=col.alpha(2,0.5) )
        lines( c(i,i) , c(0,sum(positions[1:n]==i)) , lwd=8 , col=ifelse(i==positions[n],2,1) )
    }
}

ani.record(reset=TRUE)
par(mfrow=c(1,2))
par(bg="light blue")
p <- positions
for ( f in 5000:5200 ) {

    # draw previous position
    drawbg()
    xy <- p2c( aseq[ p[f-1] ] )
    points( xy[1] , xy[2] , cex=3 , pch=16 , col=2 )
    drawhist(f-1)
    ani.record()

    # draw line to next position
    drawbg()
    xy2 <- p2c( aseq[ p[f] ] )
    points( xy[1] , xy[2] , cex=3 , pch=16 , col=2 )
    lines( c(xy[1],xy2[1]) , c(xy[2],xy2[2]) , lwd=3 , col=2 )
    drawhist(f-1)
    ani.record()

    # draw next position
    drawbg()
    points( xy[1] , xy[2] , cex=3 , pch=16 , col=2 )
    lines( c(xy[1],xy2[1]) , c(xy[2],xy2[2]) , lwd=3 , col=2 )
    points( xy2[1] , xy2[2] , cex=3 , pch=16 , col=2 )
    drawhist(f)
    ani.record()

    # update histogram
    ani.record()

}#f

oopts = ani.options(interval = 0.05)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 5 -loop 0 frame*.png king_markov.gif
# convert -delay 5 king_markov.gif king_markov.gif

# plot chain

blank(bty="n",w=3)

plot( positions[1:2000] , type="l" , lwd=4 , col=2 , xlab="time" , ylab="island" )

# animate chain

ani.record(reset=TRUE)
xmax <- 300
start <- 1
by <- 50
for ( f in 200:1000 ) {
    end <- f
    if ( end > xmax ) start <- start + 1
    plot( positions[ start:end ] , type="l" , lwd=4 , col=2 , xlab="" , ylab="island" , xaxt="n" , ylim=c(1,10) , yaxt="n" , xlim=c(0,xmax) )
    labs <- seq(from=start-1,to=max(end,xmax) , by=by )
    axis( 1 , at=seq(from=0,to=xmax , by=by ) , labels=labs )
    axis( 2 , at=1:10 , labels=1:10 )
    points( min(end,xmax) , positions[end] , cex=1.5 , lwd=3 , col=2 )
    ani.record()
}

oopts = ani.options(interval = 0.01)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 9 -loop 0 frame*.png king_markov_chain_200.gif

# convert -delay 15 king_markov_chain_50.gif king_markov_chain_50.gif

##############
# HMC visualizations

##################
# a 2D example: mu and sigma
# Uses HMC2 function from rethinking
#
# U needs to return neg-log-probability
# y ~ normal( mu , exp(log_sigma) )
# mu ~ normal( a, b )
# log_sigma ~ normal( k, d )
myU2 <- function( q , a=0 , b=1 , k=0 , d=0.5 ) {
    s <- exp(q[2]) # sigma on log latent scale
    mu <- q[1]
    U <- sum( dnorm(y,mu,s,log=TRUE) ) + dnorm(mu,a,b,log=TRUE) + dnorm(q[2],k,d,log=TRUE)
    return( -U )
}

# gradient function
# need vector of partial derivatives of U with respect to vector q
myU_grad2 <- function( q , a=0 , b=1 , k=0 , d=0.5 ) {
    mu <- q[1]
    s <- exp(q[2])
    G1 <- sum( y - mu ) * exp(-2*q[2]) + (a - mu)/b^2 #dU/dmu
    G2 <- sum( (y - mu)^2 ) * exp(-2*q[2]) - length(y) + (k-q[2])/d^2 #dU/ds
    return( c( -G1 , -G2 ) ) # negative bc energy is neg-log-prob
}

# test data
set.seed(7)
y <- abs(rnorm(50))
y <- c( y , -y ) # ensure mean is zero

###########
# example paths
library(shape) # for good arrow heads
# blank(bty="n")

# priors
priors <- list()
priors$a <- 0
priors$b <- 1
priors$k <- 0
priors$d <- 0.3

#ss <- ss + 1
set.seed(42) # seed 9 for examples

# init
n_samples <- 4
Q <- list()
Q$q <- c(-0.4,0.2)
xr <- c(-0.6,0.6)
yr <- c(-0.25,0.4)

step <- 0.02
L <- 12 # 0.02/12 okay sampling --- 0.02/20 is good for showing u-turns
xpos <- c(4,2,1,2) # for L=20
#xpos <- c(2,3,2,1) # for L=55
path_col <- col.alpha("black",0.5)

draw_bg <- function() {
    plot( NULL , ylab="log_sigma" , xlab="mu" , xlim=xr , ylim=yr )
    # draw contour of log-prob
    cb <- 0.2
    mu_seq <- seq(from=xr[1]-cb,to=xr[2]+cb,length.out=50) 
    logsigma_seq <- seq(from=yr[1]-cb,to=yr[2]+cb,length.out=50)
    z <- matrix(NA,length(mu_seq),length(logsigma_seq))
    for ( i in 1:length(mu_seq) )
        for ( j in 1:length(logsigma_seq) )
            z[i,j] <- myU2( c( mu_seq[i] , logsigma_seq[j] ) , a=priors$a , b=priors$b , k=priors$k , d=priors$d )
    cl <- contourLines( mu_seq , logsigma_seq , z , nlevels=30 )
    for ( i in 1:length(cl) ) lines( cl[[i]]$x , cl[[i]]$y , col=col.alpha("black",0.5) , lwd=1 )
}

step <- 0.01
L <- 12 # 0.02/12 okay sampling --- 0.02/20 is good for showing u-turns
xpos <- c(4,2,1,2) # for L=20
#xpos <- c(2,3,2,1) # for L=55
path_col <- col.alpha("black",0.5)

# animation loop
ani.record(reset=TRUE)

set.seed(12)

Q <- list()
Q$q <- c(-0.4,0.2) # start point
xr <- c(-0.4,0.4) # x range in plot
yr <- c(-0.25,0.3) # y range in plot

draw_bg()

n_samples <- 10
# points( Q$q[1] , Q$q[2] , pch=4 , col="black" )
pts <- matrix(NA,nrow=n_samples,ncol=3)

for ( i in 1:n_samples ) {

    Q <- HMC2( myU2 , myU_grad2 , step , L , Q$q , a=priors$a , b=priors$b , k=priors$k , d=priors$d )

    draw_bg()

    # draw previous points
    if ( i > 1 ) {
        for ( j in 1:(i-1) ) {
            V <- 0.9
            points( pts[j,1] , pts[j,2] , pch=ifelse( pts[j,3]==1 , 1 , 16 ) , col=grau(V) , lwd=2 )
        }
    }

    # draw trajectory
    for ( l in 1:L ) {
        lines( Q$traj[l:(l+1),1] , Q$traj[l:(l+1),2] , col="white" , lwd=8 )
        lines( Q$traj[l:(l+1),1] , Q$traj[l:(l+1),2] , col=4 , lwd=5 )
        ani.record()
    }
    #points( Q$traj[2:L+1,] , pch=16 , col="white" , cex=0.3 )

    # draw new point
    pts[i,1:2] <- Q$traj[L+1,]
    pts[i,3] <- Q$accept

    #Arrows( Q$traj[L,1] , Q$traj[L,2] , Q$traj[L+1,1] , Q$traj[L+1,2] , arr.length=0.3 , arr.adj = 0.7 , col=4 )
    #text( Q$traj[L+1,1] , Q$traj[L+1,2] , i , cex=1.2 , pos=1 , offset=0.4 )
    
    points( Q$traj[L+1,1] , Q$traj[L+1,2] , pch=ifelse( Q$accept==1 , 1 , 16 ) , col=ifelse( Q$accept==1 , 4 , 2 ) , lwd=2 )

    invisible( replicate( 3 , ani.record() ) )

}

#points( 0 , log(sd(y)) , pch=3 ) # mark MAP point

# mtext( concat("Leapfrog steps = ",L) )

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 15 -loop 0 frame*.png HMC_2D_1.gif
# convert -delay 5 HMC_2D_1.gif HMC_2D_1.gif

############
# Metropolis to compare


METRO <- function( U , step , q_current , ... ) {
  d <- length(q_current)
  # generate proposals
  p <- q_current + rnorm(d,0,step)
  # compte U at current and proposals
  U_current <- -U( q_current , ... )
  U_proposed <- -U( p , ... )
  # accept?
  if ( exp(U_proposed - U_current) > runif(1) ) {
    q <- p
    accept <- 1
  } else {
    q <- q_current
    accept <- 0
  }
  return( list( q=q , p=p , accept=accept ) )
}

# test data
set.seed(7)
y <- abs(rnorm(50))
y <- c( y , -y ) # ensure mean is zero

# animation loop
ani.record(reset=TRUE)

set.seed(12)

# priors
priors <- list()
priors$a <- 0
priors$b <- 1
priors$k <- 0
priors$d <- 0.3

Q <- list()
Q$q <- c(-0.1,0.1) # start point
xr <- c(-0.4,0.4) # x range in plot
yr <- c(-0.25,0.3) # y range in plot

draw_bg()

n_samples <- 100
# points( Q$q[1] , Q$q[2] , pch=4 , col="black" )
pts <- matrix(NA,nrow=n_samples,ncol=3)

step <- 0.03

for ( i in 1:n_samples ) {

    #Q <- HMC2( myU2 , myU_grad2 , step , L , Q$q , a=priors$a , b=priors$b , k=priors$k , d=priors$d )
    prev_q <- Q$q
    Q <- METRO( myU2 , step , Q$q , a=priors$a , b=priors$b , k=priors$k , d=priors$d )

    draw_bg()

    # draw previous points
    if ( i > 1 ) {
        for ( j in 1:(i-1) ) {
            V <- 0.9
            if ( pts[j,3]==1 )
                points( pts[j,1] , pts[j,2] , pch=ifelse( pts[j,3]==1 , 1 , 16 ) , col=grau(V) , lwd=2 )
        }
    }

    # draw proposal
    #lines( Q$traj[l:(l+1),1] , Q$traj[l:(l+1),2] , col="white" , lwd=8 )
    lines( c(prev_q[1],Q$p[1]) , c(prev_q[2],Q$p[2]) , col=4 , lwd=5 )
    ani.record()
    #points( Q$traj[2:L+1,] , pch=16 , col="white" , cex=0.3 )

    # draw new point
    pts[i,1:2] <- Q$q
    pts[i,3] <- Q$accept

    #Arrows( Q$traj[L,1] , Q$traj[L,2] , Q$traj[L+1,1] , Q$traj[L+1,2] , arr.length=0.3 , arr.adj = 0.7 , col=4 )
    #text( Q$traj[L+1,1] , Q$traj[L+1,2] , i , cex=1.2 , pos=1 , offset=0.4 )
    
    points( Q$p[1] , Q$p[2] , pch=ifelse( Q$accept==1 , 1 , 16 ) , col=ifelse( Q$accept==1 , 4 , 2 ) , lwd=2 )

    invisible( replicate( 3 , ani.record() ) )

}

#points( 0 , log(sd(y)) , pch=3 ) # mark MAP point

# mtext( concat("Leapfrog steps = ",L) )

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 15 -loop 0 frame*.png Metro_2D.gif
# convert -delay 5 Metro_2D.gif Metro_2D.gif



#################
# now 2D Gaussian alpha,beta , can control correlation

# model:
# y ~ normal( mu , 1 )
# mu <- a1 + a2
# c(a1,a2) ~ normal(0,1)

myU3 <- function( q , a=0 , b=1 ) {
    mu <- q[1] + q[2]
    U <- sum( dnorm(y,mu,1,log=TRUE) ) + dnorm(q[1],a,b,log=TRUE) + dnorm(q[2],a,b,log=TRUE)
    return( -U )
}

# gradient function
# need vector of partial derivatives of U with respect to vector q
myU_grad3 <- function( q , a=0 , b=1 ) {
    mu <- q[1] + q[2]
    G1 <- sum( y - mu ) / 1 + (a - q[1])/b^2 #dU/da1
    G2 <- sum( y - mu ) / 1 + (a - q[2])/b^2 #dU/da2
    return( c( -G1 , -G2 ) ) # negative bc energy is neg-log-prob
}


# test data
set.seed(7)
y <- abs(rnorm(20))
y <- c( y , -y ) # ensure mean is zero

# what is posterior correlation?
m0 <- quap(
  alist(
    y ~ dnorm(mu,1),
    mu <- a1+a2,
    c(a1,a2) ~ dnorm(0,0.5)
  ), data=list(y=y) )
post <- extract.samples(m0)
cor( post )


# priors
priors <- list()
priors$a <- 0 # mean
priors$b <- 0.5 # sd

#ss <- ss + 1
set.seed(42) # seed 9 for examples

# init
Q <- list()
Q$q <- c(-0.4,0.2)
xr <- c(-1.5,1.5)
yr <- c(-1.5,1.5)

step <- 0.03
L <- 12 # 0.02/12 okay sampling --- 0.02/20 is good for showing u-turns
xpos <- c(4,2,1,2) # for L=20
#xpos <- c(2,3,2,1) # for L=55
path_col <- col.alpha("black",0.5)

draw_bg <- function() {
    plot( NULL , ylab="a[1]" , xlab="a[2]" , xlim=xr , ylim=yr , xaxt="n" , yaxt="n" )
    at <- c(-1.5,0,1.5)
    axis( 1 , at=at , labels=at )
    axis( 2 , at=at , labels=at )
    # draw contour of log-prob
    cb <- 0.2
    a1_seq <- seq(from=xr[1]-cb,to=xr[2]+cb,length.out=50) 
    a2_seq <- seq(from=yr[1]-cb,to=yr[2]+cb,length.out=50)
    z <- matrix(NA,length(a1_seq),length(a2_seq))
    for ( i in 1:length(a1_seq) )
        for ( j in 1:length(a2_seq) )
            z[i,j] <- myU3( c( a1_seq[i] , a2_seq[j] ) , a=priors$a , b=priors$b )
    cl <- contourLines( a1_seq , a2_seq , z , nlevels=30 )
    for ( i in 1:length(cl) ) lines( cl[[i]]$x , cl[[i]]$y , col=col.alpha("black",0.5) , lwd=1 )
}

step <- 0.15
L <- 15 # 0.02/12 okay sampling --- 0.02/20 is good for showing u-turns
xpos <- c(4,2,1,2) # for L=20
#xpos <- c(2,3,2,1) # for L=55
path_col <- col.alpha("black",0.5)

# animation loop
ani.record(reset=TRUE)

# set.seed(12)

Q <- list()
Q$q <- c(-0.4,-0.4) # start point
xr <- c(-1.5,1.5)
yr <- c(-1.5,1.5)

draw_bg()

n_samples <- 3
# points( Q$q[1] , Q$q[2] , pch=4 , col="black" )
pts <- matrix(NA,nrow=n_samples,ncol=3)

stop_traj_draw <- 20
for ( i in 1:n_samples ) {

    Q <- HMC2( myU3 , myU_grad3 , step , L , Q$q , a=priors$a , b=priors$b )

    draw_bg()

    # draw previous points
    if ( i > 1 ) {
        for ( j in 1:(i-1) ) {
            V <- 0.9
            points( pts[j,1] , pts[j,2] , pch=ifelse( pts[j,3]==1 , 1 , 16 ) , col=ifelse( pts[j,3]==1 , grau(V) , 2 ) , lwd=2 )
        }
    }

    # draw trajectory
    if ( i < stop_traj_draw )
    for ( l in 1:L ) {
        lines( Q$traj[l:(l+1),1] , Q$traj[l:(l+1),2] , col="white" , lwd=8 )
        lines( Q$traj[l:(l+1),1] , Q$traj[l:(l+1),2] , col=4 , lwd=5 )
        ani.record()
    }
    #points( Q$traj[2:L+1,] , pch=16 , col="white" , cex=0.3 )

    # draw new point
    pts[i,1:2] <- Q$traj[L+1,]
    pts[i,3] <- Q$accept

    #Arrows( Q$traj[L,1] , Q$traj[L,2] , Q$traj[L+1,1] , Q$traj[L+1,2] , arr.length=0.3 , arr.adj = 0.7 , col=4 )
    #text( Q$traj[L+1,1] , Q$traj[L+1,2] , i , cex=1.2 , pos=1 , offset=0.4 )
    
    points( Q$traj[L+1,1] , Q$traj[L+1,2] , pch=ifelse( Q$accept==1 , 1 , 16 ) , col=ifelse( Q$accept==1 , 4 , 2 ) , lwd=2 )

    ani.record()
    if ( i < stop_traj_draw )
        invisible( replicate( 2 , ani.record() ) )

}

#points( 0 , log(sd(y)) , pch=3 ) # mark MAP point

# mtext( concat("Leapfrog steps = ",L) )

oopts = ani.options(interval = 0.15)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 15 -loop 0 frame*.png HMC_div.gif
# convert -delay 5 HMC_2D_1.gif HMC_2D_1.gif

#########
# new workflow examples

library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

dat <- list(
    D = standardize(d$Divorce),
    M = standardize(d$Marriage),
    A = standardize(d$MedianAgeMarriage)
)

f <- alist(
        D ~ dnorm(mu,sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0,0.2),
        bM ~ dnorm(0,0.5),
        bA ~ dnorm(0,0.5),
        sigma ~ dexp(1)
    )

mq <- quap( f , data=dat )

library(cmdstanr)
mHMC <- ulam( f , data=dat )

mHMC <- ulam( f , data=dat , chains=4 , cores=4 )

blank(bty="n",ex=2)

traceplot(mHMC,n_cols=2,lwd=1)

trankplot(mHMC,n_cols=2,lwd=2)

stancode(mHMC)

mHMC_stan <- cstan( file="08_mHMC.stan" , data=dat , chains=4 , cores=4 )

post <- extract.samples(mHMC_stan)

# illustrate R-hat

post <- extract(mHMC@stanfit, permuted = FALSE, inc_warmup = TRUE)

fRhat <- function(t) {
    W <- sapply( 1:4 , function(m) var( post[1:t,m,1] ) )
    chain_means <- sapply( 1:4 , function(m) mean( post[1:t,m,1] ) )
    B <- var(chain_means)
    return(c(mean(W),B))
}

r <- sapply( 1:1000 , fRhat )
plot( r[1,] , col=4 , lwd=4 , type="l" , xlab="sample" , ylab="variance" )
points( 1:1000 , r[2,] , col=2 , lwd=4 , type="l" )

# autocorrelation

x <- acf( post[,1,1] , lwd=7 , col=2 )

# bad chains

y <- c(-1,1)
set.seed(11)
m9.2 <- ulam(
    alist(
        y ~ dnorm( mu , sigma ),
        mu <- alpha,
        alpha ~ dnorm( 0 , 1000 ),
        sigma ~ dexp( 0.0001 )
    ) , data=list(y=y) , chains=3 , cores=3 )

traceplot(m9.2,n_cols=2,lwd=1)

trankplot(m9.2,n_cols=2,lwd=2)

m9.3 <- ulam(
    alist(
        y ~ dnorm( mu , sigma ),
        mu <- alpha,
        alpha ~ dnorm( 1 , 10 ),
        sigma ~ dexp( 1 )
    ) , data=list(y=y) , chains=3 , cores=3 )

traceplot(m9.3,n_cols=2,lwd=1)

trankplot(m9.3,n_cols=2,lwd=2)
