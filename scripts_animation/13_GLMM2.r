# lecture 12 - intro to multilevel models
# random intercepts and partial pooling

library(rethinking)
library(animation)
library(ellipse)

####################################
# UCBadmit as varying effect model

data(UCBadmit)
d <- UCBadmit

dat <- list(
    A = d$admit,
    N = d$applications,
    D = as.integer(d$dept),
    G = ifelse(d$applicant.gender=="female",1,2)
)

m0 <- ulam(
    alist(
        A ~ binomial(N,p),
        logit(p) <- a[G,D],
        matrix[G,D]:a ~ normal( 0 , 1.5 )
    ) , data=dat , chains=4 )

# pool all estimates, ignoring department structure
m1 <- ulam(
    alist(
        A ~ binomial(N,p),
        logit(p) <- a[G,D],
        matrix[G,D]:a ~ normal( a_bar , sigma ),
        a_bar ~ normal(0,1),
        sigma ~ exponential(1)
    ) , data=dat , chains=4 )

# pool while noticing department structure
# this is a "random slopes" model, without correlation
dat$G1 <- ifelse(dat$G==1,1,0)
m2 <- ulam(
    alist(
        A ~ binomial(N,p),
        logit(p) <- a[D] + b[D]*G1,
        a[D] ~ normal( a_bar , sigma ),
        b[D] ~ normal( 0 , tau ),
        a_bar ~ normal(0,1),
        sigma ~ exponential(1),
        tau ~ exponential(1)
    ) , data=dat , chains=4 )

####################################
# chimpanzees

data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition
dat <- list(
    P = d$pulled_left,
    A = d$actor,
    B = d$block,
    T = d$treatment ) 

# block interactions
mBT <- ulam(
    alist(
        P ~ bernoulli( p ) ,
        logit(p) <- b[T,B] + a[A],
      ## adaptive priors
        matrix[T,B]:b ~ dnorm( 0 , sigma_B ),
        a[A] ~ dnorm( a_bar , sigma_A ),
      ## hyper-priors
        a_bar ~ dnorm( 0 , 1.5 ),
        sigma_A ~ dexp(1),
        sigma_B ~ dexp(1)
    ) , data=dat , chains=4 , cores=4 )

mBTnc <- ulam(
    alist(
        P ~ bernoulli( p ) ,
        logit(p) <- a_bar + z_a[A]*sigma_A + z_b[T,B]*sigma_B ,
      ## adaptive priors
        matrix[T,B]:z_b ~ dnorm( 0 , 1 ),
        z_a[A] ~ dnorm( 0 , 1 ),
      ## hyper-priors
        a_bar ~ dnorm( 0 , 1.5 ),
        sigma_A ~ dexp(1),
        sigma_B ~ dexp(1),
        gq> vector[A]:a <<- a_bar + z_a*sigma_A,
        gq> matrix[T,B]:b <<- z_b*sigma_B
    ) , data=dat , chains=4 , cores=4 )

# plot treatments across blocks
post <- extract.samples(mBTnc)
plot( NULL , xlim=c(0.5,4.5) , ylim=c(-2,2) , xlab="treatment" , ylab="log-odds" , xaxt="n" )
axis( 1 , at=1:4 , labels=c("R/N","L/N","R/P","L/P") )
abline( h=0 , lty=3 , lwd=2 )
for( i in 1:6 ) points( jitter(1:4) , apply(post$b[,,i],2,mean) , lwd=4 , col=2 )

# plot actor effects
plot( NULL , xlim=c(0.5,7.5) , ylim=c(0,1) , xlab="actor" , ylab="probability pull left" , xaxt="n" )
axis( 1 , at=1:7 , labels=1:7 )
abline( h=0.5 , lty=3 , lwd=2 )
for ( i in 1:7 ) lines( c(i,i) , PI(inv_logit(post$a[,i])) , lwd=8 , col=col.alpha(4,0.5) )
points( 1:7 , apply(inv_logit(post$a),2,mean) , lwd=4 , col=4 )


blank(bty="n",h=0.5)

dens( post$sigma_a , xlab="sigma_A" , lwd=4 , col=4 , xlim=c(0,5) )

dens( post$sigma_b , xlab="sigma_B" , lwd=4 , col=2 , xlim=c(0,5) )

# show n_eff

the_pars <- c("a_bar","sigma_A","sigma_B","b","a")
precis_c <- precis( mBT , depth=3 , pars=the_pars )
precis_nc <- precis( mBTnc , depth=3 , pars=the_pars )
neff_table <- cbind( precis_c[,"n_eff"] , precis_nc[,"n_eff"] )
cols <- ifelse( neff_table[,1] > neff_table[,2] , 2 , 4 )
plot( neff_table , xlim=range(neff_table) , ylim=range(neff_table) ,
    xlab="effective sample size (centered)" , ylab="effective sample size (non-centered)" , lwd=3 , col=cols )
abline( a=0 , b=1 , lty=2 )

trankplot(mBTnc,pars=c("sigma_A","sigma_B"),n_cols=1)




####################################
# trolley model with S and U
data(Trolley)
d <- Trolley
dat <- list(
    R = d$response,
    A = d$action,
    I = d$intention,
    C = d$contact
)
# index variables for story and id
dat$S <- as.integer( d$story )
dat$U <- as.integer( d$id )

mRXSU <- ulam(
    alist(
        R ~ ordered_logistic( phi , alpha ),
        phi <- bA*A + bI*I + bC*C +
               as[S] + au[U],
        alpha ~ normal( 0 , 1 ),
        bA ~ normal( 0 , 0.5 ),
        bI ~ normal( 0 , 0.5 ),
        bC ~ normal( 0 , 0.5 ),
        # varying intercepts
        as[S] ~ normal(0,sigmaS),
        au[U] ~ normal(0,sigmaU),
        sigmaS ~ exponential(1),
        sigmaU ~ exponential(1)
    ), data=dat , chains=4 , cores=4 )

mRXSUnc <- ulam(
    alist(
        R ~ ordered_logistic( phi , alpha ),
        phi <- bA*A + bI*I + bC*C +
               zs[S]*sigmaS + au[U],
        alpha ~ normal( 0 , 1 ),
        bA ~ normal( 0 , 0.5 ),
        bI ~ normal( 0 , 0.5 ),
        bC ~ normal( 0 , 0.5 ),
        # varying intercepts
        zs[S] ~ normal(0,1),
        au[U] ~ normal(0,sigmaU),
        sigmaS ~ half_normal(0,1),
        sigmaU ~ half_normal(0,1)
    ), data=dat , chains=4 , cores=4 , threads=2 )

####################################
# trolley model mit allem
data(Trolley)
d <- Trolley
dat <- list(
    R = d$response,
    A = d$action,
    I = d$intention,
    C = d$contact
)
dat$G <- ifelse(d$male==1,2L,1L)
edu_levels <- c( 6 , 1 , 8 , 4 , 7 , 2 , 5 , 3 )
edu_new <- edu_levels[ d$edu ]
dat$E <- edu_new
dat$a <- rep(2,7) # dirichlet prior
dat$Y <- standardize(d$age)
#dat$Y <- d$age
dat$G1 <- ifelse(dat$G==1,1,0)
dat$G2 <- ifelse(dat$G==2,1,0)

dat$S <- as.integer( d$story )
dat$U <- as.integer( d$id )

mRXma <- ulam(
    alist(
        R ~ ordered_logistic( phi , alpha ),
        phi <- G1*bE[G]*sum( deltaF_j[1:E] ) + 
               G2*bE[G]*sum( deltaM_j[1:E] ) + 
               bA[G]*A + bI[G]*I + bC[G]*C +
               bY[G]*Y +
               as[S] + au[U],
        alpha ~ normal( 0 , 1 ),
        bA[G] ~ normal( 0 , 0.5 ),
        bI[G] ~ normal( 0 , 0.5 ),
        bC[G] ~ normal( 0 , 0.5 ),
        bE[G] ~ normal( 0 , 0.5 ),
        bY[G] ~ normal( 0 , 0.5 ),
        vector[8]: deltaF_j <<- append_row( 0 , deltaF ),
        vector[8]: deltaM_j <<- append_row( 0 , deltaM ),
        simplex[7]: deltaF ~ dirichlet( a ),
        simplex[7]: deltaM ~ dirichlet( a ),
        # varying intercepts
        as[S] ~ normal(0,sigmaS),
        au[U] ~ normal(0,sigmaU),
        sigmaS ~ exponential(1),
        sigmaU ~ exponential(1)
    ), data=dat , chains=4 , cores=4 )

precis(mRXma)

# non-center
mRXma2 <- ulam(
    alist(
        R ~ ordered_logistic( phi , alpha ),
        phi <- G1*bE[G]*sum( deltaF_j[1:E] ) + 
               G2*bE[G]*sum( deltaM_j[1:E] ) + 
               bA[G]*A + bI[G]*I + bC[G]*C +
               bY[G]*Y +
               zs[S]*sigmaS + zu[U]*sigmaU,
        alpha ~ normal( 0 , 1 ),
        bA[G] ~ normal( 0 , 0.5 ),
        bI[G] ~ normal( 0 , 0.5 ),
        bC[G] ~ normal( 0 , 0.5 ),
        bE[G] ~ normal( 0 , 0.5 ),
        bY[G] ~ normal( 0 , 0.5 ),
        vector[8]: deltaF_j <<- append_row( 0 , deltaF ),
        vector[8]: deltaM_j <<- append_row( 0 , deltaM ),
        simplex[7]: deltaF ~ dirichlet( a ),
        simplex[7]: deltaM ~ dirichlet( a ),
        # varying intercepts
        zs[S] ~ normal(0,1),
        zu[U] ~ normal(0,1),
        sigmaS ~ half_normal(0,1),
        sigmaU ~ half_normal(0,1)
    ), data=dat , chains=4 , cores=4 , threads=2 )

precis(mRXma2,2)

######
# prediction example

library(rethinking)
data(reedfrogs)
d <- reedfrogs

dat <- list(
    S = d$surv,
    D = d$density,
    T = 1:nrow(d),
    P = ifelse( d$pred=="no" , 1L , 2L ),
    G = ifelse( d$size=="small" , 1L , 2L )
)

mSPG <- ulam(
    alist(
        S ~ binomial( D , p ),
        logit(p) <- a[T] + b[P,G],
        a[T] ~ normal( 0 , sigma ),
        matrix[P,G]:b ~ normal( 0 , 1 ),
        sigma ~ exponential( 1 )
    ), data=dat , chains=4 , cores=4 )

# simulate causal effect of changing size distribution
# target population high density 35, 50% groups have predation
# 25% large. intervention makes 75% large. what is change in survival?

post <- extract.samples(mSPG)
# sim first under status quo
n_groups <- 1000
n_samples <- 2000
S1 <- matrix(0,nrow=n_samples,ncol=n_groups)
for ( s in 1:n_groups ) {
    # sim a tank from posterior population
    aT <- rnorm(n_samples,0,post$sigma)
    # sample P and G for this group
    P <- sample( 1:2 , size=1 , prob=c(0.5,0.5) ) # 50% pred
    G <- sample( 1:2 , size=1 , prob=c(0.75,0.25) ) # 25% large
    # sim survival
    p <- inv_logit( aT + post$b[,P,G] )
    S1[,s] <- rbinom(n_samples,35,p)
}

# intervention - 50% large
S2 <- matrix(0,nrow=n_samples,ncol=n_groups)
for ( s in 1:n_groups ) {
    # sim a tank from posterior population
    aT <- rnorm(n_samples,0,post$sigma)
    # sample P and G for this group
    P <- sample( 1:2 , size=1 , prob=c(0.5,0.5) ) # 50% pred
    G <- sample( 1:2 , size=1 , prob=c(0.25,0.75) ) # 75% large
    # sim survival
    p <- inv_logit( aT + post$b[,P,G] )
    S2[,s] <- rbinom(n_samples,35,p)
}

simplehist( as.numeric(S1) , lwd=4 , col=2 , xlab="number surviving" )
simplehist( as.numeric(S2) , lwd=4 , col=2 , xlab="number surviving" )

x <- table( as.numeric( S2 - S1 ) )
x <- x/sum(x)

plot( NULL , xlim=c(-35,35) , ylim=c(0,max(x)) , xlab="change in survival" , ylab="proportion of simulation" )
for ( i in 1:length(x) ) {
    xvals <- as.numeric(names(x))
    lines( rep(xvals[i],2) , c(0,x[i]) , lwd=4 , col=ifelse(xvals[i]<0,2,4) )
}

abline(v=mean(as.numeric( S2 - S1 )),lty=3,lwd=3)

##################################################
# devil's funnel

m13.7 <- ulam(
    alist(
        v ~ normal(0,3),
        x ~ normal(0,exp(v))
    ), data=list(N=1) , chains=4 )
precis( m13.7 )

traceplot(m13.7,n_cols=1,lwd=3)

m13.7nc <- ulam(
    alist(
        v ~ normal(0,3),
        z ~ normal(0,1),
        gq> real[1]:x <<- z*exp(v)
    ), data=list(N=1) , chains=4 )
precis( m13.7nc )

post <- extract.samples(m13.7nc)
dens(post$x,post$v)

traceplot(m13.7nc,n_cols=1,lwd=3)

# log-prob and gradient functions for animation

U_funnel <- function( q , s=3 ) {
    v <- q[2]
    x <- q[1]
    U <- sum( dnorm(x,0,exp(v),log=TRUE) ) + dnorm(v,0,s,log=TRUE)
    return( -U )
}

U_funnel_gradient <- function( q , s=3 ) {
    v <- q[2]
    x <- q[1]
    Gv <- (-v)/s^2 - length(x) + exp(-2*v)*sum( x^2 ) #dU/dv
    Gx <- -exp(-2*v)*x #dU/dx
    return( c( -Gx , -Gv ) ) # negative bc energy is neg-log-prob
}


# priors
priors <- list()
priors$s <- 3

#ss <- ss + 1
set.seed(42) # seed 9 for examples

# init
n_samples <- 4
Q <- list()
Q$q <- c(-0.4,0.2)
xr <- c(-5,5)
yr <- c(-4,4)

step <- 0.02
L <- 12 # 0.02/12 okay sampling --- 0.02/20 is good for showing u-turns
xpos <- c(4,2,1,2) # for L=20
#xpos <- c(2,3,2,1) # for L=55
path_col <- col.alpha("black",0.5)

draw_bg <- function() {
    plot( NULL , ylab="v" , xlab="x" , xlim=xr , ylim=yr )
    # draw contour of log-prob
    cb <- 1
    xlen <- 200
    xseq <- seq(from=xr[1]-cb,to=xr[2]+cb,length.out=xlen) 
    vseq <- seq(from=yr[1]-cb,to=yr[2]+cb,length.out=xlen)
    z <- matrix(NA,length(xseq),length(vseq))
    for ( i in 1:length(xseq) )
        for ( j in 1:length(vseq) )
            z[i,j] <- U_funnel( c( xseq[i] , vseq[j] ) , s=priors$s )
    
    #cl <- contourLines( xseq , vseq , z , nlevels=100 )

    lz5 <- log(z)
    lvls <- exp( pretty( range(lz5), 30 ) )
    cl <- contourLines( xseq , vseq , z , level=lvls/12 )

    depth <- 2
    for ( i in 1:length(cl) ) lines( cl[[i]]$x , cl[[i]]$y , col=col.alpha("black",0.5) , lwd=1 + depth/cl[[i]]$level )
}

draw_bg()

########################################
# animate transition through priors

ani.record(reset=TRUE)
s_seq <- seq(from=0.5,to=3,len=30)
for ( s in s_seq ) {
    priors$s <- s
    draw_bg()
    mtext(concat("v ~ normal(0,",round(s,2),")"), cex=1.6)
    ani.record()
}

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 15 -loop 0 frame*.png HMC_devil_prior_smooth.gif
# convert -delay 5 HMC_devil_funnel.gif HMC_devil_funnel.gif

########################################
# below is HMC loop

xpos <- c(4,2,1,2) # for L=20
#xpos <- c(2,3,2,1) # for L=55
path_col <- col.alpha("black",0.5)

# animation loop
set.seed(11)
Q <- list()
Q$q <- c(-0.4,0.2) # start point
xr <- c(-5,5)
yr <- c(-4,4)

step <- 0.025
L <- 30
priors$s <- 3

draw_bg()

n_samples <- 35
# points( Q$q[1] , Q$q[2] , pch=4 , col="black" )
pts <- matrix(NA,nrow=n_samples,ncol=3)
ani.record(reset=TRUE)
for ( i in 1:n_samples ) {

    Q <- HMC2( U_funnel , U_funnel_gradient , step , L , Q$q , s=priors$s )

    draw_bg()

    # draw previous points
    if ( i > 1 ) {
        for ( j in 1:(i-1) ) {
            V <- 0.9
            points( pts[j,1] , pts[j,2] , pch=ifelse( pts[j,3]==1 , 1 , 16 ) , col=ifelse( pts[j,3]==1 , 1 , 2 ) , lwd=2 )
        }
    }

    # draw trajectory
    for ( l in 1:L ) {
        if ( step > 0.05 )
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

# points( 0 , log(sd(y)) , pch=3 ) # mark MAP point

# mtext( concat("Leapfrog steps = ",L) )

oopts = ani.options(interval = 0.01)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 5 -loop 0 frame*.png HMC_devil_funnel.gif
# convert -delay 5 HMC_devil_funnel.gif HMC_devil_funnel.gif

################################################
# non-centered log-prob and gradient
# v ~ normal(0,3)
# z ~ normal(0,1)
# x = z*exp(v)

U_funnel_NC <- function( q , s=3 ) {
    v <- q[2]
    z <- q[1]
    U <- sum( dnorm(z,0,1,log=TRUE) ) + dnorm(v,0,s,log=TRUE)
    return( -U )
}

U_funnel_NC_gradient <- function( q , s=3 ) {
    v <- q[2]
    z <- q[1]
    Gv <- (-v)/s^2 #dU/dv
    Gz <- (-z) #dU/dz
    return( c( -Gz , -Gv ) ) # negative bc energy is neg-log-prob
}

draw_bg <- function() {
    plot( NULL , ylab="v" , xlab="z" , xlim=xr , ylim=yr )
    # draw contour of log-prob
    cb <- 1
    xlen <- 200
    xseq <- seq(from=xr[1]-cb,to=xr[2]+cb,length.out=xlen) 
    vseq <- seq(from=yr[1]-cb,to=yr[2]+cb,length.out=xlen)
    z <- matrix(NA,length(xseq),length(vseq))
    for ( i in 1:length(xseq) )
        for ( j in 1:length(vseq) )
            z[i,j] <- U_funnel_NC( c( xseq[i] , vseq[j] ) , s=priors$s )
    
    cl <- contourLines( xseq , vseq , z , nlevels=19 )

    #lz5 <- log(z)
    #lvls <- exp( pretty( range(lz5), 30 ) )
    #cl <- contourLines( xseq , vseq , z , level=lvls/12 )

    depth <- 1
    for ( i in 1:length(cl) ) lines( cl[[i]]$x , cl[[i]]$y , col=col.alpha("black",0.5) , lwd=1 + depth/cl[[i]]$level )
}

draw_bg()

# animation loop
set.seed(12)
Q <- list()
xr <- c(-4,4)
yr <- c(-7,7)

n_chains <- 4
for ( k in 1:n_chains ) Q[[k]] <- list()
for ( k in 1:n_chains )
    Q[[k]]$q <- c(runif(1,xr[1],xr[2]),runif(1,yr[1],yr[2])) # start point

chain_cols <- c(4,6,3,5,2)

step <- 0.3
L <- 12
priors$s <- 3

draw_bg()

n_samples <- 20
pts <- list()
for ( k in 1:n_chains ) pts[[k]] <- matrix(NA,nrow=n_samples,ncol=3)

ani.record(reset=TRUE)
for ( i in 1:n_samples ) {

    for ( j in 1:n_chains )
        Q[[j]] <- HMC2( U_funnel_NC , U_funnel_NC_gradient , step , L , Q[[j]]$q , s=priors$s )

    draw_bg()

    # draw previous points
    for ( k in 1:n_chains ) {
        if ( i > 1 ) {
            for ( j in 1:(i-1) ) {
                points( pts[[k]][j,1] , pts[[k]][j,2] , pch=ifelse( pts[[k]][j,3]==1 , 1 , 16 ) , col=ifelse( pts[[k]][j,3]==1 , chain_cols[k] , 2 ) , lwd=2 )
            }#j
        }
    }#k

    # draw trajectories
    for ( l in 1:L ) {
        for ( k in 1:n_chains ) {
            if ( step > 0.05 )
                lines( Q[[k]]$traj[l:(l+1),1] , Q[[k]]$traj[l:(l+1),2] , col="white" , lwd=8 )
            lines( Q[[k]]$traj[l:(l+1),1] , Q[[k]]$traj[l:(l+1),2] , col=chain_cols[k] , lwd=5 )
            ani.record()
        }#k
    }#l
    #points( Q$traj[2:L+1,] , pch=16 , col="white" , cex=0.3 )

    for ( k in 1:n_chains ) {
        # draw new points
        pts[[k]][i,1:2] <- Q[[k]]$traj[L+1,]
        pts[[k]][i,3] <- Q[[k]]$accept

        points( Q[[k]]$traj[L+1,1] , Q[[k]]$traj[L+1,2] , pch=ifelse( Q[[k]]$accept==1 , 1 , 16 ) , col=ifelse( Q[[k]]$accept==1 , chain_cols[k] , 2 ) , lwd=2 )
    }#k

    invisible( replicate( 3 , ani.record() ) )

}

# points( 0 , log(sd(y)) , pch=3 ) # mark MAP point

# mtext( concat("Leapfrog steps = ",L) )

oopts = ani.options(interval = 0.01)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 5 -loop 0 frame*.png HMC_devil_funnel.gif
# convert -delay 5 HMC_devil_funnel.gif HMC_devil_funnel.gif

# convert HMC_devil_funnel_z.gif -coalesce -fuzz 2% +dither -layers Optimize +map output.gif
