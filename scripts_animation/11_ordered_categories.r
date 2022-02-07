# lecture 11 - ordered categories

library(rethinking)
library(animation)
library(ellipse)

data(Trolley)
d <- Trolley

###################
# ordered logit link construction

blank(bty="n")

y <- table( d$response )

# normal freq plot
plot( NULL , xlim=c(0.5,7.5) , ylim=c(0,max(y)+1000) , xlab="outcome" , ylab="frequency" )
for ( i in 1:7 ) lines( c(i,i) , c(0,y[i]) , lwd=8 , col=2 )

# cumulative plot
ycum <- cumsum(y)

ani.record(reset=TRUE)
nf <- 30
p <- seq(from=0,to=1,len=nf)
for ( f in 1:nf ) {

    ymax <- max(max(y),p[f]*max(ycum))+1000
    plot( NULL , xlim=c(0.5,7.5) , ylim=c(0,ymax) , xlab="outcome" , ylab="cumulative frequency" )

    lines( c(1,1) , c(0,ycum[1]) , lwd=8 , col=2 )

    for ( i in 2:7 ) {
        xboost <- p[f]*ycum[i-1]
        lines( c(i,i) , c(0,y[i])+xboost , lwd=8 , col=2 )
    }

    ani.record()
}

# black supports
for ( i in 1:7 ) lines( c(i,i) , c(0,ycum[i]) , lwd=6 )
lines( c(1,1) , c(0,ycum[1]) , lwd=8 , col=2 )
for ( i in 2:7 ) lines( c(i,i) , c(ycum[i-1],ycum[i]) , lwd=8 , col=2 )

ani.record()

# segments
for ( i in 1:6 ) lines( c(i,i+1) , c(ycum[i],ycum[i]) , lty=3 , lwd=2 , col=2 )

ani.record()

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 5 -loop 0 frame*.png a_out.gif
# convert -delay 3 a_out.gif a_out.gif

###################
# show cumulative proportions
p <- y/sum(y)
pcum <- ycum/max(ycum)

ani.record(reset=TRUE)
for ( l in 1:7 ) {

    plot( NULL , xlim=c(0.5,7.5) , ylim=c(0,1) , xlab="outcome" , ylab="cumulative proportion" )

    for ( i in 1:7 ) lines( c(i,i) , c(0,pcum[i]) , lwd=6 )
    lines( c(1,1) , c(0,pcum[1]) , lwd=8 , col=2 )
    for ( i in 2:7 ) {
        xboost <- pcum[i-1]
        lines( c(i,i) , c(0,p[i])+xboost , lwd=8 , col=2 )
    }

    # lines
    for ( i in 1:l )
        lines( c(-2,i) , c(pcum[i],pcum[i]) , lwd=2 , lty=3 , col=2 )

    ani.record()

}

oopts = ani.options(interval = 0.2)
ani.replay()

###################
# cumulative proportions to cumulative log-odds

nsteps <- 30
phi_seq <- seq( from=0 , to=2 , len=nsteps )
phi_seq <- c(phi_seq,rep(2,10))
phi_seq <- c(phi_seq, seq( from=2 , to=-2 , len=nsteps*2 ) )
phi_seq <- c(phi_seq,rep(-2,10))
phi_seq <- c(phi_seq, seq( from=-2 , to=0 , len=nsteps ) )

blank(bty="n",w=2.5)
par(mfrow=c(1,2))

ani.record(reset=TRUE)
par(mfrow=c(1,2))
for ( phi in phi_seq) {

    # left: cumulative prob vs cumulative log odds

    lpc <- logit(pcum) + phi
    ilpc <- inv_logit(lpc)

    xmax <- max(lpc[1:6]) + 1.5
    xmin <- min(lpc[1:6]) - 1
    xmax <- 6
    xmin <- -5
    plot( lpc[1:6] , ilpc[1:6] , lwd=4 , col=4 , pch=1 , xlab="cumulative log-odds" , ylab="cumulative proportion" , xlim=c(xmin,xmax) , ylim=c(0,1) , xaxt="n" )
    points( xmax , 1 , lwd=4 , col=grau() )

    at <- c(-4,-2,0,2,4)
    axis(1,at=c(at,xmax),labels=c(at,"Inf"))
    
    for ( i in 1:7 ) {
        lines( c(xmin-1,lpc[i]) , c(ilpc[i],ilpc[i]) , lwd=2 , lty=3 , col=2 )
        lines( c(lpc[i],lpc[i]) , c(0,ilpc[i]) , lwd=2 , lty=3 , col=4 )
    }
    # 7
    lines( c(xmin-1,xmax) , c(1,1) , lwd=2, lty=3 , col=2 )

    # point on x-axis for value of phi
    par( xpd=TRUE )
    points( phi , -0.04 , cex=2 , col=1 , pch=16 )
    par( xpd=FALSE )

    # right: implied histogram

    p <- sapply( 2:7 , function(i) ilpc[i] - ilpc[i-1] )
    p <- c( ilpc[1] , p )
    plot( NULL , xlim=c(0.5,7.5) , ylim=c(0,0.6) , xlab="observed value" , ylab="probability" )
    for ( i in 1:7 ) {
        lines( c(i,i) , c(0,p[i]) , lwd=12 , col=col.alpha(2,0.4) )
        lines( c(i,i) , c(0,p[i]) , lwd=8 , col=2 )
    }

    # record frame
    ani.record()

}

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 5 -loop 0 frame*.png a_out.gif
# convert -delay 3 a_out.gif a_out.gif

########################
# now do some modeling

data(Trolley)
d <- Trolley
dat <- list(
    R = d$response,
    A = d$action,
    I = d$intention,
    C = d$contact
)

mRX <- ulam(
    alist(
        R ~ dordlogit(phi,alpha),
        phi <- bA*A + bI*I + bC*C,
        c(bA,bI,bC) ~ normal(0,0.5),
        alpha ~ normal(0,1)
    ) , data=dat , chains=4 , cores=4 )

precis(mRX,2)

# plot predictive distributions for each treatment

vals <- c(0,1,1) # A,I,C
Rsim <- mcreplicate( 100 , sim(mRX,data=list(A=vals[1],I=vals[2],C=vals[3])) , mc.cores=6 )
simplehist(as.vector(Rsim),lwd=8,col=2,xlab="Response")
mtext(concat("A=",vals[1],", I=",vals[2],", C=",vals[3]))

# total effect of gender
dat$G <- ifelse(d$male==1,2,1)
mRXG <- ulam(
    alist(
        R ~ dordlogit(phi,alpha),
        phi <- bA[G]*A + bI[G]*I + bC[G]*C,
        bA[G] ~ normal(0,0.5),
        bI[G] ~ normal(0,0.5),
        bC[G] ~ normal(0,0.5),
        alpha ~ normal(0,1)
    ) , data=dat , chains=4 , cores=4 )

precis(mRXG,2)

vals <- c(0,1,1,2) # A,I,C,G
Rsim <- mcreplicate( 100 , sim(mRXG,data=list(A=vals[1],I=vals[2],C=vals[3],G=vals[4])) , mc.cores=6 )
simplehist(as.vector(Rsim),lwd=8,col=2,xlab="Response")
mtext(concat("A=",vals[1],", I=",vals[2],", C=",vals[3],", G=",vals[4]))

# distributions of education and age

edu_levels <- c( 6 , 1 , 8 , 4 , 7 , 2 , 5 , 3 )
edu_new <- edu_levels[ d$edu ]

simplehist(edu_new,xlab="education level (ordered)",lwd=8,col=2)

simplehist(d$age,xlab="age (years)",lwd=6,col=4)

# education model

library(gtools) # provides rdirichlet

nf <- 40
a <- 10
delta <- rdirichlet( nf , a=rep(a,7) )
delta <- rdirichlet( nf , a=1:7 )
ani.record(reset=TRUE)
for ( f in 1:nf) {

    plot( NULL , xlim=c(1,7) , ylim=c(0,0.4) , xlab="index" , ylab="probability" )

    if ( f > 1 ) {
        start <- max(f-3,1)
        for ( i in start:(f-1) )
            lines( 1:7 , delta[i,] , type="l" , lwd=4 , col=grau(0.25*i/f) )
    }
    
    lines( 1:7 , delta[f,] , type="b" , lwd=4 , col=2 )

    ani.record()

}

oopts = ani.options(interval = 0.3)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 20 -loop 0 frame*.png a_out.gif
# convert -delay 3 a_out.gif a_out.gif

edu_levels <- c( 6 , 1 , 8 , 4 , 7 , 2 , 5 , 3 )
edu_new <- edu_levels[ d$edu ]

dat$E <- edu_new
dat$a <- rep(2,7) # dirichlet prior

mRXE <- ulam(
    alist(
        R ~ ordered_logistic( phi , alpha ),
        phi <- bE*sum( delta_j[1:E] ) + bA*A + bI*I + bC*C,
        alpha ~ normal( 0 , 1 ),
        c(bA,bI,bC,bE) ~ normal( 0 , 0.5 ),
        vector[8]: delta_j <<- append_row( 0 , delta ),
        simplex[7]: delta ~ dirichlet( a )
    ), data=dat , chains=4 , cores=4 )

precis(mRXE,2)



# BIG MODEL

dat$Y <- standardize(d$age)

# single-threaded version
mRXEYG <- ulam(
    alist(
        R ~ ordered_logistic( phi , alpha ),
        phi <- bE[G]*sum( delta_j[1:E] ) + 
               bA[G]*A + bI[G]*I + bC[G]*C +
               bY[G]*Y,
        alpha ~ normal( 0 , 1 ),
        bA[G] ~ normal( 0 , 0.5 ),
        bI[G] ~ normal( 0 , 0.5 ),
        bC[G] ~ normal( 0 , 0.5 ),
        bE[G] ~ normal( 0 , 0.5 ),
        bY[G] ~ normal( 0 , 0.5 ),
        vector[8]: delta_j <<- append_row( 0 , delta ),
        simplex[7]: delta ~ dirichlet( a )
    ), data=dat , chains=4 , cores=4 )

# multi-threaded version
mRXEYGt <- ulam(
    alist(
        R ~ ordered_logistic( phi , alpha ),
        phi <- bE[G]*sum( delta_j[1:E] ) + 
               bA[G]*A + bI[G]*I + bC[G]*C +
               bY[G]*Y,
        alpha ~ normal( 0 , 1 ),
        bA[G] ~ normal( 0 , 0.5 ),
        bI[G] ~ normal( 0 , 0.5 ),
        bC[G] ~ normal( 0 , 0.5 ),
        bE[G] ~ normal( 0 , 0.5 ),
        bY[G] ~ normal( 0 , 0.5 ),
        vector[8]: delta_j <<- append_row( 0 , delta ),
        simplex[7]: delta ~ dirichlet( a )
    ), data=dat , chains=4 , cores=4 , threads=2 )

precis(mRXEYGt,2)

# version with different delta for each G
dat$G1 <- ifelse(dat$G==1,1,0)
dat$G2 <- ifelse(dat$G==2,1,0)
mRXEYG2t <- ulam(
    alist(
        R ~ ordered_logistic( phi , alpha ),
        phi <- G1*bE[G]*sum( deltaF_j[1:E] ) + 
               G2*bE[G]*sum( deltaM_j[1:E] ) + 
               bA[G]*A + bI[G]*I + bC[G]*C +
               bY[G]*Y,
        alpha ~ normal( 0 , 1 ),
        bA[G] ~ normal( 0 , 0.5 ),
        bI[G] ~ normal( 0 , 0.5 ),
        bC[G] ~ normal( 0 , 0.5 ),
        bE[G] ~ normal( 0 , 0.5 ),
        bY[G] ~ normal( 0 , 0.5 ),
        vector[8]: deltaF_j <<- append_row( 0 , deltaF ),
        vector[8]: deltaM_j <<- append_row( 0 , deltaM ),
        simplex[7]: deltaF ~ dirichlet( a ),
        simplex[7]: deltaM ~ dirichlet( a )
    ), data=dat , chains=4 , cores=4 , threads=2 )

precis(mRXEYG2t,2)


