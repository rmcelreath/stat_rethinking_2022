# animate prior predictive simulation

library(rethinking)
library(animation)
library(ellipse)

# spline basis functions animation

data(cherry_blossoms)
d <- cherry_blossoms

library(splines)

blank(w=2,h=0.8,bty="n")

# spline on blossom doy
d3 <- d[ complete.cases(d$doy) , ] # complete cases on doy
idx <- 1:nrow(d3)
dd <- d3[ idx , ]
dd <- dd[ order(dd$year) , ]

num_knots <- 3
knot_list <- seq( from=min(dd$year) , to=max(dd$year) , length.out=num_knots )
B3 <- t(bs(dd$year, knots=knot_list , degree=2, intercept = FALSE))

w <- rep(c(1,-1),length=nrow(B3))
mu <- as.vector( w %*% B3 )

plot( dd$year, B3[1,]*w[1] , type="l" , lwd=3 , col=2 , xlab="year" , ylab="basis function" , ylim=c(-1,1) )
for ( i in 2:(nrow(B3)-1) ) lines( dd$year , B3[i,]*w[i] , lwd=3 , col=i+1 )

lines_w( dd$year , mu , lwd=6 )

# cycle through each function and animate range
ani.record(reset=TRUE)
wseq <- seq( from=1 , to=-1 , length=20 )
ww <- w
for ( i in 1:nrow(B3) ) {
    for ( wx in wseq ) {
        ww[i] <- w[i]*wx
        mu <- as.vector( ww %*% B3 )
        plot( dd$year, B3[1,]*ww[1] , type="l" , lwd=3 , col=2 , xlab="" , xaxt="n" , ylab=" weighted basis function" , ylim=c(-1,1) )
        for ( k in 2:(nrow(B3)-1) ) lines( dd$year , B3[k,]*ww[k] , lwd=3 , col=k+1 )
        lines_w( dd$year , mu , lwd=6 )
        ani.record()
    }
}

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 10 -loop 0 frame*.png spline_basis.gif

#######################################
# spline prior predictive animation

num_knots <- 20
knot_list <- seq( from=min(dd$year) , to=max(dd$year) , length.out=num_knots )
B3 <- t(bs(dd$year, knots=knot_list , degree=3, intercept = FALSE))

tau <- 10
al <- alist(
    Y ~ dnorm( mu , exp(log_sigma) ) ,
    mu <- a0 + as.vector( a %*% B ),
    a0 ~ dnorm(100,1),
    a ~ dnorm(0,tau),
    log_sigma ~ dnorm(0,0.5)
)

q <- quap( al , data=list( Y=dd$doy , B=B3 ) , start=list(a=rep(0,nrow(B3))) , dofit=FALSE )
pp <- extract.prior(q)

the_pf <- function(samp,step) {
    xlims <- c(850,2000)
    plot( NULL , cex=1.2 , ylab="doy first bloom" , xlim=xlims , bty="n" , axes=TRUE , xlab="year" , ylim=range(d3$doy) )
    mtext(concat("Prior, a ~ N(0,",tau,")"))
}

# anim_prior_predictive( prior=pp , linkf=a_link , n=5 , n_to_show=3 , n_steps=NSTEPS , ylim=ylims , ylab="outcome" , xlim=c(850,2000) , pf=the_pf , do_reset=FALSE , accel=-0.5 , vel0=1 )

l <- link(q,post=pp)
xseq <- dd$year

# function to draw line with white outline
lines_w <- function( x , y , lwd=2 , owd=3 , ... ) {
    lines( x , y , lwd=lwd+owd , col="white" )
    lines( x , y , lwd=lwd , ... )
}

ani.record(reset=TRUE)
ns <- 20
for ( i in 1:ns ) {
    the_pf()
    lines( xseq , l[i,] , lwd=4 , col=1 )
    lines_w( xseq , l[i+ns,] , lwd=4 , col=4 )
    #lines_w( xseq , l[i+2*ns,] , lwd=4 , col=5 )
    ani.record()
}

oopts = ani.options(interval = 0.5)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 30 -loop 0 frame*.png spline_prior.gif
# convert -delay 30 spline_prior.gif spline_prior.gif
# mogrify -layers 'optimize' -fuzz 7% ols_learn.gif

# posterior samples animation of same model

blossom_col <- sapply( d3$doy , function(y) hsv(1, rbeta2(1, inv_logit(logit(0.1)+0.02*y) ,10) ,1,0.8) )

the_pf <- function(samp,step) {
    xlims <- c(850,2000)
    plot( NULL , cex=1.2 , ylab="doy first bloom" , xlim=xlims , bty="n" , axes=TRUE , xlab="year" , ylim=range(d3$doy) )
    points( dd$year , dd$doy , col=blossom_col ,  pch=8  , cex=1.2 , lwd=2 )
    mtext("Posterior")
}

q2 <- quap( al , data=list( Y=dd$doy , B=B3 ) , start=list(a=rep(0,nrow(B3))) , dofit=TRUE )
pp <- extract.samples(q2)

l <- link(q2)
xseq <- dd$year

ani.record(reset=TRUE)
ns <- 20
for ( i in 1:ns ) {
    the_pf()
    lines_w( xseq , l[i,] , lwd=4 , col=1 )
    lines_w( xseq , l[i+ns,] , lwd=4 , col=4 )
    #lines_w( xseq , l[i+2*ns,] , lwd=4 , col=5 )
    ani.record()
}

oopts = ani.options(interval = 0.3)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 30 -loop 0 frame*.png spline_post.gif

#################
# now using Howell height ~ age data

library(rethinking)
library(splines)
library(animation)

data(Howell1)
d <- Howell1

blank(bty="n",w=1.6)
plot( d$age , d$height , col=2 , lwd=2 , xlab="age (years)" , ylab="height (cm)" )

num_knots <- 20
knot_list <- seq( from=min(d$age) , to=max(d$age) , length.out=num_knots )
B3 <- t(bs(d$age, knots=knot_list , degree=3, intercept = FALSE))

tau <- 25
al <- alist(
    Y ~ dnorm( mu , exp(log_sigma) ) ,
    mu <- a0 + as.vector( a %*% B ),
    a0 ~ dnorm(120,1),
    a ~ dnorm(0,tau),
    log_sigma ~ dnorm(0,0.5)
)

### now loop
ani.record(reset=TRUE)
npts_seq <- floor( seq( from=10 , to=nrow(d) , len=20 ) )
for ( npts in npts_seq ) {

plot( d$age[1:npts] , d$height[1:npts] , col=2 , lwd=2 , xlab="age (years)" , ylab="height (cm)" , xlim=range(d$age) , ylim=range(d$height) )

q <- quap( al , data=list( Y=d$height[1:npts] , B=B3[,1:npts] ) , start=list(a=rep(0,nrow(B3))) , dofit=TRUE )

pp <- extract.prior(q)
pp <- extract.samples(q)

the_pf <- function(samp,step) {
    xlims <- range(d$age)
    plot( NULL , cex=1.2 , ylab="height (cm)" , xlim=xlims , bty="n" , axes=TRUE , xlab="age (years)" , ylim=range(d$height) )
    points( d$age[1:npts] , d$height[1:npts] , col=2 , lwd=3 )
    mtext(concat(npts," individuals"))
}

#l <- link(q,post=pp)
# predictions for all points
l <- link(q,post=pp,data=list(Y=d$height,B=B3))

# function to draw line with white outline
lines_w <- function( x , y , lwd=2 , owd=3 , ... ) {
    lines( x , y , lwd=lwd+owd , col="white" )
    lines( x , y , lwd=lwd , ... )
}


ns <- 20
o <- order(d$age)
xseq <- d$age[o]
for ( i in 1:ns ) {
    the_pf()
    
    #lines( xseq , l[i,o] , lwd=4 , col=1 )
    #lines_w( xseq , l[i+2*ns,] , lwd=4 , col=5 )

    # basis functions
    for ( j in 1:23 ) lines_w( xseq , pp$a0[i] + pp$a[i,j]*B3[j,o] , lwd=2 , col=1 )

    lines_w( xseq , l[i,o] , lwd=4 , col=4 )

    ani.record()
}

}# frames

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=120)
# convert -alpha remove -background white -delay 10 -loop 0 frame*.png spline_howell_age.gif
# convert -delay 10 spline_howell_prior.gif spline_howell_age.gif

# convert -alpha remove -background white -delay 10 -loop 0 frame*.png spline_howell_age_basis.gif
