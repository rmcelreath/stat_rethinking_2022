# animate prior predictive simulation

library(rethinking)
library(animation)
library(ellipse)

blank(bty="n",w=1.4,h=0.8)

# OLS

al <- alist(
     y ~ dnorm(mu,1),
     mu <- a + b*x ,
     a ~ dnorm(0,1),
     b ~ dnorm(0,1)
)

q <- quap( al , data=list(y=1) , dofit=FALSE )
pp <- extract.prior(q)

a_link <- function(x) a + b*x

anim_prior_predictive( prior=pp , linkf=a_link , n=20 , n_to_show=3 , n_steps=20 , ylim=c(-2,2) , ylab="outcome" )

oopts = ani.options(interval = 0.03)
ani.replay()

# OLS bayesian updating with animation

n_points <- 10
n_points <- 0

x <- rnorm(n_points)
x <- ifelse( x > 2 , 2 , x )
x <- ifelse( x < -2 , -2 , x )
y <- rnorm(n_points, 0 + 0.7*x,0.5)
xlims <- c(-2,2)
ylims <- c(-2,2)

ellipse_hist <- list()
ani.record(reset=TRUE)
NSTEPS <- 20
for ( i in 0:n_points ) {

    xmean <- mean( x[1:i] )
    al <- alist(
        y ~ dnorm(mu,1),
        mu <- a + b*x ,
        a ~ dnorm(0,1),
        b ~ dnorm(0,1)
    )

    if ( i > 0 ) {
        dati <- list( x=x[1:i] , y=y[1:i] )
        q <- quap( al , data=dati , dofit=TRUE )
        # print( precis(q) )
        pp <- extract.samples(q)
    } else {
        # i = 0
        q <- quap( al , data=list(y=1,x=1) , dofit=FALSE )
        pp <- as.data.frame(extract.prior(q))
    }

    a_link <- function(x) a + b*x

    the_text <- concat("n = ",i)
    xseq <- seq(from=-2.2,to=2.2,length=50)
    ahat <- mean(pp$a)
    bhat <- mean(pp$b)
    if ( i > 0 ) {
        mu <- link( q , data=list(x=xseq) )
        ci <- apply(mu,2,PI)
    }

    use_cols <- c(2, 4, 5, 6, 3)

    the_pf <- function(samp,step) {
            par(mfrow=c(1,2))

            plot(NULL,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),bty="n",xlab="intercept",ylab="slope")
            pm <- colMeans(pp)
            #lines( c(pm[1],pm[1]) , c(-4,pm[2]) , lty=2 )
            #lines( c(-4,pm[1]) , c(pm[2],pm[2]) , lty=2 )
            
            # interpolate between previous and current contours
            the_contour <- list()
            for ( l in c(0.5,0.8,0.95) )
                the_contour[[ concat(l) ]] <- ellipse(cov(pp),centre=pm,level=l)
            if ( length(ellipse_hist) > 0 & samp==1 ) {
                mixf <- (NSTEPS-step)/NSTEPS
                for ( l in c(0.5,0.8,0.95) ) {
                    #lines(ellipse_hist[[i]],col=col.alpha(1,0.3),lwd=1)
                    the_contour[[ concat(l) ]] <- (1-mixf)*the_contour[[concat(l)]] + mixf*ellipse_hist[[concat(l)]]
                }
            }

            # draw posterior contours
            thelwd <- 2
            if ( samp==1 ) thelwd <- 2+3*(NSTEPS-step)/NSTEPS
            for ( l in c(0.5,0.8,0.95) )
                lines(the_contour[[concat(l)]],col=1,lwd=thelwd )

            # draw point for each line on righthand plot
            pb <- get("prior_blend",envir=parent.frame(n = 1))
            points( pb[[1]][,step] , pb[[2]][,step] , lwd=2 , col=use_cols[1:3] , pch=16 , cex=1.6 )

            # x/y plot with moving lines
            plot(NULL,xlim=xlims,ylim=ylims,xlab="x value",ylab="y value",bty="n")
            if ( i > 0 ) {
                shade( ci , xseq )
                abline( a=ahat , b=bhat , lwd=3 , col="gray" , lty=2 )
                points( x[1:i] , y[1:i] , pch=1 , lwd=3 , cex=1.2 )
            }
            mtext(the_text)
        }

    if ( i==0 ) r <- NULL
    n_prior_samples <- 20
    r <- anim_prior_predictive( prior=pp , linkf=a_link , n=n_prior_samples , n_to_show=3 , n_steps=NSTEPS , ylim=ylims , ylab="outcome" , xlim=xlims , pf=the_pf , do_reset=FALSE , start_from=r , accel=-0.5 , vel0=1 )

    for ( l in c(0.5,0.8,0.95) )
        ellipse_hist[[ concat(l) ]] <- ellipse(cov(pp),centre=colMeans(pp),level=l)

}

oopts = ani.options(interval = 0.02)
ani.replay()

# ani.saveqz(dpi=130)
# convert -alpha remove -background white -delay 3 -loop 0 frame*.png ols_learn.gif
# convert -delay 3 ols_learn.gif ols_learn.gif
# mogrify -layers 'optimize' -fuzz 7% ols_learn.gif

############
# OLS parabolic model

set.seed(12)

n_points <- 10
x <- rnorm(n_points)
x <- ifelse( x > 2 , 2 , x )
x <- ifelse( x < -2 , -2 , x )
y <- rnorm( n_points , 0 + 0.7*x - 1*x^2 , 0.5 )
x[9] <- 3
y[9] <- -1

xlims <- c(-4,5)
ylims <- c(-7,2)

ellipse_hist <- list()
ani.record(reset=TRUE)
NSTEPS <- 14
for ( i in 0:n_points ) {

    al <- alist(
        y ~ dnorm(mu,1),
        mu <- a + b1*x + b2*x^2 ,
        a ~ dnorm(0,1),
        b1 ~ dnorm(0,1),
        b2 ~ dnorm(0,1)
    )

    if ( i > 0 ) {
        dati <- list( x=x[1:i] , y=y[1:i] )
        q <- quap( al , data=dati , dofit=TRUE )
        pp <- extract.samples(q)
    } else {
        # i = 0
        q <- quap( al , data=list(y=1,x=1) , dofit=FALSE )
        pp <- as.data.frame(extract.prior(q))
    }

    a_link <- function(x) a + b1*x + b2*x^2

    the_text <- concat("n = ",i)
    xseq <- seq(from=xlims[1],to=xlims[2],length=50)
    ahat <- mean(pp$a)
    b1hat <- mean(pp$b1)
    b2hat <- mean(pp$b2)
    if ( i > 0 ) {
        mu <- link( q , data=list(x=xseq) )
        ci <- apply(mu,2,PI)
    }

    the_pf <- function(samp,step) {
            par(mfrow=c(1,2))

            plot(NULL,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),bty="n",xlab="b1",ylab="b2")
            mtext(the_text)
            pm <- colMeans(pp)
            #lines( c(pm[1],pm[1]) , c(-4,pm[2]) , lty=2 )
            #lines( c(-4,pm[1]) , c(pm[2],pm[2]) , lty=2 )
            
            # interpolate between previous and current contours
            the_contour <- list()
            for ( l in c(0.5,0.8,0.95) )
                the_contour[[ concat(l) ]] <- ellipse(cov(pp)[2:3,2:3],centre=pm[2:3],level=l)
            if ( length(ellipse_hist) > 0 & samp==1 ) {
                mixf <- (NSTEPS-step)/NSTEPS
                for ( l in c(0.5,0.8,0.95) ) {
                    #lines(ellipse_hist[[i]],col=col.alpha(1,0.3),lwd=1)
                    the_contour[[ concat(l) ]] <- (1-mixf)*the_contour[[concat(l)]] + mixf*ellipse_hist[[concat(l)]]
                }
            }

            # draw posterior contours
            thelwd <- 2
            if ( samp==1 ) thelwd <- 2+3*(NSTEPS-step)/NSTEPS
            for ( l in c(0.5,0.8,0.95) )
                lines(the_contour[[concat(l)]],col=2,lwd=thelwd )

            plot(NULL,xlim=xlims,ylim=ylims,xlab="x value",ylab="y value",bty="n")
            if ( i > 0 ) {
                shade( ci , xseq )
                points( x[1:i] , y[1:i] , pch=1 , lwd=3 , cex=1.2 )
            }
            mtext( "y = a + b1*x + b2*x^2" )
            
        }

    if ( i==0 ) r <- NULL
    r <- anim_prior_predictive( prior=pp , linkf=a_link , n=5 , n_to_show=3 , n_steps=NSTEPS , ylim=ylims , ylab="outcome" , xlim=xlims , pf=the_pf , do_reset=FALSE , start_from=r , accel=-0.5 , vel0=1 )

    for ( l in c(0.5,0.8,0.95) )
        ellipse_hist[[ concat(l) ]] <- ellipse(cov(pp)[2:3,2:3],centre=colMeans(pp)[2:3],level=l)

}

oopts = ani.options(interval = 0.02)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 5 -loop 0 frame*.png quad_learn.gif
# convert -delay 3 ols_learn.gif quad_learn.gif
# mogrify -layers 'optimize' -fuzz 7% quad_learn.gif

# ffmpeg -i das_model.mp3 -ignore_loop 0 -i quad_learn.gif -shortest -strict -2 -c:v libx264 -threads 4 -c:a aac -b:a 192k -pix_fmt yuv420p -shortest out.mp4

# ffmpeg -i das_model.mp3 0 -i quad_learn.gif -ignore_loop -shortest -strict -2 -shortest out.mp4


############
# OLS cubic model

blank(bty="n",w=1.4,h=0.8)

set.seed(13)

n_points <- 10
x <- rnorm(n_points)
x <- ifelse( x > 2 , 2 , x )
x <- ifelse( x < -2 , -2 , x )
y <- rnorm( n_points , 0 + 0.7*x - 2*x^2 + 3*x^3 , 0.5 )
# plot(x,y)

xlims <- range(x) + c(-1,1)*2
ylims <- range(y) + c(-1,1)*2

ellipse_hist <- list()
ani.record(reset=TRUE)
NSTEPS <- 14
for ( i in 0:n_points ) {

    al <- alist(
        y ~ dnorm(mu,1),
        mu <- a + b1*x + b2*x^2 + b3*x^3,
        a ~ dnorm(0,1),
        b1 ~ dnorm(0,1),
        b2 ~ dnorm(0,1),
        b3 ~ dnorm(0,1)
    )

    if ( i > 0 ) {
        dati <- list( x=x[1:i] , y=y[1:i] )
        q <- quap( al , data=dati , dofit=TRUE )
        pp <- extract.samples(q)
    } else {
        # i = 0
        q <- quap( al , data=list(y=1,x=1) , dofit=FALSE )
        pp <- as.data.frame(extract.prior(q))
    }

    a_link <- function(x) a + b1*x + b2*x^2 + b3*x^3

    the_text <- concat("n = ",i)
    xseq <- seq(from=xlims[1],to=xlims[2],length=50)
    ahat <- mean(pp$a)
    b1hat <- mean(pp$b1)
    b2hat <- mean(pp$b2)
    if ( i > 0 ) {
        mu <- link( q , data=list(x=xseq) )
        ci <- apply(mu,2,PI)
    }

    the_pf <- function(samp,step) {
            par(mfrow=c(1,2))

            plot(NULL,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),bty="n",xlab="b2",ylab="b3")
            mtext(the_text)
            pm <- colMeans(pp)
            #lines( c(pm[1],pm[1]) , c(-4,pm[2]) , lty=2 )
            #lines( c(-4,pm[1]) , c(pm[2],pm[2]) , lty=2 )
            
            # interpolate between previous and current contours
            the_contour <- list()
            for ( l in c(0.5,0.8,0.95) )
                the_contour[[ concat(l) ]] <- ellipse(cov(pp)[3:4,3:4],centre=pm[3:4],level=l)
            if ( length(ellipse_hist) > 0 & samp==1 ) {
                mixf <- (NSTEPS-step)/NSTEPS
                for ( l in c(0.5,0.8,0.95) ) {
                    #lines(ellipse_hist[[i]],col=col.alpha(1,0.3),lwd=1)
                    the_contour[[ concat(l) ]] <- (1-mixf)*the_contour[[concat(l)]] + mixf*ellipse_hist[[concat(l)]]
                }
            }

            # draw posterior contours
            thelwd <- 2
            if ( samp==1 ) thelwd <- 2+3*(NSTEPS-step)/NSTEPS
            for ( l in c(0.5,0.8,0.95) )
                lines(the_contour[[concat(l)]],col=2,lwd=thelwd )

            plot(NULL,xlim=xlims,ylim=ylims,xlab="x value",ylab="y value",bty="n")
            if ( i > 0 ) {
                shade( ci , xseq )
                points( x[1:i] , y[1:i] , pch=1 , lwd=3 , cex=1.2 )
            }
            mtext( "y = a + b1*x + b2*x^2 + b3*x^3" )
            
        }

    if ( i==0 ) r <- NULL
    r <- anim_prior_predictive( prior=pp , linkf=a_link , n=5 , n_to_show=3 , n_steps=NSTEPS , ylim=ylims , ylab="outcome" , xlim=xlims , pf=the_pf , do_reset=FALSE , start_from=r , accel=-0.5 , vel0=1 )

    for ( l in c(0.5,0.8,0.95) )
        ellipse_hist[[ concat(l) ]] <- ellipse(cov(pp)[3:4,3:4],centre=colMeans(pp)[3:4],level=l)

}

oopts = ani.options(interval = 0.03)
ani.replay()

# ani.saveqz(dpi=130)
# convert -alpha remove -background white -delay 6 -loop 0 frame*.png cubic_learn.gif
# convert -delay 3 ols_learn.gif quad_learn.gif
# mogrify -layers 'optimize' -fuzz 7% quad_learn.gif
