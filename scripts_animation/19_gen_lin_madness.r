# lecture 19
# generalized linear madness
# principled scientific causal models

library(rethinking)
library(animation)
library(ellipse)

#######################################################
###### geometric weight model

## R code 16.1
library(rethinking)
data(Howell1)
d <- Howell1

# scale observed variables
d$w <- d$weight / mean(d$weight)
d$h <- d$height / mean(d$height)

blank(bty="n")

plot( d$height , d$weight , xlim=c(0,max(d$height)) , ylim=c(0,max(d$weight)) , col=col.alpha(2,0.7) ,
    lwd=3 , xlab="height (cm)" , ylab="weight (kg)" )
mw <- mean(d$weight)
mh <- mean(d$height)
lines( c(mh,mh) , c(0,mw) , lty=3 , lwd=2 )
lines( c(0,mh) , c(mw,mw) , lty=3 , lwd=2 )

plot( d$h , d$w , xlim=c(0,max(d$h)) , ylim=c(0,max(d$w)) , col=col.alpha(2,0.7) ,
    lwd=3 , xlab="height (scaled)" , ylab="weight (scaled)" )
lines( c(1,1) , c(0,1) , lty=3 , lwd=2 )
lines( c(0,1) , c(1,1) , lty=3 , lwd=2 )

# prior sim
n <- 30
p <- rbeta(1e4,25,50)
k <- rexp(1e4,0.5)
sigma <- rexp(n,1)
xseq <- seq(from=0,to=1.3,len=100)
plot(NULL,xlim=c(0,1.3),ylim=c(0,1.5),xlab="height (scaled)",ylab="weight (scaled)")
for ( i in 1:n ) {
    mu <- log( pi * k[i] * p[i]^2 * xseq^3 )
    lines( xseq , exp(mu + sigma[i]^2/2) , lwd=3 , col=col.alpha(2,runif(1,0.4,0.8)) )
}
lines( c(1,1) , c(0,1) , lty=3 , lwd=2 )
lines( c(0,1) , c(1,1) , lty=3 , lwd=2 )

curve( dbeta(x,25,50) , from=0, to=1 , xlim=c(0,1) , lwd=3 , col=2 , xlab="p")

curve( dexp(x,0.5) , from=0, to=5 , xlim=c(0,5) , lwd=3 , col=2 , xlab="k")

## R code 16.2
dat <- list(W=d$w,H=d$h)
m16.1 <- ulam(
    alist(
        W ~ dlnorm( mu , sigma ),
        exp(mu) <- 3.141593 * k * p^2 * H^3,
        p ~ beta( 25 , 50 ),
        k ~ exponential( 0.5 ),
        sigma ~ exponential( 1 )
    ), data=dat , chains=4 , cores=4 )

mWH1 <- ulam(
    alist(
        W ~ dlnorm( mu , sigma ),
        exp(mu) <- 3.141593 * theta * H^3,
        theta ~ exponential( 1 ),
        sigma ~ exponential( 1 )
    ), data=dat , chains=4 , cores=4 )

mWH2 <- ulam(
    alist(
        W ~ dlnorm( mu , sigma ),
        exp(mu) <- H^3 + 0*theta, # 0*theta is hack to get link/sim to work
        theta ~ exponential( 1 ),
        sigma ~ exponential( 1 )
    ), data=dat , chains=4 , cores=4 )

## R code 16.3
h_seq <- seq( from=0 , to=max(d$h) , length.out=100 )
w_sim <- sim( mWH2 , data=list(H=h_seq) )
mu_mean <- apply( w_sim , 2 , mean )
w_CI <- apply( w_sim , 2 , PI )
plot( d$h , d$w , xlim=c(0,max(d$h)) , ylim=c(0,max(d$w)) , col=2 ,
    lwd=2 , xlab="height (scaled)" , ylab="weight (scaled)" )
shade( w_CI , h_seq , col=col.alpha(2,0.5) )
lines( h_seq , mu_mean , lwd=3 )

post <- extract.samples(m16.1)
plot( post$p^2 , post$k , lwd=3 , col=2 , xlab="p^2" , ylab="k" )
curve( 1/(x*pi) , add=TRUE , lwd=3 )

mWH3 <- ulam(
    alist(
        W ~ dlnorm( mu , sigma ),
        exp(mu) <- 3.141593 * exp(k + ka*H) * exp(p + (pa*H))^2 * H^3,
        b ~ exponential(1),
        p ~ beta( 2 , 18 ),
        k ~ exponential( 0.5 ),
        c(ka,pa) ~ normal(0,0.5),
        sigma ~ exponential( 1 )
    ), data=dat , chains=4 , cores=4 )

h_seq <- seq( from=0 , to=max(d$h) , length.out=100 )
w_sim <- sim( mWH2 , data=list(H=h_seq) )
mu_mean <- apply( w_sim , 2 , mean )
w_CI <- apply( w_sim , 2 , PI )
plot( d$h , d$w , xlim=c(0,max(d$h)) , ylim=c(0,max(d$w)) , col=2 ,
    lwd=2 , xlab="height (scaled)" , ylab="weight (scaled)" )
shade( w_CI , h_seq , col=col.alpha(2,0.5) )
lines( h_seq , mu_mean , lwd=3 )

#######################################
# state-based example

## R code 16.4
library(rethinking)
data(Boxes)
precis(Boxes)

## R code 16.5
table( Boxes$y ) / length( Boxes$y )

## R code 16.6
set.seed(7)

N <- 100 # number of children

# half choose random color
# sample from 1,2,3 at random for each
y1 <- sample( 1:3 , size=N/2 , replace=TRUE )

# half follow majority
y2 <- rep( 2 , N/2 )

# combine and shuffle y1 and y2
y <- sample( c(y1,y2) )

# count the 2s
sum(y==2)/N

plot( table(y) , lwd=12 , col=c(4,0,7) , xlim=c(0.5,3.5) , xlab="" , ylab="frequency" , xaxt="n" )
axis(1,at=1:3,labels=c("unchosen","majority","minority"))
lines( c(2,2) , c(0,sum(y1==2)-2) , lwd=12 , col=2 )
lines( c(2,2) , c(sum(y1==2)+2,sum(y==2)) , lwd=12 , col=2 )

## R code 16.7
data(Boxes_model)
cat(Boxes_model)

## R code 16.8
# prep data
dat_list <- list(
    N = nrow(Boxes),
    y = Boxes$y,
    majority_first = Boxes$majority_first )

# run the sampler
m16.2 <- cstan( model_code=Boxes_model , data=dat_list , chains=4 , cores=4 )

# show marginal posterior for p
p_labels <- c("1 Majority","2 Minority","3 Maverick","4 Color",
    "5 First")
plot( precis(m16.2,2) , labels=p_labels )

pairs( m16.2 )

# show animation of random draws from p

post <- extract.samples(m16.2)

library(animation)
blank(bty="n",w=2)

ani.record(reset=TRUE)
nf <- 30
tf <- 20 # transition frames
ns <- 3 # number of samples to show at once
x <- tf*nf
cols <- c(2,4,5)
p_old <- list()
for ( i in 1:ns ) p_old[[i]] <- rep(1/5,5)
for ( f in 1:nf ) {

    for ( i in 0:tf ) {

        plot( NULL , xlim=c(0.75,5.25) , ylim=c(0,0.5) , xlab="" , ylab="posterior probability" , xaxt="n" )
        axis(1,at=1:5,labels=c("majority","minority","maverick","color","first"))

        v <- i/tf
        for ( s in 1:ns ) {
            p <- v*post$p[f+x*(s-1),] + (1-v)*p_old[[s]]
            lines( 1:5 , p , lwd=4 , type="b" , col=cols[s] )
        }#s

        ani.record()

    }#i

    for ( s in 1:ns ) p_old[[s]] <- post$p[f+x*(s-1),]

}#f

oopts = ani.options(interval = 0.01)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 4 -loop 0 frame*.png box_p.gif
# convert -delay 2 box_p.gif box_p.gif
# convert box_p.gif -coalesce -fuzz 2% +dither -layers Optimize +map box_p_z.gif



############################################
# ODEs

## R code 16.13
library(rethinking)
data(Lynx_Hare)

blank(bty="n",w=1.6)

plot( 1:21 , Lynx_Hare[,3] , ylim=c(0,90) , xlab="year" ,
    ylab="thousands of pelts" , xaxt="n" , type="l" , lwd=3  )
at <- c(1,11,21)
axis( 1 , at=at , labels=Lynx_Hare$Year[at] )
lines_w( 1:21 , Lynx_Hare[,2] , lwd=3 , col=2 )
points( 1:21 , Lynx_Hare[,3] , bg="black" , col="white" , pch=21 , cex=1.4 )
points( 1:21 , Lynx_Hare[,2] , bg=2 , col="white" , pch=21 , cex=1.4 )
text( 17 , 80 , "Lepus" , pos=2 )
text( 19 , 50 , "Lynx" , pos=2 , col=2 )

## R code 16.14
sim_lynx_hare <- function( n_steps , init , theta , dt=0.002 ) {
    L <- rep(NA,n_steps)
    H <- rep(NA,n_steps)
    L[1] <- init[1]
    H[1] <- init[2]
    for ( i in 2:n_steps ) {
        H[i] <- H[i-1] + dt*H[i-1]*( theta[1] - theta[2]*L[i-1] )
        L[i] <- L[i-1] + dt*L[i-1]*( theta[3]*H[i-1] - theta[4] )
    }
    return( cbind(L,H) )
}

## R code 16.15
theta <- c( 0.5 , 0.05 , 0.025 , 0.5 )
z <- sim_lynx_hare( 1e4 , as.numeric(Lynx_Hare[1,2:3]) , theta )

plot( z[,2] , type="l" , ylim=c(0,max(z[,2])) , lwd=4 , xaxt="n" ,
    ylab="number (thousands)" , xlab="" )
lines( z[,1] , col=2 , lwd=4 )
mtext( "time" , 1 )

##############################
### animate simulation
ani.record(reset=TRUE)
# sample from priors
theta <- rmvnorm2( 5000 , c(0.5,0.05,0.025,0.5) , rep(0.01,4) )
nf <- 30
tf <- 20 # transition frames
ns <- 1 # number of samples to show at once
x <- tf*nf
maxt <- 1e4
cols <- c(2,4,5)
p_old <- list()
for ( i in 1:ns ) p_old[[i]] <- theta[1,]
for ( f in 2:nf ) {

    for ( i in 0:tf ) {

        plot( NULL , xlim=c(0,maxt) , type="l" , ylim=c(0,60) , lwd=4 , xaxt="n" , ylab="number (thousands)" , xlab="" )

        v <- i/tf
        for ( s in 1:ns ) {
            theta_t <- v*theta[f+x*(s-1),] + (1-v)*p_old[[s]]
            z <- sim_lynx_hare( maxt , as.numeric(Lynx_Hare[1,2:3]) , theta_t )
            lines( 1:maxt , z[,2] , col=1 , lwd=4 )
            lines_w( 1:maxt , z[,1] , col=2 , lwd=4 )
        }#s

        ani.record()

    }#i

    for ( s in 1:ns ) p_old[[s]] <- theta[f+x*(s-1),]

}#f

oopts = ani.options(interval = 0.01)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 4 -loop 0 frame*.png lynx_sim.gif
# convert -delay 2 lynx_sim.gif lynx_sim.gif

##############################

## R code 16.16
N <- 1e4
Ht <- 1e4
p <- rbeta(N,2,18)
h <- rbinom( N , size=Ht , prob=p )
h <- round( h/1000 , 2 )
dens( h , xlab="thousand of pelts" , lwd=2 )

## R code 16.17
data(Lynx_Hare_model)
cat(Lynx_Hare_model)

## R code 16.18
dat_list <- list(
    N = nrow(Lynx_Hare),
    pelts = Lynx_Hare[,2:3] )

m16.5 <- cstan( model_code=Lynx_Hare_model , data=dat_list , chains=3 ,
    cores=3 , control=list( adapt_delta=0.95 ) )

## R code 16.19
post <- extract.samples(m16.5)
pelts <- dat_list$pelts
plot( 1:21 , pelts[,2] , pch=16 , ylim=c(0,140) , xlab="year" ,
    ylab="thousands of pelts" , xaxt="n" )
at <- c(1,11,21)
axis( 1 , at=at , labels=Lynx_Hare$Year[at] )
# 21 time series from posterior
for ( s in 1:21 ) {
    lines( 1:21 , post$pelts_pred[s,,2] , col=col.alpha("black",0.2) , lwd=3 )
    lines( 1:21 , post$pelts_pred[s,,1] , col=col.alpha(2,0.3) , lwd=3 )
}
# points
points( 1:21 , pelts[,2] , col="white" , pch=16 , cex=1.5 )
points( 1:21 , pelts[,2] , col=1 , pch=16 )
points( 1:21 , pelts[,1] , col="white" , pch=16 , cex=1.5 )
points( 1:21 , pelts[,1] , col=2 , pch=16 )
# text labels
text( 17 , 110 , "Lepus" , pos=2 )
text( 19 , 50 , "Lynx" , pos=2 , col=2 )

#### animated version
ani.record(reset=TRUE)
nf <- 30
tf <- 20 # transition frames
ns <- 1 # number of samples to show at once
x <- tf*nf
maxt <- 1e4
cols <- c(2,4,5)
obserr <- 0.185
p_old <- list()
for ( i in 1:ns ) p_old[[i]] <- post$pop[1,,]
for ( f in 2:nf ) {

    for ( i in 0:tf ) {

        plot( 1:21 , pelts[,2] , pch=16 , ylim=c(0,85) , xlab="year" , ylab="thousands of pelts" , xaxt="n" , col=0 )
        at <- c(1,11,21)
        axis( 1 , at=at , labels=Lynx_Hare$Year[at] )

        v <- i/tf
        for ( s in 1:ns ) {
            posts <- v*post$pop[f+x*(s-1),,] + (1-v)*p_old[[s]]
            lines( 1:21 , posts[,2] * obserr , col=1 , lwd=4 )
            lines_w( 1:21 , posts[,1] * obserr , col=2 , lwd=4 )
        }#s

        # points
        points( 1:21 , pelts[,2] , col="white" , pch=16 , cex=1.5 )
        points( 1:21 , pelts[,2] , col=1 , pch=16 )
        points( 1:21 , pelts[,1] , col="white" , pch=16 , cex=1.5 )
        points( 1:21 , pelts[,1] , col=2 , pch=16 )

        ani.record()

    }#i

    for ( s in 1:ns ) p_old[[s]] <- post$pop[f+x*(s-1),,]

}#f

oopts = ani.options(interval = 0.01)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 4 -loop 0 frame*.png lynx_post.gif
# convert -delay 2 lynx_post.gif lynx_post.gif

############


