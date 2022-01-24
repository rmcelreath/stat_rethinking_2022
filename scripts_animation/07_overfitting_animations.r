# overfitting animations etc for chapter 7

# animate cross-validation

library(rethinking)
library(animation)

# function to draw line with white outline
lines_w <- function( x , y , lwd=2 , owd=3 , ... ) {
    lines( x , y , lwd=lwd+owd , col="white" )
    lines( x , y , lwd=lwd , ... )
}

sppnames <- c( "afarensis","africanus","habilis","boisei",
    "rudolfensis","ergaster","sapiens")
brainvolcc <- c( 438 , 452 , 612, 521, 752, 871, 1350 )
masskg <- c( 37.0 , 35.5 , 34.5 , 41.5 , 55.5 , 61.0 , 53.5 )
d <- data.frame( species=sppnames , brain=brainvolcc , mass=masskg )
d <- d[ order(d$mass) , ]
d$B <- standardize(d$brain)
d$M <- standardize(d$mass)

# simulated larger dataset with same variables
NN <- 20
set.seed(13)

Msim <- sort( runif( NN , 30 , 65 ) )
# curve( 430 + 12*(x-30) , from=30 , to=65 )
Bsim <- rnorm( NN , 430 + 12*(Msim-30) , Msim*2 )
plot( Msim , Bsim )

d <- data.frame( brain=Bsim , mass=Msim )
d <- d[ order(d$mass) , ]
d$B <- standardize(d$brain)
d$M <- standardize(d$mass)


blank(bty="n")

waic_list <- matrix(NA,6,4)
psis_list <- matrix(NA,6,4)

# linear - 619 249 / 318 154.4
nterms <- 1
ylims <- range(d$brain)
flist <- alist(
            B ~ dnorm( mu , exp(log_sigma) ),
            mu <- a + b*M,
            a ~ dnorm( 0.5 , 1 ),
            b ~ dnorm( 0 , 10 ),
            log_sigma ~ dnorm( 0 , 1 ) )

# quadratic - 865 164.2 / 289 418.2
nterms <- 2
ylims <- range(d$brain)
flist <- alist(
            B ~ dnorm( mu , exp(log_sigma) ),
            mu <- a + b[1]*M + b[2]*M^2,
            a ~ dnorm( 0.5 , 1 ),
            b ~ dnorm( 0 , 10 ),
            log_sigma ~ dnorm( 0 , 1 ) )

# cubic - 12 538 298 / 201 039.5
nterms <- 3
ylims <- c(200,1500)
flist <- alist(
            B ~ dnorm( mu , exp(log_sigma) ),
            mu <- a + b[1]*M + b[2]*M^2 + b[3]*M^3,
            a ~ dnorm( 0.5 , 1 ),
            b ~ dnorm( 0 , 10 ),
            log_sigma ~ dnorm( 0 , 1 ) )

# 4 - 25 530 529 / 120 782.4
nterms <- 4
ylims <- c(200,1500)
flist <- alist(
            B ~ dnorm( mu , exp(log_sigma) ),
            mu <- a + b[1]*M + b[2]*M^2 + b[3]*M^3 + b[4]*M^4,
            a ~ dnorm( 0.5 , 1 ),
            b ~ dnorm( 0 , 10 ),
            log_sigma ~ dnorm( 0 , 1 ) )

# 5 - 293 840 335 / 7138.931
nterms <- 5
ylims <- c(-500,2500)
flist <- alist(
            B ~ dnorm( mu , exp(log_sigma) ),
            mu <- a + b[1]*M + b[2]*M^2 + b[3]*M^3 + b[4]*M^4 + b[5]*M^5,
            a ~ dnorm( 0.5 , 1 ),
            b ~ dnorm( 0 , 10 ),
            log_sigma ~ dnorm( 0 , 1 ) )

# 6 - 8 442 220 / 0.0003483536
nterms <- 6
ylims <- c(-500,2500)
flist <- alist(
            B ~ dnorm( mu , exp(log_sigma) ),
            mu <- a + b[1]*M + b[2]*M^2 + b[3]*M^3 + b[4]*M^4 + b[5]*M^5 + b[6]*M^6,
            a ~ dnorm( 0.5 , 1 ),
            b ~ dnorm( 0 , 10 ),
            log_sigma ~ dnorm( 0 , 1 ) )

# regularize
flist[[4]] <- b ~ dnorm(0,0.1)

#####################################
# ANIMATION LOOP
if ( nrow(d) > 7 ) ylims <- range(d$brain)
the_seed <- 1
ani.record(reset=TRUE)
yout <- rep(NA,nrow(d))
xseq <- seq(from=30,to=69,len=100)
xseq_std <- ( xseq - mean(d$mass) ) / sd(d$mass)
yprep <- matrix(NA,nrow=nrow(d),ncol=length(xseq))
transition_frames <- 10
for ( f in 1:nrow(d) ) {

    cols <- rep(1,nrow(d))
    cols[f] <- 2

    #m <- lm( brain ~ mass , d[-f,] )
    dat <- list(B=d$B[-f],M=d$M[-f])
    set.seed(the_seed)
    m <- quap(
        flist , data=dat , start=list(b=rep(0,nterms)) )
    
    set.seed(the_seed)
    mu <- link(m,data=list(M=xseq_std))
    yprep[f,] <- apply(mu,2,mean) * sd(d$brain) + mean(d$brain)

    ntf <- transition_frames
    #if ( f==1 ) ntf <- 1
    for ( tf in 1:ntf ) {

        plot( d$mass , d$brain , lwd=3 , col=cols , pch=1 , xlab="mass (kg)" , ylab="brain volume (cc)" , cex=1.2 , ylim=ylims )

        # draw previous curves
        if ( f > 1 )
            for ( g in 1:(f-1) ) {
                lines( xseq , yprep[g,] , lwd=3 , col=col.alpha(1,0.3) )
            }
        # draw current LOO curve
        if ( f==1 ) {
            # first point transitions from intercept
            ydraw <- (tf/ntf)*yprep[f,] + (ntf-tf)/ntf*rep(870,ncol(yprep))
        } else {
            # interpolate
            ydraw <- (tf/ntf)*yprep[f,] + (ntf-tf)/ntf*yprep[f-1,]
        }
        lines( xseq , ydraw , lwd=10 , col=col.alpha(2,0.5) )
        lines( xseq , ydraw , lwd=3 , col=2 )

        # draw LOO point prediction and error line
        set.seed(the_seed)
        yout[f] <- mean( link(m,data=list(M=d$M[f])) ) * sd(d$brain) + mean(d$brain)
        if ( f > 1 )
            for ( g in 1:(f-1) ) {
                lines( rep(d$mass[g],2) , c(yout[g],d$brain[g]) , lwd=3 , col=2 , lty=3 )
                points( d$mass[g] , yout[g] , pch=16 , col=2 , cex=1.4 )
            }

        ani.record()

    }#tf

    lines( rep(d$mass[f],2) , c(yout[f],d$brain[f]) , lwd=3 , col=2 , lty=3 )
    points( d$mass[f] , yout[f] , pch=16 , col=2 , cex=1.4 )

    invisible( replicate( 10 , ani.record() ) )
    
}#f

# final frame
nff <- 20
# compute in sample curve
dat <- list(B=d$B,M=d$M)
set.seed(the_seed)
m <- quap( flist , data=dat , start=list(b=rep(0,nterms)) )
set.seed(the_seed)
mu <- link(m,data=list(M=xseq_std))
ypfull <- apply(mu,2,mean)*sd(d$brain) + mean(d$brain)
for ( ff in 1:nff ) {

    plot( d$mass , d$brain , lwd=3 , col=4 , pch=1 , xlab="mass (kg)" , ylab="brain volume (cc)" , cex=1.2 , ylim=ylims )
    for ( g in 1:f ) {
        lines( xseq , yprep[g,] , lwd=3 , col=col.alpha(1,0.3) )
        lines( rep(d$mass[g],2) , c(yout[g],d$brain[g]) , lwd=3 , col=2 , lty=3 )
        points( d$mass[g] , yout[g] , pch=16 , col=2 , cex=1.4 )
    }
    # add in-sample curve
    V <- ff/nff
    lines( xseq , ypfull , lwd=11 , col=col.alpha("white",V) )
    lines( xseq , ypfull , lwd=8 , col=col.alpha(4,V) )

    # record
    ani.record()

}
# freeze frame
invisible( replicate( 25 , ani.record() ) )
#####################################

# prediction error
sum( (yout-d$brain)^2 ) / 1000
# error in sample
mu <- link(m)
yhat <- apply(mu,2,mean)*sd(d$brain) + mean(d$brain)
sum( (yhat-d$brain)^2 ) / 1000

waic_list[nterms,] <- as.numeric( WAIC(m) )
colnames(waic_list) <- names(WAIC(m))
psis_list[nterms,] <- as.numeric( PSIS(m) )
colnames(psis_list) <- names(PSIS(m))

# show
oopts = ani.options(interval = 0.01)
ani.replay()
# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 5 -loop 0 frame*.png overfit_sim_06.gif

# graph fit in and out

# lppd of single point
ll <- sim(m, ll = TRUE)
dens(ll[,2],lwd=4,col=2,xlab="log posterior prob of observation")


# original data
x <- c(
619249 , 318154.4 ,
865164.2 , 289418.2 ,
12538298 , 201039.5 ,
25530529 , 120782.4 ,
293840335 , 7138.931
)

r <- t( matrix( x , nrow=2 , ncol=5 ) )
colnames(r) <- c("out","in")

plot( 1:5 , r[,1]/sum(r[,1]) , lwd=4 , col=2 , xlab="polynomial terms" , ylab="relative error" , cex=1.5 , type="b" )
points( 1:5 , r[,2]/sum(r[,2]) , lwd=4 , col=4 , cex=1.5 , type="b" )

# synthetic scores: degree, out , in
x <- c(
258.4789 , 212.3881 ,
214.4751 , 160.2875 ,
235.6202 , 142.0885 ,
249.3383 , 127.7217 ,
326.4922 , 124.8697 ,
1315.365 , 124.6889
)

r <- t( matrix( x , nrow=2 , ncol=6 ) )
colnames(r) <- c("out","in")

pr <- 1:6
plot( pr , r[pr,1]/sum(r[pr,1]) , lwd=4 , col=2 , xlab="polynomial terms" , ylab="relative error" , cex=1.5 , type="b" , ylim=c(0.05,0.3) )
points( pr , r[pr,2]/sum(r[pr,2]) , lwd=4 , col=4 , cex=1.5 , type="b" )

pr <- 1:5
plot( pr , r[pr,1] , lwd=4 , col=2 , xlab="polynomial terms" , ylab="prediction error" , cex=1.5 , type="b" , ylim=range(r[pr,]) )
points( pr , r[pr,2] , lwd=4 , col=4 , cex=1.5 , type="b" )

# waic / psis

pr <- 1:5
plot( pr , waic_list[pr,1] , lwd=4 , col=2 , xlab="polynomial terms" , ylab="log pointwise predictive density" , cex=1.5 , type="b" , ylim=c(30,60) )
points( pr , -2*waic_list[pr,2] , lwd=4 , col=4 , cex=1.5 , type="b" )

points( pr , psis_list[pr,1] , lwd=4 , col=6 , cex=1.5 , type="b" )
#points( pr , -2*(psis_list[pr,2]+psis_list[pr,3]) , lwd=4 , col=4 , cex=1.5 , type="b" )



# plot distance out-in

penalty <- sapply( 1:5 , function(i) r[i,1] - r[i,2] )
plot( 1:5 , penalty , xlab="polynomial terms" , ylab="out-of-sample penalty" , lwd=4 , col=1 , cex=1.5 , ylim=c(0,200) , type="b" )

# a ~ dnorm(0,1)
# out , in
x1 <- c(
258.6889 , 212.491 , # x
212.8562 , 160.3505 , # x
226.5399 , 143.32 , # x
227.1352 , 130.074 , # x
307.727 , 124.8697 , # x
476.5134 , 129.2839 # x
)

r1 <- t( matrix( x1 , nrow=2 , ncol=6 ) )
colnames(r) <- c("out","in")

pr <- 1:5
plot( pr , r[pr,1] , lwd=4 , col=2 , xlab="polynomial terms" , ylab="prediction error" , cex=1.5 , type="b" , ylim=range(r[pr,]) )
points( pr , r[pr,2] , lwd=4 , col=4 , cex=1.5 , type="b" )

points( pr , r1[pr,1] , lwd=3 , col=2 , cex=1.5 , type="b" )
points( pr , r1[pr,2] , lwd=3 , col=4 , cex=1.5 , type="b" )

# a ~ dnorm(0,0.5)
x2 <- c(
260.6711 , 214.1892 , # x
209.0693 , 161.0249 , # x
218.6809 , 152.4985 , # x
216.3307 , 144.2392 , # x
277.3484 , 144.6256 , # x
350.5387 , 140.6424 # x
)

r2 <- t( matrix( x2 , nrow=2 , ncol=6 ) )
colnames(r) <- c("out","in")

points( pr , r2[pr,1] , lwd=2 , col=2 , cex=1.5 , type="b" )
points( pr , r2[pr,2] , lwd=2 , col=4 , cex=1.5 , type="b" )

# a ~ dnorm(0,0.1)
x3 <- c(
402.673 , 357.5877 , # x
340.6418 , 292.4061 , # x
263.9933 , 224.6419 , # x
263.7808 , 220.5813 , # x
272.8136 , 228.147 , # x
279.0085 , 226.4686 # x
)

r3 <- t( matrix( x3 , nrow=2 , ncol=6 ) )
colnames(r) <- c("out","in")

pr <- 1:5
plot( pr , r[pr,1] , lwd=4 , col=2 , xlab="polynomial terms" , ylab="prediction error" , cex=1.5 , type="b" , ylim=c(min(r[pr,]),max(r3[pr,])) )
points( pr , r[pr,2] , lwd=4 , col=4 , cex=1.5 , type="b" )

points( pr , r1[pr,1] , lwd=3 , col=2 , cex=1.5 , type="b" )
points( pr , r1[pr,2] , lwd=3 , col=4 , cex=1.5 , type="b" )

points( pr , r2[pr,1] , lwd=2 , col=2 , cex=1.5 , type="b" )
points( pr , r2[pr,2] , lwd=2 , col=4 , cex=1.5 , type="b" )

points( pr , r3[pr,1] , lwd=1 , col=2 , cex=1.5 , type="b" )
points( pr , r3[pr,2] , lwd=1 , col=4 , cex=1.5 , type="b" )

###############
# draw leverage example

y <- c(-0.8,0.18,2.5,-0.13)

idx <- c(1,2,3,4)

m00 <- quap(
    alist(
        y ~ dnorm(mu,exp(log_sigma)),
        mu ~ dnorm(0,1),
        log_sigma ~ dnorm(0,1)
    ) , data=list(y=y[idx]) )

post <- extract.samples(m00)
s <- rnorm(1e4,post$mu,exp(post$log_sigma))
dens( s , lwd=4 , col=1 , xlab="" , xlim=c(-5,5) , ylab="" , yaxt="n" )

dens( ss , add=TRUE , lty=3 , lwd=3 )

for ( i in idx ) abline(v=y[i],lwd=3,col=2)

abline( v=y[-idx] , lwd=3 , col="gray" )

curve( dnorm( y[1] , x , 1.6 ) , from=-4 , to=4 )

# importance sampling

sd( rnorm(1e4,post$mu,exp(post$log_sigma)) )

ll <- sim(m00, ll = TRUE )
r <- 1/exp(ll[,3])

yy <- sample( ss[1:1000] , size=1e3 , replace=TRUE , prob=r )

#################
# plant example again, this time with model comparison

set.seed(71)
# number of plants
N <- 100
# simulate initial heights
h0 <- rnorm(N,10,2)
# assign treatments and simulate fungus and growth
treatment <- rep( 0:1 , each=N/2 )
fungus <- rbinom( N , size=1 , prob=0.5 - treatment*0.4 )
h1 <- h0 + rnorm(N, 5 - 3*fungus)
# compose a clean data frame
d <- data.frame( h0=h0 , h1=h1 , treatment=treatment , fungus=fungus )

plot( jitter(d$fungus) , d$h1/d$h0 , col=ifelse(d$treatment==1,2,4) , lwd=3 , xaxt="n" , xlab="" , ylab="growth" )
axis( 1 , at=c(0,1) , labels=c("no fungus","yo fungus") )

plot( jitter(d$treatment) , d$h1/d$h0 , col=ifelse(d$fungus==1,2,4) , lwd=3 , xaxt="n" , xlab="" , ylab="growth" )
axis( 1 , at=c(0,1) , labels=c("control","treatment") )

abline(lm( I(d$h1/d$h0)[d$fungus==1] ~ treatment[d$fungus==1] ))
abline(lm( I(d$h1/d$h0)[d$fungus==0] ~ treatment[d$fungus==0] ))

m6.6 <- quap(
    alist(
        h1 ~ dnorm( mu , sigma ),
        mu <- h0*p,
        p ~ dlnorm( 0 , 0.25 ),
        sigma ~ dexp( 1 )
    ), data=d )

m6.7 <- quap(
    alist(
        h1 ~ dnorm( mu , sigma ),
        mu <- h0 * p,
        p <- a + bt*treatment + bf*fungus,
        a ~ dlnorm( 0 , 0.2 ) ,
        bt ~ dnorm( 0 , 0.5 ),
        bf ~ dnorm( 0 , 0.5 ),
        sigma ~ dexp( 1 )
    ), data=d )

m6.8 <- quap(
    alist(
        h1 ~ dnorm( mu , sigma ),
        mu <- h0 * p,
        p <- a + bt*treatment,
        a ~ dlnorm( 0 , 0.2 ),
        bt ~ dnorm( 0 , 0.5 ),
        sigma ~ dexp( 1 )
    ), data=d )

plot(precis(m6.7,pars=c("bt","bf")))

plot(precis(m6.8,pars=c("bt")))

post <- extract.samples(m6.7)
dens( post$bt , lwd=4 , xlab="effect of treatment (posterior)" , xlim=c(-0.15,0.2) )
post <- extract.samples(m6.8)
dens( post$bt , lwd=4 , col=2 , add=TRUE )
abline( v=0 , lwd=2 , lty=3 )

compare( m6.7 , m6.8 , func=PSIS )

plot(compare( m6.7 , m6.8 , func=PSIS ))

###############
# outlier example

library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d$A <- standardize( d$MedianAgeMarriage )
d$D <- standardize( d$Divorce )
d$M <- standardize( d$Marriage )

m5.1 <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bA * A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
) , data = d )
m5.2 <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM * M ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = d )
m5.3 <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
) , data = d )

compare( m5.1 , m5.2 , m5.3 , func=PSIS )

PSIS_m5.3 <- PSIS(m5.3,pointwise=TRUE)
set.seed(24071847)
WAIC_m5.3 <- WAIC(m5.3,pointwise=TRUE)

plot( PSIS_m5.3$k , WAIC_m5.3$penalty , xlab="PSIS Pareto k" ,
    ylab="WAIC penalty" , col=2 , lwd=3 , cex=1.5 )

identify( PSIS_m5.3$k , WAIC_m5.3$penalty , d$Location )

plot( d$A , d$D , col=ifelse(d$Loc %in% c("ME","ID"),2,grau(0.7)) , lwd=3 , cex=1.5 , xlab="Age at marriage (std)" , ylab="Divorce rate (std)" )

identify( d$A , d$D , d$Location )

# sim student t

curve( dnorm( x , 0 , 0.5 ) , from=-7 , to=7 , lwd=4 , xlab="value" , ylab="density" , col=grau(0.3) )
for ( i in 1:10 ) curve( dnorm( x , 0 , runif(1,0.5,2) ) , add=TRUE , lwd=4 , col=grau(0.3) )

x <- rnorm( 1e4 , 0 , runif(1e4,0.5,2) )
dens(x)
curve(dnorm(x,0,1),add=TRUE,lwd=3,col=2)

curve( dnorm( x , 0 , 1) , from=-7 , to=7 , lwd=4 , xlab="value" , ylab="density" )
curve( dstudent( x , 2 , 0 , 1 ) , add=TRUE , col=2 , lwd=4 )

y <- rstudent( 1e4 , 2 , 0 , 1 )
dens( y , xlim=c(-7,7) , lwd=4 )

# model

library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
dat <- list(
    A = standardize( d$MedianAgeMarriage ),
    D = standardize( d$Divorce ),
    M = standardize( d$Marriage )
)

m5.3 <- quap(
    alist(
        D ~ dnorm( mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
) , data = dat )

m5.3t <- quap(
    alist(
        D ~ dstudent( 2 , mu , sigma ) ,
        mu <- a + bM*M + bA*A ,
        a ~ dnorm( 0 , 0.2 ) ,
        bM ~ dnorm( 0 , 0.5 ) ,
        bA ~ dnorm( 0 , 0.5 ) ,
        sigma ~ dexp( 1 )
    ) , data = dat )

post <- extract.samples(m5.3)
dens( post$bA , xlab="bA (effect of age of marriage)" , lwd=4 , ylim=c(0,3) )
post <- extract.samples(m5.3t)
dens( post$bA , add=TRUE , lwd=4 , col=2 )


library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d$A <- standardize( d$MedianAgeMarriage )
d$D <- standardize( d$Divorce )
d$M <- standardize( d$Marriage )


muG <- link(m5.3)
muT <- link(m5.3t)
muG_mean <- apply( muG , 2, mean )
muT_mean <- apply( muT , 2, mean )

plot( muG_mean , muT_mean , col=ifelse(d$Loc %in% c("ME","ID"),2,grau(0.7)) , lwd=3 , cex=1.5 , xlab="Predicted divorce (Gaussian)" , ylab="Predicted divorce (Student-t)" , xlim=c(-2,2) , ylim=c(-2,2) )
abline(a=0,b=1,lty=3,lwd=2)

identify( d$A , d$D , d$Location )

