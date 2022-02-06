# lecture 9 - binomial GLMs

library(rethinking)
library(animation)
library(ellipse)

####
# UCBadmit is motivating example

# priors for logistic/binomial models

curve( inv_logit(x) , from=-6 , to=6 , lwd=4 , col=2 , xlab="logit(p)" , ylab="p" , yaxt="n" )
axis(2,at=c(0,0.5,1),labels=c(0,0.5,1))
abline(h=0.5,lty=3)

# intercept

dens( rnorm(1e4,0,1) , lwd=4 , col=4 , xlab="alpha" , xlim=c(-30,30) )

dens( inv_logit( rnorm(1e4,0,10) ) , lwd=4 , col=2 , xlab="probability of event" )

# slope

a <- rnorm(1e4,0,1.5)
b <- rnorm(1e4,0,0.5)
xseq <- seq(from=-3,to=3,len=100)
p <- sapply( xseq , function(x) inv_logit(a+b*x) )
plot( NULL , xlim=c(-2.5,2.5) , ylim=c(0,1) , xlab="x value" , ylab="probability" )
for ( i in 1:10 ) lines( xseq , p[i,] , lwd=3 , col=2 )

# animate priors and updating

n_points <- 20
a_true <- (-1)
b_true <- 1

x <- runif(n_points,-2,2)
y <- rbern(n_points, inv_logit(a_true + b_true*x) )
xlims <- c(-2,2)
ylims <- c(0,1)

ellipse_hist <- list()
ani.record(reset=TRUE)
NSTEPS <- 10
for ( i in 0:n_points ) {

    xmean <- mean( x[1:i] )
    al <- alist(
        y ~ dbern(p),
        logit(p) <- a + b*x ,
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

    a_link <- function(x) inv_logit( a + b*x )

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

            plot(NULL,xlim=c(-2.5,2.5),ylim=c(-2.5,2.5),bty="n",xlab="alpha",ylab="beta")
            points( a_true , b_true , pch=3 , col=2 , lwd=3 )
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
            if ( samp==1 ) thelwd <- thelwd + 4*(NSTEPS-step)/NSTEPS
            for ( l in c(0.5,0.8,0.95) )
                lines(the_contour[[concat(l)]],col=1,lwd=thelwd )

            # draw point for each line on righthand plot
            pb <- get("prior_blend",envir=parent.frame(n = 1))
            points( pb[[1]][,step] , pb[[2]][,step] , lwd=2 , col=use_cols[1:3] , pch=16 , cex=1.6 )

            # x/y plot with moving lines
            plot(NULL,xlim=xlims,ylim=ylims,xlab="x value",ylab="probability of event",bty="n" ,yaxt="n")
            axis(2,at=c(0,0.5,1),labels=c(0,0.5,1))
            abline(h=0.5,lty=3)
            if ( i > 0 ) {
                #shade( ci , xseq )
                #abline( a=ahat , b=bhat , lwd=3 , col="gray" , lty=2 )
                points( x[1:i] , y[1:i] , pch=1 , lwd=3 , cex=1.2 )
            }
            mtext(the_text)
        }

    if ( i==0 ) r <- NULL
    n_prior_samples <- 5
    if ( i==0 | i==n_points ) n_prior_samples <- 20
    r <- anim_prior_predictive( prior=pp , linkf=a_link , n=n_prior_samples , n_to_show=3 , n_steps=NSTEPS , ylim=ylims , ylab="probability of event" , xlim=xlims , pf=the_pf , do_reset=FALSE , start_from=r , accel=-0.5 , vel0=1 )

    for ( l in c(0.5,0.8,0.95) )
        ellipse_hist[[ concat(l) ]] <- ellipse(cov(pp),centre=colMeans(pp),level=l)

}

oopts = ani.options(interval = 0.02)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 5 -loop 0 frame*.png logistic_updating.gif
# convert -delay 3 ols_learn.gif ols_learn.gif

######
# UCBadmit priors

a <- rnorm(1e4,0,1.5) # gender means
delta <- rnorm(1e4,0,1.5) # dept means

xseq <- seq(from=-3,to=3,len=100)
p <- sapply( xseq , function(x) inv_logit(a) )
plot( NULL , xlim=c(-2.5,2.5) , ylim=c(0,1) , xlab="x value" , ylab="probability" )
for ( i in 1:10 ) lines( xseq , p[i,] , lwd=3 , col=2 )

# generative model, basic mediator scenario

N <- 1000 # number of applicants
# even gender distribution
G <- sample( 1:2 , size=N , replace=TRUE )
# gender 1 tends to apply to department 1, 2 to 2
D <- rbern( N , ifelse( G==1 , 0.3 , 0.8 ) ) + 1
# matrix of acceptance rates [dept,gender]
accept_rate <- matrix( c(0.1,0.3,0.1,0.3) , nrow=2 )
accept_rate <- matrix( c(0.05,0.2,0.1,0.3) , nrow=2 )
# simulate acceptance
p <- sapply( 1:N , function(i) accept_rate[D[i],G[i]] )
A <- rbern( N , p )

table(G,D)
table(G,A)
table( G , A , D )

# total effect gender
dat_sim <- list( A=A , D=D , G=G )

m1 <- ulam(
    alist(
        A ~ bernoulli(p),
        logit(p) <- a[G],
        a[G] ~ normal(0,1)
    ), data=dat_sim , chains=4 , cores=4 )

precis(m1,depth=2)

# direct effects
m2 <- ulam(
    alist(
        A ~ bernoulli(p),
        logit(p) <- a[G,D],
        matrix[G,D]:a ~ normal(0,1)
    ), data=dat_sim , chains=4 , cores=4 )

precis(m2,depth=3)

# contrasts

post1 <- extract.samples(m1)
G_contrast <- post1$a[,1] - post1$a[,2]
dens(G_contrast,lwd=4,col=2,xlab="Gender contrast")

post2 <- extract.samples(m2)
G_contrast_D1 <- post2$a[,1,1] - post2$a[,2,1]
G_contrast_D2 <- post2$a[,1,2] - post2$a[,2,2]
dens(G_contrast_D1,lwd=4,col=2,xlab="Gender contrast",ylim=c(0,2))
dens(G_contrast_D2,lwd=4,col=4,add=TRUE)

# on probability scale

diff_prob <- inv_logit( post1$a[,1] ) - inv_logit( post1$a[,2] )
dens(diff_prob,lwd=4,col=2,xlab="Gender contrast (probability)")

diff_prob_D1 <- inv_logit(post2$a[,1,1]) - inv_logit(post2$a[,2,1])
diff_prob_D2 <- inv_logit(post2$a[,1,2]) - inv_logit(post2$a[,2,2])
dens(diff_prob_D1,lwd=4,col=2,xlab="Gender contrast (probability)",xlim=c(-0.15,0.15))
dens(diff_prob_D2,lwd=4,col=4,add=TRUE)

# aggregate the dat_sim

x <- as.data.frame(cbind( A=dat_sim$A , G=dat_sim$G , D=dat_sim$D  ))
head(x,20)

dat_sim2 <- aggregate( A ~ G + D , dat_sim , sum )
dat_sim2$N <- aggregate( A ~ G + D , dat_sim , length )$A

m2_bin <- ulam(
    alist(
        A ~ binomial(N,p),
        logit(p) <- a[G,D],
        matrix[G,D]:a ~ normal(0,1)
    ), data=dat_sim2 , chains=4 , cores=4 )

precis(m2_bin,3)

# model real data - binomial data structure

data(UCBadmit)
d <- UCBadmit

dat <- list( 
    A = d$admit,
    N = d$applications,
    G = ifelse(d$applicant.gender=="female",1,2),
    D = as.integer(d$dept)
)

# total effect gender
mG <- ulam(
    alist(
        A ~ binomial(N,p),
        logit(p) <- a[G],
        a[G] ~ normal(0,1)
    ), data=dat , chains=4 , cores=4 )

precis(mG,2)

# direct effects
mGD <- ulam(
    alist(
        A ~ binomial(N,p),
        logit(p) <- a[G,D],
        matrix[G,D]:a ~ normal(0,1)
    ), data=dat , chains=4 , cores=4 )

precis(mGD,3)

# check chains

traceplot(mGD)
trankplot(mGD)

# contrasts
# on probability scale

post1 <- extract.samples(mG)
PrA_G1 <- inv_logit( post1$a[,1] )
PrA_G2 <- inv_logit( post1$a[,2] )
diff_prob <- PrA_G1 - PrA_G2
dens(diff_prob,lwd=4,col=2,xlab="Gender contrast (probability)")

post2 <- extract.samples(mGD)
PrA <- inv_logit( post2$a ) 
diff_prob_D_ <- sapply( 1:6 , function(i) PrA[,1,i] - PrA[,2,i] )
plot(NULL,xlim=c(-0.2,0.3),ylim=c(0,25),xlab="Gender contrast (probability)",ylab="Density")
for ( i in 1:6 ) dens( diff_prob_D_[,i] , lwd=4 , col=1+i , add=TRUE )
abline(v=0,lty=3)

# marginal effect of gender perception (direct effect)

# compute department weights via simulation
# we can just compute predictions as if all applications had been perceived as men
# and then again as if all had been perceived as women
# difference is marginal effect of perception, beause does not change department assignments (G -> A only, no G -> D)

# OLD WRONG CODE!
#p_G1 <- link( mGD , data=list(N=dat$N,D=dat$D,G=rep(1,12)) )
#p_G2 <- link( mGD , data=list(N=dat$N,D=dat$D,G=rep(2,12)) )

# NEW CORRECT CODE

# number of applicatons to simulate
total_apps <- sum(dat$N)

# number of applications per department
apps_per_dept <- sapply( 1:6 , function(i) sum(dat$N[dat$D==i]) )

# simulate as if all apps from women
p_G1 <- link(mGD,data=list(
    D=rep(1:6,times=apps_per_dept),
    N=rep(1,total_apps),
    G=rep(1,total_apps)))

# simulate as if all apps from men
p_G2 <- link(mGD,data=list(
    D=rep(1:6,times=apps_per_dept),
    N=rep(1,total_apps),
    G=rep(2,total_apps)))

# summarize
dens( p_G1 - p_G2 , lwd=4 , col=2 , xlab="effect of gender perception" )
abline(v=0,lty=3)

# show each dept density with weight as in population
w <- xtabs( dat$N ~ dat$D ) / sum(dat$N)
w <- w/max(w)
plot(NULL,xlim=c(-0.2,0.3),ylim=c(0,25),xlab="Gender contrast (probability)",ylab="Density")
for ( i in 1:6 ) dens( diff_prob_D_[,i] , lwd=2+8*w[i]^3 , col=1+i , add=TRUE )
abline(v=0,lty=3)

