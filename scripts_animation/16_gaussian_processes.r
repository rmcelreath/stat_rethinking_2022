# lecture 16
# gaussian processes

library(rethinking)
library(animation)
library(ellipse)

# animate samples from gaussian process
# build kernel and then draw y vector from multi_normal given x vector
# https://mc-stan.org/docs/2_29/stan-users-guide/simulating-from-a-gaussian-process.html

# start with example with fixed kernel

# precompile model
ma <- cmdstan_model(stan_file="16_gp_fast_f.stan", compile = TRUE )

#data {
#  int<lower=1> N1;
#  array[N1] real x1;
#  vector[N1] y1;
#  int<lower=1> N2;
#  array[N2] real x2;
#}

cols <- c(2,4,5)
ani.record(reset=TRUE)
post_prev <- NULL
n_frames <- 20
tf <- 10 # transition frames
nfunc <- 3 # number of draws to animate simultaneously

N1 <- 4
x1 <- c(0.5,0.2,0.85,0.3)[1:N1]
y1 <- c(0,1,-2,-1)[1:N1]
N2 <- 100
x2 <- seq(from=-0.2,to=1.2,len=N2)
# alpha=1, rho=0.2, sigma=0.25
datx <- list( N1=N1 , x1=x1 , y1=y1 , N2=N2 , x2=x2 , 
    alpha=sqrt(2) , rho=0.35 , sigma=0.35 )
max <- ma$sample(data = datx, chains = 4, 
            parallel_chains = 4, iter_sampling = 1000, iter_warmup = 1000, 
            adapt_delta = 0.9, max_treedepth = 15, 
            save_warmup = FALSE , refresh=0 )
maxx <- rstan::read_stan_csv(max$output_files())
post <- extract.samples(maxx)

# the loop
par(mfrow=c(1,2))
for( f in 1:n_frames ) {

    for ( g in 1:tf ) {

        plot( NULL , xlim=c(0,1) , ylim=c(-3,3) , xlab="x variable" , ylab="y variable" )

        for ( i in 1:nfunc ) {

            o <- (i-1)*n_frames
            if ( f > 1 ) {
                f2a <- post$f2[f-1+o,]
            } else {
                if ( !is.null(post_prev) ) {
                    f2a <- post_prev$f2[f+o,]
                } else {
                    # just start from same spot
                    f2a <- post$f2[f+o,]
                }
            }
            f2b <- post$f2[f+o,]
            f2 <- (g/tf)*f2b + (1-g/tf)*f2a
            lines_w( x2 , f2 , lwd=4 , col=cols[i] )

        }#i
        # observed points
        if ( datx$sigma > 0.1 )
            for ( i in 1:N1 )
                lines( c(x1[i],x1[i]) , y1[i] + c(-1,1)*2*datx$sigma , lwd=12 , col=grau(0.3) )
        points( x1 , y1 , pch=1 , lwd=4 , cex=1.5 )
        text( x1 , y1 , 1:N1 , pos=3 )

        # now draw kernel on right
        if ( !is.null(post_prev) && (f==1) ) {
            a0 <- datx_prev$alpha
            r0 <- datx_prev$rho
        } else {
            a0 <- datx$alpha
            r0 <- datx$rho
        }
        a <- (g/tf)*datx$alpha + (1-g/tf)*a0
        r <- (g/tf)*datx$rho + (1-g/tf)*r0
        curve( a^2 * exp( -1/(2*r^2)*x^2 ) , from=0, to=1 , lwd=5 , col=1 , xlab="distance" , ylab="covariance" , ylim=c(0,2.1) )

        # distance points on kernel
        ox <- 0.05
        if ( N1 > 1 ) {
            for ( i in 1:(N1-1) )
                for ( j in i:N1 ) {
                    if ( i !=j ) {
                        dij <- abs( x1[i] - x1[j] )
                        kij <- a^2 * exp( -1/(2*r^2)*dij^2 )
                        points( dij , kij , pch=16 , col="white" )
                        points( dij , kij , lwd=4 , cex=1.5 )
                        text( dij + ox , kij , concat(i,",",j) , pos=3 )
                    }
                }#ij
        }

        # record
        ani.record()
    }#g - transition frames
}#f
post_prev <- post # for smooth transitions btw simulations
datx_prev <- datx # store prev kernel shape

oopts = ani.options(interval = 0.01)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 4 -loop 0 frame*.png gp_demo.gif
# convert -delay 4 gp_demo.gif gp_demo.gif
# convert gp_demo.gif -coalesce -fuzz 2% +dither -layers Optimize +map gp_demo_z.gif

####
# GP animation - learning the kernel version

# precompile model
ma <- cmdstan_model(stan_file="16_gp_fast_l.stan" , compile = TRUE )

# blank(bty="n",w=1.6)

#data {
#  int<lower=1> N1;
#  array[N1] real x1;
#  vector[N1] y1;
#  int<lower=1> N2;
#  array[N2] real x2;
#}

cols <- c(2,4,5)
ani.record(reset=TRUE)
post_prev <- NULL
n_frames <- 20
tf <- 10 # transition frames
nfunc <- 3 # number of draws to animate simultaneously

N1 <- 4
# initital examples
x1 <- c(0.5,0.2,0.85,0.3,0.6,0.7)[1:N1]
y1 <- c(0,1,-2,-1,-1,2)[1:N1]
# regularization example
x1 <- c(0.2,0.25,0.8,0.75,0.8)[1:N1]
y1 <- c(-1,0,-2,2,-1.5)[1:N1]
N2 <- 100
x2 <- seq(from=-0.2,to=1.2,len=N2)
# alpha=1, rho=0.2, sigma=0.25
datx <- list( N1=N1 , x1=x1 , y1=y1 , N2=N2 , x2=x2 , 
alphal=0 , alphah=1.7 , ra=2 , rb=2 )
max <- ma$sample(data = datx, chains = 4, 
            parallel_chains = 4, iter_sampling = 1000, iter_warmup = 1000, 
            adapt_delta = 0.9, max_treedepth = 15, 
            save_warmup = FALSE , refresh=0 )
maxx <- rstan::read_stan_csv(max$output_files())
post <- extract.samples(maxx)

# the loop
par(mfrow=c(1,2))
for( f in 1:n_frames ) {

    for ( g in 1:tf ) {

        plot( NULL , xlim=c(0,1) , ylim=c(-3,3) , xlab="x variable" , ylab="y variable" )

        for ( i in 1:nfunc ) {

            o <- (i-1)*n_frames
            if ( f > 1 ) {
                f2a <- post$f2[f-1+o,]
            } else {
                if ( !is.null(post_prev) ) {
                    f2a <- post_prev$f2[f+o,]
                } else {
                    # just start from same spot
                    f2a <- post$f2[f+o,]
                }
            }
            f2b <- post$f2[f+o,]
            f2 <- (g/tf)*f2b + (1-g/tf)*f2a
            lines_w( x2 , f2 , lwd=4 , col=cols[i] )

        }#i
        # observed points
        for ( i in 1:N1 ) {
            # intervals
            for ( j in 1:nfunc ) {
                o <- (j-1)*n_frames
                if ( f > 1 ) {
                    s <- post$sigma[f-1+o]
                } else {
                    if ( !is.null(post_prev) ) {
                        s <- post_prev$sigma[f+o]
                    } else {
                        # just start from same spot
                        s <- post$sigma[f+o]
                    }
                }
                s1 <- post$sigma[f+o]
                s <- (g/tf)*s1 + (1-g/tf)*s
                ox <- 0
                od <- 0.025
                if ( j==1 ) ox <- (-od)
                if ( j==3 ) ox <- od
                if ( nfunc==1 ) ox <- 0
                lines( c(x1[i],x1[i])+ox , y1[i] + c(-1,1)*2*s , lwd=16/nfunc , col=col.alpha(cols[j],0.6) )
            }#j
        }#i
        points( x1 , y1 , pch=1 , lwd=4 , cex=1.5 )

        if ( TRUE ) {
        plot( NULL , xlim=c(0,1) , ylim=c(0,3) , xlab="distance" , ylab="covariance" )
        # now draw kernel on right
        for ( i in 1:nfunc ) {

            o <- (i-1)*n_frames
            if ( f > 1 ) {
                a0 <- post$alpha[f-1+o]
                r0 <- post$rho[f-1+o]
            } else {
                if ( !is.null(post_prev) ) {
                    a0 <- post_prev$alpha[f+o]
                    r0 <- post_prev$rho[f+o]
                } else {
                    # just start from same spot
                    a0 <- post$alpha[f+o]
                    r0 <- post$rho[f+o]
                }
            }
            a1 <- post$alpha[f+o]
            r1 <- post$rho[f+o]
            a <- (g/tf)*a1 + (1-g/tf)*a0
            r <- (g/tf)*r1 + (1-g/tf)*r0
            curve( a^2 * exp( -1/(2*r^2)*x^2 ) , from=0, to=1 , lwd=5 , col=cols[i] , add=TRUE )

        }#i
        }#draw kernel

        # record
        ani.record()
    }#g - transition frames
}#f
post_prev <- post # for smooth transitions btw simulations

oopts = ani.options(interval = 0.01)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 4 -loop 0 frame*.png gp_demo.gif
# convert -delay 4 gp_demo.gif gp_demo.gif
# convert gp_demo.gif -coalesce -fuzz 5% +dither -layers Optimize +map gp_demo_z.gif
# convert gp_demo.gif -coalesce -layers Optimize +map gp_demo_z.gif


# example kernels

# L2
curve( exp( -(x^2)/0.45^2 ) , from=0 , to=1 , xlab="|x1 - x2|" , lwd=4 , col=2 , ylim=c(0,1) , ylab="covariance" )

# L1
curve( exp( -(x)/0.45^2 ) , from=0 , to=1 , xlab="|x1 - x2|" , lwd=4 , col=2 , ylim=c(0,1) , ylab="covariance" )

# sin^2
curve( exp( -(2/0.7)*sin(x/2)^2 ) , from=0 , to=2*pi , xlab="|x1 - x2|" , lwd=4 , col=2 , ylim=c(0,1) , ylab="covariance" )


############################
# kline distance model

# sim priors for distance model
n <- 30
etasq <- rexp(n,2)
rhosq <- rexp(n,1)
plot( NULL , xlim=c(0,7) , ylim=c(0,2) , xlab="distance (thousand km)" , ylab="covariance" )
for ( i in 1:n )
    curve( etasq[i]*exp(-rhosq[i]*x^2) , add=TRUE , lwd=4 , col=col.alpha(2,0.5) )


data(Kline2) # load the ordinary data, now with coordinates
d <- Kline2
data(islandsDistMatrix)
# display (measured in thousands of km)
Dmat <- islandsDistMatrix
colnames(Dmat) <- c("Ml","Ti","SC","Ya","Fi","Tr","Ch","Mn","To","Ha")
round(Dmat,1)

dat_list <- list(
    T = d$total_tools,
    P = d$population,
    S = 1:10,
    D = islandsDistMatrix )

mTdist <- ulam(
    alist(
        T ~ dpois(lambda),
        log(lambda) <- abar + a[S],
        vector[10]:a ~ multi_normal( 0 , K ),
        transpars> matrix[10,10]:K <- cov_GPL2(D,etasq,rhosq,0.01),
        abar ~ normal(3,0.5),
        etasq ~ dexp( 2 ),
        rhosq ~ dexp( 0.5 )
    ), data=dat_list , chains=4 , cores=4 , iter=4000 , log_lik=TRUE )

precis(mTdist,2)

mTdist_nc <- ulam(
    alist(
        T ~ dpois(lambda),
        log(lambda) <- abar + a[S],
        # non-centered Gaussian Process prior
        transpars> vector[10]: a <<- L_K * z,
        vector[10]: z ~ normal( 0 , 1 ),
        transpars> matrix[10,10]: L_K <<- cholesky_decompose( K ),
        transpars> matrix[10,10]: K <- cov_GPL2(D,etasq,rhosq,0.01),
        abar ~ normal(3,0.5),
        etasq ~ dexp( 2 ),
        rhosq ~ dexp( 0.5 )
    ), data=dat_list , chains=4 , cores=4 , iter=2000 , log_lik=TRUE )

precis(mTdist_nc,2)

# plot posterior kernel

# repeat priors
n <- 30
etasq <- rexp(n,2)
rhosq <- rexp(n,1)
plot( NULL , xlim=c(0,7) , ylim=c(0,2) , xlab="distance (thousand km)" , ylab="covariance" )
for ( i in 1:n )
    curve( etasq[i]*exp(-rhosq[i]*x^2) , add=TRUE , lwd=4 , col=col.alpha(1,0.5) )
# post
n <- 30
post <- extract.samples(mTdist_nc)
for ( i in 1:n )
    curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE , lwd=4 , col=col.alpha(2,0.5) )

# show spatial associations

# compute posterior median covariance among societies
K <- extract.samples(mTdist_nc)$K
# K <- extract.samples(mTDP)$K

# scale point size to logpop
psize <- d$logpop / max(d$logpop)
psize <- exp(psize*1.5)-1

# plot raw data and labels
plot( d$lon2 , d$lat , xlab="longitude" , ylab="latitude" ,
    col=2 , cex=psize , pch=16 , xlim=c(-50,30) , ylim=c(-25,25) )
# overlay lines shaded by cov
Kmean <- apply(K,2:3,mean)
Kmean <- Kmean / max(Kmean)
for( i in 1:10 )
    for ( j in 1:10 )
        if ( i < j )
            lines( c( d$lon2[i],d$lon2[j] ) , c( d$lat[i],d$lat[j] ) ,
                lwd=4 , col=col.alpha("black",Kmean[i,j]) )

points( d$lon2 , d$lat , col=2 , cex=psize , pch=16 )
labels <- as.character(d$culture)
text( d$lon2 , d$lat , labels=labels , cex=0.7 , pos=c(2,4,3,3,4,1,3,2,4,2) )

# now full model

dat_list <- list(
    T = d$total_tools,
    P = d$population,
    S = 1:10,
    D = islandsDistMatrix )

mTDP <- ulam(
    alist(
        T ~ dpois(lambda),
        lambda <- (abar*P^b/g)*exp(a[S]),
        vector[10]:a ~ multi_normal( 0 , K ),
        transpars> matrix[10,10]:K <- cov_GPL2(D,etasq,rhosq,0.01),
        c(abar,b,g) ~ dexp( 1 ),
        etasq ~ dexp( 2 ),
        rhosq ~ dexp( 0.5 )
    ), data=dat_list , chains=4 , cores=4 , iter=4000 , log_lik=TRUE )

precis( mTDP , depth=2 )

post <- extract.samples(mTDP)

# repeat previous model
n <- 30
post0 <- extract.samples(mTdist_nc)
plot( NULL , xlim=c(0,7) , ylim=c(0,2) , xlab="distance (thousand km)" , ylab="covariance" )
for ( i in 1:n )
    curve( etasq[i]*exp(-rhosq[i]*x^2) , add=TRUE , lwd=4 , col=col.alpha(1,0.5) )
# post
n <- 30
post <- extract.samples(mTDP)
for ( i in 1:n )
    curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE , lwd=4 , col=col.alpha(2,0.5) )


# plot relationship with population
post <- extract.samples(mTDP)

plot( log(dat_list$P) , dat_list$T , xlab="log population" , ylab="tools" , col=2 , cex=psize , pch=16 )

xseq <- seq( from=6 , to=14 , len=100 )
l <- sapply( xseq , function(x) with(post, (abar * exp(x)^b/g) ) )
lines_w( xseq , apply(l,2,mean) , lwd=4 , col=4 )
shade( apply(l,2,PI) , xseq , col=col.alpha(4,0.25) )

# overlay lines shaded by cov
K <- post$K
Kmean <- apply(K,2:3,mean)
Kmean <- Kmean / max(Kmean)
x <- log(dat_list$P)
y <- dat_list$T
for( i in 1:10 )
    for ( j in 1:10 )
        if ( i < j )
            lines( c( x[i],x[j] ) , c( y[i],y[j] ) ,
                lwd=4 , col=col.alpha("black",Kmean[i,j]) )
points( x , y , col=2 , cex=psize , pch=16 )



##################################################################
# primates phylogeny

## R code 14.47
library(rethinking)
data(Primates301)
data(Primates301_nex)
d <- Primates301

# plot it using ape package - install.packages('ape') if needed
library(ape)
d <- Primates301
l <- ladderize(Primates301_nex)

cols <- rep(1,301)
cols[1:30] <- 4 # NWM
cols[31:106] <- 2 # lemurs
cols[107:173] <- 3
cols[174:198] <- 5 # apes
cols[199:208] <- 3
cols[209:286] <- 4
cols[287:301] <- 7
cols[164:173] <- 7

par(bg="black")
plot( l , type="fan" , font=1 , no.margin=TRUE ,
    label.offset=1 , cex=0.55 , edge.width=1.6 , edge.col="white" , tip.col=cols , rotate.tree=37 )

## R code 14.48
d <- Primates301
d$name <- as.character(d$name)
dstan <- d[ complete.cases( d$group_size , d$body , d$brain ) , ]
spp_obs <- dstan$name

# tree again for complete cases only
lcc <- keep.tip( l , spp_obs )
idx <- sapply( lcc$tip.label , function(spp) which( l$tip.label == spp ) )
colscc <- cols[ idx ]
par(bg="black")
plot( lcc , type="fan" , font=1 , no.margin=TRUE ,
    label.offset=1 , cex=0.8 , edge.width=1.6 , edge.col="white" , tip.col=colscc , rotate.tree=30 )

# tree again highlighting complete cases
colsh <- rep(grau(),301)
idx2 <- which( l$tip.label %in% spp_obs )
colsh <- ifelse( l$tip.label %in% spp_obs , cols , grau() )
par(bg="black")
plot( l , type="fan" , font=1 , no.margin=TRUE ,
    label.offset=1 , cex=0.55 , edge.width=1.6 , edge.col="white" , tip.col=colsh , rotate.tree=37 )

# show distributions of focal variables
blank(bty="n",ex=0.7)
par(bg="black",col.axis="white",fg="white")
plot( log(dstan$body) , log(dstan$brain) , lwd=3 , col=2 )
plot( log(dstan$group_size) , log(dstan$brain) , lwd=3 , col=2 )
plot( log(dstan$body) , log(dstan$group_size) , lwd=3 , col=2 )

# traits on tree
par(bg="black")

xx <- plot( lcc , type="fan" , font=1 , no.margin=TRUE ,
    label.offset=1 , cex=0.8 , edge.width=1.6 , edge.col="white" , tip.col=colscc , rotate.tree=30 , show.tip.label=FALSE )

idx <- sapply( lcc$tip.label , function(spp) which( dstan$name == spp ) )
B <- (dstan$brain)[idx]
M <- log(dstan$body)[idx]
G <- log(dstan$group_size)[idx]
tiplabels( pch=16 , col=colscc , cex=normalize(B)+0.75 )
tiplabels( pch=1 , lwd=2 , col=colscc , cex=normalize(M)+0.5 , offset=2.2 )
tiplabels( pch=2 , lwd=2 , col=colscc , cex=normalize(G)+0.5 , offset=4.5 )

# plot key
par(bg="black")
plot(NULL,xlim=c(-1,1),ylim=c(-1,1),bty="n")
points( rep(0,3) , c(0.25,0,-0.25) , pch=c(16,1,2) , col=5 , lwd=4 , cex=3 )

## R code 14.49
data(Primates301)
d <- Primates301
# complete case analysis
dstan <- d[complete.cases(d$group_size,d$body,d$brain),]
dat_list <- list(
    N_spp = nrow(dstan),
    M = standardize(log(dstan$body)),
    B = standardize(log(dstan$brain)),
    G = standardize(log(dstan$group_size)),
    Imat = diag(nrow(dstan)) )

# classical regression form
mBMG0 <- ulam(
    alist(
        B ~ normal( mu , sigma ),
        mu <- a + bM*M + bG*G,
        a ~ normal( 0 , 1 ),
        c(bM,bG) ~ normal( 0 , 0.5 ),
        sigma ~ exponential( 1 )
    ), data=dat_list , chains=4 , cores=4 )

# multivariate form
mBMG <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + bM*M + bG*G,
        matrix[N_spp,N_spp]: K <- Imat*(sigma^2),
        a ~ normal( 0 , 1 ),
        c(bM,bG) ~ normal( 0 , 0.5 ),
        sigma ~ exponential( 1 )
    ), data=dat_list , chains=4 , cores=4 )

precis( mBMG )
precis( mBMG0 )

## simulaton
par(bg="black")

s <- rphylo(151,3,0.01)
xx <- plot( s , type="fan" , font=1 , no.margin=TRUE ,
    label.offset=1 , cex=0.8 , edge.width=1.6 , edge.col="white" , tip.col=colscc , rotate.tree=30 , show.tip.label=FALSE )

# s <- lcc

nf <- 20
tf <- 10
CA <- col.alpha
ani.record(reset=TRUE)
for ( f in 1:nf ) {
    # sim traits under brownian motion
    st <- rTraitDisc( s , model=matrix(c(0, 0.1, 0.1, 0), 2))
    for ( i in 1:tf ) {
        a <- min( 1, (tf-i+5)/tf )
        cols <- c( CA(2,a) , CA(4,a) )
        ew <- 0.5 * (1+4*a)
        plot( s , type="fan" , font=1 , no.margin=TRUE ,
        label.offset=1 , cex=0.8 , edge.width=ew , edge.col="white" , tip.col=colscc , rotate.tree=30 , show.tip.label=FALSE )
        tiplabels( pch=16 , col=ifelse(st=="A",cols[1],cols[2]) , cex=2 , offset=0 )
        ani.record()
    }
}

oopts = ani.options(interval = 0.02)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 11 -loop 0 frame*.png tree_sim.gif
# convert -delay 11 tree_sim.gif tree_sim.gif
# convert tree_sim.gif -coalesce -fuzz 2% +dither -layers Optimize +map tree_sim_z.gif
# convert tree_sim.gif -coalesce -layers Optimize +map tree_sim_z.gif



## R code 14.50
library(ape)
tree_trimmed <- keep.tip( Primates301_nex, spp_obs )
Rbm <- corBrownian( phy=tree_trimmed )
V <- vcv(Rbm)
Dmat <- cophenetic( tree_trimmed )

plot( Dmat , V , xlab="phylogenetic distance" , ylab="covariance" , lwd=3 , col=2 , type="l" )
curve( max(V)*exp(-(1/28)*x) , add=TRUE , lwd=4 , col=4 )

## R code 14.51
# put species in right order
dat_list$V <- V[ spp_obs , spp_obs ]
# convert to correlation matrix
dat_list$R <- dat_list$V / max(V)

# Brownian motion model
mBMG_brownian <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + bM*M + bG*G,
        matrix[N_spp,N_spp]: K <- R * sigma_sq,
        a ~ normal( 0 , 1 ),
        c(bM,bG) ~ normal( 0 , 0.5 ),
        sigma_sq ~ exponential( 1 )
    ), data=dat_list , chains=4 , cores=4 )

precis( mBMG_brownian )

# Ornstein-Uhlenbeck (L1 gaussian process)
# add scaled and reordered distance matrix
dat_list$Dmat <- Dmat[ spp_obs , spp_obs ] / max(Dmat)
mB_OU <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + 0*M,
        matrix[N_spp,N_spp]:K <- cov_GPL1(Dmat,etasq,rho,0.01),
        a ~ normal(0,1),
        etasq ~ half_normal(1,0.25),
        rho ~ half_normal(3,0.25)
    ), data=dat_list , chains=4 , cores=4 )

precis( mB_OU )

mBM_OU <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + bM*M,
        matrix[N_spp,N_spp]:K <- cov_GPL1(Dmat,etasq,rho,0.01),
        a ~ normal(0,1),
        bM ~ normal(0,0.5),
        etasq ~ half_normal(1,0.25),
        rho ~ half_normal(3,0.25)
    ), data=dat_list , chains=4 , cores=4 )

precis( mBM_OU )

mBMG_OU <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + bM*M + bG*G,
        matrix[N_spp,N_spp]:K <- cov_GPL1(Dmat,etasq,rho,0.01),
        a ~ normal(0,1),
        c(bM,bG) ~ normal(0,0.5),
        etasq ~ half_normal(1,0.25),
        rho ~ half_normal(3,0.25)
    ), data=dat_list , chains=4 , cores=4 )

precis( mBMG_OU )

## R code 14.53
post1 <- extract.samples(mBMG_OU)
post0 <- extract.samples(mB_OU)

plot( NULL , xlim=c(0,max(dat_list$Dmat)) , ylim=c(0,1.5) ,
    xlab="phylogenetic distance" , ylab="covariance" )

# prior mean and 89% interval
n <- 30
eta <- abs(rnorm(n,1,0.25))
rho <- abs(rnorm(n,3,0.25))
for ( i in 1:n ) curve( eta[i]*exp(-rho[i]*x) , add=TRUE , col=grau(0.4) , lwd=2 )

# posterior
for ( i in 1:n )
    curve( post0$etasq[i]*exp(-post$rho[i]*x) , add=TRUE , col=4 , lwd=2 )

for ( i in 1:n )
    curve( post1$etasq[i]*exp(-post$rho[i]*x) , add=TRUE , col=2 , lwd=2 )



# compare bG in mBMG and mBMG_OU

post1 <- extract.samples(mBMG_OU)
post0 <- extract.samples(mBMG)

dens( post0$bG , lwd=4 , col=1 , xlab="effect of group size" , xlim=c(-0.05,0.22) )
dens( post1$bG , lwd=4 , col=2 , add=TRUE )
abline(v=0,lty=3)
