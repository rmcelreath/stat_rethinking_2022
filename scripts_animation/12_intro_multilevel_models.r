# lecture 12 - intro to multilevel models
# random intercepts and partial pooling

library(rethinking)
library(animation)
library(ellipse)

# introduce with story & individual example from previous lecture

library(rethinking)
data(Trolley)
d <- Trolley

blank(bty="n")

# stories
S <- as.integer( d$story )
mu_story <- sapply( 1:max(S) , function(i) mean(d$response[S==i]) )

plot( 1:12 , mu_story , lwd=4 , col=2 , ylim=c(1,7) , xlab="Story" , ylab="Response" )
for ( i in 1:12 ) {
    pi <- PI(d$response[S==i],0.5)
    lines( c(i,i) , pi , lwd=12 , col=col.alpha(2,0.5) )
}

# individuals

U <- as.integer( d$id )
mu_u <- sapply( 1:max(U) , function(i) mean(d$response[U==i]) )

plot( 1:max(U) , mu_u , lwd=4 , col=2 , ylim=c(1,7) , xlab="Participant" , ylab="Response" , xlim=c(166,330) )

for ( i in 1:max(U) ) {
    pi <- PI(d$response[U==i],0.5)
    lines( c(i,i) , pi , lwd=5 , col=col.alpha(2,0.5) )
}

########################################
# coffee golem
# visit cafes, record waiting time, update posterior

# sim data first
set.seed(7)
n_cafes <- 5
w <- rpois(n_cafes,5)

n_visits <- 20
cafe <- sample(1:n_cafes,size=n_visits,replace=TRUE)
y <- rpois(n_visits,w[cafe])

mc <- "
data{
    int N;
    int n_cafes;
    int y[N];
    int cafe[N];
}
parameters{
    vector[n_cafes] z;
    real a_bar;
    real<lower=0> sigma;
}
transformed parameters{
    vector[n_cafes] a;
    a = a_bar + z*sigma;
}
model{
    vector[N] lambda;
    sigma ~ exponential( 10 );
    a_bar ~ normal( 2 , 0.5 );
    z ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        lambda[i] = a[cafe[i]];
        lambda[i] = exp(lambda[i]);
    }
    y ~ poisson( lambda );
}
"

# blank(bty="n",w=1.7)

# precompile model
file <- cmdstanr_model_write(mc)
ma <- cmdstan_model(file, compile = TRUE )

n_cafes <- 5
w <- c(3,15,3,1,3)

set.seed(3)
n_visits <- 20
cafe <- c( rep(1,1) , rep(2,1) , 1 , rep(c(1,3,4,5),times=4) , 2 )
y <- rpois(n_visits,w[cafe]) + 1
y[2] <- 18

# loop
ani.record(reset=TRUE)
par(mfrow=c(2,3))
xlims <- c(0,15)
ylims <- c(0,0.5)
for ( f in 0:n_visits ) {

    if ( f > 0 ) {
        m <- f
        dat <- list( N=m , n_cafes=n_cafes , y=y[1:m] , cafe=cafe[1:m] )
        max <- ma$sample(data = dat, chains = 4, 
                    parallel_chains = 4, iter_sampling = 1000, iter_warmup = 1000, 
                    adapt_delta = 0.9, max_treedepth = 15, 
                    save_warmup = FALSE , refresh=0 )
        maxx <- rstan::read_stan_csv(max$output_files())
        post <- extract.samples(maxx)
    } else {
        post$a_bar <- rnorm(4000,2,0.5)
        post$sigma <- rexp(4000,1)
        post$a <- replicate(n_cafes,rnorm(4000,post$a_bar,post$sigma))
    }

    # population distribution
    if ( f == 0  )
        curve( dlnorm(x,mean(post$a_bar),mean(post$sigma)) , from=0, to=20 , lwd=4, col=4 , xlab="waiting time" , ylab="" , ylim=ylims )
    else {
        curve( dlnorm(x,mean(post$a_bar),mean(post$sigma)) , from=0, to=20 , lwd=4, col="white" , xlab="waiting time" , ylab="" , ylim=ylims )
        curve( dlnorm(x,mean(post_old$a_bar),mean(post_old$sigma)) , from=0, to=20 , lwd=4, col=grau() , add=TRUE )
        curve( dlnorm(x,mean(post$a_bar),mean(post$sigma)) , from=0, to=20 , lwd=8, col="white" , add=TRUE )
        curve( dlnorm(x,mean(post$a_bar),mean(post$sigma)) , from=0, to=20 , lwd=4, col=4 , add=TRUE )
    }

    mtext("Population of cafes")

    for ( j in 1:n_cafes ) {
        xlwd <- 4
        if ( f==0 )
            curve( dlnorm(x,mean(post$a[,j]),sd(post$a[,j])) , from=0, to=20 , lwd=4, col=2 , xlab="waiting time" , ylab="" , ylim=ylims )
        else {
            if ( j==cafe[f] ) xlwd <- 8
            curve( dlnorm(x,mean(post$a[,j]),sd(post$a[,j])) , from=0, to=20 , lwd=4, col="white" , xlab="waiting time" , ylab="" , ylim=ylims )
            curve( dlnorm(x,mean(post_old$a[,j]),sd(post_old$a[,j])) , from=0, to=20 , lwd=4, col=grau() , add=TRUE )
            curve( dlnorm(x,mean(post$a[,j]),sd(post$a[,j])) , from=0, to=20 , lwd=8, col="white" , add=TRUE )
            curve( dlnorm(x,mean(post$a[,j]),sd(post$a[,j])) , from=0, to=20 , lwd=4, col=2 , add=TRUE )
            if ( j==cafe[f] ) lines( c(y[f],y[f]) , c(0,10) , lwd=4 , col=1 )
        }
        mtext(concat( sum(cafe[1:f]==j) , " visits" ))
        
    }

    post_old <- post

    ani.record()
}

oopts = ani.options(interval = 1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 8 -loop 0 frame*.png a_out.gif
# convert -delay 10 a_out.gif a_out.gif


########################################
# tadpoles 

library(rethinking) 
data(reedfrogs)
d <- reedfrogs 

# just data

# display raw proportions surviving in each tank
plot( d$propsurv , ylim=c(0,1) , pch=16 , xaxt="n" ,
    xlab="tank" , ylab="proportion survival" , col=1 , yaxt="n" )
axis( 1 , at=c(1,16,32,48) , labels=c(1,16,32,48) )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
abline( h=mean(d$propsurv) , lwd=3 , col=1 , lty=3 )
# draw vertical dividers between tank densities
abline( v=16.5 , lwd=0.5 )
abline( v=32.5 , lwd=0.5 )
text( 8 , 0 , "small tanks" )
text( 16+8 , 0 , "medium tanks" )
text( 32+8 , 0 , "large tanks" )

# show priors in multilevel model
curve( dnorm(x,0,1.5) , from=-6 , to=6 , lwd=4 , col=2 , xlab="a_bar" , ylab="Density" )

curve( dexp(x,1) , from=0 , to=5 , lwd=4 , col=2 , xlab="sigma" , ylab="Density" )

abar <- rnorm(1e5,0,1.5)
s <- rexp(1e5,1)
a <- rnorm(1e5,abar,s)
dens( a , xlim=c(-6,6) , lwd=4 , col=2 , xlab="a_j" , ylab="Density" , adj=1 )

library(rethinking) 
data(reedfrogs)
d <- reedfrogs 
d$tank <- 1:nrow(d)
dat <- list(
    S = d$surv,
    D = d$density,
    T = d$tank )

mST <- ulam( 
    alist(
        S ~ dbinom( D , p ) ,
        logit(p) <- a[T] ,
        a[T] ~ dnorm( a_bar , sigma ) , 
        a_bar ~ dnorm( 0 , 1.5 ) ,
        sigma ~ dexp( 1 )
    ), data=dat , chains=4 , log_lik=TRUE )

mSTnomem <- ulam( 
    alist(
        S ~ dbinom( D , p ) ,
        logit(p) <- a[T] ,
        a[T] ~ dnorm( a_bar , 1 ) , 
        a_bar ~ dnorm( 0 , 1.5 )
    ), data=dat , chains=4 , log_lik=TRUE )

compare( mST , mSTnomem , func=WAIC )

# display partial pooling estimates colored by size/predation
plot( d$propsurv , ylim=c(0,1) , pch=16 , xaxt="n" ,
    xlab="tank" , ylab="proportion survival" , col=1 , yaxt="n" , cex=1.2 )
axis( 1 , at=c(1,16,32,48) , labels=c(1,16,32,48) )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
#abline( h=mean(d$propsurv) , lwd=3 , col=1 , lty=3 )
# draw vertical dividers between tank densities
abline( v=16.5 , lwd=0.5 )
abline( v=32.5 , lwd=0.5 )
text( 8 , 0 , "small tanks" )
text( 16+8 , 0 , "medium tanks" )
text( 32+8 , 0 , "large tanks" )

post <- extract.samples(mST)
a_est <- apply(inv_logit(post$a),2,mean)

# cols <- rep(2,48)
cols <- ifelse(d$pred=="pred",2,4)

for ( i in 1:48 ) {
    pi <- PI( inv_logit(post$a[,i]) , 0.89 )
    lines( c(i,i) , pi , lwd=8 , col=col.alpha(cols[i],0.5) )
}
points( 1:48 , a_est , lwd=3 , col=cols )

# pred model

dat$P <- ifelse(d$pred=="pred",1,0)
mSTP <- ulam( 
    alist(
        S ~ dbinom( D , p ) ,
        logit(p) <- a[T] + bP*P ,
        bP ~ dnorm( 0 , 0.5 ),
        a[T] ~ dnorm( a_bar , sigma ) , 
        a_bar ~ dnorm( 0 , 1.5 ) ,
        sigma ~ dexp( 1 )
    ), data=dat , chains=4 , log_lik=TRUE )

post <- extract.samples(mSTP)
dens( post$bP , lwd=4 , col=2 , xlab="bP (effect of predators)" )


# compare estimates from the two models
p1 <- link(mST)
p2 <- link(mSTP)

# cols <- rep(2,48)
cols <- ifelse(d$pred=="pred",2,4)

plot( apply(p1,2,mean) , apply(p2,2,mean) , lwd=4 , col=cols , xlab="prob survive (model without predators)" , ylab="prob survive (model with predators)" )
abline(a=0,b=1,lty=3)

for ( i in 1:48 ) {
    pi <- PI( p[,i] , 0.89 )
    lines( c(i,i) , pi , lwd=8 , col=col.alpha(cols[i],0.5) )
}
points( 1:48 , apply(p,2,mean) , lwd=3 , col=cols )

# sigmas
post1 <- extract.samples(mST)
post2 <- extract.samples(mSTP)
dens( post1$sigma , xlab="sigma" , lwd=4, col=4 , xlim=c(0,2.5) , ylim=c(0,2.4) )
dens( post2$sigma , xlab="sigma" , lwd=4, col=2 , add=TRUE )

compare( mST , mSTP , func=PSIS )

###################################
# try to animate the tank plot
# illustrate overfitting perspective on varying effects

dat$sigma <- 1.5
mx <- ulam(
    alist(
        S ~ dbinom( N , p ) , 
        logit(p) <- a[tank] , 
        a[tank] ~ dnorm( a_bar , sigma ),
        a_bar ~ normal(0,1.5)
    ), data=dat , chains=4 , cores=4 , sample=FALSE , log_lik=TRUE )
mc <- stancode(mx)
dat$sigma <- 0.5
# precompile model
file <- cmdstanr_model_write(mc)
ma <- cmdstan_model(file, compile = TRUE )


# loop #######
ani.record(reset=TRUE)
par(mfrow=c(2,1))

# sigma_seq <- seq(from=log(5),to=log(0.1),len=20)
sigma_seq <- seq(from=log(0.1),to=log(5),len=20)
psis_list <- matrix(NA,nrow=length(sigma_seq),ncol=4)
for ( f in 1:length(sigma_seq) ) {

    dat$sigma <- exp(sigma_seq[f])
    max <- ma$sample(data = dat, chains = 4, 
            parallel_chains = 4, iter_sampling = 1000, iter_warmup = 500, 
            adapt_delta = 0.9, max_treedepth = 15, 
            save_warmup = FALSE , refresh=0 )

    maxx <- rstan::read_stan_csv(max$output_files())

    post <- extract.samples(maxx)

    if ( sigma_seq[1] < sigma_seq[2] )
        psis_list[f,1:4] <- as.numeric(PSIS(post))

    # compute mean intercept for each tank
    # also transform to probability with logistic
    d$propsurv.est <- inv_logit( apply( post$a , 2 , mean ) )

    # display raw proportions surviving in each tank
    plot( d$propsurv , ylim=c(0,1) , pch=16 , xaxt="n" ,
        xlab="tank" , ylab="proportion survival" , col=1 , yaxt="n" )
    axis( 1 , at=c(1,16,32,48) , labels=c(1,16,32,48) )
    axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
    abline( h=mean(d$propsurv) , lwd=3 , col=1 , lty=3 )

    # overlay posterior means
    points( d$propsurv.est , lwd=3 , col=2 , cex=1.2 )
        pi <- apply( post$p , 2 , PI , 0.89 )
        for ( j in 1:48 ) lines( c(j,j) , pi[,j] , lwd=12, col=col.alpha(2,0.5) )
        pi <- apply( post$p , 2 , PI , 0.5 )
        for ( j in 1:48 ) lines( c(j,j) , pi[,j] , lwd=6, col=2 )

    # data again
    points( 1:48 , d$propsurv , pch=16 , col=1 )

    # mark posterior mean probability across tanks
    abline( h=mean(inv_logit(post$a_bar)) , lty=3 , col=2 , lwd=3 )

    # draw vertical dividers between tank densities
    abline( v=16.5 , lwd=0.5 )
    abline( v=32.5 , lwd=0.5 )
    text( 8 , 0 , "small tanks" )
    text( 16+8 , 0 , "medium tanks" )
    text( 32+8 , 0 , "large tanks" )

    mtext(concat("sigma=",round( exp(sigma_seq[f]) ,2)))

    # posterior random effects
    if ( FALSE ) {
    ymax <- dnorm(mean(post$a_bar),mean(post$a_bar),dat$sigma)
    curve( dnorm(x,mean(post$a_bar),dat$sigma) , from=-5 , to=7 , lwd=4 , col=grau() , xlab="log-odds survival" , n=200 , ylab="Density" , xlim=c(-4.8,6.8) , ylim=c(0,ymax) )
    for ( j in 1:nrow(d) ) {
        mlo <- mean(post$a[,j])
        yy <- dnorm(mlo,mean(post$a_bar),dat$sigma)
        #pi <- apply( post$p , 2 , PI , 0.5 )
        #w <- abs(pi[1,j] - pi[2,j])
        lines( c(mlo,mlo) , c(0, yy ) , lwd=5 , col=2)
        #points( mlo , ymax/2 ,lwd=3 , col=2)
    }
    curve( dnorm(x,mean(post$a_bar),dat$sigma) , add=TRUE , lwd=10 , col="white" , n=200)
    lines( rep(mean(post$a_bar),2) , c(0,ymax) , col="white" , lty=1 , lwd=4 )
    lines( rep(mean(post$a_bar),2) , c(0,ymax) , col=2 , lty=3 , lwd=3 )
    curve( dnorm(x,mean(post$a_bar),dat$sigma) , add=TRUE , lwd=6 , col=1 , n=200)
    }

    # PSIS graph
    if ( TRUE ) {
        rev <- FALSE
        if ( sigma_seq[1] > sigma_seq[2] ) rev <- TRUE
        nn <- sum( !is.na(psis_list[,1]) )

        xseq <- 1:nn
        if ( rev==TRUE ) xseq <- nn:1
        plot( sigma_seq[xseq] , psis_list[xseq,1] , type="b" , lwd=4 , col=2 , xlim=range(sigma_seq) , xaxt="n" , xlab="sigma" , ylab="PSIS score" , ylim=c(190,650) )
        axis(1,at=sigma_seq,labels=round(exp(sigma_seq),2))

        # se of psis
        for ( j in 1:nn ) {
            jj <- j
            if ( rev==TRUE ) jj <- nrow(psis_list) - j
            lines( rep(sigma_seq[j],2) , psis_list[jj,1]+c(-1,1)*psis_list[jj,4] , lwd=10 , col=col.alpha(2,0.5) )
        }

        # highlight current sigma
        points( sigma_seq[f] , psis_list[f,1] , lwd=12 , col=col.alpha(2,0.5) )
        points( sigma_seq[f] , psis_list[f,1] , lwd=8 , col=2 )
    }

    # record frame
    ani.record()

    if ( f==1 | f==length(sigma_seq) ) for ( j in 1:10 ) ani.record()

}

oopts = ani.options(interval = 0.05)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 8 -loop 0 frame*.png a_out.gif
# convert -delay 10 a_out.gif a_out.gif

plot( sigma_seq , psis_list[,1] , type="b" , lwd=4 , col=2 , xlim=range(sigma_seq) , xaxt="n" , xlab="sigma" , ylab="PSIS score" , ylim=c(190,650) )
axis(1,at=sigma_seq,labels=round(exp(sigma_seq),2))
# overlay sigma density from m13.2
post132 <- extract.samples(m13.2)
dplot <- density( log(post132$sigma) , adj=0.5 )
lines( dplot$x , normalize(dplot$y)*400 + 200 , lwd=4 )

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
    ), data=dat , chains=4 , cores=4 )

precis(mRXma2,2)
