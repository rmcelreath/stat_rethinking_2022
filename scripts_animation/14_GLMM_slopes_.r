# lecture 14
# varying slopes, non-centered covariance, group centering and fixed effects

library(rethinking)
library(animation)
library(ellipse)

########################################
# UCBadmit learning covariance

data(UCBadmit)
d <- UCBadmit

dat <- list( 
    A=d$admit , 
    N=d$applications,
    D=as.integer(d$dept) , 
    F=ifelse(d$applicant.gender=="female",1,0) )

m0 <- ulam(
    alist(
        A ~ binomial(N,p),
        logit(p) <- abar + a[D] + (bbar + b[D])*F,
        c(a,b)[D] ~ multi_normal(0,Rho,Sigma),
        abar ~ normal(0,1),
        bbar ~ normal(0,0.5),
        Rho ~ lkj_corr(4),
        Sigma ~ half_normal(0,1)
    ) , data=dat , chains=4 , cores=4 , sample=TRUE )

post <- extract.samples(m0)
plot( apply(post$a,2,mean) , apply(post$b,2,mean) , lwd=3 , col=2 , pch=levels(d$dept) )
abline(a=0,b=1,lty=3)
abline(h=0,lty=3)
abline(v=0,lty=3)

mc <- "
data{
    int n;
    int A[n];
    int N[n];
    int F[n];
    int D[n];
    int onlypriors;
}
parameters{
    vector[6] b;
    vector[6] a;
    real abar;
    real bbar;
    corr_matrix[2] Rho;
    vector<lower=0>[2] Sigma;
}
model{
    vector[n] p;
    Sigma ~ normal( 0 , 1 );
    Rho ~ lkj_corr( 4 );
    bbar ~ normal( 0 , 0.5 );
    abar ~ normal( 0 , 1 );
    {
    vector[2] YY[6];
    for ( j in 1:6 ) YY[j] = [ a[j] , b[j] ]';
    YY ~ multi_normal( rep_vector(0,2) , quad_form_diag(Rho , Sigma) );
    }
    if ( onlypriors==0 ) { 
        for ( i in 1:n ) {
            p[i] = abar + a[D[i]] + (bbar + b[D[i]]) * F[i];
            p[i] = inv_logit(p[i]);
        }
        A ~ binomial( N , p );
    }
}
"

# precompile model
file <- cmdstanr_model_write(mc)
ma <- cmdstan_model(file, compile = TRUE )

blank(bty="n",w=1.6)

# simulated dataset with exaggerated covariance in population
set.seed(1)
#a <- seq(from=1.4,to=-1.4,len=5) # log-odds m admit
#b <- seq(from=2,to=-2,len=5) # diff log-odds f admit
a <- c(1.4,0.7,0,-1,-2)
b <- c(2,1,0,-1.5,-1)
D <- rep(1:5,each=2)
F <- rep(0:1,times=5)
N <- rep(c(100,100),each=5)
A <- rbinom(10,size=N,p=inv_logit(a[D]+b[D]*F))

dat <- list( A=A , N=N , D=D , F=F )

# loop
ani.record(reset=TRUE)
par(mfrow=c(2,3))
xlims <- c(0,15)
ylims <- c(0,0.5)
xlist <- c(1,2,9,10,3,4,7,8,11,12) # for real data
xlist <- 1:10
Sigma_prev <- list()
Mu_prev <- list()
for ( f in 0:length(xlist) ) {

    if ( f > 0 ) {

        xx <- xlist[1:f]
        datx <- list( n=length(xx) , A=as.array(dat$A[xx]) , N=as.array(dat$N[xx]) , D=as.array(dat$D[xx]) , F=as.array(dat$F[xx]) , onlypriors=0 )
        max <- ma$sample(data = datx, chains = 4, 
                    parallel_chains = 4, iter_sampling = 1000, iter_warmup = 1000, 
                    adapt_delta = 0.9, max_treedepth = 15, 
                    save_warmup = FALSE , refresh=0 )
        maxx <- rstan::read_stan_csv(max$output_files())
        post <- extract.samples(maxx)

    } else {
        
        xx <- xlist[1:2]
        datx <- list( n=length(xx) , A=dat$A[xx] , N=dat$N[xx] , D=dat$D[xx] , F=dat$F[xx] , onlypriors=1 )
        max <- ma$sample(data = datx, chains = 4, 
                    parallel_chains = 4, iter_sampling = 1000, iter_warmup = 1000, 
                    adapt_delta = 0.9, max_treedepth = 15, 
                    save_warmup = FALSE , refresh=0 )
        maxx <- rstan::read_stan_csv(max$output_files())
        post <- extract.samples(maxx)

    }

    # population distribution
    
    Mu <- c( mean(post$abar) , mean(post$bbar) )
    popMu <- Mu
    rho_est <- mean( post$Rho[,1,2] )
    sa_est <- mean( post$Sigma[,1] )
    sb_est <- mean( post$Sigma[,2] )
    cov_ab <- sa_est*sb_est*rho_est
    Sigma <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )
    plot( NULL , xlim=c(-3,3) , ylim=c(-2,2) , xlab="avg log-odds admit" , ylab="avg diff F-M" )
    #plot( NULL , xlim=c(0,1) , ylim=c(0,1) , xlab="prob admit m" , ylab="prob admit w" , xaxt="n" )
    #axis(1,at=c(0,0.5,1),labels=c(0,0.5,1))
    abline( h=0 , lty=3 )
    abline( v=0 , lty=3 )
    # draw shadow
    if ( f > 0 ) {
        for ( l in c(0.25,0.5,0.75,0.95) )
            lines( ( ellipse(Sigma_old,centre=Mu_old,level=l,npoints=500) ), col=grau(),lwd=2)
    }
    # draw post
    for ( l in c(0.25,0.5,0.75,0.95) ) {
        if ( f > 0 ) lines( (ellipse(Sigma,centre=Mu,level=l,npoints=500)) ,col="white",lwd=6)
        lines( ( ellipse(Sigma,centre=Mu,level=l,npoints=500) ), col=4,lwd=3)
    }
    # store for drawing shadow next time
    Sigma_old <- Sigma
    Mu_old <- Mu
    

    mtext("Population of departments")

    # 5 departments
    ds <- c(1,2,4,5,6) # for real data
    ds <- 1:5
    for ( j in ds ) {

        xlwd <- 4
        
        aj <- post$abar+post$a[,j]
        bj <- post$bbar+post$b[,j]
        Mu <- c( mean(aj) , mean(bj) )
        sa_est <- sd( aj )
        sb_est <- sd( bj )
        cov_ab <- cov( aj , bj )
        Sigma <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )
        plot( NULL , xlim=c(-3,3) , ylim=c(-2,2) , xlab="a" , ylab="b" )
        #plot( NULL , xlim=c(0,1) , ylim=c(0,1) , xlab="prob admit m" , ylab="diff prob admit w" , xaxt="n" )
        #axis(1,at=c(0,0.5,1),labels=c(0,0.5,1))
        #abline( h=0 , lty=3 )
        #abline( v=0 , lty=3 )
        abline( v=popMu[1] , lty=1 , col=4 , lwd=2 )
        abline( h=popMu[2] , lty=1 , col=4 , lwd=2 )

        if ( f > 0 ) {
            # shadow
            for ( l in c(0.25,0.5,0.75,0.95) )
                lines( (ellipse(Sigma_prev[[j]],centre=Mu_prev[[j]],level=l,npoints=500)) ,col=grau(),lwd=1)
        }

        for ( l in c(0.25,0.5,0.75,0.95) ) {
            if ( f > 0 ) lines( (ellipse(Sigma,centre=Mu,level=l,npoints=500)) ,col="white",lwd=4)
            lines( (ellipse(Sigma,centre=Mu,level=l,npoints=500)) ,col=2,lwd=2)
        }
        
        mtext(concat("Department ",j))

        Sigma_prev[[j]] <- Sigma
        Mu_prev[[j]] <- Mu
        
    }#j

    post_old <- post

    ani.record()
}

oopts = ani.options(interval = 1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 8 -loop 0 frame*.png a_out.gif
# convert -delay 10 a_out.gif a_out.gif

# now show evolution of shrinkage as sample size grows

data(UCBadmit)
d <- UCBadmit

dat <- list( 
    A=d$admit , 
    N=d$applications,
    D=as.integer(d$dept) , 
    F=ifelse(d$applicant.gender=="female",1,0),
    R=d$reject )


# simulated dataset with exaggerated covariance in population
set.seed(1)
#a <- seq(from=1.4,to=-1.4,len=5) # log-odds m admit
#b <- seq(from=2,to=-2,len=5) # diff log-odds f admit
a <- c(1.4,0.7,0,-1,-2)
b <- c(2,1,0,-1.5,-1)
D <- rep(1:5,each=2)
F <- rep(0:1,times=5)
N <- rep(c(100,100),each=5)
A <- rbinom(10,size=N,p=inv_logit(a[D]+b[D]*F))

dat <- list( A=A , N=N , D=D , F=F , R=N-A )


# partial pooling model
mc1 <- "
data{
    int n;
    int N[n];
    int A[n];
    int F[n];
    int D[n];
    int onlypriors;
}
parameters{
    vector[6] b;
    vector[6] a;
    real abar;
    real bbar;
    corr_matrix[2] Rho;
    vector<lower=0>[2] Sigma;
}
model{
    vector[n] p;
    Sigma ~ normal( 0 , 1 );
    Rho ~ lkj_corr( 4 );
    bbar ~ normal( 0 , 0.5 );
    abar ~ normal( 0 , 1 );
    {
    vector[2] YY[6];
    for ( j in 1:6 ) YY[j] = [ a[j] , b[j] ]';
    YY ~ multi_normal( rep_vector(0,2) , quad_form_diag(Rho , Sigma) );
    }
    if ( onlypriors==0 ) { 
        for ( i in 1:n ) {
            p[i] = abar + a[D[i]] + (bbar + b[D[i]]) * F[i];
            p[i] = inv_logit(p[i]);
        }
        A ~ binomial( N , p );
    }
}
"

# precompile model
file <- cmdstanr_model_write(mc1)
ma1 <- cmdstan_model(file, compile = TRUE )

# no pooling model
mc2 <- "
data{
    int n;
    int N[n];
    int A[n];
    int F[n];
    int D[n];
    int onlypriors;
}
parameters{
    vector[6] b;
    vector[6] a;
    real abar;
    real bbar;
}
model{
    vector[n] p;
    bbar ~ normal( 0 , 0.5 );
    abar ~ normal( 0 , 1 );
    a ~ normal( 0 , 1 );
    b ~ normal( 0 , 1 );
    if ( onlypriors==0 ) { 
        for ( i in 1:n ) {
            p[i] = abar + a[D[i]] + (bbar + b[D[i]]) * F[i];
            p[i] = inv_logit(p[i]);
        }
        A ~ binomial( N , p );
    }
}
"

# precompile model
file2 <- cmdstanr_model_write(mc2)
ma2 <- cmdstan_model(file2, compile = TRUE )


ani.record(reset=TRUE)
# seq for proportion of data
pseq <- seq(from=0.1,to=1,len=20)
for ( ps in pseq ) {

    # sub-sample each row
    A <- ceiling(dat$A * ps)
    R <- ceiling(dat$R * ps)
    N <- A + R
    datx <- list( n=10 , A=A , N=N , D=dat$D , F=dat$F , onlypriors=0 )

    # run models
    max1 <- ma1$sample(data = datx, chains = 4, 
                parallel_chains = 4, iter_sampling = 1000, iter_warmup = 1000, 
                adapt_delta = 0.9, max_treedepth = 15, 
                save_warmup = FALSE , refresh=0 )
    maxx1 <- rstan::read_stan_csv(max1$output_files())
    post1 <- extract.samples(maxx1)

    max2 <- ma2$sample(data = datx, chains = 4, 
                parallel_chains = 4, iter_sampling = 1000, iter_warmup = 1000, 
                adapt_delta = 0.9, max_treedepth = 15, 
                save_warmup = FALSE , refresh=0 )
    maxx2 <- rstan::read_stan_csv(max2$output_files())
    post2 <- extract.samples(maxx2)

    # 2 is no pooling, 1 is partial pooling
    a2 <- sapply( 1:5 , function(j) mean( post2$abar + post2$a[,j] ) )
    b2 <- sapply( 1:5 , function(j) mean( post2$bbar + post2$b[,j] ) )
    a1 <- sapply( 1:5 , function(j) mean( post1$abar + post1$a[,j] ) )
    b1 <- sapply( 1:5 , function(j) mean( post1$bbar + post1$b[,j] ) )
    plot( a2 , b2 , lwd=3 , col=1 , xlim=c(-2,2) , ylim=c(-1.5,1.5) , xlab="a" , ylab="b" )

    Mu <- c( mean(post1$abar) , mean(post1$bbar) )
    rho_est <- mean( post1$Rho[,1,2] )
    sa_est <- mean( post1$Sigma[,1] )
    sb_est <- mean( post1$Sigma[,2] )
    cov_ab <- sa_est*sb_est*rho_est
    Sigma <- matrix( c(sa_est^2,cov_ab,cov_ab,sb_est^2) , ncol=2 )
    for ( l in c(0.25,0.5,0.75,0.95) )
        lines( ( ellipse(Sigma,centre=Mu,level=l,npoints=500) ), col=4,lwd=1)

    points( a1 , b1 , lwd=3 , col=2 )
    for ( i in 1:6 ) lines( c(a1[i],a2[i]) , c(b1[i],b2[i]) , lwd=2 , col=2 )

    mtext( concat( "number of applications: ",sum(datx$N) ) )

    ani.record()

}#ps

oopts = ani.options(interval = 0.2)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 8 -loop 0 frame*.png a_out.gif
# convert -delay 10 a_out.gif a_out.gif

#######
# chimpanzees

library(rethinking)
data(chimpanzees)
d <- chimpanzees
dat <- list(
    P = d$pulled_left,
    A = as.integer(d$actor),
    B = as.integer(d$block),
    T = 1L + d$prosoc_left + 2L*d$condition)

m14.2 <- ulam(
alist(
    P ~ bernoulli(p),
    logit(p) <- abar[A] + a[A,T] + bbar[B] + b[B,T],
    
    # adaptive priors
    vector[4]:a[A] ~ multi_normal(0,Rho_A,sigma_A),
    vector[4]:b[B] ~ multi_normal(0,Rho_B,sigma_B),
    abar[A] ~ normal(0,tau_A),
    bbar[B] ~ normal(0,tau_B),

    # fixed priors
    c(tau_A,tau_B) ~ exponential(1),
    sigma_A ~ exponential(1),
    Rho_A ~ dlkjcorr(4),
    sigma_B ~ exponential(1),
    Rho_B ~ dlkjcorr(4)
) , data=dat , chains=4 , cores=4 )

dashboard(m14.2)
precis(m14.2,3)

m14.3 <- ulam(
alist(
    P ~ bernoulli(p),
    logit(p) <- abar[A] + a[A,T] + bbar[B] + b[B,T],

    # adaptive priors - non-centered
    transpars> matrix[A,4]:a <- 
            compose_noncentered( sigma_A , L_Rho_A , zA ),
    transpars> matrix[B,4]:b <- 
            compose_noncentered( sigma_B , L_Rho_B , zB ),
    matrix[4,A]:zA ~ normal( 0 , 1 ),
    matrix[4,B]:zB ~ normal( 0 , 1 ),
    zAbar[A] ~ normal(0,1),
    zBbar[B] ~ normal(0,1),
    transpars> vector[A]:abar <<- zAbar*tau_A,
    transpars> vector[B]:bbar <<- zBbar*tau_B,
    
    # fixed priors
    c(tau_A,tau_B) ~ exponential(1),
    vector[4]:sigma_A ~ exponential(1),
    cholesky_factor_corr[4]:L_Rho_A ~ lkj_corr_cholesky(4),
    vector[4]:sigma_B ~ exponential(1),
    cholesky_factor_corr[4]:L_Rho_B ~ lkj_corr_cholesky(4),

    # compute ordinary correlation matrixes
    gq> matrix[4,4]:Rho_A <<- Chol_to_Corr(L_Rho_A),
    gq> matrix[4,4]:Rho_B <<- Chol_to_Corr(L_Rho_B)

) , data=dat , chains=4 , cores=4 , log_lik=TRUE )

precis(m14.3,3)

precis(m14.3,3,pars="Rho_A")
precis(m14.3,3,pars="Rho_B")

the_pars <- c("abar","a","b","sigma_A","sigma_B","Rho_A","Rho_B")
prc <- precis(m14.2,depth=3,pars=the_pars)
prnc <- precis(m14.3,depth=3,pars=the_pars)

cols <- ifelse( prc[,"n_eff"] > prnc[,"n_eff"] , 2, 4 )
plot( prc[,"n_eff"] , prnc[,"n_eff"] , lwd=3 , col=cols , xlim=c(0,2500) , ylim=c(0,2500) , xlab="effective samples (centered model)" , ylab="effective samples (non-centered model)" )
abline(a=0,b=1,lty=3)

# plot effects

post <- extract.samples(m14.3)

plot( NULL , xlim=c(1,4*7+6) , ylim=c(0,1) , xlab="actor-treatment" , ylab="probability pull left" , xaxt="n" )
for ( i in 1:7 ) for ( j in 1:4 ) {
    aij <- inv_logit( post$abar[,i] + post$a[,i,j] )
    x <- (i-1)*5+j
    lines( c(x,x) , PI(aij) , lwd=8 , col=col.alpha(4,0.5) )
    points( x , mean(aij) , lwd=3 , col=4 )
}

plot( NULL , xlim=c(1,4*6+5) , ylim=c(0,1) , xlab="block-treatment" , ylab="probability pull left" , xaxt="n" )
for ( i in 1:6 ) for ( j in 1:4 ) {
    bij <- inv_logit( post$bbar[,i] + post$b[,i,j] )
    x <- (i-1)*5+j
    lines( c(x,x) , PI(bij) , lwd=8 , col=col.alpha(2,0.5) )
    points( x , mean(bij) , lwd=3 , col=2 )
}

dens( post$tau_A , lwd=4 , col=4 , xlab="standard dev among actors" , xlim=c(0,5) , ylim=c(0,2) )
for ( i in 1:4 ) dens( post$sigma_A[,i] , lwd=4 , col=col.alpha(6,0.5) , add=TRUE )

dens( post$tau_B , lwd=4 , col=2 , xlab="standard dev among blocks" , xlim=c(0,5) , ylim=c(0,3.2) )
for ( i in 1:4 ) dens( post$sigma_B[,i] , lwd=4 , col=col.alpha(6,0.5) , add=TRUE )

# cholesky factor example

# define 2D Gaussian with correlation 0.6
N <- 1e4
sigma1 <- 2
sigma2 <- 0.5
rho <- 0.6

# independent z-scores
z1 <- rnorm( N )
z2 <- rnorm( N )

# use Cholesky to blend in correlation
a1 <- z1 * sigma1
a2 <- ( rho*z1 + sqrt( 1-rho^2 )*z2 )*sigma2

# random correlation matrix
R <- rlkjcorr(1,4)
# Cholesky factor
L <- chol(R)

