# week 9

################################################
# 1
# analyze success as function of age

library(rethinking)
data(Achehunting)
d <- Achehunting

dat <- list(
    S = ifelse(d$kg.meat>0,1,0),
    A = standardize(log(d$age))
)

# log age model
m1a <- ulam(
    alist(
        S ~ bernoulli(p),
        logit(p) <- a + bA*A,
        a ~ normal(0,1),
        bA ~ normal(0,0.5)
    ) , data=dat , chains=4 , cores=4 , log_lik=TRUE )

# round age to nearest decade
Ar <- round( d$age / 10 ) * 10
Aseq <- c(10,20,30,40,50,60,70,80)
SA <- sapply( Aseq , function(a) mean(dat$S[Ar==a]) )
SAsd <- sapply( Aseq , function(a) sqrt( var(dat$S[Ar==a])/sum(dat$S[Ar==a]) ) )
plot( Aseq , SA , ylim=c(0,1) , xlim=c(0,90) , type="b" , lwd=3 , col=2 , xlab="age (years)" , ylab="proportion successful" )
for ( i in 1:length(Aseq) ) lines( rep(Aseq[i],2) , SA[i] + SAsd[i]*c(-1,1) , col=2 )

post <- extract.samples(m1a)
for ( i in 1:20 ) curve( inv_logit( post$a[i] + post$bA[i]*( log(x)-3.79 )/0.335 ) , add=TRUE , lwd=2 , col=grau() , from=1 )

# try something else
# normalize age
dat$A2 <- d$age / 80

curve( 0.6*(1-exp(-3*x))*exp(-0.5*x) , from=0 , to=1 , ylim=c(0,1) )

m1b <- ulam(
    alist(
        S ~ bernoulli(p),
        p <- a*(1-exp(-b1*A2))*exp(-b2*A2),
        a ~ beta(4,4),
        c(b1,b2) ~ exponential(2)
    ) , data=dat , chains=4 , cores=4 , log_lik=TRUE )

# round age to nearest decade
Ar <- round( d$age / 10 ) * 10
Aseq <- c(10,20,30,40,50,60,70,80)
SA <- sapply( Aseq , function(a) mean(dat$S[Ar==a]) )
plot( Aseq , SA , ylim=c(0,1) , xlim=c(0,90) , type="b" , lwd=3 , col=2 , xlab="age (years)" , ylab="proportion successful" )

post <- extract.samples(m1b)
for ( i in 1:20 ) with( post , 
    curve( a[i]*(1-exp(-b1[i]*x/80))*exp(-b2[i]*x/80) , add=TRUE , lwd=2 , col=grau() , from=1 ) )

# now with elasticity
curve( 0.6*exp(-0.5*x)*(1-exp(-3*x))^2 , from=0 , to=1 , ylim=c(0,1) )

m1c <- ulam(
    alist(
        S ~ bernoulli(p),
        p <- a*exp(-b2*A2)*(1-exp(-b1*A2))^g,
        a ~ beta(4,4),
        g ~ exponential(0.5),
        c(b1,b2) ~ exponential(2)
    ) , data=dat , chains=4 , cores=4 , log_lik=TRUE )

post <- extract.samples(m1c)
for ( i in 1:20 ) with( post , 
    curve( a[i]*exp(-b2[i]*x/80)*(1-exp(-b1[i]*x/80))^g[i] , add=TRUE , lwd=2 , col=grau() , from=1 ) )

compare( m1a , m1b , m1c , func=PSIS )

################################################
# 2
# add individual hunter varying effects

dat$H <- as.integer(as.factor(d$id))
dat$NH <- max(dat$H)

check_index(dat$H)

m2 <- ulam(
    alist(
        S ~ bernoulli(p),
        p <- a*exp(-b2H[H]*A2)*(1-exp(-b1H[H]*A2))^g,
        # centered varying effects
        transpars> vector[NH]:b1H <<- exp(b1+V[1:NH,1]),
        transpars> vector[NH]:b2H <<- exp(b2+V[1:NH,2]),
        # non-centered varying effects
        transpars> matrix[NH,2]:V <-
            compose_noncentered( sigma_H , L_Rho_H , Z ),
        matrix[2,NH]:Z ~ normal( 0 , 1 ),
        cholesky_factor_corr[2]:L_Rho_H ~ lkj_corr_cholesky( 4 ),
        vector[2]:sigma_H ~ exponential(1),
        # fixed priors
        a ~ beta(4,4),
        g ~ exponential(0.5),
        c(b1,b2) ~ normal(0,0.5),
        gq> matrix[2,2]:Rho_H <<- Chol_to_Corr( L_Rho_H )
    ) , data=dat , chains=4 , cores=4 , iter=4000 )

precis(m2,3,pars=c("a","g","b1","b2","sigma_H","Rho_H"))
post <- extract.samples(m2)

blank(bty="n",ex=2)

# round age to nearest decade
Ar <- round( d$age / 10 ) * 10
Aseq <- c(10,20,30,40,50,60,70,80)
SA <- sapply( Aseq , function(a) mean(dat$S[Ar==a]) )

par(mfrow=c(4,4))
for (k in 1:16 ) {

    plot( Aseq , SA , ylim=c(0,1) , xlim=c(0,90) , type="b" , lwd=3 , col=0 , xlab="age (years)" , ylab="proportion successful" )

    # plot nj random hunters
    cols <- c(2,4,5,6)
    nj <- 2
    hseq <- sample(1:dat$NH,size=nj)
    for ( j in hseq ) 
    for ( i in 1:10 ) with( post , 
        curve( a[i]*exp(-b2H[i,j]*x/80)*(1-exp(-b1H[i,j]*x/80))^g[i] , add=TRUE , lwd=2 , from=1 , col=cols[which(hseq==j)] ) )
    points( Aseq , SA , type='b' , col="white" , lwd=6 )
    points( Aseq , SA , type='b' , col=1 , lwd=3 )

}#k

# post <- extract.samples(m3)
# post <- extract.samples(m2L)
# a <- 1
# compute variation across age for each hunter
vHA <- rep(NA,dat$NH)
mHA <- rep(NA,dat$NH) # mean success
Aseq <- 10:80
for ( i in 1:dat$NH ) {
    # compute variation across age for each sample from posterior
    # then average across samples
    # v has margins [samples,ages]
    v <- sapply( Aseq , function(x) 
        with( post , 
            a*exp(-b2H[,i]*x/80)*(1-exp(-b1H[,i]*x/80))^g ) )
    # now average variation across ages
    vHA[i] <- mean( apply( v , 1 , var ) )
    mHA[i] <- mean( apply( v , 1 , mean ) )
}#i

# compute variation across individuals averaged by age
# variation across all hunters at each age, then average over ages
vAH <- rep(NA,length(Aseq))
for ( j in 1:length(Aseq) ) {
    # variation at age j across all hunters
    # v has margins [samples,hunters]
    v <- sapply( 1:dat$NH , function(i) 
        with( post , 
            a*exp(-b2H[,i]*Aseq[j]/80)*(1-exp(-b1H[,i]*Aseq[j]/80))^g ) )
    # average variance across individuals
    vAH[j] <- mean( apply( v , 1 , var ) )
}#j

# average across hunters (of variation across age)
mean(vHA)
# average across age (of variation across hunters)
mean(vAH)

plot( Aseq , vAH , xlab="age (years)" , ylab="variation across hunters" , type="l" , lwd=3 , col=2 )

plot( sort(vHA) , lwd=3 , col=2 , xlab="individual hunter" , ylab="variation across age" )

plot( mHA , vHA , lwd=3 , col=2 , xlab="mean success across age" , ylab="variation across age" )

################################################
# 3
# duration and impute missing values

f <- function(x,a=0.9,b1=8,b2=0.8,g=6) a*exp(-b2*x)*(1-exp(-b1*x))^g

# use poisson prob > 0, 1-exp(-lambda), and make lambda = L^d * f(age)
curve( 1-exp(-f(x)) , from=0 , to=1 , ylim=c(0,1) )
for ( L in c(0.1,0.5,1,10,100) )
curve( 1-exp(-(L^0.5)*f(x)) , add=TRUE , col=2 )

L = d$hours / max(d$hours,na.rm=TRUE)
dat$log_L <- log(L)

ccidx <- which(!is.na(dat$log_L))
datcc <- list(
    S = dat$S[ccidx],
    A2 = dat$A2[ccidx],
    H = dat$H[ccidx],
    NH = dat$NH,
    log_L = dat$log_L[ccidx] 
)

flist3 <- alist(
    S ~ bernoulli(p),
    p <- 1-exp( -exp(lambda*log_L) * f ),
    f <- exp(-b2H[H]*A2)*(1-exp(-b1H[H]*A2))^g,
    # centered varying effects
    transpars> vector[NH]:b1H <<- exp(b1+V[1:NH,1]),
    transpars> vector[NH]:b2H <<- exp(b2+V[1:NH,2]),
    # non-centered varying effects
    transpars> matrix[NH,2]:V <-
        compose_noncentered( sigma_H , L_Rho_H , Z ),
    matrix[2,NH]:Z ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_H ~ lkj_corr_cholesky( 4 ),
    vector[2]:sigma_H ~ exponential(1),
    # duration prior
    log_L ~ normal(muL,sigmaL),
    muL ~ normal(-1,0.25),
    sigmaL ~ exponential(2),
    # fixed priors
    lambda ~ exponential(1),
    g ~ exponential(0.5),
    c(b1,b2) ~ normal(0,0.5),
    gq> matrix[2,2]:Rho_H <<- Chol_to_Corr( L_Rho_H )
)

m3cc <- ulam( flist3 , data=datcc , chains=4 , cores=4 , iter=1000 )

m3 <- ulam( flist3 , data=dat , chains=4 , cores=4 , warmup=1000 , iter=4000 )

precis(m3,3,pars=c("b1","b2","g","lambda","sigma_H"))
precis(m3cc,3,pars=c("b1","b2","g","lambda","sigma_H"))

post <- extract.samples(m3)

f <- function(x,a=1,b1=8,b2=0.8,g=6) a*exp(-b2*x)*(1-exp(-b1*x))^g

plot( NULL , xlim=c(-4,0) , ylim=c(0,1) , xlab="log trip duration" , ylab="probability success" )

for ( i in 1:20 )
with( post , 
curve( 1-exp( -exp(x)^lambda[i] * f(0.5,b1=exp(b1[i]),b2=exp(b2[i]),g=g[i]) ) , add=TRUE , lwd=3 , col=2 )
)

# version of model 2 that uses same success function, for comparison
dat0 <- dat
dat0$log_L <- NULL
m2L <- ulam(
    alist(
        S ~ bernoulli(p),
        p <- 1 - exp( -exp(-b2H[H]*A2)*(1-exp(-b1H[H]*A2))^g ),
        # centered varying effects
        transpars> vector[NH]:b1H <<- exp(b1+V[1:NH,1]),
        transpars> vector[NH]:b2H <<- exp(b2+V[1:NH,2]),
        # non-centered varying effects
        transpars> matrix[NH,2]:V <-
            compose_noncentered( sigma_H , L_Rho_H , Z ),
        matrix[2,NH]:Z ~ normal( 0 , 1 ),
        cholesky_factor_corr[2]:L_Rho_H ~ lkj_corr_cholesky( 4 ),
        vector[2]:sigma_H ~ exponential(1),
        # fixed priors
        a ~ beta(4,4),
        g ~ exponential(0.5),
        c(b1,b2) ~ normal(0,0.5),
        gq> matrix[2,2]:Rho_H <<- Chol_to_Corr( L_Rho_H )
    ) , data=dat0 , chains=4 , cores=4 , warmup=1000 , iter=4000 )

################################################
# 4
# model harvest size

dat$M <- d$kg.meat / mean(d$kg.meat[dat$S==1])
datcc$M <- dat$M[ccidx]

dens( log(dat$M[dat$S==1]) , lwd=3 , col=2 , xlab="log(M)|M>0" )

flist4 <- alist(
    M|S==1 ~ dlnorm(mu,tau),
    mu <- a + lambda*log_L + log(f),
    f <- exp(-b2H[H]*A2)*(1-exp(-b1H[H]*A2))^g,
    tau ~ exponential(1),
    # centered varying effects
    transpars> vector[NH]:b1H <<- exp(b1+V[1:NH,1]),
    transpars> vector[NH]:b2H <<- exp(b2+V[1:NH,2]),
    # non-centered varying effects
    transpars> matrix[NH,2]:V <-
        compose_noncentered( sigma_H , L_Rho_H , Z ),
    matrix[2,NH]:Z ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_H ~ lkj_corr_cholesky( 4 ),
    vector[2]:sigma_H ~ exponential(1),
    # duration prior
    log_L ~ normal(muL,sigmaL),
    muL ~ normal(0,0.25),
    sigmaL ~ exponential(2),
    # fixed priors
    a ~ normal(0,0.5),
    lambda ~ exponential(1),
    g ~ exponential(0.5),
    c(b1,b2) ~ normal(0,0.5),
    gq> matrix[2,2]:Rho_H <<- Chol_to_Corr( L_Rho_H )
)

m4cc <- ulam( flist4 , data=datcc , chains=4 , cores=4 , iter=1000 )

m4 <- ulam( flist4 , data=dat , chains=4 , cores=4 , warmup=1000 , iter=4000 )

precis(m4,3,pars=c("sigma_H"))

# plot hunters against empirical means at each age
blank(bty="n",ex=2)

# round age to nearest decade
Ar <- round( d$age / 10 ) * 10
Aseq <- c(10,20,30,40,50,60,70,80)
MA <- sapply( Aseq , function(a) mean(dat$M[dat$S==1 & Ar==a]) )

post <- extract.samples(m4)
LL <- log(0.45)
par(mfrow=c(3,3))
for (k in 1:9 ) {

    plot( Aseq , MA , ylim=c(0,1.3) , xlim=c(0,90) , type="b" , lwd=3 , col=0 , xlab="age (years)" , ylab="expected harvest" )

    # plot nj random hunters
    cols <- c(2,4,5,6)
    nj <- 2
    hseq <- sample(1:dat$NH,size=nj)
    for ( j in hseq ) 
    for ( i in 1:10 ) with( post , 
        curve( 
                exp( 
                    a[i] + lambda[i]*LL +   
                    (-b2H[i,j]*x/80) + g[i]*log(1-exp(-b1H[i,j]*x/80))
                    + tau[i]^2/2 
                ), 
            add=TRUE , lwd=2 , from=1 , col=cols[which(hseq==j)] ) )
    points( Aseq , SA , type='b' , col="white" , lwd=6 )
    points( Aseq , SA , type='b' , col=1 , lwd=3 )

}#k
