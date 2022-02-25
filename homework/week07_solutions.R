# week 7
# varying effects, clusters and features, non-centering

library(rethinking)

# 1
# simple varying intercepts model

data(bangladesh)
d <- bangladesh

dat <- list(
    C = d$use.contraception,
    D = as.integer(d$district) )

m1 <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D],
        vector[61]:a ~ normal(abar,sigma),
        abar ~ normal(0,1),
        sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

# plot estimates
post <- extract.samples(m1)
p <- inv_logit(post$a)
plot( apply(p,2,mean) , xlab="district" , lwd=3 , col=2 , ylim=c(0,1) , ylab="prob use contraception" )
for ( i in 1:61 ) lines( c(i,i) , PI(p[,i]) , lwd=8 , col=col.alpha(2,0.5) )

# show raw proportions - have to skip 54
n <- table(dat$D)
Cn <- xtabs(dat$C ~ dat$D)
pC <- as.numeric( Cn/n )
pC <- c( pC[1:53] , NA , pC[54:60] )
points( pC , lwd=2 )

# add sample size labels
y <- rep(c(0.8,0.85,0.9),len=61)
n <- as.numeric(n)
n <- c( n[1:53] , 0 , n[54:60] )
text( 1:61 , y , n , cex=0.8 )

# 2
# draw DAG for system, design strategy to estimate influence of urban on contraceptive use

# 3
# estimate influence of U on C
# need to adjust for D

dat <- list(
    C = d$use.contraception,
    D = as.integer(d$district),
    U = d$urban,
    A = standardize(d$age.centered),
    K = d$living.children )

# total
m3a <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + bU*U,
        bU ~ normal(0,0.5),
        vector[61]:a ~ normal(abar,sigma),
        abar ~ normal(0,1),
        sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

precis(m3a)

# direct
dat$Kprior <- rep(2,3)
m3b <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D] + bU*U + bA*A + 
                    bK*sum( delta_j[1:K] ),
        c(bU,bA,bK) ~ normal(0,0.5),
        vector[4]: delta_j <<- append_row( 0 , delta ),
        simplex[3]: delta ~ dirichlet( Kprior ),
        vector[61]:a ~ normal(abar,sigma),
        abar ~ normal(0,1),
        sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

precis(m3b)

# varying slopes models

# total
m3c <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D,1] + a[D,2]*U,
        matrix[61,2]:a ~ multi_normal(abar,Rho,sigma),
        vector[2]:abar ~ normal(0,1),
        corr_matrix[2]:Rho ~ lkj_corr(4),
        vector[2]:sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

# non-centered with age
m3cnc <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- abar[1] + a[D,1] + 
                    (abar[2] + a[D,2])*U + 
                    (abar[3] + a[D,3])*A,
        transpars> matrix[61,3]:a <-
            compose_noncentered( sigma , L_Rho , Z ),
        matrix[3,61]:Z ~ normal( 0 , 1 ),
        vector[3]:abar ~ normal(0,1),
        cholesky_factor_corr[3]:L_Rho ~ lkj_corr_cholesky( 4 ),
        vector[3]:sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

# direct
dat$Kprior <- rep(2,3)
m3d <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- a[D,1] + a[D,2]*U + a[D,3]*A +
                    a[D,4]*sum( delta_j[1:K] ),
        matrix[61,4]:a ~ multi_normal(abar,Rho,sigma),
        vector[4]: delta_j <<- append_row( 0 , delta ),
        simplex[3]: delta ~ dirichlet( Kprior ),
        vector[4]:abar ~ normal(0,1),
        corr_matrix[4]:Rho ~ lkj_corr(4),
        vector[4]:sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

# non-centered version

m3dnc <- ulam(
    alist(
        C ~ bernoulli(p),
        logit(p) <- abar[1] + a[D,1] + 
                    (abar[2] + a[D,2])*U + 
                    (abar[3] + a[D,3])*A +
                    (abar[4] + a[D,4])*sum( delta_j[1:K] ),
        transpars> matrix[61,4]:a <-
            compose_noncentered( sigma , L_Rho , Z ),
        matrix[4,61]:Z ~ normal( 0 , 1 ),
        vector[4]: delta_j <<- append_row( 0 , delta ),
        simplex[3]: delta ~ dirichlet( Kprior ),
        vector[4]:abar ~ normal(0,1),
        cholesky_factor_corr[4]:L_Rho ~ lkj_corr_cholesky( 4 ),
        vector[4]:sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

precis(m3dnc,3,pars="abar")

# 4
# Can you also go beyond the parameter estimates from problem 3 and compute a marginal causal effect of urban living for each district, using a standard age distribution? Use any population age distribution you like. The important thing is to project the estimates from the sample of each district (which does not have representative age distributions) to the population. If you think the different districts should have different age distributions, that would be even more interesting.

dat$Ks <- standardize(dat$K)
m4 <- ulam(
    alist(
        # C model
        C ~ bernoulli(p),
        logit(p) <- abar[1] + a[D,1] + 
                    (abar[2] + a[D,2])*U + 
                    (abar[3] + a[D,3])*A +
                    (abar[4] + a[D,4])*Ks,

        # K model
        Ks ~ normal( mu , tau ),
        mu <- aK + bAK*A + bUK*U,
        c(aK,bAK,bUK) ~ normal(0,1),
        tau ~ exponential(1),

        # guts of the machine below
        transpars> matrix[61,4]:a <-
            compose_noncentered( sigma , L_Rho , Z ),
        matrix[4,61]:Z ~ normal( 0 , 1 ),
        vector[4]:abar ~ normal(0,1),
        cholesky_factor_corr[4]:L_Rho ~ lkj_corr_cholesky( 4 ),
        vector[4]:sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4 )

# age distribution - total sample
Ax <- round( d$age.centered - min(d$age.centered) )
plot(table(Ax),xlab="age (centered & shifted)",lwd=4,col=2,ylab="frequency")

# sample from our idealized age distribution
f <- seq( from=80 , to=20 , len=34 )
lines( 0:33 , f , lwd=4 , col=1 )

n <- 1000
Asim <- sample(0:33,size=n,replace=TRUE,prob=f)
Asim <- standardize( Asim + min(d$age.centered) )

pU0 <- sapply( 1:61 , 
    function(dist) 
        sim(m4,vars=c("Ks","C"),data=list(A=Asim,U=rep(0,n),D=rep(dist,n)))$C 
    )

pU1 <- sapply( 1:61 , 
    function(dist) 
        sim(m4,vars=c("Ks","C"),data=list(A=Asim,U=rep(1,n),D=rep(dist,n)))$C 
    )

plot( NULL , xlim=c(1,61) , ylim=c(-0.5,0.5) , ylab="P(C|do(U=1)) - P(C|do(U=0))" , xlab="District" )
for ( i in 1:61 ) points( i , mean(pU1[,i] - pU0[,i]) , lwd=4 , col=2 )
abline(h=0,lty=3,lwd=1)

ate <- mean(pU1-pU0)
abline(h=ate,lwd=2)
