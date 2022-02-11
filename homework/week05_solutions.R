# week 5

# 1

library(rethinking)
data(NWOGrants)
d <- NWOGrants
dat <- list(
    A = as.integer(d$awards),
    N = as.integer(d$applications),
    G = ifelse( d$gender=="f" , 1L , 2L ) ,
    D = as.integer(d$discipline) 
)

# for total effect, just G, no D
m1 <- ulam(
    alist(
        A ~ binomial( N , p ),
        logit(p) <- a[G],
        a[G] ~ normal(-1,1)
    ), data=dat , chains=4 , cores=4 )

precis(m1,2)

post1 <- extract.samples(m1)
post1$pG1 <- inv_logit( post1$a[,1] )
post1$pG2 <- inv_logit( post1$a[,2] )
post1$G_contrast <- post1$pG1 - post1$pG2

dens( post1$G_contrast , lwd=4 , col=2 , xlab="F-M contrast (total)" )
abline( v=0 , lty=3 )

# 2

# now with discipline

m2 <- ulam(
    alist(
        A ~ binomial( N , p ),
        logit(p) <- a[G,D],
        matrix[G,D]:a ~ normal(-1,1)
    ), data=dat , chains=4 , cores=4 )

precis(m2,3)

total_apps <- sum(dat$N)
apps_per_disc <- sapply( 1:9 , function(i) sum(dat$N[dat$D==i]) )

pG1 <- link(m2,data=list(
    D=rep(1:9,times=apps_per_disc),
    N=rep(1,total_apps),
    G=rep(1,total_apps)))

pG2 <- link(m2,data=list(
    D=rep(1:9,times=apps_per_disc),
    N=rep(1,total_apps),
    G=rep(2,total_apps)))

dens( pG1 - pG2 , lwd=4 , col=2 , xlab="F-M contrast (marginal)" , xlim=c(-0.3,0.3) )
abline( v=0 , lty=3 )

# show contrasts for each discipline
plot( NULL , xlim=c(-0.4,0.4) , ylim=c(0,18) , xlab="F-M contrast for each discipline" , ylab="Density" )
abline( v=0 , lty=3 )
dat$disc <- as.character(d$discipline)
disc <- dat$disc[order(dat$D)]
for ( i in 1:9 ) {
    pG1Di <- link(m2,data=list(D=i,N=1,G=1))
    pG2Di <- link(m2,data=list(D=i,N=1,G=2))
    Gcont <- pG1Di - pG2Di
    dens( Gcont , add=TRUE , lwd=3 , col=i )
    xloc <- ifelse( mean(Gcont) < 0 , -0.35 , 0.35 )
    xpos <- ifelse( mean(Gcont) < 0 , 4 , 2 )
    text( xloc + 0.5*mean(Gcont) , 18-i , disc[2*i] , col=i , pos=xpos , font=2 )
}

# 3

total_f <- sum(d$applications[d$gender=="f"])
pDf <- d$applications[d$gender=="f"] / total_f

total_m <- sum(d$applications[d$gender=="m"])
pDm <- d$applications[d$gender=="m"] / total_m

# overall award rate in each discipline
n_apps <- xtabs( dat$N ~ dat$D )
n_awards <- xtabs( dat$A ~ dat$D )
p_award <- n_awards / n_apps

# f/m award rate in each discipline
post <- extract.samples(m2)
pF <- apply( inv_logit(post$a[,1,]) , 2 , mean )
pM <- apply( inv_logit(post$a[,2,]) , 2 , mean )

plot( pDf , pDm , lwd=3 , col=ifelse(pDf>pDm,2,4) , pch=ifelse(pF>pM,"F","M") )
abline(a=0,b=1, lty=3)

identify( pDf , pDm , round(p_award,2) )

# 4

library(rethinking)
data(UFClefties)
d <- UFClefties
dat <- list(
    Y = d$fighter1.win,
    L1 = d$fighter1.lefty,
    L2 = d$fighter2.lefty
)

m4a <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- b*(L1 - L2),
        b ~ normal(0,1)
    ) , data=dat , chains=4 , cores=4  )

post <- extract.samples(m4a)
dens( inv_logit(post$b) , xlab="prob lefty win over righty" , lwd=4 , col=2 , xlim=c(0,1) )
abline(v=0.5,lty=3)

# model that tries to estimate fighter ability
# not mentioned in solution guide, but this is what it would look like
# sort(unique(c(d$fighter1,d$fighter2)))
dat$id1 <- d$fighter1
dat$id2 <- d$fighter2
m4b <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- (a[id1]*sigma + b*L1) - (a[id2]*sigma + b*L2) ,
        vector[244]:a ~ normal(0,1),
        b ~ normal(0,0.5),
        sigma ~ exponential(1)
    ) , data=dat , chains=4 , cores=4  )
precis(m4b,2)

# SIMULATION
# collider explanation example
N <- 5000
L <- rbern(N,0.1)
A <- rnorm(N)
Q <- rep(1,N)

# qualify if A large enough or Lefty
Q <- ifelse( A > 2 | (A > 1.25 & L==1) , 1 , 0 )

# summarize
table(Q,L)
sum(Q==1 & L==1)/sum(Q==1) # prop lefties
mean( A[L==1 & Q==1] )
mean( A[L==0 & Q==1] )

# now sim fights among those qualifying
k <- 2 # importance of ability differences
b <- 0.5 # lefty advantage
l <- L[Q==1]
a <- A[Q==1]
M <- sum(Q==1)
W <- rep(NA,M/2) # matches
id1 <- rep(NA,M/2)
id2 <- rep(NA,M/2)
L1 <- id1
L2 <- L1
for ( i in 1:(M/2) ) {
    a1 <- a[i] + b*l[i]
    a2 <- a[M/2+i] + b*l[M/2+i]
    p <- inv_logit( k*(a1 - a2) )
    f1win <- rbern(1,p)
    W[i] <- f1win
    id1[i] <- i
    id2[i] <- M/2+i
    L1[i] <- l[i]
    L2[i] <- l[M/2+i]
}

# now estimate left hand advantage

datx <- list( Y=W , id1=id1, id2=id2 , L1=L1 , L2=L2 , nid=max(id2) )

m4x <- ulam(
    alist(
        Y ~ bernoulli(p),
        logit(p) <- b*(L1 - L2),
        b ~ normal(0,0.5)
    ) , data=datx , chains=4 , cores=4  )

precis(m4x,1)

post <- extract.samples(m4x)
dens( inv_logit(post$b) , xlab="prob lefty win over righty" , lwd=4 , col=2 , xlim=c(0,1) )
abline(v=0.5,lty=3)
