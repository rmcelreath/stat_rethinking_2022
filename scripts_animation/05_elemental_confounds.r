# new DAG intro scripts

# FORK example

n <- 1000
Z <- rbern( n , 0.5 )
X <- rbern( n , (1-Z)*0.1 + Z*0.9 )
Y <- rbern( n , (1-Z)*0.1 + Z*0.9 )

table(X,Y)
table(X,Y,Z)
cor(X,Y)
cor(X[Z==0],Y[Z==0])
cor(X[Z==1],Y[Z==1])

round(table(X[Z==1],Y[Z==1])/sum(Z==1),2)


# PIPE example

n <- 1000
X <- rbern( n , 0.5)
Z <- rbern( n , (1-X)*0.1 + X*0.9 )
Y <- rbern( n , (1-Z)*0.1 + Z*0.9 )

table(X,Y)
table(X,Y,Z)
cor(X,Y)
cor(X[Z==0],Y[Z==0])
cor(X[Z==1],Y[Z==1])

# COLLIDER example

n <- 1000
X <- rbern( n , 0.5 )
Y <- rbern( n , 0.5 )
Z <- rbern( n , ifelse(X+Y>0,0.9,0.2) )

table(X,Y)

table(X,Y,Z)

cor(X,Y)
cor(X[Z==0],Y[Z==0])
cor(X[Z==1],Y[Z==1])


######
# DIVORCE EXAMPLE

library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

blank(bty="n")


d$WHpm <- (d$WaffleHouses/d$Population)

plot( d$WHpm , d$Divorce , col=2 , lwd=3 , xlab="Waffle Houses per million" , ylab="Divorce rate" )
identify( d$WHpm , d$Divorce , d$Loc , cex=0.8 )

m <- lm( Divorce ~ WHpm , d )
post <- extract.samples(m)
xseq <- seq(from=-5,to=45,len=30)
mu <- sapply( xseq , function(x) post$Intercept + post$WHpm*x )
shade( apply(mu,2,PI) , xseq )

plot( d$Marriage , d$Divorce , col=ifelse(d$South==1,2,4) , lwd=3 , xlab="Marriage rate" , ylab="Divorce rate" )
identify( d$Marriage , d$Divorce , d$Loc , cex=0.8 )

plot( d$MedianAgeMarriage , d$Divorce , col=ifelse(d$South==1,2,4) , lwd=3 , xlab="Median age of marriage" , ylab="Divorce rate" )
identify( d$MedianAgeMarriage , d$Divorce , d$Loc , cex=0.8 )

plot( d$MedianAgeMarriage , d$Marriage , col=ifelse(d$South==1,2,4) , lwd=3 , xlab="Median age of marriage" , ylab="Marriage rate" )
identify( d$MedianAgeMarriage , d$Marriage , d$Loc , cex=0.8 )

# explain stratifying by continuous variable
plot( sort( d$MedianAgeMarriage ) , 1:50 , yaxt="n" , lwd=3 , col=2 , ylab="" , xlab="Median age of marriage" )



div <- 2
divs <- seq(from=0,to=1,len=div+1)
for ( i in 1:(length(divs)) ) abline( v=quantile(d$MedianAgeMarriage,divs[i]) , lwd=2 )

# visualize standardizing
plot( standardize(d$MedianAgeMarriage) , standardize(d$Divorce) , col=ifelse(d$South==1,2,4) , lwd=3 , xlab="Median age of marriage (standardized)" , ylab="Divorce rate (standardized)" )

# prior predictive simulation
n <- 20
a <- rnorm(n,0,10)
bM <- rnorm(n,0,10)
bA <- rnorm(n,0,10)

plot( NULL , xlim=c(-2,2) , ylim=c(-2,2) , xlab="Median age of marriage (standardized)" , ylab="Divorce rate (standardized)" )
Aseq <- seq(from=-3,to=3,len=30)
for ( i in 1:n ) {
    mu <- a[i] + bA[i]*Aseq
    lines( Aseq , mu , lwd=2 , col=2 )
}

# better priors
n <- 20
a <- rnorm(n,0,0.2)
bM <- rnorm(n,0,0.5)
bA <- rnorm(n,0,0.5)


# model
dat <- list(
    D = standardize(d$Divorce),
    M = standardize(d$Marriage),
    A = standardize(d$MedianAgeMarriage)
)

m_DMA <- quap(
    alist(
        D ~ dnorm(mu,sigma),
        mu <- a + bM*M + bA*A,
        a ~ dnorm(0,0.2),
        bM ~ dnorm(0,0.5),
        bA ~ dnorm(0,0.5),
        sigma ~ dexp(1)
    ) , data=dat )

plot(precis(m_DMA))

# posterior predictions


##########
# pipe example - plant growth

# as in book

#########
# collider examples

# real selection bias - publication

 
set.seed(1914)
N <- 200 # num grant proposals
p <- 0.1 # proportion to select
# uncorrelated newsworthiness and trustworthiness 
nw <- rnorm(N)
tw <- rnorm(N)
# select top 10% of combined scores
s<-nw+tw #totalscore
q <- quantile( s , 1-p ) # top 10% threshold 
selected <- ifelse( s >= q , TRUE , FALSE )
cor( tw[selected] , nw[selected] )

blank(bty="n")

plot( nw , tw , xlab="Newsworthiness" , ylab="Trustworthiness" , col=2 , lwd=3 )

plot( nw , tw , xlab="Newsworthiness" , ylab="Trustworthiness" , lwd=3 , col=ifelse(selected==TRUE,2,"gray") )
abline( lm(tw[selected]~nw[selected]), lwd=3  )

# happiness sim

fyears <- 100
start_years <- 100
library(animation)
ani.record(reset = TRUE)

for ( i in 1:fyears ) {

x <- sim_happiness( seed=1977 , N_years=start_years+i )

plot( x$age , x$happiness , pch=ifelse(x$married==1,16,1) , lwd=3 , col=ifelse(x$married==1,2,"gray") , xlab="Age (years)" , ylab="Happiness (standardized)" )

ani.record()

}

oopts = ani.options(interval = 0.05)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 10 -loop 0 frame*.png happiness.gif

######
# descendant example

n <- 1000
X <- rbern( n , 0.5 )
Z <- rbern( n , (1-X)*0.1 + X*0.9 )
Y <- rbern( n , (1-Z)*0.1 + Z*0.9 )
A <- rbern( n , (1-Z)*0.1 + Z*0.9 )

table(X,Y)

table(X,Y,A)

cor(X,Y)
cor(X[A==0],Y[A==0])
cor(X[A==1],Y[A==1])

############################
# d-separation plots
a <- 0.7
cols <- c( col.alpha(4,a) , col.alpha(2,a) )

# pipe

N <- 300
X <- rnorm(N)
Z <- rbern(N,inv_logit(X))
Y <- rnorm(N,(2*Z-1))

plot( X , Y , col=cols[Z+1] , lwd=3 )

abline(lm(Y[Z==1]~X[Z==1]),col="white",lwd=5)
abline(lm(Y[Z==1]~X[Z==1]),col=2,lwd=3)

abline(lm(Y[Z==0]~X[Z==0]),col="white",lwd=5)
abline(lm(Y[Z==0]~X[Z==0]),col=4,lwd=3)

abline(lm(Y~X),lwd=5,col="white")
abline(lm(Y~X),lwd=3)

# fork

cols <- c(4,2)

N <- 300
Z <- rbern(N)
X <- rnorm(N,2*Z-1)
Y <- rnorm(N,(2*Z-1))

plot( X , Y , col=cols[Z+1] , lwd=3 )

abline(lm(Y[Z==1]~X[Z==1]),col="white",lwd=5)
abline(lm(Y[Z==1]~X[Z==1]),col=2,lwd=3)

abline(lm(Y[Z==0]~X[Z==0]),col="white",lwd=5)
abline(lm(Y[Z==0]~X[Z==0]),col=4,lwd=3)

abline(lm(Y~X),lwd=5,col="white")
abline(lm(Y~X),lwd=3)

# collider

N <- 300
X <- rnorm(N)
Y <- rnorm(N)
Z <- rbern(N,inv_logit(2*X+2*Y-2))

plot( X , Y , col=cols[Z+1] , lwd=3 )

abline(lm(Y[Z==1]~X[Z==1]),col="white",lwd=5)
abline(lm(Y[Z==1]~X[Z==1]),col=2,lwd=3)

abline(lm(Y[Z==0]~X[Z==0]),col="white",lwd=5)
abline(lm(Y[Z==0]~X[Z==0]),col=4,lwd=3)

abline(lm(Y~X),lwd=5,col="white")
abline(lm(Y~X),lwd=3)
