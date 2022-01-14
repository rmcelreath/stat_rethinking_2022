# Statistical Rethinking 2022
# Homework for week 1

## Problem 1

# define grid
p_grid <- seq( from=0 , to=1 , length.out=1000 )
# define prior
prior <- rep( 1 , 1000 )
# compute likelihood at each value in grid
likelihood <- dbinom( 4 , size=15 , prob=p_grid )
# compute product of likelihood and prior
unstd.posterior <- likelihood * prior
# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="l" , xlab="proportion water" , ylab="posterior probability" )

## Problem 2

# define grid
p_grid <- seq( from=0 , to=1 , length.out=1000 )
# define prior
prior <- c( rep( 0 , 500 ) , rep( 2 , 500 ) )
# compute likelihood at each value in grid
likelihood <- dbinom( 4 , size=6 , prob=p_grid )
# compute product of likelihood and prior
unstd.posterior <- likelihood * prior
# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

plot( p_grid , posterior , type="l" , xlab="proportion water" , ylab="posterior probability" )

## Problem 3

samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )
plot( samples , ylim=c(0,1) , xlab="samples" , ylab="proportion water" )
PI( samples , 0.89 )
HPDI( samples , 0.89)

## Problem 4 - optional challenge

# Pr(W|W) = 0.8
# Pr(W|L) = 0.2
# Pr(W) = 0.7*0.8

library(rethinking)

# sim as double sampling process
N <- 1e5
trueW <- rbinom(N,size=20,prob=0.7)
obsW <- rbinom(N,size=trueW,prob=0.8)

# or as single sample
W <- rbinom(N,size=20,prob=0.7*0.8)

mean(obsW/20)
mean(W/20)

# now analyze
# Pr(p|W,N) = Pr(W|p,N)Pr(p) / Z
# Pr(W|N,p) = Pr(W)Pr(W|W)

W <- rbinom(1,size=20,prob=0.7*0.8)
grid_p <- seq(from=0,to=1,len=100)
pr_p <- dbeta(grid_p,1,1)
prW <- dbinom(W,20,grid_p*0.8)
post <- prW*pr_p

post_bad <- dbinom(W,20,grid_p)

plot(grid_p,post,type="l",lwd=4,xlab="proportion water",ylab="plausibility")
lines(grid_p,post_bad,col=2,lwd=4)
