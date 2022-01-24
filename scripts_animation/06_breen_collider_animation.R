library(rethinking)

N <- 200  # number of grandparent-parent-child triads
b_GP <- 1 # direct effect of G on P
b_GC <- 0 # direct effect of G on C
b_PC <- 1 # direct effect of P on C
b_U <- 2 #direct effect of U on P and C

set.seed(1)
U <- 2*rbern( N , 0.5 ) - 1
G <- rnorm( N )
P <- rnorm( N , b_GP*G + b_U*U )
C <- rnorm( N , b_PC*P + b_GC*G + b_U*U )
d <- data.frame( C=C , P=P , G=G , U=U )

m6.11 <- quap(
    alist(
        C ~ dnorm( mu , sigma ),
        mu <- a + b_PC*P + b_GC*G,
        a ~ dnorm( 0 , 1 ),
        c(b_PC,b_GC) ~ dnorm( 0 , 1 ),
        sigma ~ dexp( 1 )
    ), data=d )

plot(precis(m6.11,pars=c("b_GC","b_PC")))

points( b_GC , 1 , col=2 , pch=3 , lwd=4 )
points( b_PC , 2 , col=2 , pch=3 , lwd=4 )

# blank2()
# blank(bty="n")
plot( C ~ G , pch=1 , col=ifelse( U > 0 , 2 , 1 ) , cex=exp(P/6) , xlab="Grandparent education" , ylab="child education" , lwd=3 )

# highlight band of parental incomes
# and draw regression line through them
# band given in percentiles

library(animation)


Pband_seq <- seq( from=0 , to=1 , len=10 )

ani.record(reset = TRUE)
for ( i in 1:(length(Pband_seq)-1) ) {

Pband <- Pband_seq[i:(i+1)]

b <- quantile( P , Pband )
in_band <- ifelse( P>b[1] & P<b[2] , 1 , 0 )
plot( C ~ G , col=ifelse( U > 0 , 2 , 1 ) , pch=ifelse( in_band , 16 , 1 ) , cex=0.5*in_band+1 , xlab="Grandparent education" , ylab="Child education" , lwd=2 )
idx <- which( in_band==1 )
abline( lm( C[idx] ~ G[idx] ) , lwd=5 , col="white" )
abline( lm( C[idx] ~ G[idx] ) , lwd=3 )
text( -1 , max(C) , "good neighborhoods" , col="red" , cex=0.8 )
text( 1 , min(C) , "bad neighborhoods" , col="black" , cex=0.8 )
mtext( concat("Parents between ", round(Pband[1]*100) ,"th and ", round(Pband[2]*100) ,"th centiles") )

ani.record()

}

oopts = ani.options(interval = 0.5)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 20 -loop 0 frame*.png breen_collider.gif
# convert -delay 100 breen_collider.gif breen_collider.gif
