
library(rethinking)
library(animation)
library(ellipse)

# 

data(Howell1)
d <- Howell1[Howell1$age>=18,]

blank(bty="n")
plot( d$weight , d$height , col=2 , lwd=2 , xlab="weight (kg)" , ylab="height (cm)" )

# show height ~ weight priors progressively

blank(bty="n")

set.seed(2971)
n <- 100
a <- rnorm(n,178,20)
b <- rnorm(n,0,10)
b <- rlnorm(n,0,1)

plot( NULL , xlim=c(30,70) , ylim=c(-100,400) , xlab="weight (kg)" , ylab="height (cm)" , col=2 , lwd=2 )
abline(h=272 , lty=2)
abline(h=0,lty=2)

ani.record(reset = TRUE)

xseq <- seq(from=30,to=80,length=10)
for ( i in 1:n ) {
    y <- a[i] + b[i]*(xseq-mean(d$weight))
    lines( xseq , y , lwd=2 , col=col.alpha(2,0.5) )
    ani.record()
}

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 10 -loop 0 frame*.png howell_height_weight_prior2.gif
# convert -delay 20 gaussian_sim.gif gaussian_sim.gif
