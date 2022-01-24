library(rethinking)
library(animation)

ani.saveqz <- function (list , dpi=150) {   
    if (missing(list)) 
        list = animation:::.ani.env$.images
    lapply(1:length(list), function(x) {
        dev.hold()
        replayPlot( list[[x]] )
        ani.pause()
        quartz.save( concat("frames/frame_",1000+x,".png") , type = "png", device = dev.cur(), dpi = dpi )
    })
    invisible(NULL)
}

# posterior and prior predictive simulation
# uses globe tossing example
# sample: W L W W W L W L W (6 W, 3 L)

blank(w=2,h=0.7)
par(mfrow=c(1,3))

# left: posterior (prior)
# middle: sampling distribution and highlighted sample from a specific value of [left]
# right: accumulated samples from posterior (prior) predictive distribution

W <- 6
L <- 3
N <- W + L
a <- 1 + W #alpha of posterior
b <- 1 + L #beta of posterior

n_frames <- 500
pp <- rep(0,10) # 0 to 9 water observed
p <- NA
ani.record(reset = TRUE)
set.seed(8675)
for ( f in 1:n_frames ) {

    p_old <- p
    p <- rbeta(1,a,b)
    w <- rbinom(1,size=9,prob=p)
    pp[w+1] <- pp[w+1] + 1

    # flash cycle
    gr <- 1:3
    if ( f > 0 ) gr <- 4 # fast after first 50 frames
    for ( g in gr ) {

        par(mfrow=c(1,3))
        par(cex.lab=1.3, cex.axis=1.3)

        # sample from posterior
        curve( dbeta(x,a,b) , from=0 , to=1 , xlab="proportion water" , ylab="" , bty="n" , lwd=3 , yaxt="n" , col=1 )
        mtext("Posterior distribution")
        if ( g==1 | g>3 )
            lines( c(p,p) , c(0,dbeta(p,a,b)) , lwd=10 , col=col.alpha(2,0.4) )
        lines( c(p,p) , c(0,dbeta(p,a,b)) , lwd=3 , col=2 )

        # sampling distribution of binomial(size=9) for p
        plot(NULL,xlim=c(0,9),ylim=c(0,dbinom(round(p*9),9,p)),bty="n",yaxt="n",xlab="number of water samples",ylab="",xaxt="n")
        axis(1,at=c(0,3,6,9),labels=c(0,3,6,9))
        mtext("Predictive distribution for p")
        for ( i in 0:9 ) {
            if ( i==w & g>1 ) lines( c(i,i) , c(0,dbinom(i,9,p)) , lwd=10 , col=col.alpha(2,0.4) )
            lines( c(i,i) , c(0,dbinom(i,9,p)) , lwd=3 , col=ifelse(i==w & g>1,2,1) )
        }

        # posterior predictive (accumulated)
        plot(NULL,xlim=c(0,9),ylim=c(0,max(105,pp)),bty="n",yaxt="n",xlab="number of water samples",ylab="",xaxt="n")
        axis(1,at=c(0,3,6,9),labels=c(0,3,6,9))
        mtext("Posterior predictive")
        for ( i in 0:9 ) {
            if ( i==w & g>2 ) lines( c(i,i) , c(0,pp[i+1]) , lwd=10 , col=col.alpha(2,0.4) )
            lines( c(i,i) , c(0,pp[i+1]) , lwd=3 , col=ifelse(i==w & g>2,2,1) )
        }

        ani.record()

    }

}

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=150)

# convert -alpha remove -background white -delay 10 -loop 0 frame*.png post_predict_sim.gif

# convert -delay 10 post_predict_sim.gif post_predict_sim.gif

# convert post_predict_sim.gif \( -clone [1..150] -set delay 50 \) -swap [1..150],-[1..150] +delete post_predict_sim.gif
