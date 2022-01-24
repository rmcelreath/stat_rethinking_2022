# cherry blossom cover plot

library(rethinking)
library(animation)

data(cherry_blossoms)
d <- cherry_blossoms

library(splines)
num_knots <- 30

# spline on blossom doy
d3 <- d[ complete.cases(d$doy) , ] # complete cases on doy

# subsample

blossom_col <- sapply( d3$doy , function(y) hsv(1, rbeta2(1, inv_logit(logit(0.1)+0.02*y) ,10) ,1,0.8) )

blank2(w=2,h=0.8)

#set.seed(1375)
ani.record(reset = TRUE)  # clear history before recording
nf <- 1
ntwinkle <- 50
head_frames <- 30
twinkle_frames <- 0
tot_frames <- head_frames + twinkle_frames
ynorm <- normalize(d3$year)
yseq <- seq( from=0 , to=1 , length.out=twinkle_frames )
do_save <- TRUE

knothist <- rep(30,tot_frames)
knothist[c(10,20,28)] <- 50
knothist[25] <- 10

for ( i in 1:tot_frames ) {

    print(concat("frame ",i))

    if ( i > head_frames ) {
        # weight twinkle set by year
        pvec <- 1 - exp( -10*(ynorm - yseq[i-head_frames])^2 )
        idx <- sample( 1:nrow(d3) , size=nrow(d3)-ntwinkle , prob=pvec )
    } else {
        # stable head still so all points
        idx <- 1:nrow(d3)
    }

    dd <- d3[ idx , ]
    dd <- dd[ order(dd$year) , ]

    num_knots <- knothist[i]
    knot_list <- seq( from=min(dd$year) , to=max(dd$year) , length.out=num_knots )
    B3 <- t(bs(dd$year, knots=knot_list , degree=3, intercept = FALSE))

    m2 <- quap(
        alist(
            Y ~ dnorm( mu , exp(log_sigma) ) ,
            mu <- a0 + as.vector( a %*% B ),
            a0 ~ dnorm(100,10),
            a ~ dnorm(0,10),
            log_sigma ~ dnorm(0,0.5)
        ),
        data=list( Y=dd$doy , B=B3 ) , start=list(a=rep(0,nrow(B3))) )

    # PLOT

    par( mfrow=c(1,1) , mgp = c(1.25, 0.25, 0), mar = c(0.75, 2.5, 0.75, 0.75) + 0.1, 
            tck = -0.02, cex.axis = 0.8 )

    xcex <- 1.2
    xpch <- 16
    xcol1 <- col.alpha(rangi2,0.3)
    col_spline <- col.alpha("black",0.4)
    xlims <- c(850,2000)

    y <- d3$doy
    y <- y - min(y)
    y <- y/max(y)

    #blossom_col <- sapply( d3$doy , function(y) hsv(1, rbeta2(1, inv_logit(logit(0.1)+0.02*y) ,10) ,1,0.8) )

    plot( NULL , cex=xcex , ylab="" , xlim=xlims , bty="n" , axes=FALSE , xlab="" , ylim=range(d3$doy) )
    l <- link( m2 )
    li <- apply(l,2,PI,0.9)
    points( dd$year , dd$doy , col=blossom_col ,  pch=8  , cex=xcex , lwd=2 )

    shade( li , dd$year , col=grau(0.3) )

    ani.record()  # record the current frame
    if ( do_save==TRUE ) {
        quartz.save( concat("frames/frame_",1000+nf,".png") , type = "png", device = dev.cur(), dpi = 200 )
        nf <- nf + 1
    }

}

oopts = ani.options(interval = 0.1)
ani.replay()

# convert -alpha remove -background white -delay 10 -loop 0 frame*.png title_twinkle.gif
