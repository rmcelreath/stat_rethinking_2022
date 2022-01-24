library(rethinking)
library(animation)

# football field gaussian

blank(w=2)

n_tosses <- 100
n_coins <- 1000

pos <- rep(0,n_coins)
xcol <- "white"
probs <- c(1,1)
ani.record(reset = TRUE)
hmax <- 20
fmax <- 20
seq_view <- FALSE
sequences <- matrix(NA,nrow=n_coins,ncol=n_tosses)
set.seed(384)
for ( f in 1:n_tosses ) {
    par( mfrow=c(1,2) , bg = 'black')
    par( xpd=FALSE )

    # field
    pos_old <- pos
    if ( f > 1 ) {
        # toss coins and update positions
        toss <- sample( c(-1,1) , size=n_coins , replace=TRUE , prob=probs )
        #toss <- runif(n_coins,-1,1)
        pos <- pos + toss
        xcol <- ifelse( pos > 0 , 4 , 2 )
        #xcol <- ifelse( pos==0 , "white" , xcol )
        #text( rep( -3 , n_coins ) , 1:n_coins , ifelse( toss==1 , "+" , "–" ) , col=xcol )
        sequences[,f] <- toss
    }
    if ( seq_view==FALSE ) {
        fmax <- max(fmax,abs(pos))
        plot(NULL, xlim=c(-fmax,fmax) ,ylim=c(1,n_coins),bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
        abline( v=0 , lty=2 , col="white" )
        #for ( v in c(1,-1,2,-2,3,-3) ) abline( v=v , lty=2 , col="white" )
        par( xpd=NA )
        #xcol2 <- ifelse( pos_old > 0 , col.alpha(4,0.2) , col.alpha(2,0.2) )
        #points( pos_old , 1:n_coins , pch=16 , cex=1.5 , col=xcol2 )
        points( pos , 1:n_coins , pch=16 , cex=0.5 , col=xcol )
    } else {
        # sequence view
        fmax <- max(fmax,f)
        plot(NULL, xlim=c(1,n_tosses) ,ylim=c(1,n_coins),bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
        # sort by sum
        #seq_sort <- order( pos )
        #seq_copy <- sequences[ seq_sort , ]
        for ( y in 1:n_coins ) {
            xcol <- ifelse( sequences[ y , 1:f ] > 0 , 4 , 2 )
            points( 1:f , rep(y,f) , col=xcol , pch=16 , cex=0.5 )
        }
    }

    # histogram
    par( xpd=FALSE )
    #plot.new()
    hmax <- max(hmax,abs(pos))
    tab <- table( round( pos ) )
    plot( NULL , yaxt="n" , bty="n" , col="white" , xlim=c(-hmax,hmax) , ylim=c(0,max(10,tab)) )
    abline( v=0 , lty=2 , col="white" )
    for ( i in 1:length(tab) ) {
        x <- as.numeric( names(tab)[i] )
        y <- as.numeric( tab[i] )
        xcol <- ifelse( x > 0 , 4 , 2 )
        #if ( x==0 ) xcol <- "white"
        lines( c(x,x) , c(0,y) , lwd=6 , col=xcol )
    }
    ani.record()
}

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 10 -loop 0 frame*.png gaussian_sim.gif
# convert -delay 20 gaussian_sim.gif gaussian_sim.gif

# ways to get x in n tosses
n <- 50
z <- sapply( 0:n , function(x) choose(n,x) )
plot( z , col="white" )


# growth model gaussian

blank(w=2)

n_tosses <- 100
n_coins <- 1000

pos <- rep(0,n_coins)
xcol <- "white"
probs <- c(1,1)
ani.record(reset = TRUE)
hmax <- 20
fmax <- 20
seq_view <- FALSE
sequences <- matrix(NA,nrow=n_coins,ncol=n_tosses)
gf <- rnorm(n_tosses,0,0.1)
set.seed(384)
for ( f in 1:n_tosses ) {
    par( mfrow=c(1,2) , bg = 'black')
    par( xpd=FALSE )

    # field
    pos_old <- pos
    if ( f > 1 ) {
        # toss coins and update positions
        #toss <- sample( c(-1,1) , size=n_coins , replace=TRUE , prob=probs )
        toss <- runif(n_coins,0,1+gf)
        pos <- pos + toss
        xcol <- ifelse( pos > mean(pos) , 4 , 2 )
        #xcol <- ifelse( pos==0 , "white" , xcol )
        #text( rep( -3 , n_coins ) , 1:n_coins , ifelse( toss==1 , "+" , "–" ) , col=xcol )
        sequences[,f] <- toss
    }
    if ( seq_view==FALSE ) {
        fmax <- max(fmax,abs(pos))
        plot(NULL, xlim=c(0,fmax) ,ylim=c(1,n_coins),bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
        abline( v=mean(pos) , lty=2 , col="white" )
        #for ( v in c(1,-1,2,-2,3,-3) ) abline( v=v , lty=2 , col="white" )
        par( xpd=NA )
        #xcol2 <- ifelse( pos_old > 0 , col.alpha(4,0.2) , col.alpha(2,0.2) )
        #points( pos_old , 1:n_coins , pch=16 , cex=1.5 , col=xcol2 )
        points( pos , 1:n_coins , pch=16 , cex=0.5 , col=xcol )
    } else {
        # sequence view
        fmax <- max(fmax,f)
        plot(NULL, xlim=c(1,n_tosses) ,ylim=c(1,n_coins),bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
        # sort by sum
        #seq_sort <- order( pos )
        #seq_copy <- sequences[ seq_sort , ]
        for ( y in 1:n_coins ) {
            xcol <- ifelse( sequences[ y , 1:f ] > 0 , 4 , 2 )
            points( 1:f , rep(y,f) , col=xcol , pch=16 , cex=0.5 )
        }
    }

    # histogram
    par( xpd=FALSE )
    #plot.new()
    hmax <- max(hmax,abs(pos))
    tab <- table( round( pos ) )
    plot( NULL , yaxt="n" , bty="n" , col="white" , xlim=c(0,hmax) , ylim=c(0,max(10,tab)) )
    abline( v=0 , lty=2 , col="white" )
    for ( i in 1:length(tab) ) {
        x <- as.numeric( names(tab)[i] )
        y <- as.numeric( tab[i] )
        xcol <- ifelse( x > mean(pos) , 4 , 2 )
        #if ( x==0 ) xcol <- "white"
        lines( c(x,x) , c(0,y) , lwd=6 , col=xcol )
    }
    ani.record()
}

oopts = ani.options(interval = 0.1)
ani.replay()
