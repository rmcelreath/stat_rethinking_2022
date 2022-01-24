# simulation of simple Kopernican model showing epicycle

library(rethinking)
library(animation)

my_circle <- function(x=0,y=0,r=1,angle=0,...) {
    a <- seq(angle, angle + 2 * pi, length = 360)
    lines( r*cos(a)+x , r*sin(a)+y , ... )
}
# plot(NULL,xlim=c(-1,1),ylim=c(-1,1)); my_circle(0,0,0.5,lty=2)
p2c <- function(a) c(cos(a) , sin(a))

# blank()

par(bg = 'black')
plot(NULL,xlim=c(-1,1),ylim=c(-1,1),bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

draw_kop_epi <- function(tx,history=FALSE,r1=1,r2=0.2,rays=FALSE,epi=FALSE,sun_x=-0.2) {

    if ( history==TRUE ) {
        xs <- rep(NA,tx)
        ys <- rep(NA,tx)
        for ( ttx in 1:tx ) {
            xy <- p2c( k1[ttx] )*r1
            xy2 <- xy + p2c( k2[ttx] )*r2
            xs[ttx] <- xy2[1]
            ys[ttx] <- xy2[2]
        }
        lines( xs , ys , lty=1 , lwd=2 , col=col.alpha(4,0.5) )
    }

    if ( epi==TRUE )
    my_circle(0,0,r=r1,lty=2,lwd=0.5,col="white")

    xy <- p2c( k1[tx] )*r1
    if ( epi==TRUE ) {
        my_circle( xy[1] , xy[2] , angle=k1[tx] , r=r2 , lty=2 , lwd=0.5 , col="white" )
        points(0,0,pch=3,col="white",lwd=2)
    }

    xy2 <- xy + p2c( k2[tx] )*r2
    points( xy2[1] , xy2[2] , pch=16 , col=4 , cex=2 ) # earth
    points( sun_x , 0 , pch=16 , col="yellow" , cex=4 ) # sun

    if ( rays==TRUE ) {
        # lines from sun to earth
        for ( j in 1:tx )
            if ( j/25 == floor(j/25) ) {
                xy <- p2c( k1[j] )*r1
                xy2 <- xy + p2c( k2[j] )*r2
                lines( c(sun_x,xy2[1]) , c(0,xy2[2]) , lwd=0.5 , col="yellow" , lty=2 )
            }
    }
}

mpts <- 300
k1 <- seq( 0 , 2 * pi , length = mpts)
k2 <- seq( 0 , 2 * 2 * pi , length = mpts )

ani.record(reset = TRUE)  # clear history before recording

# first orbit without epicycles shown
for ( i in 1:length(k1) ) {
    par(bg = 'black')
    par(xpd=NA)
    plot(c(-1, 1), c(-1, 1), type = "n", asp = 1, axes = FALSE, 
            xlab = "", ylab = "")
    draw_kop_epi( i , history=TRUE , r1=1 , r2=0.1 , rays=FALSE , epi=FALSE )
    ani.record()  # record the current frame
    par(bg = 'white')
}
# second orbit with epicycles
for ( i in 1:length(k1) ) {
    par(bg = 'black')
    par(xpd=NA)
    plot(c(-1, 1), c(-1, 1), type = "n", asp = 1, axes = FALSE, 
            xlab = "", ylab = "")
    draw_kop_epi( i , history=TRUE , r1=1 , r2=0.1 , rays=FALSE , epi=TRUE )
    ani.record()  # record the current frame
    par(bg = 'white')
}
# third orbit with epicycles and sun segments
for ( i in 1:length(k1) ) {
    par(bg = 'black')
    par(xpd=NA)
    plot(c(-1, 1), c(-1, 1), type = "n", asp = 1, axes = FALSE, 
            xlab = "", ylab = "")
    draw_kop_epi( i , history=TRUE , r1=1 , r2=0.1 , rays=TRUE , epi=TRUE )
    ani.record()  # record the current frame
    par(bg = 'white')
}

oopts = ani.options(interval = 0.02)
ani.replay()

# ani.saveqz(dpi=150)

# convert -alpha remove -background white -delay 5 -loop 0 frame*.png kopernicus_epi.gif

# convert -delay 5 kopernicus_epi.gif kopernicus_epi.gif
