# simulation of simple ptolemaic epicycles model to show retrograde motion

library(rethinking)
library(animation)

# two body version
# earth at center
# mars orbits an orbit around earth

my_circle <- function(x=0,y=0,r=1,angle=0,...) {
    a <- seq(angle, angle + 2 * pi, length = 360)
    lines( r*cos(a)+x , r*sin(a)+y , ... )
}
# plot(NULL,xlim=c(-1,1),ylim=c(-1,1)); my_circle(0,0,0.5,lty=2)
p2c <- function(a) c(cos(a) , sin(a))

# need updating for radial position of (1) mars epicycle center (2) mars itself


blank()

par(bg = 'black')
plot(NULL,xlim=c(-1,1),ylim=c(-1,1),bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

draw_ptolemy <- function(tx,history=FALSE,r1=1,r2=0.2) {

    if ( history==TRUE ) {
        xs <- rep(NA,tx)
        ys <- rep(NA,tx)
        for ( ttx in 1:tx ) {
            xy <- p2c( k1[ttx] )*r1
            xy2 <- xy + p2c( k2[ttx] )*r2
            xs[ttx] <- xy2[1]
            ys[ttx] <- xy2[2]
        }
        lines( xs , ys , lty=1 , lwd=2 , col=2 )
    }

    my_circle(0,0,r=r1,lty=2,lwd=0.5,col="white")

    xy <- p2c( k1[tx] )*r1
    my_circle( xy[1] , xy[2] , angle=k1[tx] , r=r2 , lty=2 , lwd=0.5 , col="white" )

    xy2 <- xy + p2c( k2[tx] )*r2
    points( xy2[1] , xy2[2] , pch=16 , col=2 , cex=2 ) # mars
    points( 0 , 0 , pch=16 , col=4 , cex=4 ) # earth
}

mpts <- 300
k1 <- seq( 0 , 2 * pi , length = mpts)
mu <- 2
k2 <- seq( 0 , mu * 2 * pi , length = mpts )

ani.record(reset = TRUE)  # clear history before recording
for ( i in 1:length(k1) ) {
    par(bg = 'black')
    par(xpd=NA)
    plot(c(-1, 1), c(-1, 1), type = "n", asp = 1, axes = FALSE, 
            xlab = "", ylab = "")
    draw_ptolemy( i , history=TRUE , r1=0.77 , r2=0.5 )
    ani.record()  # record the current frame
    par(bg = 'white')
}

oopts = ani.options(interval = 0.01)
ani.replay()

# ani.saveqz(dpi=150)

# convert -alpha remove -background white -delay 5 -loop 0 frame*.png geocentric.gif

# convert -delay 5 geocentric.gif geocentric.gif

############################################
# now kopernican with earth and mars

blank(w=3)

par(bg = 'black')
plot(NULL,xlim=c(-3,1),ylim=c(-1,1),bty="n",xaxt="n",yaxt="n",xlab="",ylab="")


kop_drawtime <- function(tx,history=FALSE,r1=0.5,r2=1,r3=7) {

    my_circle(0,0,r=r1,lty=2,lwd=0.5,col="white")
    my_circle(0,0,r=r2,lty=2,lwd=0.5,col="white")

    # earth
    xy <- p2c( k1[tx] )*r1
    
    # mars
    xy2 <- p2c( k2[tx] )*r2
    
    # line between earth-mars
    lines( c(xy[1],xy2[1]) , c(xy[2],xy2[2]) , lwd=2 , col="white" )

    points( xy2[1] , xy2[2] , pch=16 , col=2 , cex=2 ) # mars
    points( xy[1] , xy[2] , pch=16 , col=4 , cex=4 ) # earth

    # sky
    my_circle(0,0,r=r3,lty=2,lwd=0.5,col="white")

    # now need line projecting out to sky "orbit" with same slope as line btw earth-mars
    # equation for the sky is x^2 + y^2 = r3^2
    # line has slope m = (xy2[2]-xy[2])/(xy2[1]-xy[1])
    # y1 - y2 = m*(x1 - x2)
    # let (X,Y) be our solution points, then:
    # Y - xy[2] = m*(X - xy[1]) , X^2 + Y^2 = r^2
    # will have quadratic form
    # X = -((m (-m x + y) + Sqrt[(1 + m^2) r^2 - (-m x + y)^2])/(1 + m^2))
    # Y = (y - m (x + Sqrt[(1 + m^2) r^2 - (-m x + y)^2]))/(1 + m^2)
    m <- (xy2[2]-xy[2])/(xy2[1]-xy[1])
    x <- xy2[1]
    y <- xy2[2]
    X <- -((m*(-m*x + y) + sqrt((1 + m^2)*r3^2 - (-m*x + y)^2))/(1 + m^2))
    Y <- (y - m*(x + sqrt((1 + m^2)*r3^2 - (-m*x + y)^2)))/(1 + m^2)
    if ( x < 0 ) lines( c(X,x) , c(Y,y) , lty=2 , col="white" )

    # ghost mars on sky
    if ( x < 0 ) points( X , Y , pch=1 , lwd=3 , col=2 , cex=2 ) # mars

}

mpts <- 300
k2 <- seq( pi/2 , 2 * pi + pi/2 , length = mpts)
mu <- 2
k1 <- seq( 0 , mu * 2 * pi , length = mpts )

ani.record(reset = TRUE)  # clear history before recording
for ( i in 1:length(k1) ) {
    par(bg = 'black')
    par(xpd=NA)
    plot(c(-7, 2), c(-2, 2), type = "n", asp = 1, axes = FALSE, 
            xlab = "", ylab = "")
    kop_drawtime( i , history=TRUE , r1=0.8 , r2=1.2 , r3=7.5 )
    ani.record()  # record the current frame
    par(bg = 'white')
}

oopts = ani.options(interval = 0.03)
ani.replay()

# ani.saveqz(dpi=120)

# convert -alpha remove -background white -delay 5 -loop 0 frame*.png kopernicus.gif
