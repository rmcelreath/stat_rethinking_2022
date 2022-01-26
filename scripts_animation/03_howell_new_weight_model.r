# new Howell script
# model weight as function of height (rather than reverse)

library(rethinking)
data(Howell1)
d <- Howell1

# draw data

blank(bty="n")

col2 <- col.alpha(2,0.8)
plot( d$height , d$weight , col=col2 , lwd=3 , cex=1.2 , xlab="height (cm)" , ylab="weight (kg)" )

d <- d[ d$age>=18 , ]
plot( d$height , d$weight , col=col2 , lwd=3 , cex=1.2 , xlab="height (cm)" , ylab="weight (kg)" )

# sample from example prior

n_samples <- 10
alpha <- rnorm(n_samples,0,1)
beta <- rnorm(n_samples,0,1)
plot(NULL,xlim=c(-2,2),ylim=c(-2,2),xlab="x",ylab="y")
for ( i in 1:n_samples )
    abline(alpha[i],beta[i],lwd=4,col=2)

# sample from generative model of weight

alpha <- 0
beta <- 0.5
sigma <- 5
n_individuals <- 100
H <- runif(n_individuals,130,170)
mu <- alpha + beta*H
W <- rnorm(n_individuals,mu,sigma)

plot(H,W,xlab="height (cm)",ylab="weight (kg)",col=2,lwd=3)

# sample prior for weight model

n <- 10
alpha <- rnorm(n,60,10)
beta <- rnorm(n,0,10)
Hbar <- 150
Hseq <- seq(from=130,to=170,len=30)
plot(NULL,xlim=c(130,170),ylim=c(10,100),xlab="height (cm)",ylab="weight (kg)")
for ( i in 1:n ) 
    lines( Hseq , alpha[i] + beta[i]*(Hseq-Hbar) , lwd=3 , col=2 )

# sample again with better priors

n <- 10
alpha <- rnorm(n,60,5)
beta <- rlnorm(n,0,1)
Hbar <- 150
Hseq <- seq(from=130,to=170,len=30)
plot(NULL,xlim=c(130,170),ylim=c(10,100),xlab="height (cm)",ylab="weight (kg)")
for ( i in 1:n ) 
    lines( Hseq , alpha[i] + beta[i]*(Hseq-Hbar) , lwd=3 , col=2 )

# validate statistical model with simulated data

alpha <- 70
beta <- 0.5
sigma <- 5
n_individuals <- 100
H <- runif(n_individuals,130,170)
mu <- alpha + beta*(H-mean(H))
W <- rnorm(n_individuals,mu,sigma)

dat <- list( H=H , W=W , Hbar=mean(H) )

m_validate <- quap(
    alist(
        W ~ dnorm(mu,sigma),
        mu <- a + b*(H-Hbar),
        a ~ dnorm(60,10),
        b ~ dlnorm(0,1),
        sigma ~ dunif(0,10)
    ), data=dat )

precis(m_validate)

# now repeat 1000 times, sampling inits from prior
fvalidate <- function( n=200 ) {
    alpha <- rnorm(1,60,10)
    beta <- rlnorm(1,0,1)
    sigma <- runif(1,0,10)
    n_individuals <- n
    H <- runif(n_individuals,130,170)
    mu <- alpha + beta*(H-mean(H))
    W <- rnorm(n_individuals,mu,sigma)
    dat <- list( H=H , W=W , Hbar=mean(H) )
    m_validate <- quap(
        alist(
            W ~ dnorm(mu,sigma),
            mu <- a + exp(logb)*(H-Hbar),
            a ~ dnorm(60,10),
            logb ~ dnorm(0,1),
            sigma ~ dunif(0,10)
        ), data=dat , start=list(a=60,logb=-1,sigma=5) )
    return(coef(m_validate))
}

rm(r)
r <- replicate( 100 , fvalidate() )


# statistical of weight

data(Howell1)
d <- Howell1
d <- d[ d$age>=18 , ]

dat <- list(
    W = d$weight,
    H = d$height,
    Hbar = mean(d$height) )

m_adults <- quap(
    alist(
        W ~ dnorm(mu,sigma),
        mu <- a + b*(H-Hbar),
        a ~ dnorm(60,10),
        b ~ dlnorm(0,1),
        sigma ~ dunif(0,10)
    ), data=dat )

# plot sample
col2 <- col.alpha(2,0.8)
plot( d$height , d$weight , col=col2 , lwd=3 , cex=1.2 , xlab="height (cm)" , ylab="weight (kg)" )

# expectation with 99% compatibility interval
xseq <- seq(from=130,to=190,len=50)
mu <- link(m_adults,data=list(H=xseq,Hbar=mean(d$height)))
lines( xseq , apply(mu,2,mean) , lwd=4 )
shade( apply(mu,2,PI,prob=0.99) , xseq , col=col.alpha(2,0.5) )

# 89% prediction interval
W_sim <- sim(m_adults,data=list(H=xseq,Hbar=mean(d$height)))
shade( apply(W_sim,2,PI,prob=0.89) , xseq , col=col.alpha(1,0.3) )


#########
# Howell adult weight updating animation

blank(w=1.6,bty="n")

data(Howell1)
d <- Howell1
d <- d[ d$age>=18 , ]

set.seed(17)
j <- sample(1:nrow(d))

dat <- list(
    W = d$weight[j],
    H = d$height[j],
    Hbar = mean(d$height) )

n_points_per <- 1
n_points <- ceiling( length(dat$W)/n_points_per )
n_points <- 30

xlims <- c(130,180)
ylims <- c(30,70)

ellipse_hist <- list()
ani.record(reset=TRUE)
NSTEPS <- 10
for ( i in 0:n_points ) {

    al <- alist(
        W ~ dnorm(mu,1),
        mu <- a + b*(H-Hbar) ,
        a ~ dnorm(60,10),
        b ~ dlnorm(0,1)
    )

    if ( i > 0 ) {
        iend <- min( n_points_per*i , length(dat$W) )
        dati <- list( W=dat$W[1:iend] , H=dat$H[1:iend] , Hbar=dat$Hbar )
        q <- quap( al , data=dati , dofit=TRUE )
        # print( precis(q) )
        pp <- extract.samples(q)
    } else {
        # i = 0
        q <- quap( al , data=list(y=1,x=1) , dofit=FALSE )
        pp <- as.data.frame(extract.prior(q))
    }

    a_link <- function(x) a + b*(x - 155)

    the_text <- concat("n = ",i)
    xseq <- seq(from=130,to=180,length=50)
    ahat <- mean(pp$a)
    bhat <- mean(pp$b)
    if ( i > 0 ) {
        mu <- link( q , data=list(H=xseq,Hbar=155) )
        ci <- apply(mu,2,PI)
    }

    use_cols <- c(2, 4, 5, 6, 3)

    if ( i > 0 & i==5 ) {
        qprecis <- precis(q,prob=0.999)
        post_xlims <- as.numeric(qprecis[1,3:4]) * c(1,1)
        post_ylims <- as.numeric(qprecis[2,3:4]) * c(0.9,1.2)
    }
    if ( i==0 ) {
        post_xlims <- c(60-30,60+30)
        post_ylims <- c(-2,5)
    }

    the_pf <- function(samp,step) {
            par(mfrow=c(1,2))

            plot(NULL,xlim=post_xlims,ylim=post_ylims,bty="n",xlab="intercept",ylab="slope")
            pm <- colMeans(pp)
            #lines( c(pm[1],pm[1]) , c(-4,pm[2]) , lty=2 )
            #lines( c(-4,pm[1]) , c(pm[2],pm[2]) , lty=2 )
            
            # interpolate between previous and current contours
            the_contour <- list()
            for ( l in c(0.5,0.8,0.95) )
                the_contour[[ concat(l) ]] <- ellipse(cov(pp),centre=pm,level=l)
            if ( length(ellipse_hist) > 0 & samp==1 ) {
                mixf <- (NSTEPS-step)/NSTEPS
                for ( l in c(0.5,0.8,0.95) ) {
                    #lines(ellipse_hist[[i]],col=col.alpha(1,0.3),lwd=1)
                    the_contour[[ concat(l) ]] <- (1-mixf)*the_contour[[concat(l)]] + mixf*ellipse_hist[[concat(l)]]
                }
            }

            # draw posterior contours
            thelwd <- 2
            if ( samp==1 ) thelwd <- 2+3*(NSTEPS-step)/NSTEPS
            for ( l in c(0.5,0.8,0.95) )
                lines(the_contour[[concat(l)]],col=1,lwd=thelwd )

            # draw point for each line on righthand plot
            pb <- get("prior_blend",envir=parent.frame(n = 1))
            points( pb[[1]][,step] , pb[[2]][,step] , lwd=2 , col=use_cols[1:3] , pch=16 , cex=1.6 )

            # x/y plot with moving lines
            plot(NULL,xlim=xlims,ylim=ylims,xlab="height (cm)",ylab="weight (kg)",bty="n")
            if ( i > 0 ) {
                #shade( ci , xseq )
                #abline( a=ahat , b=bhat , lwd=3 , col="gray" , lty=2 )
                points( dat$H[1:iend] , dat$W[1:iend] , pch=1 , lwd=3 , cex=1.2 )
            }
            mtext(the_text)
        }

    if ( i==0 ) r <- NULL
    n_prior_samples <- 5
    r <- anim_prior_predictive( prior=pp , linkf=a_link , n=n_prior_samples , n_to_show=3 , n_steps=NSTEPS , ylim=ylims , ylab="outcome" , xlim=xlims , pf=the_pf , do_reset=FALSE , start_from=r , accel=-0.5 , vel0=1 )

    for ( l in c(0.5,0.8,0.95) )
        ellipse_hist[[ concat(l) ]] <- ellipse(cov(pp),centre=colMeans(pp),level=l)

}

oopts = ani.options(interval = 0.02)
ani.replay()

# ani.saveqz(dpi=130)
# convert -alpha remove -background white -delay 3 -loop 0 frame*.png ols_howell_learn.gif
# convert -delay 6 ols_howell_learn.gif ols_howell_learn.gif

######
# weight model with categories for male/female

data(Howell1)
d <- Howell1
d <- d[ d$age>=18 , ]

dat <- list(
    W = d$weight,
    H = d$height,
    Hbar = mean(d$height),
    S = d$male + 1 ) # S=1 female, S=2 male

dd <- data.frame(H=dat$H,W=dat$W,S=dat$S)

# plot distributions of weight and height for each sex

blank(ex=2,bty="n")
par(mfrow=c(2,2))
dens(dat$H[dat$S==1],adj=0.2,lwd=4,col=2,xlim=c(130,180),xlab="height (cm)")
dens(dat$W[dat$S==1],adj=0.2,lwd=4,col=2,xlim=c(30,65),xlab="weight (kg)")
dens(dat$H[dat$S==2],adj=0.2,lwd=4,col=4,xlim=c(130,180),xlab="height (cm)")
dens(dat$W[dat$S==2],adj=0.2,lwd=4,col=4,xlim=c(30,65),xlab="weight (kg)")

# model
m_adults2 <- quap(
    alist(
        W ~ dnorm(mu,sigma),
        mu <- a[S] + b[S]*(H-Hbar),
        a[S] ~ dnorm(60,10),
        b[S] ~ dlnorm(0,1),
        sigma ~ dunif(0,10)
    ), data=dat )

# compute contrasts

post <- extract.samples(m_adults2)
a_contrast <- post$a[,1] - post$a[,2]
b_contrast <- post$b[,1] - post$b[,2]

blank(w=1.8,bty="n")
par(mfrow=c(1,2))
dens(a_contrast,show.zero=TRUE,lwd=4,col=4,xlab="intercept (a) contrast")
dens(b_contrast,show.zero=TRUE,lwd=4,col=2,xlab="slope (b) contrast")

# contrast of expectation at each height

H_seq <- seq(from=130,to=180,len=30)
Hbar <- 155
mu_F <- sapply( H_seq , function(H) with( post , a[,1] + b[,1]*(H-Hbar) ) )
mu_M <- sapply( H_seq , function(H) with( post , a[,2] + b[,2]*(H-Hbar) ) )
plot( H_seq , apply(mu_F,2,mean) )
lines( H_seq , apply(mu_M,2,mean) )

mu_contrast <- mu_F - mu_M


# plot sample
col2 <- col.alpha(2,0.6)
col4 <- col.alpha(4,0.6)
plot( d$height , d$weight , col=ifelse(d$male==1,col4,col2) , lwd=3 , cex=1.2 , xlab="height (cm)" , ylab="weight (kg)" )

xseq <- seq(from=130,to=190,len=50)

muF <- link(m_adults2,data=list(S=rep(1,50),H=xseq,Hbar=mean(d$height)))
lines( xseq , apply(muF,2,mean) , lwd=3 , col=2 )

muM <- link(m_adults2,data=list(S=rep(2,50),H=xseq,Hbar=mean(d$height)))
lines( xseq , apply(muM,2,mean) , lwd=3 , col=4 )

mu_contrast <- muF - muM
plot( NULL , xlim=range(xseq) , ylim=c(-6,8) , xlab="height (cm)" , ylab="weight contrast (F–M)" )
for ( p in c(0.5,0.6,0.7,0.8,0.9,0.99) )
    shade( apply(mu_contrast,2,PI,prob=p) , xseq )
abline(h=0,lty=2)

# negative correlation example for contrast

blank(bty="n")

x <- rmvnorm2(500,c(0,1),sigma=c(0.5,0.5),Rho=matrix(c(1,0.95,0.95,1),2,2))

col2 <- col.alpha(2,0.5)
plot(x,xlab="parameter 1",ylab="parameter 2",col=col2,lwd=3)

dens( x[,1] , lwd=4 , xlab="value" , col=2 , xlim=c(-3,3) , ylim=c(0,max(density(x[,1]-x[,2])$y)) )
dens( x[,2] , lwd=4 , col=4 , add=TRUE )
dens( x[,1]-x[,2], lwd=4 , lty=1 , add=TRUE )

# coded example of vector parameter

color <- c("C","M","Y","K")

#####################
# polynomials


data(Howell1)
d <- Howell1

dat <- list(
    W = d$weight,
    H = standardize(d$height),
    S = d$male + 1 ) # S=1 female, S=2 male

m_parabola <- quap(
    alist(
        W ~ dnorm(mu,sigma),
        mu <- a + b1*H + b2*H^2,
        a ~ dnorm(60,10),
        b1 ~ dnorm(0,10),
        b2 ~ dnorm(0,10),
        sigma ~ dunif(0,10)
    ), data=dat )

m_cubic <- quap(
    alist(
        W ~ dnorm(mu,sigma),
        mu <- a + b1*H + b2*H^2 + b3*H^3 + b4*H^4,
        a ~ dnorm(60,10),
        b1 ~ dnorm(0,10),
        b2 ~ dnorm(0,10),
        b3 ~ dnorm(0,10),
        b4 ~ dnorm(0,10),
        sigma ~ dunif(0,10)
    ), data=dat )

m_crazy <- quap(
    alist(
        W ~ dnorm(mu,sigma),
        mu <- a + b[1]*H + b[2]*H^2 + b[3]*H^3 + b[4]*H^4 + b[5]*H^5 + b[6]*H^6 + b[7]*H^7 + b[8]*H^8 + b[9]*H^9,
        a ~ dnorm(60,10),
        b ~ dnorm(0,10),
        sigma ~ dunif(0,10)
    ), data=dat , start=list(b=rep(0,9)) )

# plot sample
col2 <- col.alpha(2,0.6)
col4 <- col.alpha(4,0.6)
plot( d$height , d$weight , col=col2 , lwd=3 , cex=1.2 , xlab="height (cm)" , ylab="weight (kg)" , xlim=c(30,200) , ylim=c(0,75) )

mplot <- m_cubic
#mplot <- m_parabola
mplot <- m_crazy

xseq <- seq(from=-4,to=4,len=150)

mu <- link( mplot , data=list(H=xseq) )
Hseq <- xseq*sd(d$height) + mean(d$height)
lines( Hseq , apply(mu,2,mean) , lwd=4 )
shade( apply(mu,2,PI,0.99) , Hseq )

wsim <- sim( mplot , data=list(H=xseq) )
lines( Hseq , apply(wsim,2,PI,0.89)[1,] , lwd=2 , lty=2 )
lines( Hseq , apply(wsim,2,PI,0.89)[2,] , lwd=2 , lty=2 )

#################
# log(W) model


data(Howell1)
d <- Howell1

dat <- list(
    logW = log(d$weight),
    H = d$height,
    Hbar=mean(d$height) ) # S=1 female, S=2 male

m_logW <- quap(
    alist(
        logW ~ dnorm(mu,sigma),
        mu <- a + b*(H-Hbar),
        a ~ dnorm(60,10),
        b ~ dnorm(0,10),
        sigma ~ dunif(0,10)
    ), data=dat )

# plot sample
col2 <- col.alpha(2,0.6)
col4 <- col.alpha(4,0.6)
plot( d$height , d$weight , col=col2 , lwd=3 , cex=1.2 , xlab="height (cm)" , ylab="weight (kg)"  )

Hseq <- seq(from=50,to=180,len=30)
mu <- link(m_logW,data=list(H=Hseq,Hbar=mean(d$height)))
lines( Hseq , apply(exp(mu),2,mean) , lwd=3 )


#######
# height and age as piecewise linear

data(Howell1)
d <- Howell1

dat <- list(
    H = d$height,
    A = d$age )

m <- quap(
    alist(
        H ~ dnorm(mu,sigma),
        mu <- a1*(1-exp(-b1*A)) + a2*(1-exp(-b2*(A))
    )
)

