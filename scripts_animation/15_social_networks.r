# lecture 15
# social relations model
# instruments

library(rethinking)
library(animation)
library(ellipse)

data(KosterLeckie)
d1 <- kl_dyads
d2 <- kl_households

CA <- col.alpha
plot( d1$giftsAB , d1$giftsBA , lwd=3 , col=CA(2,0.7) , xlim=c(0,110) , ylim=c(0,110) , xlab="A gives B" , ylab="B gives A" )
abline( a=0 , b=1 , lty=3 )
precis(lm(giftsAB~giftsBA,d1))

####
# generative model - just dyad effects

# N households
N <- 25
dyads <- t(combn(N,2))
N_dyads <- nrow(dyads)

# simulate "friendships" in which ties are shared
f <- rbern(N_dyads,0.1) # 10% of dyads are friends

# now simulate directed ties for all individuals
# there can be ties that are not reciprocal
alpha <- (-3) # base rate of ties; -3 ~= 0.05
y <- matrix(NA,N,N) # matrix of ties
for ( i in 1:N ) for ( j in 1:N ) {
    if ( i != j ) {
        # directed tie from i to j
        ids <- sort( c(i,j) )
        the_dyad <- which( dyads[,1]==ids[1] & dyads[,2]==ids[2] )
        p_tie <- f[the_dyad] + (1-f[the_dyad])*inv_logit( alpha )
        y[i,j] <- rbern( 1 , p_tie )
    }
}#ij

# now simulate gifts
giftsAB <- rep(NA,N_dyads)
giftsBA <- rep(NA,N_dyads)
lambda <- log(c(0.5,2)) # rates of giving for y=0,y=1
for ( i in 1:N_dyads ) {
    A <- dyads[i,1]
    B <- dyads[i,2]
    giftsAB[i] <- rpois( 1 , exp( lambda[1+y[A,B]] ) )
    giftsBA[i] <- rpois( 1 , exp( lambda[1+y[B,A]] ) )
}

if ( FALSE ) {
# draw network
library(igraph)
sng <- graph_from_adjacency_matrix(y)
lx <- layout_nicely(sng)
vcol <- "#DE536B"
plot(sng , layout=lx , vertex.size=8 , edge.arrow.size=0.75 , edge.width=2 , edge.curved=0.35 , vertex.color=vcol , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA )

}

# analyze synthetic data
sim_data <- list(
    N_dyads = N_dyads,
    N_households = N,
    D = 1:N_dyads,
    HA = dyads[,1],
    HB = dyads[,2],
    GAB = giftsAB,
    GBA = giftsBA )

# dyad model
f_dyad <- alist(
    GAB ~ poisson( lambdaAB ),
    GBA ~ poisson( lambdaBA ),
    log(lambdaAB) <- a + T[D,1] ,
    log(lambdaBA) <- a + T[D,2] ,
    a ~ normal(0,1),

    ## dyad effects - non-centered
    transpars> matrix[N_dyads,2]:T <-
            compose_noncentered( rep_vector(sigma_T,2) , L_Rho_T , Z ),
    matrix[2,N_dyads]:Z ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_T ~ lkj_corr_cholesky( 2 ),
    sigma_T ~ exponential(1),

    ## compute correlation matrix for dyads
    gq> matrix[2,2]:Rho_T <<- Chol_to_Corr( L_Rho_T )
)

mGD <- ulam( f_dyad , data=sim_data , chains=4 , cores=4 , iter=2000 )

precis( mGD , depth=3 , pars=c("a","Rho_T","sigma_T") )

post <- extract.samples(mGD)
T_est <- apply(post$T,2:3,mean)
# convert to adjacency matrix
y_est <- y
for ( i in 1:N_dyads ) {
    y_est[ dyads[i,1] , dyads[i,2] ] <- T_est[i,1]
    y_est[ dyads[i,2] , dyads[i,1] ] <- T_est[i,2]
}#i

# show discrimination as densities for each true tie state
dens( y_est[y==0] , xlim=c(-0.5,1.5) , lwd=4 , col=2 , xlab="posterior mean T" )
dens( y_est[y==1] , add=TRUE , lwd=4 , col=4 )

# show correlation by true friend state
plot( T_est[,1] , T_est[,2] , lwd=3 , col=ifelse(f==1,6,1) , xlab="Household A" , ylab="Household B" )

# show reciprocity
dens( post$Rho_T[,1,2] , lwd=4 , col=2 , xlab="correlation within dyads" , xlim=c(-1,1) )

# analyze sample
kl_data <- list(
    N_dyads = nrow(kl_dyads),
    N_households = max(kl_dyads$hidB),
    D = 1:nrow(kl_dyads),
    HA = kl_dyads$hidA,
    HB = kl_dyads$hidB,
    GAB = kl_dyads$giftsAB,
    GBA = kl_dyads$giftsBA )

mGDkl <- ulam( f_dyad , data=kl_data , chains=4 , cores=4 , iter=4000 )

precis( mGDkl , depth=3 , pars=c("a","Rho_T","sigma_T") )

post <- extract.samples(mGDkl)
dens( post$Rho_T[,1,2] , lwd=4 , col=2 , xlab="correlation within dyads" , xlim=c(-1,1) )

################
# generative model - full model with generalized giving/receiving

N <- 25
dyads <- t(combn(N,2))
N_dyads <- nrow(dyads)
f <- rbern(N_dyads,0.1) # 10% of dyads are friends

# now simulate directed ties for all individuals
alpha <- (-3) # base rate of ties
y <- matrix(NA,N,N) # network matrix
for ( i in 1:N ) for ( j in 1:N ) {
    if ( i != j ) {
        # directed tie from i to j
        ids <- sort( c(i,j) )
        the_dyad <- which( dyads[,1]==ids[1] & dyads[,2]==ids[2] )
        p_tie <- f[the_dyad] + (1-f[the_dyad])*inv_logit( alpha )
        y[i,j] <- rbern( 1 , p_tie )
    }
}#ij

# simulate wealth
W <- rnorm(N) # standardized relative wealth in community
bWG <- 0.5 # effect of wealth on giving - rich give more
bWR <- (-1) # effect of wealth on receiving - rich get less / poor get more

# now simulate gifts
giftsAB <- rep(NA,N_dyads)
giftsBA <- rep(NA,N_dyads)
lambda <- log(c(0.5,2)) # rates of giving for y=0,y=1
for ( i in 1:N_dyads ) {
    A <- dyads[i,1]
    B <- dyads[i,2]
    giftsAB[i] <- rpois( 1 , exp( lambda[1+y[A,B]] + bWG*W[A] + bWR*W[B] ) )
    giftsBA[i] <- rpois( 1 , exp( lambda[1+y[B,A]] + bWG*W[B] + bWR*W[A] ) )
}

if ( FALSE ) {
library(igraph)
sng <- graph_from_adjacency_matrix(y)
lx <- layout_nicely(sng)
vcol <- "#DE536B"
plot( sng , layout=lx , vertex.size=8 , edge.arrow.size=0.75 , edge.width=2 , edge.curved=0.35 , vertex.color=vcol , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA  )
}

plot( jitter(giftsAB) , jitter(giftsBA) )

sim_data <- list(
    N_dyads = N_dyads,
    N_households = N,
    D = 1:N_dyads,
    HA = dyads[,1],
    HB = dyads[,2],
    GAB = giftsAB,
    GBA = giftsBA )

# general model
f_general <- alist(
    GAB ~ poisson( lambdaAB ),
    GBA ~ poisson( lambdaBA ),
    log(lambdaAB) <- a + T[D,1] + gr[HA,1] + gr[HB,2],
    log(lambdaBA) <- a + T[D,2] + gr[HB,1] + gr[HA,2],
    a ~ normal(0,1),

    ## dyad effects - non-centered
    transpars> matrix[N_dyads,2]:T <-
            compose_noncentered( rep_vector(sigma_T,2) , L_Rho_T , Z ),
    matrix[2,N_dyads]:Z ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_T ~ lkj_corr_cholesky( 2 ),
    sigma_T ~ exponential(1),

   ## gr matrix of varying effects
    transpars> matrix[N_households,2]:gr <-
            compose_noncentered( sigma_gr , L_Rho_gr , Zgr ),
    matrix[2,N_households]:Zgr ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_gr ~ lkj_corr_cholesky( 2 ),
    vector[2]:sigma_gr ~ exponential(1),

   ## compute correlation matrixes
    gq> matrix[2,2]:Rho_T <<- Chol_to_Corr( L_Rho_T ),
    gq> matrix[2,2]:Rho_gr <<- Chol_to_Corr( L_Rho_gr )
)

mGDGR <- ulam( f_general , data=sim_data , chains=4 , cores=4 , iter=2000 )

precis( mGDGR, depth=3 , pars=c("a","Rho_gr","sigma_gr","Rho_T","sigma_T") )

## R code 14.33
post <- extract.samples( mGDGR )
g <- sapply( 1:kl_data$N_households , function(i) post$a + post$gr[,i,1] )
r <- sapply( 1:kl_data$N_households , function(i) post$a + post$gr[,i,2] )
Eg_mu <- apply( exp(g) , 2 , mean )
Er_mu <- apply( exp(r) , 2 , mean )

## R code 14.34
xymax <- max(c(Eg_mu,Er_mu)) + 0.5
plot( NULL , xlim=c(0,xymax) , ylim=c(0,xymax) , xlab="generalized giving" ,
    ylab="generalized receiving" , lwd=2 , col=2 )
abline(a=0,b=1,lty=2)

# ellipses
library(ellipse)
for ( i in 1:kl_data$N_households ) {
    Sigma <- cov( cbind( g[,i] , r[,i] ) )
    Mu <- c( mean(g[,i]) , mean(r[,i]) )
    for ( l in c(0.5) ) {
        el <- ellipse( Sigma , centre=Mu , level=l )
        lines( exp(el) , col=grau() )
    }
}
# household means
points( Eg_mu , Er_mu , pch=21 , bg="white" , lwd=2 , col=2 )

T_est <- apply(post$T,2:3,mean)
# convert to adjacency matrix
y_est <- y
for ( i in 1:N_dyads ) {
    y_est[ dyads[i,1] , dyads[i,2] ] <- T_est[i,1]
    y_est[ dyads[i,2] , dyads[i,1] ] <- T_est[i,2]
}#i

# show discrimination as densities for each true tie state
dens( y_est[y==0] , xlim=c(-1,2) , lwd=4 , col=2 , xlab="posterior mean T" )
dens( y_est[y==1] , add=TRUE , lwd=4 , col=4 )

# show correlation by true friend state
plot( T_est[,1] , T_est[,2] , lwd=3 , col=ifelse(f==1,6,1) , xlab="Household A" , ylab="Household B" )

# show reciprocity
dens( post$Rho_T[,1,2] , lwd=4 , col=2 , xlab="correlation within dyads" , xlim=c(-1,1) )

# show correlation give-receive
dens( post$Rho_gr[,1,2] , lwd=4 , col=4 , xlab="correlation giving-receiving" , xlim=c(-1,1) )

####
# analyze real sample
kl_data <- list(
    N_dyads = nrow(kl_dyads),
    N_households = max(kl_dyads$hidB),
    D = 1:nrow(kl_dyads),
    HA = kl_dyads$hidA,
    HB = kl_dyads$hidB,
    GAB = kl_dyads$giftsAB,
    GBA = kl_dyads$giftsBA )

mGDGR <- ulam( f_general , data=kl_data , chains=4 , cores=4 , iter=2000 )

precis( mGDGR, depth=3 , pars=c("a","Rho_gr","sigma_gr","Rho_T","sigma_T") )


## R code 14.33
post <- extract.samples( mGDGR )
g <- sapply( 1:kl_data$N_households , function(i) post$a + post$gr[,i,1] )
r <- sapply( 1:kl_data$N_households , function(i) post$a + post$gr[,i,2] )
Eg_mu <- apply( exp(g) , 2 , mean )
Er_mu <- apply( exp(r) , 2 , mean )

## R code 14.34
xymax <- max(c(Eg_mu,Er_mu)) + 0.5
plot( NULL , xlim=c(0,xymax) , ylim=c(0,xymax) , xlab="generalized giving" ,
    ylab="generalized receiving" , lwd=2 , col=2 )
abline(a=0,b=1,lty=2)

# ellipses
library(ellipse)
for ( i in 1:kl_data$N_households ) {
    Sigma <- cov( cbind( g[,i] , r[,i] ) )
    Mu <- c( mean(g[,i]) , mean(r[,i]) )
    for ( l in c(0.5) ) {
        el <- ellipse( Sigma , centre=Mu , level=l )
        lines( exp(el) , col=grau() )
    }
}
# household means
points( Eg_mu , Er_mu , pch=21 , bg="white" , lwd=2 , col=2 )

T_est <- apply(post$T,2:3,mean)
# convert to adjacency matrix
y_est <- y
for ( i in 1:N_dyads ) {
    y_est[ dyads[i,1] , dyads[i,2] ] <- T_est[i,1]
    y_est[ dyads[i,2] , dyads[i,1] ] <- T_est[i,2]
}#i

# show correlation by true friend state
plot( T_est[,1] , T_est[,2] , lwd=3 , col=2 , xlab="Household A" , ylab="Household B" )

# show reciprocity
dens( post$Rho_T[,1,2] , lwd=4 , col=2 , xlab="correlation within dyads" , xlim=c(-1,1) )

# show correlation give-receive
dens( post$Rho_gr[,1,2] , lwd=4 , col=4 , xlab="correlation giving-receiving" , xlim=c(-1,1) )

# posterior mean social network
library(igraph)
gr <- graph_from_adjacency_matrix(y)
lx <- layout_nicely(gr)
vcol <- "#DE536B"


# animate network layout?
# start with posterior mean and get stable node positions
T_est <- apply(post$T,2:3,mean)
# convert to adjacency matrix
y_est <- y
for ( i in 1:N_dyads ) {
    y_est[ dyads[i,1] , dyads[i,2] ] <- T_est[i,1]
    y_est[ dyads[i,2] , dyads[i,1] ] <- T_est[i,2]
}#i
library(igraph)
sng <- graph_from_adjacency_matrix(y_est)
lx <- layout_nicely(sng)
vcol <- "#DE536B"
plot(sng , layout=lx , vertex.size=8 , edge.arrow.size=1 , edge.width=2 , edge.curved=0.35 , vertex.color=vcol , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA )

# draw 100 networks from posterior
ani.record(reset=TRUE)
T_est <- post$T
for ( f in 1:100 ) {

    y_est <- y
    s <- f
    for ( i in 1:N_dyads ) {
        y_est[ dyads[i,1] , dyads[i,2] ] <- T_est[s,i,1]
        y_est[ dyads[i,2] , dyads[i,1] ] <- T_est[s,i,2]
    }#i
    grp <- graph_from_adjacency_matrix(y_est)
    plot(grp , layout=lx , vertex.size=8 , edge.arrow.size=1 , edge.width=2 , edge.curved=0.35 , vertex.color=vcol , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA )

    ani.record()
}

oopts = ani.options(interval = 0.1)
ani.replay()

# ani.saveqz(dpi=150)
# convert -alpha remove -background white -delay 15 -loop 0 frame*.png kl_networks.gif
# convert -delay 15 kl_networks.gif kl_networks.gif
# convert gp_demo.gif -coalesce -fuzz 5% +dither -layers Optimize +map gp_demo_z.gif
# convert kl_networks.gif -coalesce -layers Optimize +map kl_networks_z.gif

#######
# kl_households features
# add household wealth to model

kl_data <- list(
    N_dyads = nrow(kl_dyads),
    N_households = max(kl_dyads$hidB),
    D = 1:nrow(kl_dyads),
    HA = kl_dyads$hidA,
    HB = kl_dyads$hidB,
    GAB = kl_dyads$giftsAB,
    GBA = kl_dyads$giftsBA,
    W = standardize(kl_households$hwealth), # wealth
    A = standardize(kl_dyads$dass) ) # production

# general model with features
f_houses <- alist(
    GAB ~ poisson( lambdaAB ),
    GBA ~ poisson( lambdaBA ),
    # A to B
    log(lambdaAB) <- a + TAB + GA + RB,
    TAB <- T[D,1] + bA*A,
    GA <- gr[HA,1] + bW[1]*W[HA] ,
    RB <- gr[HB,2] + bW[2]*W[HB] ,
    # B to A
    log(lambdaBA) <- a + TBA + GB + RA,
    TBA <- T[D,2] + bA*A,
    GB <- gr[HB,1] + bW[1]*W[HB] ,
    RA <- gr[HA,2] + bW[2]*W[HA] ,
    # priors
    a ~ normal(0,1),
    vector[2]:bW ~ normal(0,1),
    bA ~ normal(0,1),

    ## dyad effects - non-centered
    transpars> matrix[N_dyads,2]:T <-
            compose_noncentered( rep_vector(sigma_T,2) , L_Rho_T , Z ),
    matrix[2,N_dyads]:Z ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_T ~ lkj_corr_cholesky( 2 ),
    sigma_T ~ exponential(1),

   ## gr matrix of varying effects
    transpars> matrix[N_households,2]:gr <-
            compose_noncentered( sigma_gr , L_Rho_gr , Zgr ),
    matrix[2,N_households]:Zgr ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_gr ~ lkj_corr_cholesky( 2 ),
    vector[2]:sigma_gr ~ exponential(1),

   ## compute correlation matrixes
    gq> matrix[2,2]:Rho_T <<- Chol_to_Corr( L_Rho_T ),
    gq> matrix[2,2]:Rho_gr <<- Chol_to_Corr( L_Rho_gr )
)

mGDGRW <- ulam( f_houses , data=kl_data , chains=4 , cores=4 , iter=4000 )

precis( mGDGRW, depth=3 , pars=c("a","Rho_gr","sigma_gr","Rho_T","sigma_T","bW","bA") )

precis( mGDGR , depth=3 , pars=c("a","Rho_gr","sigma_gr","Rho_T","sigma_T") )

post1 <- extract.samples(mGDGR)
post2 <- extract.samples(mGDGRW)

dens( post1$sigma_T , lwd=4 , col=2 , xlab="standard deviation ties" , xlim=c(0.5,1.5) , ylim=c(0,8) )
dens( post2$sigma_T , lwd=4 , col=4 , add=TRUE )

dens( post2$bW[,1] , lwd=4 , col=2 , xlab="effect of wealth" , ylim=c(0,5) )
dens( post2$bW[,2] , lwd=4 , col=4 , add=TRUE )
abline(v=0,lty=3)

####
# plot raw network and posterior mean

post <- extract.samples(mGDGR)
T_est <- apply(post$T,2:3,mean)
# convert to adjacency matrix
y_est <- y
for ( i in 1:N_dyads ) {
    y_est[ dyads[i,1] , dyads[i,2] ] <- T_est[i,1]
    y_est[ dyads[i,2] , dyads[i,1] ] <- T_est[i,2]
}#i

library(igraph)
sng <- graph_from_adjacency_matrix(y_est)
lx <- layout_nicely(sng)
vcol <- "#DE536B" # 2
plot(sng , layout=lx , vertex.size=8 , edge.arrow.size=1 , edge.width=2 , edge.curved=0.35 , vertex.color=vcol , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA )

y_raw <- y_est
for ( i in 1:N_dyads ) {
    y_raw[ dyads[i,1] , dyads[i,2] ] <- kl_dyads$giftsAB[i]
    y_raw[ dyads[i,2] , dyads[i,1] ] <- kl_dyads$giftsBA[i]
}#i
y_raw <- y_raw / 6
sng2 <- graph_from_adjacency_matrix(y_raw)
vcol <- "#2297E6" # 4
plot(sng2 , layout=lx , vertex.size=8 , edge.arrow.size=1 , edge.width=2 , edge.curved=0.35 , vertex.color=vcol , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA )

