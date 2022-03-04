# week 8 social networks and gaussian processes

library(rethinking)
dat <- read.csv("week08_Monks.csv")

# 1

# dyad model
f_dyad <- alist(
    likeAB ~ binomial( 3 , pAB ),
    likeBA ~ binomial( 3 , pBA ),
    logit(pAB) <- a + T[D,1] ,
    logit(pBA) <- a + T[D,2] ,
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

n_dyads <- nrow(dat)
monk_dat <- list(
    N_dyads = nrow(dat),
    D = dat$dyad_id,
    likeAB = dat$like_AB,
    likeBA = dat$like_BA )

m1 <- ulam( f_dyad , data=monk_dat , chains=4 , cores=4 , iter=2000 )

# 2

# dyad model with dislike
f_dyad2 <- alist(
    likeAB ~ binomial( 3 , pAB ),
    likeBA ~ binomial( 3 , pBA ),
    logit(pAB) <- a[1] + T[D,1] ,
    logit(pBA) <- a[1] + T[D,2] ,

    dislikeAB ~ binomial( 3 , qAB ),
    dislikeBA ~ binomial( 3 , qBA ),
    logit(qAB) <- a[2] + T2[D,1] ,
    logit(qBA) <- a[2] + T2[D,2] ,
    vector[2]:a ~ normal(0,1),

    ## like dyad effects - non-centered
    transpars> matrix[N_dyads,2]:T <-
            compose_noncentered( rep_vector(sigma_T,2) , L_Rho_T , Z ),
    matrix[2,N_dyads]:Z ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_T ~ lkj_corr_cholesky( 2 ),
    sigma_T ~ exponential(1),

    ## dislike dyad effects - non-centered
    transpars> matrix[N_dyads,2]:T2 <-
            compose_noncentered( rep_vector(sigma_T2,2) , L_Rho_T2 , Z2 ),
    matrix[2,N_dyads]:Z2 ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_T2 ~ lkj_corr_cholesky( 2 ),
    sigma_T2 ~ exponential(1),

    ## compute correlation matrix for dyads
    gq> matrix[2,2]:Rho_T <<- Chol_to_Corr( L_Rho_T ),
    gq> matrix[2,2]:Rho_T2 <<- Chol_to_Corr( L_Rho_T2 )
)

n_dyads <- nrow(dat)
monk_dat <- list(
    N_dyads = nrow(dat),
    D = dat$dyad_id,
    likeAB = dat$like_AB,
    likeBA = dat$like_BA,
    dislikeAB = dat$dislike_AB,
    dislikeBA = dat$dislike_BA )

m2 <- ulam( f_dyad2 , data=monk_dat , chains=4 , cores=4 , iter=2000 )

precis( m2 , depth=3 , pars=c("Rho_T","Rho_T2") )

post <- extract.samples(m2)
quantile( post$Rho_T[,1,2] - post$Rho_T2[,1,2] )

dens( post$Rho_T[,1,2] - post$Rho_T2[,1,2] )

T_est <- apply(post$T,2:3,mean)
T2_est <- apply(post$T2,2:3,mean)

N <- 18
n_dyads <- choose(18,2)
dyads <- t(combn(18,2))
# convert to adjacency matrix
y_est <- matrix(NA,18,18)
y2_est <- matrix(NA,18,18)
for ( i in 1:n_dyads ) {
    y_est[ dyads[i,1] , dyads[i,2] ] <- T_est[i,1]
    y_est[ dyads[i,2] , dyads[i,1] ] <- T_est[i,2]
    y2_est[ dyads[i,1] , dyads[i,2] ] <- T2_est[i,1]
    y2_est[ dyads[i,2] , dyads[i,1] ] <- T2_est[i,2]
}#i

library(igraph)
# like
sng <- graph_from_adjacency_matrix(y_est)
lx <- layout_nicely(sng)
plot(sng , layout=lx , vertex.size=8 , edge.arrow.size=0.7 , edge.width=2 , edge.curved=0.35 , vertex.color=2 , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA )

# dislike
sng2 <- graph_from_adjacency_matrix(y2_est)
plot(sng2 , layout=lx , vertex.size=8 , edge.arrow.size=0.7 , edge.width=2 , edge.curved=0.35 , vertex.color=6 , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA )


# 3

# general receiving model
f_dyad3 <- alist(
    likeAB ~ binomial( 3 , pAB ),
    likeBA ~ binomial( 3 , pBA ),
    logit(pAB) <- a[1] + T[D,1] + R[B,1],
    logit(pBA) <- a[1] + T[D,2] + R[A,1],

    dislikeAB ~ binomial( 3 , qAB ),
    dislikeBA ~ binomial( 3 , qBA ),
    logit(qAB) <- a[2] + T2[D,1] + R[B,2],
    logit(qBA) <- a[2] + T2[D,2] + R[A,2],
    vector[2]:a ~ normal(0,1),

    ## dyad effects - non-centered
    transpars> matrix[N_dyads,2]:T <-
            compose_noncentered( rep_vector(sigma_T,2) , L_Rho_T , Z ),
    matrix[2,N_dyads]:Z ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_T ~ lkj_corr_cholesky( 2 ),
    sigma_T ~ exponential(1),

    ## dyad effects - non-centered
    transpars> matrix[N_dyads,2]:T2 <-
            compose_noncentered( rep_vector(sigma_T2,2) , L_Rho_T2 , Z2 ),
    matrix[2,N_dyads]:Z2 ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_T2 ~ lkj_corr_cholesky( 2 ),
    sigma_T2 ~ exponential(1),

    ## R matrix of receiving effects
    transpars> matrix[18,2]:R <-
            compose_noncentered( sigma_R , L_Rho_R , ZR ),
    matrix[2,18]:ZR ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_R ~ lkj_corr_cholesky( 2 ),
    vector[2]:sigma_R ~ exponential(1),

    ## compute correlation matrix for dyads
    gq> matrix[2,2]:Rho_T <<- Chol_to_Corr( L_Rho_T ),
    gq> matrix[2,2]:Rho_T2 <<- Chol_to_Corr( L_Rho_T2 ),
    gq> matrix[2,2]:Rho_R <<- Chol_to_Corr( L_Rho_R )
)

n_dyads <- nrow(dat)
monk_dat <- list(
    N_dyads = nrow(dat),
    D = dat$dyad_id,
    likeAB = dat$like_AB,
    likeBA = dat$like_BA,
    dislikeAB = dat$dislike_AB,
    dislikeBA = dat$dislike_BA,
    A = dat$A,
    B = dat$B )

m3 <- ulam( f_dyad3 , data=monk_dat , chains=4 , cores=4 , iter=2000 )

precis( m3 , depth=3 , pars=c("Rho_T","Rho_T2","sigma_R","Rho_R") )

plot(precis(m3 , depth=3 , pars="R"))

plot_precis <- 
function (x, y, pars, col=1 , col.ci = "black", xlab = "Value", add = FALSE, xlim = NULL, labels = rownames(x)[1:n], iwd=3  , ...) 
{
    if (!missing(pars)) {
        x <- x[pars, ]
    }
    n <- nrow(x)
    mu <- x[n:1, 1]
    left <- x[[3]][n:1]
    right <- x[[4]][n:1]
    set_nice_margins()
    labels <- labels[n:1]
    if (is.null(xlim)) 
        xlim <- c(min(left), max(right))
    if (add == FALSE) 
        dotchart(mu, labels = labels, xlab = xlab, xlim = xlim, lwd=iwd , 
            ...)
    else points(mu[n:1], n:1, ...)
    col <- rep(col,len=length(mu))
    for (i in 1:length(mu)) lines(c(left[i], right[i]), c(i, 
        i), lwd = iwd, col = col[i])
    if (add == FALSE) 
        abline(v = 0, lty = 1, col = col.alpha("black", 0.15))
}

plot_precis(precis(m3 , depth=3 , pars="R"),col=2:1,iwd=3)

# 4

# add factions
factions <- c( 1 ,1, 3, 2, 2, 2, 1, 4, 2, 4, 2, 1, 4, 1, 1, 1, 3, 3 )

monk_dat <- list(
    N_dyads = nrow(dat),
    D = dat$dyad_id,
    likeAB = dat$like_AB,
    likeBA = dat$like_BA,
    dislikeAB = dat$dislike_AB,
    dislikeBA = dat$dislike_BA,
    A = dat$A,
    B = dat$B,
    F = factions )

# block model

f_block <- alist(
    likeAB ~ binomial( 3 , pAB ),
    likeBA ~ binomial( 3 , pBA ),
    logit(pAB) <- a[1] + T[D,1] + bl[F[A],F[B]],
    logit(pBA) <- a[1] + T[D,2] + bl[F[B],F[A]],

    dislikeAB ~ binomial( 3 , qAB ),
    dislikeBA ~ binomial( 3 , qBA ),
    logit(qAB) <- a[2] + T2[D,1] + bd[F[A],F[B]],
    logit(qBA) <- a[2] + T2[D,2] + bd[F[B],F[A]],

    vector[2]:a ~ normal(0,1),
    matrix[4,4]:Zl ~ normal(0,1),
    matrix[4,4]:Zd ~ normal(0,1),
    c(tau_l,tau_d) ~ half_normal(0,1),
    transpars> matrix[4,4]:bl <<- Zl*tau_l,
    transpars> matrix[4,4]:bd <<- Zd*tau_d,

    ## like dyad effects - non-centered
    transpars> matrix[N_dyads,2]:T <-
            compose_noncentered( rep_vector(sigma_T,2) , L_Rho_T , Z ),
    matrix[2,N_dyads]:Z ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_T ~ lkj_corr_cholesky( 2 ),
    sigma_T ~ exponential(1),

    ## dislike dyad effects - non-centered
    transpars> matrix[N_dyads,2]:T2 <-
            compose_noncentered( rep_vector(sigma_T2,2) , L_Rho_T2 , Z2 ),
    matrix[2,N_dyads]:Z2 ~ normal( 0 , 1 ),
    cholesky_factor_corr[2]:L_Rho_T2 ~ lkj_corr_cholesky( 2 ),
    sigma_T2 ~ exponential(1),

    ## compute correlation matrix for dyads
    gq> matrix[2,2]:Rho_T <<- Chol_to_Corr( L_Rho_T ),
    gq> matrix[2,2]:Rho_T2 <<- Chol_to_Corr( L_Rho_T2 )
)

m4 <- ulam( f_block , data=monk_dat , chains=4 , cores=4 , iter=2000 )

precis( m4 , depth=3 , pars=c("a","Rho_T","sigma_T","Rho_T2","sigma_T2","bl","bd","tau_l","tau_d") )

post <- extract.samples(m4)
( BL <- round(apply(post$bl,2:3,mean),2) )
( BD <- round(apply(post$bd,2:3,mean),2) )

N <- 18
n_dyads <- choose(18,2)
dyads <- t(combn(18,2))
# convert to adjacency matrix
y_est <- matrix(NA,18,18)
F <- factions
for ( i in 1:n_dyads ) {
    A <- dyads[i,1]
    B <- dyads[i,2]
    y_est[ A , B ] <- mean( inv_logit( with( post , a[,1]+T[,i,1]+bl[,F[A],F[B]] ) ) )
    y_est[ B , A ] <- mean( inv_logit( with( post , a[,1]+T[,i,2]+bl[,F[B],F[A]] ) ) )
}#i

library(igraph)
# like
sng <- graph_from_adjacency_matrix(3*y_est/max(y_est,na.rm=TRUE))
lx <- layout_nicely(sng)
plot(sng , layout=lx , vertex.size=8 , edge.arrow.size=0.7 , edge.width=2 , edge.curved=0.35 , vertex.color=factions , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA )

# dislike
sng2 <- graph_from_adjacency_matrix(y2_est)
plot(sng2 , layout=lx , vertex.size=8 , edge.arrow.size=0.7 , edge.width=2 , edge.curved=0.35 , vertex.color=factions , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA )
