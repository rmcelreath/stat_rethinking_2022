# lecture 18
# missing data

library(rethinking)
library(animation)
library(ellipse)

##################################################################
# Dogs

# Dog eats random homework
# aka missing completely at random
N <- 100
S <- rnorm(N)
H <- rnorm(N,0.5*S)
# dog eats 50% of homework at random
D <- rbern(N,0.5) 
Hstar <- H
Hstar[D==1] <- NA

plot( S , H , col=grau(0.8) , lwd=2 )
points( S , Hstar , col=2 , lwd=3 )

abline( lm(H~S) , lwd=3 , col=grau(0.8) )
abline( lm(Hstar~S) , lwd=3 , col=2 )

# Dog eats homework of students who study too much
# aka missing at random
N <- 100
S <- rnorm(N)
H <- rnorm(N,0.5*S)
# dog eats 80% of homework where S>0
D <- rbern(N, ifelse(S>0,0.8,0) ) 
Hstar <- H
Hstar[D==1] <- NA

plot( S , H , col=grau(0.8) , lwd=2 )
points( S , Hstar , col=2 , lwd=3 )

abline( lm(H~S) , lwd=3 , col=grau(0.8) )
abline( lm(Hstar~S) , lwd=3 , col=2 )

# Dog eats homework of students who study too much
# BUT NOW NONLINEAR WITH CEILING EFFECT
N <- 100
S <- rnorm(N)
H <- rnorm(N,(1-exp(-0.7*S)))
# dog eats 100% of homework where S>0
D <- rbern(N, ifelse(S>0,1,0) ) 
Hstar <- H
Hstar[D==1] <- NA

plot( S , H , col=grau(0.8) , lwd=2 )
points( S , Hstar , col=2 , lwd=3 )

abline( lm(H~S) , lwd=3 , col=grau(0.8) )
abline( lm(Hstar~S) , lwd=3 , col=2 )

# Dog eats bad homework
# aka missing not at random
N <- 100
S <- rnorm(N)
H <- rnorm(N,0.5*S)
# dog eats 90% of homework where H<0
D <- rbern(N, ifelse(H<0,0.9,0) ) 
Hstar <- H
Hstar[D==1] <- NA

plot( S , H , col=grau(0.8) , lwd=2 )
points( S , Hstar , col=2 , lwd=3 )

abline( lm(H~S) , lwd=3 , col=grau(0.8) )
abline( lm(Hstar~S) , lwd=3 , col=2 )

##################################################################
# primates phylogeny - now with imputation

## R code 14.47
library(rethinking)
data(Primates301)
data(Primates301_nex)
d <- Primates301

# plot it using ape package - install.packages('ape') if needed
library(ape)
d <- Primates301
l <- ladderize(Primates301_nex)

cols <- rep(1,301)
cols[1:30] <- 4 # NWM
cols[31:106] <- 2 # lemurs
cols[107:173] <- 3
cols[174:198] <- 5 # apes
cols[199:208] <- 3
cols[209:286] <- 4
cols[287:301] <- 7
cols[164:173] <- 7

par(bg="black")
plot( l , type="fan" , font=1 , no.margin=TRUE ,
    label.offset=1 , cex=0.55 , edge.width=1.6 , edge.col="white" , tip.col=cols , rotate.tree=37 )

# missingness on tree
par(bg="black")
xx <- plot( l , type="fan" , font=1 , no.margin=TRUE ,
    label.offset=1 , cex=0.8 , edge.width=1.6 , edge.col="white" , tip.col=cols , rotate.tree=30 , show.tip.label=FALSE )

idx <- sapply( l$tip.label , function(spp) which( d$name == spp ) )
B <- !is.na(d$brain)[idx]
M <- !is.na(d$body)[idx]
G <- !is.na(d$group_size)[idx]
tiplabels( pch=16 , lwd=2 , col=ifelse(B==1,cols,grau(0.01)) , cex=1 )
tiplabels( pch=16 , lwd=2, col=ifelse(M==1,cols,grau(0.01)), cex=1 , offset=2.2 )
tiplabels( pch=16 , lwd=2, col=ifelse(G==1,cols,grau(0.01)), cex=1 , offset=4.5 )

# complete cases only
par(bg="black")
xx <- plot( l , type="fan" , font=1 , no.margin=TRUE ,
    label.offset=1 , cex=0.8 , edge.width=1.6 , edge.col="white" , tip.col=cols , rotate.tree=30 , show.tip.label=FALSE )
cccols <- ifelse( B==1 & M==1 & G==1 , cols , grau(0.01) )
tiplabels( pch=16 , lwd=2 , col=ifelse(B==1,cccols,grau(0.01)) , cex=1 )
tiplabels( pch=16 , lwd=2, col=ifelse(M==1,cccols,grau(0.01)), cex=1 , offset=2.2 )
tiplabels( pch=16 , lwd=2, col=ifelse(G==1,cccols,grau(0.01)), cex=1 , offset=4.5 )

## R code 14.48
d <- Primates301
d$name <- as.character(d$name)
dstan <- d[ complete.cases( d$group_size , d$body , d$brain ) , ]
spp_obs <- dstan$name

# tree again for complete cases only
lcc <- keep.tip( l , spp_obs )
idx <- sapply( lcc$tip.label , function(spp) which( l$tip.label == spp ) )
colscc <- cols[ idx ]
par(bg="black")
plot( lcc , type="fan" , font=1 , no.margin=TRUE ,
    label.offset=1 , cex=0.8 , edge.width=1.6 , edge.col="white" , tip.col=colscc , rotate.tree=30 )

# tree again highlighting complete cases
colsh <- rep(grau(),301)
idx2 <- which( l$tip.label %in% spp_obs )
colsh <- ifelse( l$tip.label %in% spp_obs , cols , grau() )
par(bg="black")
plot( l , type="fan" , font=1 , no.margin=TRUE ,
    label.offset=1 , cex=0.55 , edge.width=1.6 , edge.col="white" , tip.col=colsh , rotate.tree=37 )

## R code 14.49
data(Primates301)
d <- Primates301
d$name <- as.character(d$name)
dstan <- d[ complete.cases( d$group_size , d$body , d$brain ) , ]
spp_cc <- dstan$name

# complete case analysis
dstan <- d[complete.cases(d$group_size,d$body,d$brain),]
dat_cc <- list(
    N_spp = nrow(dstan),
    M = standardize(log(dstan$body)),
    B = standardize(log(dstan$brain)),
    G = standardize(log(dstan$group_size)),
    Imat = diag(nrow(dstan)) )

# drop just missing brain cases
dd <- d[complete.cases(d$brain),]
dat_all <- list(
    N_spp = nrow(dd),
    M = standardize(log(dd$body)),
    B = standardize(log(dd$brain)),
    G = standardize(log(dd$group_size)),
    Imat = diag(nrow(dd)) )

dd <- d[complete.cases(d$brain),]
table( M=!is.na(dd$body) , G=!is.na(dd$group_size) )

# now with phylogeny

library(ape)
spp <- as.character(dd$name)
tree_trimmed <- keep.tip( Primates301_nex, spp )
Rbm <- corBrownian( phy=tree_trimmed )
V <- vcv(Rbm)
Dmat <- cophenetic( tree_trimmed )

# distance matrix
dat_all$Dmat <- Dmat[ spp , spp ] / max(Dmat)
dat_cc$Dmat <- Dmat[ spp_obs , spp_obs ] / max(Dmat)

# imputation ignoring models of M and G
fMBG_OU <- alist(
    B ~ multi_normal( mu , K ),
    mu <- a + bM*M + bG*G,
    M ~ normal(0,1),
    G ~ normal(0,1),
    matrix[N_spp,N_spp]:K <- cov_GPL1(Dmat,etasq,rho,0.01),
    a ~ normal( 0 , 1 ),
    c(bM,bG) ~ normal( 0 , 0.5 ),
    etasq ~ half_normal(1,0.25),
    rho ~ half_normal(3,0.25)
)
mBMG_OU <- ulam( fMBG_OU , data=dat_all , chains=4 , cores=4 , sample=TRUE )

# complete case version
mBMG_OU_cc <- ulam( fMBG_OU , data=dat_cc , chains=4 , cores=4 , sample=TRUE )

precis(mBMG_OU_cc)
precis(mBMG_OU)

post <- extract.samples(mBMG_OU)

# show M imputed values against regression of B on M
plot( dat_all$M , dat_all$B , lwd=3 , col=grau(0.2) , xlab="body mass (standardized)" , ylab="brain volume (standardized)" )
points( apply(post$M_impute,2,mean) , dat_all$B[mBMG_OU@data$M_missidx] , lwd=3 , col=2 )
for ( i in 1:2 ) {
    y <- dat_all$B[mBMG_OU@data$M_missidx][i]
    lines( PI(post$M_impute[,i]) , c(y,y) , lwd=8 , col=col.alpha(2,0.7) )
}

# show relation between G estimates and M
Gest <- apply(post$G_impute,2,mean)
idx <- which(is.na(dat_all$G))
plot( dat_cc$M , dat_cc$G , lwd=2 , col=grau() , xlab="Body mass (standardized)" , ylab="Group size (standardized)" )
points( dat_all$M[idx] , Gest , lwd=3 , col=2 )

# compare posterior bG of complete case and imputation models
postcc <- extract.samples(mBMG_OU_cc)
dens( postcc$bG , lwd=3 , col=grau(0.8) , xlab="effect of G on B" , ylim=c(0,25) )
dens( post$bG , lwd=3 , col=2 , add=TRUE )
abline(v=0,lty=3)

# no phylogeny on G but have submodel M -> G
mBMG_OU_G <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + bM*M + bG*G,
        G ~ normal(nu,sigma),
        nu <- aG + bMG*M,
        M ~ normal(0,1),
        matrix[N_spp,N_spp]:K <- cov_GPL1(Dmat,etasq,rho,0.01),
        c(a,aG) ~ normal( 0 , 1 ),
        c(bM,bG,bMG) ~ normal( 0 , 0.5 ),
        c(etasq) ~ half_normal(1,0.25),
        c(rho) ~ half_normal(3,0.25),
        sigma ~ exponential(1)
    ), data=dat_all , chains=4 , cores=4 , sample=TRUE )

# phylogeny information for G imputation (but no M -> G model)
mBMG_OU2 <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + bM*M + bG*G,
        M ~ normal(0,1),
        G ~ multi_normal( 'rep_vector(0,N_spp)' ,KG),
        matrix[N_spp,N_spp]:K <- cov_GPL1(Dmat,etasq,rho,0.01),
        matrix[N_spp,N_spp]:KG <- cov_GPL1(Dmat,etasqG,rhoG,0.01),
        a ~ normal( 0 , 1 ),
        c(bM,bG) ~ normal( 0 , 0.5 ),
        c(etasq,etasqG) ~ half_normal(1,0.25),
        c(rho,rhoG) ~ half_normal(3,0.25)
    ), data=dat_all , chains=4 , cores=4 , sample=TRUE )

mBMG_OU2_cc <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + bM*M + bG*G,
        M ~ normal(0,1),
        G ~ multi_normal( 'rep_vector(0,N_spp)' ,KG),
        matrix[N_spp,N_spp]:K <- cov_GPL1(Dmat,etasq,rho,0.01),
        matrix[N_spp,N_spp]:KG <- cov_GPL1(Dmat,etasqG,rhoG,0.01),
        a ~ normal( 0 , 1 ),
        c(bM,bG) ~ normal( 0 , 0.5 ),
        c(etasq,etasqG) ~ half_normal(1,0.25),
        c(rho,rhoG) ~ half_normal(3,0.25)
    ), data=dat_cc , chains=4 , cores=4 , sample=TRUE )

precis(mBMG_OU2_cc)
precis(mBMG_OU2)
precis(mBMG_OU_G)

# show relation between G estimates and M
post <- extract.samples(mBMG_OU_G)
post2 <- extract.samples(mBMG_OU2)
Gest <- apply(post$G_impute,2,mean)
Gest2 <- apply(post2$G_impute,2,mean)
idx <- which(is.na(dat_all$G))
plot( dat_cc$M , dat_cc$G , lwd=2 , col=grau() , xlab="Body mass (standardized)" , ylab="Group size (standardized)" )
points( dat_all$M[idx] , Gest , lwd=3 , col=2 ) # just M -> G
points( dat_all$M[idx] , Gest2 , lwd=3 , col=4 ) # just phylogeny

# GP kernel for group size
plot( NULL , xlim=c(0,1) , ylim=c(0,1) ,  xlab="phylogenetic distance" , ylab="G covariance" )
for ( i in 1:30 ) curve( post2$etasqG[i]*exp(-post2$rhoG[i]*x) , add=TRUE , col=4 , lwd=2 )

# add in G sub-model
mBMG_OU3 <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + bM*M + bG*G,
        G ~ multi_normal( nu , KG ),
        nu <- aG + bMG*M,
        M ~ normal(0,1),
        matrix[N_spp,N_spp]:K <- cov_GPL1(Dmat,etasq,rho,0.01),
        matrix[N_spp,N_spp]:KG <- cov_GPL1(Dmat,etasqG,rhoG,0.01),
        c(a,aG) ~ normal( 0 , 1 ),
        c(bM,bG,bMG) ~ normal( 0 , 0.5 ),
        c(etasq,etasqG) ~ half_normal(1,0.25),
        c(rho,rhoG) ~ half_normal(3,0.25)
    ), data=dat_all , chains=4 , cores=4 , sample=TRUE )

# complete case comparison
mBMG_OU3_cc <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + bM*M + bG*G,
        G ~ multi_normal( nu , KG ),
        nu <- aG + bMG*M,
        M ~ normal(0,1),
        matrix[N_spp,N_spp]:K <- cov_GPL1(Dmat,etasq,rho,0.01),
        matrix[N_spp,N_spp]:KG <- cov_GPL1(Dmat,etasqG,rhoG,0.01),
        c(a,aG) ~ normal( 0 , 1 ),
        c(bM,bG,bMG) ~ normal( 0 , 0.5 ),
        c(etasq,etasqG) ~ half_normal(1,0.25),
        c(rho,rhoG) ~ half_normal(3,0.25)
    ), data=dat_cc , chains=4 , cores=4 , sample=TRUE )

precis(mBMG_OU3)
precis(mBMG_OU3_cc)

# show relation between G estimates and M
post <- extract.samples(mBMG_OU_G)
post2 <- extract.samples(mBMG_OU2)
post3 <- extract.samples(mBMG_OU3)
Gest <- apply(post$G_impute,2,mean)
Gest2 <- apply(post2$G_impute,2,mean)
Gest3 <- apply(post3$G_impute,2,mean)
idx <- which(is.na(dat_all$G))
plot( dat_cc$M , dat_cc$G , lwd=2 , col=grau() , xlab="Body mass (standardized)" , ylab="Group size (standardized)" )
points( dat_all$M[idx] , Gest , lwd=3 , col=2 ) # just M -> G
points( dat_all$M[idx] , Gest2 , lwd=3 , col=4 ) # just phylogeny
points( dat_all$M[idx] , Gest3 , lwd=3 , col=6 ) # M->G and phylogeny

# GP kernel for group size
plot( NULL , xlim=c(0,1) , ylim=c(0,1) )
for ( i in 1:30 ) curve( post3$etasqG[i]*exp(-post3$rhoG[i]*x) , add=TRUE , col=grau() )

precis(lm(G ~ M ,dat_all))

# posterior bG
dens( post$bG , lwd=3 , col=2 , xlim=c(-0.06,0.2) , xlab="effect of G on B (bG)" )
dens( post2$bG , lwd=3 , col=4 , add=TRUE )
dens( post3$bG , lwd=3 , col=6 , add=TRUE )

post0 <- extract.samples(mBMG0cc)
dens( post0$bG , lwd=3 , add=TRUE )

# dens( postv$bG , add=TRUE , lwd=4 )

# model thinks about half of association G ~ M is phylogeny, so imputed values don't move much
# but to get model right, we need to include *either* sub-model for M -> G *or* phylogeny on G, otherwise imputed values don't have right association with M

# now model that imputes missing B values and uses phylogeny for structure
# do this in raw Stan for control

fna <- function(x) ifelse( is.na(x) , -99 , x )
dat_voll <- list(
    N_spp = nrow(d),
    M = fna( standardize(log(d$body)) ),
    B = fna( standardize(log(d$brain)) ),
    G = fna( standardize(log(d$group_size)) ),
    N_B_miss = sum(is.na(d$brain)) ,
    N_G_miss = sum(is.na(d$group_size)) ,
    N_M_miss = sum(is.na(d$body)) ,
    B_missidx = which( is.na(d$brain) ),
    G_missidx = which( is.na(d$group_size) ),
    M_missidx = which( is.na(d$body) )
)

spp_all <- as.character(d$name)
tree_voll <- keep.tip( Primates301_nex, spp_all )
Dmat <- cophenetic( tree_voll )
dat_voll$Dmat <- Dmat[ spp_all , spp_all ] / max(Dmat)

mv <- cstan( file="18_BMG_OU.stan" , data=dat_voll , chains=4 , cores=4 )

# complete case comparison
mBMG_OU4_cc <- ulam(
    alist(
        B ~ multi_normal( mu , K ),
        mu <- a + bM*M + bG*G,
        G ~ multi_normal( nu , KG ),
        nu <- aG + bMG*M,
        M ~ multi_normal('rep_vector(0,N_spp)',KM),
        matrix[N_spp,N_spp]:K <- cov_GPL1(Dmat,etasq,rho,0.01),
        matrix[N_spp,N_spp]:KG <- cov_GPL1(Dmat,etasqG,rhoG,0.01),
        matrix[N_spp,N_spp]:KM <- cov_GPL1(Dmat,etasqM,rhoM,0.01),
        c(a,aG) ~ normal( 0 , 1 ),
        c(bM,bG,bMG) ~ normal( 0 , 0.5 ),
        c(etasq,etasqG,etasqM) ~ half_normal(1,0.25),
        c(rho,rhoG,rhoM) ~ half_normal(3,0.25)
    ), data=dat_cc , chains=4 , cores=4 , sample=TRUE )

precis(mv)
precis(mBMG_OU4_cc) # complete case comparison

postv <- extract.samples(mv)
post4cc <- extract.samples(mBMG_OU4_cc)

# posterior coefficients - compare full to cc

dens( postv$bG , lwd=4 , col=2 , xlab="effect of G on B" )
dens( post4cc$bG , lwd=4 , col=grau() , add=TRUE)
abline(v=0,lty=3)

dens( postv$bM , lwd=4 , col=2 , xlab="effect of M on B" )
dens( post4cc$bM , lwd=4 , col=grau() , add=TRUE)
abline(v=0,lty=3)

dens( postv$bMG , lwd=4 , col=2 , xlab="effect of M on G" )
dens( post4cc$bMG , lwd=4 , col=grau() , add=TRUE)
abline(v=0,lty=3)

# GP kernels

# B
plot( NULL , xlim=c(0,1) , ylim=c(0,1) , xlab="phylogenetic distance" , ylab="B covariance" )
for ( i in 1:30 ) curve( post4cc$etasq[i]*exp(-post4cc$rho[i]*x) , add=TRUE , col=grau() )
for ( i in 1:30 ) curve( postv$etasq[i]*exp(-postv$rho[i]*x) , add=TRUE , col=2 )

# M
plot( NULL , xlim=c(0,1) , ylim=c(0,1) , xlab="phylogenetic distance" , ylab="M covariance" )
for ( i in 1:30 ) curve( post4cc$etasqM[i]*exp(-post4cc$rhoM[i]*x) , add=TRUE , col=grau() )
for ( i in 1:30 ) curve( postv$etasqM[i]*exp(-postv$rhoM[i]*x) , add=TRUE , col=2 )

# G
plot( NULL , xlim=c(0,1) , ylim=c(0,1) , xlab="phylogenetic distance" , ylab="G covariance" )
for ( i in 1:30 ) curve( post4cc$etasqG[i]*exp(-post4cc$rhoG[i]*x) , add=TRUE , col=grau() )
for ( i in 1:30 ) curve( postv$etasqG[i]*exp(-postv$rhoG[i]*x) , add=TRUE , col=2 )


plot( apply(postv$M_est,2,mean) , apply(postv$B_est,2,mean) , col=ifelse(is.na(d$brain),2,1) )

plot( apply(postv$G_est,2,mean) , apply(postv$B_est,2,mean) , col=ifelse(is.na(d$brain),2,1) )

plot( apply(postv$M_est,2,mean) , apply(postv$G_est,2,mean) , col=ifelse(is.na(d$group_size),2,1) )

#################
# cat adoption example

L <- 0.25
curve( exp(-L*x) , from=0 , to=10 , lwd=4 , xlab="time" , ylab="proportion not yet adopted" , ylim=c(0,1) )

n <- 10
y <- rexp(n,L)
for ( i in 1:length(y) ) lines( c(y[i],y[i]) , c(0,exp(-L*y[i])) , lwd=4 , col=ifelse(y[i]>5,2,4) )

abline(v=5,lty=3,lwd=2)

# simulated example with analyses (wasn't in lecture)

N <- 100
L <- 0.15
y <- rexp(N,L)
k <- ifelse( y > 5 , 1 , 0 ) # censoring indicator
y_obs <- y
y_obs[k==1] <- 5

# show that ignoring censored values bad
dat0 <- list( y=y_obs[k==0] )
m0 <- ulam(
    alist(
        y ~ exponential(lambda),
        lambda ~ exponential(0.5)
    ) , data=dat0 , chains=4 , cores=4 )

# model that imputes censored values
dat <- list( y=y_obs , k=k , N=N )
mc1 <- "
data{
    int N;
    vector[N] y;
    int k[N];
}
parameters{
    real<lower=0> lambda;
    real<lower=5> y_cens[N];
}
model{
    lambda ~ exponential(0.5);
    for ( i in 1:N ) {
        if ( k[i]==0 ) {
            y[i] ~ exponential(lambda);
            y_cens[i] ~ normal(5,1);
        }
        if ( k[i]==1 )
            y_cens[i] ~ exponential(lambda);
    }
}
"
m1 <- cstan( model_code=mc1 , data=dat , chains=4 , cores=4 )

# model with marginalization
m2 <- ulam(
    alist(
        y|k==0 ~ exponential(lambda),
        y|k==1 ~ custom( exponential_lccdf( y | lambda ) ),
        lambda ~ exponential(0.5)
    ) , data=dat , chains=4 , cores=4 )

# m1 and m2 should arrive at same inference
precis(m0)
precis(m1)
precis(m2)
