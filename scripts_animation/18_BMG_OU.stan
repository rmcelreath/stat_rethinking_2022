functions{

    matrix cov_GPL1(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * x[i,j] );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }

    vector merge_missing( int[] miss_indexes , vector x_obs , vector x_miss ) {
        int N = dims(x_obs)[1];
        int N_miss = dims(x_miss)[1];
        vector[N] merged;
        merged = x_obs;
        for ( i in 1:N_miss )
            merged[ miss_indexes[i] ] = x_miss[i];
        return merged;
    }
}
data{
    int N_spp;

    vector[N_spp] B;
    int N_B_miss;
    int B_missidx[N_B_miss];

    vector[N_spp] M;
    int N_M_miss;
    int M_missidx[N_M_miss];

    vector[N_spp] G;
    int N_G_miss;
    int G_missidx[N_G_miss];

    matrix[N_spp,N_spp] Dmat;
}
parameters{
    real a;
    real aG;
    real bG;
    real bM;
    real bMG;
    real<lower=0> etasq;
    real<lower=0> rho;
    real<lower=0> etasqG;
    real<lower=0> rhoG;
    real<lower=0> etasqM;
    real<lower=0> rhoM;
    vector[N_M_miss] M_impute;
    vector[N_G_miss] G_impute;
    vector[N_B_miss] B_impute;
}
model{
    vector[N_spp] mu;
    vector[N_spp] nu;
    vector[N_spp] M_merge;
    vector[N_spp] G_merge;
    vector[N_spp] B_merge;

    matrix[N_spp,N_spp] K;
    matrix[N_spp,N_spp] KG;
    matrix[N_spp,N_spp] KM;

    rho ~ normal( 3 , 0.25 );
    etasq ~ normal( 1 , 0.25 );
    rhoG ~ normal( 3 , 0.25 );
    etasqG ~ normal( 1 , 0.25 );
    rhoM ~ normal( 3 , 0.25 );
    etasqM ~ normal( 1 , 0.25 );

    K = cov_GPL1(Dmat, etasq, rho, 0.01);
    KG = cov_GPL1(Dmat, etasqG, rhoG, 0.01);
    KM = cov_GPL1(Dmat, etasqM, rhoM, 0.01);

    bM ~ normal( 0 , 0.5 );
    bG ~ normal( 0 , 0.5 );
    bMG ~ normal( 0 , 0.5 );
    a ~ normal( 0 , 1 );
    aG ~ normal( 0 , 1 );

    G_merge = merge_missing(G_missidx, to_vector(G), G_impute);
    M_merge = merge_missing(M_missidx, to_vector(M), M_impute);
    B_merge = merge_missing(B_missidx, to_vector(B), B_impute);

    for ( i in 1:N_spp ) {
        mu[i] = a + bM * M_merge[i] + bG * G_merge[i];
        nu[i] = aG + bMG * M_merge[i];
    }
    M_merge ~ multi_normal( rep_vector(0,N_spp) , KM );
    G_merge ~ multi_normal( nu , KG );
    B_merge ~ multi_normal( mu , K );
}
generated quantities {
    vector[N_spp] M_est;
    vector[N_spp] G_est;
    vector[N_spp] B_est;
    G_est = merge_missing(G_missidx, to_vector(G), G_impute);
    M_est = merge_missing(M_missidx, to_vector(M), M_impute);
    B_est = merge_missing(B_missidx, to_vector(B), B_impute);
}
