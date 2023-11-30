data {
    int  N_spp;
    vector[N_spp] B;
    vector[N_spp] M;
    vector[N_spp] G;
    matrix[N_spp, N_spp] R;
}
parameters {
    real a;
    real bG;
    real bM;
    real<lower=0> sigma;
}
transformed parameters {
    real<lower=0> sigma_sq;
    sigma_sq = square(sigma);
    vector[N_spp] mu; 
    matrix[N_spp, N_spp] K;
    // Linear model
    mu = a + bG * G + bM * M;
    // Kernel matrix
    K = R * sigma_sq;
}
model { 
    bM ~ normal(0, 0.5);
    bG ~ normal(0, 0.5);
    a ~ normal(0,1);
    sigma ~ exponential(1);
    B ~ multi_normal(mu, K);
}

