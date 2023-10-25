data{
     int N_dyads;
    array[300] int GAB;
    array[300] int GBA;
    array[300] int D;
}
parameters{
     real a;
     matrix[2,N_dyads] Z;
     cholesky_factor_corr[2] L_Rho_T;
     real<lower=0> sigma_T;
}
transformed parameters{
     matrix[N_dyads,2] T;
    T = (diag_pre_multiply(rep_vector(sigma_T, 2), L_Rho_T) * Z)';
}
model{
     vector[300] lambdaAB;
     vector[300] lambdaBA;
    sigma_T ~ exponential( 1 );
    L_Rho_T ~ lkj_corr_cholesky( 2 );
    to_vector( Z ) ~ normal( 0 , 1 );
    a ~ normal( 0 , 1 );
    for ( i in 1:300 ) {
        lambdaBA[i] = a + T[D[i], 2];
        lambdaBA[i] = exp(lambdaBA[i]);
    }
    for ( i in 1:300 ) {
        lambdaAB[i] = a + T[D[i], 1];
        lambdaAB[i] = exp(lambdaAB[i]);
    }
    GBA ~ poisson( lambdaBA );
    GAB ~ poisson( lambdaAB );
}
generated quantities{
     matrix[2,2] Rho_T;
    Rho_T = multiply_lower_tri_self_transpose(L_Rho_T);
}
