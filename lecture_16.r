#
# Lecture 16: Gaussian processes
# 
library(rethinking)

# Kline data
data(Kline2)
d <- Kline2
data(islandsDistMatrix)
# display (measured in thousands of km)
Dmat <- islandsDistMatrix
colnames(Dmat) <- abbreviate(rownames(Dmat), 2)
round(Dmat,1)

# Prior predictive simulation 
n <- 50 
etasq <- rexp(n,2) ; rhosq <- rexp(n,1)
plot( NULL , xlim=c(0,7) , ylim=c(0,2) , xlab="distance (thousand km)" , ylab="covariance" )
for ( i in 1:n ) {
    curve(etasq[i] * exp(-rhosq[i] * x^2), add = TRUE, lwd = 4, 
    col = "steelblue")
}

# Data list
dat_list <- list(
    N = nrow(d),
    T = d$total_tools,
    P = d$population,
    S = 1:10,
    D = islandsDistMatrix )

mTdist <- ulam(
    alist(
        T ~ dpois(lambda),
        log(lambda) <- abar + a[S],
        vector[10]:a ~ multi_normal( 0 , K ),
        transpars> matrix[10,10]:K <- cov_GPL2(D,etasq,rhosq,0.01),
        abar ~ normal(3,0.5),
        etasq ~ dexp( 2 ),
        rhosq ~ dexp( 0.5 )
    ), data=dat_list , chains=4 , cores=4 , iter=4000 , log_lik=TRUE )

stancode(mTdist)

#
# Stancode
#
# Centered parametrization
path <- "~/projects/rethinking22"
file <- file.path(path, "stan", "16", "1_cnt.stan")
mdl <- cmdstanr::cmdstan_model(file, pedantic=TRUE)
fit <- mdl$sample(data=dat_list)
# Note that the alpha parameters are on the log scale!
fit$print(max_rows=200)
# Diagnostics
#
fit$sampler_diagnostics()
fit$cmdstan_diagnose() # No divergences, non-centered not needed!
fit$diagnostic_summary()

# Note test if non-centered works worse!
mTdist_nc <- ulam(
    alist(
        T ~ dpois(lambda),
        log(lambda) <- abar + a[S],
        # non-centered Gaussian Process prior
        transpars> vector[10]: a <<- L_K * z,
        vector[10]: z ~ normal( 0 , 1 ),
        transpars> matrix[10,10]: L_K <<- cholesky_decompose( K ),
        transpars> matrix[10,10]: K <- cov_GPL2(D,etasq,rhosq,0.01),
        abar ~ normal(3,0.5),
        etasq ~ dexp( 2 ),
        rhosq ~ dexp( 0.5 )
    ), data=dat_list , chains=4 , cores=4 , iter=2000 , log_lik=TRUE )

stancode(mTdist_nc)

# Non-Centered parametrization (worse than centered)
path <- "~/projects/rethinking22"
file <- file.path(path, "stan", "16", "1_nct.stan")
mdl <- cmdstanr::cmdstan_model(file, pedantic=TRUE)
fit <- mdl$sample(data=dat_list)
fit$print(max_rows=200)

# Diagnostics
#
fit$sampler_diagnostics()
fit$cmdstan_diagnose() # Divergences, centered-parametrization better!
fit$diagnostic_summary()
