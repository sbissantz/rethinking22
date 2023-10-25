#
# Lecture 16
#

# Load required packages
library("rethinking")
library("dagitty")
library("igraph")

# Data for example
data(KosterLeckie)

# Inspect the data
View(kl_dyads)

# (1) Estimand
# Tie_AB -> Giving_AB <- Tie_BA
# Adjusting for the Households
dag <- dagitty('dag {
Giving_AB [outcome,pos="0,0"]
Tie_BA [exposure,pos="0,-0.5"]
Tie_AB [exposure,pos="0,0.5"]
Household_A [adjusted,pos="-0.5,0"]
Household_B [adjusted,pos="0.5,0"]
Household_A -> Giving_AB
Household_A -> Tie_AB
Household_A -> Tie_BA
Household_B -> Giving_AB
Household_B -> Tie_AB
Household_B -> Tie_BA
Tie_AB -> Giving_AB
Tie_BA -> Giving_AB
}')
plot(dag)

# (2) Generative model (fake data)
#
# Reproducibility
set.seed(112)
#
# Number of households
N <- 25 
# Generate all combinations of the elements of x, i.e. seq_along(N) taken 2 at
# a time – and transpose to get a column vector
# Note: Each row is a dyad of households 1....25
# Note: Read "In row 1 (dyad 1) where houshold 1 and 2 are listed"
dyads <- t( combn(N,2) ) 
#
# Number of dyads
N_dyads <- nrow(dyads)
#
# Simulate reciprocal ties "friendships" 
# Note: 10% of dyads are friends
# Friendships share reciprocal ties 
f <- rbinom(N_dyads, 1, 0.1)
#
# Simulte directed ties (non reciprocal) for all individuals
# Sharing flows only in one direction
# e.g. son giving to their mothers (not vice versa)
alpha <- (-3) # plogis(-3) = 0.05 or 5% directed ties present
y <- matrix(NA, N, N) # Social network (i.e. all ties)
# Loop
for ( i in 1:N ) for ( j in 1:N ) {
   if( i != j )  {
    #directed tie from i to j
      ids <- sort( c(i,j) )
      the_dyad <- which( dyads[,1] == ids[1] & dyads[,2] == ids[2] )
      # 1. possibility: friendship; look up in the f vector; if 1; then tie 
      # 2. role the dice and give directed tie with prob. plogis(alpha) (~5%)
      p_tie <- f[the_dyad] + (1 - f[the_dyad]) * plogis(alpha)
    # "y" Is now the social network (y) that holds reciprocal and directed ties
      y[i,j] <- rbinom(1, 1, p_tie)
   }
}
# Giving from A to B and B to A
giftsAB <- rep(NA, N_dyads)
giftsBA <- rep(NA, N_dyads)
# Average poisson rate of giving 
# y=0: No friendship 0.5 (generalized giving, receiving)
# y=1: If friendship: 2  (tie specific giving, receiving)
# Note: 2 > 0.5, i.e. more giving in friendships
lambda <- log(c(0.5, 2))
# Loop
for( i in 1:N_dyads ) {
    A <- dyads[i,1] ; B <- dyads[i,2]
    giftsAB[i] <- rpois(1, exp(lambda[1+y[A,B]]))
    giftsBA[i] <- rpois(1, exp(lambda[1+y[B,A]]))
}
# Draw the simulated network
sng <- graph_from_adjacency_matrix(y)
lx <- layout_nicely(sng)
vcol <- "#DE536B"
plot(sng, layout=lx, vertex.color=vcol, vertex.size=10, vertex.label=NA,
     edge.color="black", edge.arrow.size=0.5, edge.curved=0.35, asp=0.9)

sim_data <- list( N_dyads = N_dyads, D = seq(N_dyads), GAB=giftsAB , GBA=giftsBA )

# (3) Statistical model
#

# Ulam
f_dyad <- alist(
  GAB ~ poisson( lambdaAB ),
  GBA ~ poisson( lambdaBA ),
  log(lambdaAB) <- a + T[D,1] ,
  log(lambdaBA) <- a + T[D,2] ,
  a ~ normal(0,1),
## dyad effects
  transpars> matrix[N_dyads,2]:T <-
  compose_noncentered( rep_vector(sigma_T,2) , L_Rho_T , Z ),
  matrix[2,N_dyads]:Z ~ normal( 0 , 1 ),
  cholesky_factor_corr[2]:L_Rho_T ~ lkj_corr_cholesky( 2 ),
  sigma_T ~ exponential(1),
## compute correlation matrix for dyads
  gq> matrix[2,2]:Rho_T <<- Chol_to_Corr( L_Rho_T )
)
mGD <- ulam( f_dyad , data=sim_data , chains=4 , cores=4 , iter=2000 )
precis(mGD, depth = 2)

rethinking::stancode(mGD)

# Stan
#
path <- "~/projects/rethinking22"
file <- file.path(path, "stan", "16", "1.stan")
mdl <- cmdstanr::cmdstan_model(file, pedantic=TRUE)
fit <- mdl$sample(data=sim_data)
# Note that the alpha parameters are on the log-odds scale!
fit$print(max_rows=200)
# Diagnostics
#
fit$sampler_diagnostics()
fit$cmdstan_diagnose()
fit$diagnostic_summary()

# Posterior
#
post <- fit$draws(format = "matrix") 
sigma <- fit$draws("sigma_T", format = "matrix") 
mean(sigma)

# Posteror mean T 
T <- fit$draws("T", format = "matrix") 
plot(density(colMeans( T[,1:300] )))
lines(density(colMeans( T[,301:600] )))




# Poster correlation within dyads 
Rho <- fit$draws("Rho_T", format = "matrix") 
plot(density(Rho[,2]), main = "Correlation within dyads")



# (4) Analyze sample
#






