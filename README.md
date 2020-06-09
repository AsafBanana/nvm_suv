multivariate frailty for clustered survival data with competing risks

Based on the following paper:  https://pubmed.ncbi.nlm.nih.gov/20707868/

Performs survival analysis based on N clusters M observations within each cluster and K competing events.

Each observation within each cluster is having Q covariates. 

Each cluster is having a K-dimensional frailty vector from a multivariate normal distribution

Main function is: sur_mvn(X,Z,delta)

X[N,M] event times

Z[N,Q] covarites. Must be a two dimensional array.

delta[N,M,K] observed events 

Returns:

Coef

Var cov matrix of frailty

Culmative hazared

see nvm_example.R for basic usage
