# nvm_suv
multi variate frailty survival 
Based on this artical https://pubmed.ncbi.nlm.nih.gov/20707868/

Perform survival analysis on N cluster M subject and K competing events.
Each cluzer has Q covariates and a multivariate normal distribution frailty vector


Main function is:
sur_mvn(X,Z,delta)

X[N,M] event times
Z[N,Q] covarites. Must be a two dimensional array.
delta[N,M,K] event witnesed 

Returns:
Coef
Var cov matrix of frailty
Culmative hazared




see nvm_example.R for basic usuage

Cuda support also exist but currently not in master
