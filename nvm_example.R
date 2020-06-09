

source("surv_nvm_multi.R")
source("nvm_data.R")

aa = sur_mvn(data$X,data$Z,data$delta,cores=4)

print(aa$coef)
print(aa$var_cov)
