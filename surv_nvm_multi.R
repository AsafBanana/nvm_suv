library(MASS)
library(survival)
library(abind)
library(mvtnorm)
library(Rcpp)
library(statmod)
library(doParallel)


N_P = 20
sourceCpp("comp_risk_real.cpp", verbose=TRUE)



hermite <- function (points, z) {
  p1 <- 1/pi^0.4
  p2 <- 0
  for (j in 1:points) {
    p3 <- p2
    p2 <- p1
    p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
  }
  pp <- sqrt(2 * points) * p2
  c(p1, pp)
}

gauss.hermite <- function (points, iterlim = 50) {
  x <- w <- rep(0, points)
  m <- (points + 1)/2
  for (i in 1:m) {
    z <- if (i == 1) 
      sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
    else if (i == 2) 
      z - sqrt(points)/z
    else if (i == 3 || i == 4) 
      1.9 * z - 0.9 * x[i - 2]
    else 2 * z - x[i - 2]
    for (j in 1:iterlim) {
      z1 <- z
      p <- hermite(points, z)
      z <- z1 - p[1]/p[2]
      if (abs(z - z1) <= 1e-15) 
        break
    }
    if (j == iterlim) 
      warning("iteration limit exceeded")
    x[points + 1 - i] <- -(x[i] <- z)
    w[i] <- w[points + 1 - i] <- 2/p[2]^2
  }
  r <- cbind(x * sqrt(2), w/sum(w))
  colnames(r) <- c("Points", "Weights")
  r
}


mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  
  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  
  ## rotate, scale, translate points
  eig <- eigen(sigma) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}



flatten_array <- function(arr,N_EVENTS) {
	ret = arr[,,1]
	for (i in seq(2,N_EVENTS)) {
		ret = rbind(ret,arr[,,i])
	}
	return(as.matrix(ret))
}

get_beta_Z <- function(beta_hat,Z,N,N_QUESTION,N_EVENTS) {

	ret = array(rep(0,N*N_QUESTION*N_EVENTS),c(N,N_QUESTION,N_EVENTS))
	for (i in seq(1,N_EVENTS)) {
		for (j in seq(1,N_QUESTION)) {
			ret[,j,i] = as.matrix(Z) %*% beta_hat[j,,i]
		}
	}
	return(ret)
}




cox_phase <- function(X,delta,Z,omega_hat,cores=1) {
	registerDoParallel(cores)
	N = dim(delta)[1]	
	N_EVENTS = dim(delta)[3]
	N_QUESTION = dim(delta)[2]
	if (length(dim(Z)) == 2) {
		N_COEF = dim(Z)[2]
	} else {
		N_COEF = dim(Z)[3]
	}
	stf_ret = list()
	ret_all = array(rep(0,N_QUESTION*N_EVENTS*(N_COEF+1)),dim=c(N_QUESTION,N_COEF,N_EVENTS))
	ret_hazared_all = array(rep(0,N*N_QUESTION*N_EVENTS),c(N,N_QUESTION,N_EVENTS))
	Q_res <- foreach (i=1:N_QUESTION)  %dopar% {
		ret = rep(0,N_EVENTS)
		ret_hazared = array(rep(0,N*N_EVENTS),c(N,N_EVENTS))
		stf_ret = list()
		for (j in seq(1,N_EVENTS)) {
			if (length(dim(Z)) == 2) {
				Zc  = Z
			} else {
				Zc = Z[,i]
			}
			df = data.frame(Zc)
			rv = log(omega_hat[,j])
			srv = Surv(X[,i],delta[,i,j])
			res.cox <- coxph(srv ~ . + offset(rv), data = df)
			ret[j] = res.cox$coefficients[1]
			fit.obj<-coxph.detail(res.cox)
			if (N_COEF == 1) {
				clam<-cumsum(fit.obj$hazard[]/exp(t(as.vector(res.cox$coef))%*%colMeans(Zc)))
			} else {
				clam<-cumsum(fit.obj$hazard[]/exp(res.cox$coef*mean(Zc)))
			}

			clam1<-c(0,clam)
			stf = stepfun(fit.obj$time,clam1)			
			ret_hazared[,j] = stf(X[,i])
			stf_ret[[j]] = list("x"=fit.obj$time,"y"=clam1)
		}
		return(list("coefs"=ret,"hazared_at_event"=ret_hazared,"stf"=stf_ret))
	}
	for (i in seq(1,N_QUESTION)) {
		ret_all[i,,] = as.vector(Q_res[[i]]$coefs)
		ret_hazared_all[,i,] = Q_res[[i]]$hazared_at_event
		stf_ret[[i]] = Q_res[[i]]$stf
	}
	return(list("coefs"=ret_all,"hazared_at_event"=ret_hazared_all,"stf"=stf_ret))
}



sur_mvn <- function(X,Z,delta,cores=1,max_itr = 100,con_thr=0.01) {
	N = dim(delta)[1]	
	N_EVENTS = dim(delta)[3]
	N_QUESTION = dim(delta)[2]
	if (length(dim(Z)) == 2) {
		N_COEF = dim(Z)[2]
	} else {
		N_COEF = dim(Z)[3]
	}
	delta_sums = apply(delta,3,rowSums)
	W = rep(1,N)
	answrd = matrix(rep(0,N*N_QUESTION),ncol=N_QUESTION) + 1
	skip = matrix(rep(0,N_EVENTS*N_QUESTION),ncol=N_QUESTION)
	conv = 1
	omega_hat = matrix(rep(1,N*N_EVENTS),ncol=N_EVENTS)
	beta_hat = array(rep(0,N_QUESTION*N_EVENTS*(N_COEF)),dim=c(N_QUESTION,N_COEF,N_EVENTS))
	omega_varcov_hat =  diag(N_EVENTS)

	fail =FALSE
	itrcnt = 0
	while (conv > con_thr & itrcnt < max_itr) {
		old_beta_hat = beta_hat
		old_omega_varcov_hat = omega_varcov_hat
		omega_hat_old = omega_hat
		tmp = try(expr = cox_phase(X,delta,Z,omega_hat,cores))
		if (class(tmp) == "try-error") {
			fail = TRUE		
			break
		}
		beta_hat = tmp$coefs
		hazared_at_event = tmp$hazared_at_event
		stf = tmp$stf
		beta_Z = get_beta_Z(beta_hat,Z,N,N_QUESTION,N_EVENTS)
		gh = mgauss.hermite(N_P,rep(0,N_EVENTS),omega_varcov_hat)
		omega_varcov_hat_tmp = E_phase_C(delta_sums,gh$points,gh$weights,flatten_array(hazared_at_event,N_EVENTS),flatten_array(beta_Z,N_EVENTS),answrd,W,skip)
		gh = mgauss.hermite(N_P,rep(0,N_EVENTS),omega_varcov_hat)
		omega_hat = omega_expectation_C(delta_sums,gh$points,gh$weights,flatten_array(hazared_at_event,N_EVENTS),flatten_array(beta_Z,N_EVENTS),answrd,skip)	
		omega_varcov_hat = omega_varcov_hat_tmp
		conv = sum(abs(old_beta_hat - beta_hat)) + sum(abs(old_omega_varcov_hat - omega_varcov_hat))
		print(sprintf("conv: %f",conv))
		if (is.na(conv)) {
			fail = TRUE		
			break
		}
		itrcnt = itrcnt + 1

	}
	return(list("coef"=beta_hat,"var_cov"=omega_varcov_hat,"hazared"=stf))
}







