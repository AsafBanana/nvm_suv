#include <stdio.h>
#include <stdlib.h>
/**#include <math.h>**/
#include <Rcpp.h>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif

//#include "cod_complex.h"
using namespace Rcpp;

double* R2C_mat(NumericMatrix a) {
	int n = a.nrow();
	int m = a.ncol();
	double* a_c = (double*)malloc(sizeof(double)*n*m);
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
        		a_c[i*m + j] = a(i,j);
		}
  	}
	return a_c;
}


double* R2C_vec(NumericVector a) {
	int n = a.size();
	double* a_c = (double*)malloc(sizeof(double)*n);
	for(int i = 0; i < n; ++i) {
        	a_c[i] = a[i];
  	}
	return a_c;
}



void expected_value_at_point_c(double* __restrict__ rowSums, double* __restrict__ points, double* __restrict__ hazared_at_event, double* __restrict__  beta_Z, int N, int N_Q, int N_E, int N_P, double* __restrict__ ret) {
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N_P; j++) 
			ret[i*N_P + j] = 1;	
	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < N_P; k++) {
			double in_exp = 0;
			for (int j = 0; j < N_E; j++) {
				in_exp += points[k*N_E + j]*rowSums[i*N_E + j];
				for (int z = 0; z < N_Q; z++) {
					 in_exp -= hazared_at_event[j*N*N_Q + i*N_Q + z]*exp(beta_Z[j*N*N_Q + i*N_Q + z] + points[k*N_E + j]);
				}
			}
			ret[i*N_P + k] = exp(in_exp);
		}
	}
}




// [[Rcpp::export]]
NumericVector expected_value_at_point_R(NumericMatrix rowSums, NumericMatrix P , NumericMatrix hazared_at_event , NumericMatrix beta_Z) {
	int N = rowSums.nrow();
	int N_E = rowSums.ncol();
	int N_P = P.nrow();
	int N_Q = hazared_at_event.ncol();
	double* rowSums_c = R2C_mat(rowSums);
	double* P_c = R2C_mat(P);
	double* hazared_at_event_c = R2C_mat(hazared_at_event);
	double* beta_Z_c = R2C_mat(beta_Z);
	
	double* ret = (double*)malloc(sizeof(double)*N * N_P);
	expected_value_at_point_c(rowSums_c,P_c,hazared_at_event_c,beta_Z_c,N,N_Q,N_E,N_P,ret);
	NumericMatrix out(N,N_P);

	for(int i = 0; i <N; ++i) {
		for(int j = 0; j  < N_P; ++j) {
	    		out(i,j) = ret[i*N_P + j];
		}
	}
	free(rowSums_c);
	free(P_c);
	free(ret);
	free(hazared_at_event_c);
	free(beta_Z_c);
	return out;
}


void omega_expectation(double* __restrict__ rowSums, double* __restrict__ points, double* __restrict__ weights, double* __restrict__ hazared_at_event, double* __restrict__  beta_Z, int N, int N_Q, int N_E, int N_P, double* ans, double* skip, double* __restrict__ ret) {
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N_E; j++) 
			ret[i*N_E + j] = 0;	
	double devisor[N];
	for (int i = 0; i < N; i++)
		devisor[i] = 0;
	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < N_P; k++) {
			double integrand = 0;
			for (int j = 0; j < N_E; j++) {
				double in_exp = points[k*N_E + j]*rowSums[i*N_E + j];
				for (int z = 0; z < N_Q; z++) {
					if (ans[i*N_Q + z] == 0 || skip[j*N_Q + z] == 1) 
					//if (ans[i*N_Q + z] == 0) 
						continue;
					 in_exp -= hazared_at_event[j*N*N_Q + i*N_Q + z]*exp(beta_Z[j*N*N_Q + i*N_Q + z] + points[k*N_E + j]);
				}
				integrand += in_exp;
			}
			integrand = exp(integrand)*weights[k];
			devisor[i] += integrand;
			for (int j = 0; j < N_E; j++) {
				ret[i*N_E + j] += exp(points[k*N_E + j])*integrand;
			}		
		}
	}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N_E; j++) 
			ret[i*N_E + j] = ret[i*N_E + j]/devisor[i];	
}



// [[Rcpp::export]]
NumericVector omega_expectation_C(NumericMatrix rowSums, NumericMatrix P, NumericVector weights , NumericMatrix hazared_at_event , NumericMatrix beta_Z,NumericMatrix ans, NumericMatrix skip) {
	int N = rowSums.nrow();
	int N_E = rowSums.ncol();
	int N_P = P.nrow();
	int N_Q = hazared_at_event.ncol();
	double* rowSums_c = R2C_mat(rowSums);
	double* P_c = R2C_mat(P);
	double* W_c = R2C_vec(weights);
	double* hazared_at_event_c = R2C_mat(hazared_at_event);
	double* beta_Z_c = R2C_mat(beta_Z);
	double* ans_c = R2C_mat(ans);
	//double* skip_c = NULL;	
	double* skip_c = R2C_mat(skip);	
	double* ret = (double*)malloc(sizeof(double)*N * N_E);
	memset(ret,0,sizeof(double)*N * N_E);
	omega_expectation(rowSums_c,P_c,W_c,hazared_at_event_c,beta_Z_c,N,N_Q,N_E,N_P,ans_c,skip_c,ret);
	NumericMatrix out(N,N_E);

	for(int i = 0; i <N; ++i) {
		for(int j = 0; j  < N_E; ++j) {
	    		out(i,j) = ret[i*N_E + j];
		}
	}
	free(rowSums_c);
	free(P_c);
	free(W_c);
	free(ret);
	free(hazared_at_event_c);
	free(beta_Z_c);
	free(ans_c);
	free(skip_c);
	return out;
}



void E_phase(double* __restrict__ rowSums, double* __restrict__ points, double* __restrict__ weights, double* __restrict__ hazared_at_event, double* __restrict__  beta_Z, int N, int N_Q, int N_E, int N_P, double* ans, double* W, double* skip, double* __restrict__ ret) {
	for (int i = 0; i < N_E; i++)
		for (int j = 0; j < N_E; j++) 
			ret[i*N_E + j] = 0;
	double total_W = 0;
	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		double devisor = 0;
		double PU_varcov[N_E][N_E];
		for (int j = 0; j < N_E; j++) {
			for (int m = 0; m < N_E; m++) {
				PU_varcov[j][m] = 0;
			}
		}	
		for (int k = 0; k < N_P; k++) {
			double integrand = 0;
			for (int j = 0; j < N_E; j++) {
				double in_exp = points[k*N_E + j]*rowSums[i*N_E + j];
				for (int z = 0; z < N_Q; z++) {
					if (ans[i*N_Q + z] == 0 || skip[j*N_Q + z] == 1) 
					//if (ans[i*N_Q + z] == 0) 
						continue;
					in_exp -= hazared_at_event[j*N*N_Q + i*N_Q + z]*exp(beta_Z[j*N*N_Q + i*N_Q + z] + points[k*N_E + j]);
				}
				integrand += in_exp;
			}
			integrand = exp(integrand)*weights[k];
			devisor += integrand;
			for (int j = 0; j < N_E; j++) {
				for (int m = 0; m < N_E; m++) {
					PU_varcov[j][m] += points[k*N_E + j]*points[k*N_E + m]*integrand;
				}
			}		
		}
		for (int j = 0; j < N_E; j++) {
			for (int m = 0; m < N_E; m++) {
				#pragma omp critical
				{
					ret[j*N_E + m] += W[i]*PU_varcov[j][m]/devisor;
				}
			}
		}
		total_W += W[i];		
	}
	for (int j = 0; j < N_E; j++) {
		for (int m = 0; m < N_E; m++) {
			{
				ret[j*N_E + m] = ret[j*N_E + m]/total_W;
			}
		}
	}		
}


// [[Rcpp::export]]
NumericVector E_phase_C(NumericMatrix rowSums, NumericMatrix P, NumericVector weights , NumericMatrix hazared_at_event , NumericMatrix beta_Z, NumericMatrix ans, NumericVector W,NumericMatrix skip) {
	int N = rowSums.nrow();
	int N_E = rowSums.ncol();
	int N_P = P.nrow();
	int N_Q = hazared_at_event.ncol();
	double* rowSums_c = R2C_mat(rowSums);
	double* P_c = R2C_mat(P);
	double* weights_c = R2C_vec(weights);
	double* hazared_at_event_c = R2C_mat(hazared_at_event);
	double* beta_Z_c = R2C_mat(beta_Z);
	double* ans_c = R2C_mat(ans);
	double* skip_c = R2C_mat(skip);
	//double* skip_c = NULL;
	double* W_c = R2C_vec(W);
	
	double* ret = (double*)malloc(sizeof(double)*N_E * N_E);
	memset(ret,0,sizeof(double)*N_E * N_E);
	E_phase(rowSums_c,P_c,weights_c,hazared_at_event_c,beta_Z_c,N,N_Q,N_E,N_P,ans_c,W_c,skip_c,ret);
	NumericMatrix out(N_E,N_E);

	for(int i = 0; i <N_E; ++i) {
		for(int j = 0; j  < N_E; ++j) {
	    		out(i,j) = ret[i*N_E + j];
		}
	}
	free(rowSums_c);
	free(P_c);
	free(weights_c);
	free(ret);
	free(hazared_at_event_c);
	free(beta_Z_c);
	free(ans_c);
	free(W_c);
	free(skip_c);
	return out;
}


