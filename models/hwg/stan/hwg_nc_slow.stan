//aria: compile=1
functions{
	// flatten_lower_tri: function that returns the lower-tri of a matrix, flattened to a vector
	vector flatten_lower_tri(matrix mat) {
		int n_cols = cols(mat);
		int n_uniq = (n_cols * (n_cols - 1)) %/% 2;
		vector[n_uniq] out ;
		int i = 1;
		for(c in 1:(n_cols-1)){
			for(r in (c+1):n_cols){
				out[i] = mat[r,c];
				i += 1;
			}
		}
		return(out) ;
	}
}
data{

	// k: num predictors
	int<lower=1> k ;

	// m: num entries in the observation vectors
	int<lower=k> m ;

	// y: observations
	vector[m] y ;


	// w: predictor array
	matrix[m,k] w ;

	// n: number of subj
	int<lower=1> n ;

	array[m] int<lower=1,upper=n> id ;


}
transformed data{

	// num_r: number of correlations implied by k
	int num_r = (k * (k - 1)) %/% 2 ;


}
parameters{

	// noise: measurement noise scale
	real<lower=0> noise ;

	// chol_corr: population-level correlations (on cholesky factor scale) amongst within-subject predictors
	cholesky_factor_corr[k] chol_corr ;

	// mu: mean (across subj) for each coefficient
	row_vector[k] mu ;

	// sigma: sd (across subj) for each coefficient
	vector<lower=0>[k] sigma ;

	// z_: a helper variable for implementing non-centered parameterization
	matrix[k,n] z_ ;

}
model{

	////
	// Priors
	////

	// prior on noise peaked around .8
	noise ~ weibull(2,1) ;

	// z_ must have normal(0,1) prior for non-centered parameterization
	to_vector(z_) ~ std_normal() ;

	// relatively flat prior on correlations
	chol_corr ~ lkj_corr_cholesky(2) ;

	// normal(0,1) priors on all sds
	sigma ~ std_normal() ;

	// normal(0,1) priors on all means
	mu ~ std_normal() ;

	// compute coefficients for each subject/condition
	matrix[n,k] z = (
		rep_matrix(mu,n)
		+ transpose(
			diag_pre_multiply(sigma,chol_corr)
			* z_
		)
	) ;

	// Loop over subj and conditions to compute appropriate dot product
	vector[m] z_dot_w ;
	for(i_m in 1:m){
		z_dot_w[i_m] = dot_product(z[id[i_m],], w[i_m]) ;
	}

	// Likelihood
	y ~ normal(z_dot_w, noise) ;

}
generated quantities{

	// r: lower-tri of correlation matrix flattened to a vector
	vector[num_r] r = flatten_lower_tri(multiply_lower_tri_self_transpose(chol_corr)) ;

}
