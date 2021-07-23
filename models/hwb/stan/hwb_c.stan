//aria: compile=1
//aria: syntax_ignore = "Warning: The parameter z has no priors."
//aria: run_debug=0 #because the auto-generated debug data doesn't work
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

	// w: predictor array
	matrix[k,k] uw ;

	// n: number of subj
	int<lower=1> n ;

	// m: num entries in the observation vectors
	int<lower=k> m ;

	// num_total: number of observations per subject/condition
	array[m] int<lower=1> num_total ;

	// num_successes: sum of observations per subject/condition
	array[m] int<lower=0,upper=num_total> num_successes ;

	// row: for each observation, row index in w (see TD below)
	array[m] int<lower=1,upper=n*k> row ;

}
transformed data{

	// num_r: number of correlations implied by k
	int num_r = (k * (k - 1)) %/% 2 ;

	// compute observed intercept
	real obs_intercept = logit(mean(to_vector(num_successes)./to_vector(num_total))) ;

	// w: unique predictions repeated for each subject
	array[k] matrix[n,k] w ;
	for(i_k in 1:k){
		w[i_k] = rep_matrix(uw[i_k],n);
	}

}
parameters{

	// chol_corr: population-level correlations (on cholesky factor scale) amongst within-subject predictors
	cholesky_factor_corr[k] chol_corr ;

	//for parameters below, trailing underscore denotes that they need to be un-scaled in generated quantities

	// mu: mean (across subj) for each coefficient
	row_vector[k] mu_ ;

	// sigma: sd (across subj) for each coefficient
	vector<lower=0>[k] sigma ;

	// z: by-subject coefficients
	array[n] row_vector[k] z ;

}
model{

	////
	// Priors
	////

	// relatively flat prior on correlations
	chol_corr ~ lkj_corr_cholesky(2) ;

	// normal(0,1) priors on all sds
	sigma ~ std_normal() ;

	// normal(0,1) priors on all means
	mu_ ~ std_normal() ;

	// mid-level multivariate structure
	z ~ multi_normal_cholesky(
		append_col(
			mu_[1] + obs_intercept
			, mu_[2:k]
		)
		, diag_pre_multiply(sigma, chol_corr)
	) ;

	//convert z from array of row-vectors to matrix
	// (hopefully replaceable by `rows_dot_product(to_matrix(z),...) soon,
	// see: https://github.com/stan-dev/cmdstan/issues/1015 )
	matrix[n,k] z_mat ;
	for(i_n in 1:n){
		z_mat[i_n] = z[i_n] ;
	}

	// Loop over subj and conditions to compute unique entries in design matrix
	matrix[n,k] z_dot_w ;
	for(i_k in 1:k){
		z_dot_w[,i_k] = rows_dot_product(z_mat, w[i_k]) ;
	}

	// Likelihood
	num_successes ~ binomial(
		num_total
		, inv_logit(to_vector(z_dot_w))[row]
	) ;

}
generated quantities{

	// version of mu_ with obs_intercept added to the first entry
	vector[k] mu = to_vector(mu_);
	mu[1] += obs_intercept ;

	// r: lower-tri of correlation matrix flattened to a vector
	vector[num_r] r = flatten_lower_tri(multiply_lower_tri_self_transpose(chol_corr)) ;

}
