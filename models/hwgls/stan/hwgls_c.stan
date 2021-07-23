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

	// uw: predictor array
	matrix[k,k] uw ;

	// n: number of subj
	int<lower=1> n ;

	// m: num entries in the observation vectors
	int<lower=k> m ;

	// y: observations
	vector[m] y ;

	// row: for each observation, row index in w (see TD below)
	array[m] int<lower=1,upper=n*k> row ;

}
transformed data{

	// k2: k times 2
	int k2 = k*2 ;

	// num_r: number of correlations implied by k
	int num_r = (k2 * (k2 - 1)) %/% 2 ;

	// w: unique predictions repeated for each subject
	array[k] matrix[n,k] w ;
	for(i_k in 1:k){
		w[i_k] = rep_matrix(uw[i_k],n);
	}

}
parameters{

	// chol_corr: population-level correlations (on cholesky factor scale) amongst within-subject predictors
	cholesky_factor_corr[k2] chol_corr ;

	// mu: mean (across subj) for each coefficient
	row_vector[k2] mu ;

	// sigma: sd (across subj) for each coefficient
	vector<lower=0>[k2] sigma ;

	// z: by-subject coefficients
	array[n] row_vector[k2] z ;

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
	mu ~ std_normal() ;

	// mid-level multivariate structure
	z ~ multi_normal_cholesky(mu,diag_pre_multiply(sigma, chol_corr)) ;

	//convert z from array of row-vectors to matrix
	// (hopefully replaceable by `rows_dot_product(to_matrix(z),...) soon,
	// see: https://github.com/stan-dev/cmdstan/issues/1015 )
	matrix[n,k2] z_mat ;
	for(i_n in 1:n){
		z_mat[i_n] = z[i_n] ;
	}

	// Loop over subj and conditions to compute unique entries in design matrix
	matrix[n,k] loc_z_dot_w ;
	matrix[n,k] scale_z_dot_w ;
	for(i_k in 1:k){
		loc_z_dot_w[,i_k] = rows_dot_product(z_mat[,1:k], w[i_k]) ;
		scale_z_dot_w[,i_k] = rows_dot_product(z_mat[,(k+1):k2], w[i_k]) ;
	}

	// Likelihood
	y ~ normal(
		to_vector(loc_z_dot_w)[row]
		, sqrt(exp(to_vector(to_vector(scale_z_dot_w)[row])))
	) ;

}
generated quantities{

	// r: lower-tri of correlation matrix flattened to a vector
	vector[num_r] r = flatten_lower_tri(multiply_lower_tri_self_transpose(chol_corr)) ;

}
