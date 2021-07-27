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

	// nXg: number of cols in the group-level predictor matrix
	int<lower=2> nXg ;

	// rXg: number of rows in the group-level predictor matrix
	int<lower=nXg> rXg ;

	// Xg: group-level predictor matrix
	matrix[rXg,nXg] Xg ;

	// nI: number of individuals
	int<lower=nXg> nI ;

	// iXg: which group each individual is associated with
	array[nI] int<lower=1,upper=rXg> iXg ;

	// nXq: number of cols in the condition-level predictor matrix
	int<lower=2> nXq ;

	// rXq: number of rows in the condition-level predictor matrix
	int<lower=nXq> rXq ;

	// Xq: condition-level predictor matrix
	matrix[rXq,nXq] Xq ;

	// iXq: which individual is associated with each row in Xq
	array[rXq] int<lower=1,upper=nI> iXq ;

	// nY: num entries in the observation vectors
	int<lower=nXg*nXq*nI> nY ;

	// Y: observations
	vector[nY] Y ;

	// yXq: which row in Xq is associated with each observation in Y
	array[nY] int<lower=1,upper=rXq> yXq ;

}
parameters{

	// Y_sd: magnitude of observation-level variability
	real<lower=0> Y_sd ;

	// Z: coefficients associated with predictors (group-level, condition-level, & interactions)
	matrix[nXg,nXq] Z ;

	// iZq_sd: magnitude of variability among individuals within a group
	vector<lower=0>[nXq] iZq_sd ;

	// iZq_cholcorr: cholesky-factor of correlation structure associated with variability among individuals on influence of within-individual predictors
	cholesky_factor_corr[nXq] iZq_cholcorr ;

	// iZq_: helper-variable (note underscore suffix) for non-centered parameterization of iZq
	matrix[nXq,nI] iZq_ ;

}
model{

	////
	// group-level structure
	////

	// standard-normal priors on all group-level coefficients
	to_vector(Z) ~ std_normal() ;

	// using the group predictors and coefficients, compute condition coefficients for each group
	matrix[rXg,nXq] gZq ; //n.b. transposed so that the qZs are now row_vectors
	for(this_nXq in 1:nXq){
		gZq[,this_nXq]= rows_dot_product(
			rep_matrix(to_row_vector(Z[,this_nXq]),rXg)
			, Xg
		) ;
	}

	////
	// individual-level structure
	////

	// positive-standard-normal priors on magnitude of variability among individuals within a group
	iZq_sd ~ std_normal() ;

	// flat prior on correlations
	iZq_cholcorr ~ lkj_corr_cholesky(1) ;

	// iZq_ must have a standard-normal prior for non-centered parameterization
	to_vector(iZq_) ~ std_normal() ;

	// compute values for individuals given group associated group values and the individuals'
	//   deviations therefrom, scaled by ind_sd
	matrix[nI,nXq] iZq = (
		gZq[iXg]
		+ transpose(
			diag_pre_multiply(
				iZq_sd
				, iZq_cholcorr
			)
			* iZq_
		)
	);

	// using the indivividual condition coefficients and predictors, compute
	// values for each individual
	vector[nXq] iZq_dot_Xq = rows_dot_product( iZq[iXq] , Xq ) ;

	////
	// observation-level structure
	////

	// prior peaked around .8 for magnitude of observation-level variability
	Y_sd ~ weibull(2,1) ;

	Y ~ normal(
		iZq_dot_Xq[yXq]
		, Y_sd
	) ;

}
generated quantities{

	// iZq_r_vec: lower-tri of correlation matrix flattened to a vector
	vector[(nXq*(nXq-1))%/%2] iZq_r_vec = flatten_lower_tri(multiply_lower_tri_self_transpose(iZq_cholcorr)) ;

}
