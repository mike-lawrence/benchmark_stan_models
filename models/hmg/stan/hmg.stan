//aria: compile=1
//aria: syntax_ignore = "Warning: The parameter iZq has no priors."
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
	int<lower=nI> nY ;

	// Y: observations
	vector[nY] Y ;

	// yXq: which row in Xq is associated with each observation in Y
	array[nY] int<lower=1,upper=rXq> yXq ;

	// centered: whether to monolithically-center (1) or non-center (2) mid-hierarchy parameters
	int<lower=0,upper=1> centered ;

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

	// iZq: by-individual coefficients (centered parameterization)
	array[nI*centered] row_vector[nXq] iZq ;

	// iZq_: helper-variable (note underscore suffix) for non-centered parameterization of iZq
	array[(1-centered)] matrix[nXq,nI] iZq_ ;

}
model{

	////
	// group-level structure
	////

	// standard-normal priors on all group-level coefficients
	to_vector(Z) ~ std_normal() ;

	// using the group predictors and coefficients, compute condition coefficients for each group
	// NOT SURE THIS IS THE MOST EFFICIENT WAY TO DO THIS (ESP. GIVEN NECESSARY LATER CONVERSION)
	matrix[rXg,nXq] gZq ;
	for(this_nXq in 1:nXq){
		gZq[,this_nXq]= rows_dot_product(
			rep_matrix(to_row_vector(Z[,this_nXq]),rXg)
			, Xg
		) ;
	}

	//convert gZq from matrix to array of row-vectors
	array[rXg] row_vector[nXq] gZq_arr ;
	for(this_rXg in 1:rXg){
		gZq_arr[this_rXg] = gZq[this_rXg] ;
	}


	////
	// individual-level structure
	////

	// positive-standard-normal priors on magnitude of variability among individuals within a group
	iZq_sd ~ std_normal() ;

	// flat prior on correlations
	iZq_cholcorr ~ lkj_corr_cholesky(1) ;

	matrix[nI,nXq] iZq_mat ;
	if(centered==1){
		// multi-normal structure for iZq
		iZq ~ multi_normal_cholesky(
			gZq_arr[iXg]
			, diag_pre_multiply(iZq_sd, iZq_cholcorr)
		) ;

		//convert iZq from array of row-vectors to matrix
		// (hopefully replaceable by `to_matrix(iZq)` soon,
		// see: https://github.com/stan-dev/cmdstan/issues/1015 )
		for(this_nI in 1:nI){
			iZq_mat[this_nI] = iZq[this_nI] ;
		}

	}else{
		// iZq_ must have a standard-normal prior for non-centered parameterization
		to_vector(iZq_[1]) ~ std_normal() ;

		// compute values for individuals given group associated group values and the individuals'
		//   deviations therefrom, scaled by ind_sd
		iZq_mat = (
			gZq[iXg]
			+ transpose(
				diag_pre_multiply(
					iZq_sd
					, iZq_cholcorr
				)
				* iZq_[1]
			)
		);

	}

	// using the indivividual condition coefficients and predictors, compute
	// values for each individual
	vector[rXq] iZq_dot_Xq = rows_dot_product( iZq_mat[iXq] , Xq ) ;

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
