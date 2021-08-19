//aria: compile=1
//aria: syntax_ignore = "has no priors"

////
// Glossary
////

// n.b. the following employs a mix of snake_case and camelCase that is sure to
//   vex some, but represents the author's best attempt to balance to the competing
//   aims of clarity & brevity.

// Y: observed outcome
// nY: number of observed outcomes
// X: predictor/contrast matrix
// nX: number of predictors (columns in the contrast matrix)
// rX: number of rows in the contrast matrix X
// (i)ndividual: a unit of observation within which correlated measurements may take place
// (c)ondition: a labelled set of observations within an individual that share some feature/predictor or conjunction of features/predictors
// Xc: condition-level contrast matrix
// nXc: number of predictors in the condition-level contrast matrix
// rXc: number of rows in the condition-level contrast matrix
// yXc: for each observation in y, an index indicating the associated row in Xc corresponding to that observation's individual/condition combo
// (g)roup: a collection of individuals that share some feature/predictor
// Xg: group-level contrast matrix
// nXg: number of predictors in the group-level contrast matrix
// rXg: number of rows in the group-level contrast matrix
// Z: matrix of coefficient row-vectors to be dot-product'd with a contrast matrix
// iZc: matrix of coefficient row-vectors associated with each individual

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

	// nXg: number of group-level predictors
	int<lower=2> nXg ;

	// rXg: number of rows in the group-level predictor matrix
	int<lower=nXg> rXg ;

	// Xg: group-level predictor matrix
	matrix[rXg,nXg] Xg ;

	// nI: number of individuals
	int<lower=nXg> nI ;

	// iXg: which group each individual is associated with
	array[nI] int<lower=1,upper=rXg> iXg ;

	// nXc: number of condition-level predictors
	int<lower=2> nXc ;

	// rXc: number of rows in the condition-level predictor matrix
	int<lower=nXc> rXc ;

	// Xc: condition-level predictor matrix
	matrix[rXc,nXc] Xc ;

	// iXc: which individual is associated with each row in Xc
	array[rXc] int<lower=1,upper=nI> iXc ;

	// nY: num entries in the observation vector
	int<lower=nI> nY ;

	// Y: observations
	vector[nY] Y ;

	// yXc: which row in Xc is associated with each observation in Y
	array[nY] int<lower=1,upper=rXc> yXc ;

	// centered: whether to monolithically-center (1) or non-center (2) mid-hierarchy parameters
	int<lower=0,upper=1> centered ;

}

transformed data{

	// compute indices for fast rows_dot_product computation of group-level coefficients
	array[rXg*nXc] int i_arr_Xg ; // will be rep(1:rXg,times=nXc)
	array[rXg*nXc] int i_arr_Z ; // will be rep(1:nXc,each=rXg)
	int i = 0 ;
	for(i_nXc in 1:nXc){
		for(i_rXg in 1:rXg){
			i += 1 ;
			i_arr_Xg[i] = i_rXg ;
			i_arr_Z[i] = i_nXc ;
		}
	}

	// Xg_for_dot_Z: Xg with rows repeated through index i_arr_Xg
	matrix[rXg*nXc,nXg] Xg_for_dot_Z = Xg[i_arr_Xg] ;

}

parameters{

	// Y_sd: magnitude of observation-level variability
	real<lower=0> Y_sd ;

	// Z: coefficients associated with predictors (group-level, condition-level, & interactions)
	matrix[nXc,nXg] Z ;

	// iZc_sd: magnitude of variability among individuals within a group
	vector<lower=0>[nXc] iZc_sd ;

	// iZc_cholcorr: cholesky-factor of correlation structure associated with variability among individuals on influence of within-individual predictors
	cholesky_factor_corr[nXc] iZc_cholcorr ;

	//next two parameters represent the same quantity, just with structures suitable for either
	//  centered or non-centered parameterization. The `centered` data variable will make one or the other zero-length.
	//  Hopefully Stan will soon allow matrix arguments to multi_normal_*(), which would obviate this hack.

	// iZc: by-individual coefficients (centered parameterization)
	array[nI*centered] row_vector[nXc] iZc ;

	// iZc_: helper-variable (note underscore suffix) for non-centered parameterization of iZc
	array[(1-centered)] matrix[nXc,nI] iZc_ ;

}

model{

	////
	// group-level structure
	////

	// standard-normal priors on all group-level coefficients
	to_vector(Z) ~ std_normal() ;

	// compute condition coefficients for each group using the group predictors and coefficients
	matrix[rXg,nXc] gZc = to_matrix(
		rows_dot_product( Z[i_arr_Z] , Xg_for_dot_Z )
		, rXg, nXc
	) ;

	////
	// individual-level structure
	////

	// positive-standard-normal priors on magnitude of variability among individuals within a group
	iZc_sd ~ std_normal() ;

	// flat prior on correlations
	iZc_cholcorr ~ lkj_corr_cholesky(1) ;

	matrix[nI,nXc] iZc_mat ;
	if(centered==1){

		//convert gZc from matrix to array of row-vectors
		array[rXg] row_vector[nXc] gZc_arr ;
		for(this_rXg in 1:rXg){
			gZc_arr[this_rXg] = gZc[this_rXg] ;
		}

		// multi-normal structure for iZc
		iZc ~ multi_normal_cholesky(
			gZc_arr[iXg]
			, diag_pre_multiply(iZc_sd, iZc_cholcorr)
		) ;

		//convert iZc from array of row-vectors to matrix
		// (hopefully replaceable by `to_matrix(iZc)` soon,
		// see: https://github.com/stan-dev/cmdstan/issues/1015 )
		for(this_nI in 1:nI){
			iZc_mat[this_nI] = iZc[this_nI] ;
		}

	}else{
		// iZc_ must have a standard-normal prior for non-centered parameterization
		to_vector(iZc_[1]) ~ std_normal() ;

		// compute values for individuals given group associated group values and the individuals'
		//   deviations therefrom
		iZc_mat = (
			gZc[iXg]
			+ transpose(
				diag_pre_multiply(
					iZc_sd
					, iZc_cholcorr
				)
				* (iZc_[1])
			)
		);

	}

	// using the indivividual condition coefficients and predictors, compute
	// values for each individual
	vector[rXc] iZc_dot_Xc = rows_dot_product( iZc_mat[iXc] , Xc ) ;

	////
	// observation-level structure
	////

	// prior peaked around .8 for magnitude of observation-level variability
	Y_sd ~ weibull(2,1) ;

	// observations as normal with common error variability and varying mean
	Y ~ normal(
		iZc_dot_Xc[yXc]
		, Y_sd
	) ;

}

generated quantities{

	// iZc_r_vec: lower-tri of correlation matrix flattened to a vector
	vector[(nXc*(nXc-1))%/%2] iZc_r_vec = flatten_lower_tri(multiply_lower_tri_self_transpose(iZc_cholcorr)) ;

}
