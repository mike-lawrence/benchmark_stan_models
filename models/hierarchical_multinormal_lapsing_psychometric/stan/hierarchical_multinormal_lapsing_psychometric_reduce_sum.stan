//aria: compile=1
//aria: syntax_ignore = c("has no priors","The variable partial_lupmf may not have been assigned a value before its use.")

////
// Glossary
////

// n.b. the following employs a mix of snake_case and camelCase that is sure to
//  vex some, but represents the author's best attempt to balance to the competing
//  aims of clarity & brevity.

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

	// p_quick: probability of success via the "quick" function (threshold=where performance hits halfway between chance and lapse)
	vector p_quick(
		vector x
		, vector p_chance
		, vector one_minus_p_chance
		, vector threshold
		, vector slope
		, vector p_lapse
	){
		return(
			p_chance
			+ (
				(one_minus_p_chance-p_lapse)
				.* (
					1-pow(2,-pow(x./threshold,slope))
				)
			)
		) ;
	}

	real partial_lpmf(
		int[] Y_slice
		, int start
		, int end
		, vector intensity
		, vector p_chance
		, vector one_minus_p_chance
		, vector p_lapse
		, vector threshold
		, vector slope
	){
		return(
			bernoulli_lupmf(
				Y_slice
				|
			 	p_quick(
					intensity[start:end] // data variable
					, p_chance[start:end] // data variable
					, one_minus_p_chance[start:end] // data variable
					, threshold[start:end] // threshold (parameter variable)
					, slope[start:end] //slope (parameter variable)
					, p_lapse[start:end] // p_lapse (parameter variable)
				)
			)
		) ;
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

	// Y: 1 (success) or 2 (failure) observations
	array[nY] int<lower=0,upper=1> Y ;

	// yXc: which row in Xc is associated with each observation in Y
	array[nY] int<lower=1,upper=rXc> yXc ;

	// p_chance: chance probability for each observation
	vector<lower=0,upper=1>[nY] p_chance ;

	// intensity: intensity covariate for each observation
	vector[nY] intensity ;

	// centered: whether to monolithically-center (1) or non-center (2) mid-hierarchy parameters
	int<lower=0,upper=1> centered ;

}

transformed data{

	// nXc2: nXc times two
	int nXc2 = nXc*2 ;

	// nXc21: nXc times two plus one
	int nXc21 = nXc*2+1 ;

	// one_minus_p_chance: 1-p_chance
	vector[nY] one_minus_p_chance = 1-p_chance ;

	// compute indices for fast rows_dot_product computation of group-level coefficients
	array[rXg*nXc21] int i_arr_Xg ; // will be rep(1:rXg,times=nXc21)
	array[rXg*nXc21] int i_arr_Z ; // will be rep(1:nXc21,each=rXg)
	int i = 0 ;
	for(i_nXc21 in 1:nXc21){
		for(i_rXg in 1:rXg){
			i += 1 ;
			i_arr_Xg[i] = i_rXg ;
			i_arr_Z[i] = i_nXc21 ;
		}
	}

	// Xg_for_dot_Z: Xg with rows repeated through index i_arr_Xg
	matrix[rXg*nXc21,nXg] Xg_for_dot_Z = Xg[i_arr_Xg] ;

}

parameters{

	// Z: coefficients associated with predictors (group-level, condition-level, & interactions)
	matrix[nXc21,nXg] Z ;

	real<lower=0,upper=1> p_lapse_intercept ;

	// iZc_sd: magnitude of variability among individuals within a group
	vector<lower=0>[nXc21] iZc_sd ;

	// iZc_cholcorr: cholesky-factor of correlation structure associated with variability among individuals on influence of within-individual predictors
	cholesky_factor_corr[nXc21] iZc_cholcorr ;

	//next two parameters represent the same quantity, just with structures suitable for either
	//  centered or non-centered parameterization. The `centered` data variable will make one or the other zero-length.
	//  Hopefully Stan will soon allow matrix arguments to multi_normal_*(), which would obviate this hack.

	// iZc: by-individual coefficients (centered parameterization)
	array[nI*centered] row_vector[nXc21] iZc ;

	// iZc_: helper-variable (note underscore suffix) for non-centered parameterization of iZc
	array[(1-centered)] matrix[nXc21,nI] iZc_ ;

}

model{

	////
	// group-level structure
	////
	// standard-normal priors on all group-level coefficients
	to_vector(Z) ~ std_normal() ;

	// Z_: copy of Z we can modify
	matrix[nXc21,nXg] Z_ = Z ; // would it be more efficient to not use to_vector(Z)?
	//replacing the entry in Z corresponding to the lapse intercept
	// with logit of p_lapse_intercept (which has an implicit uniform prior)
	Z_[nXc21,1] = logit(p_lapse_intercept) ;

	// using the group predictors and coefficients, compute condition coefficients for each group
	matrix[rXg,nXc21] gZc = to_matrix(
		rows_dot_product( Z_[i_arr_Z] , Xg_for_dot_Z )
		, rXg
		, nXc21
	) ;

	////
	// individual-level structure
	////

	// positive-standard-normal priors on magnitude of variability among individuals within a group
	iZc_sd ~ std_normal() ;

	// flat prior on correlations
	iZc_cholcorr ~ lkj_corr_cholesky(1) ;

	matrix[nI,nXc21] iZc_mat ;
	if(centered==1){

		//convert gZc from matrix to array of row-vectors
		array[rXg] row_vector[nXc21] gZc_arr ;
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


	////
	// observation-level structure
	////

	target += reduce_sum(
		partial_lupmf // n.b. U in name!!
		, Y
		, 1 //grainsize (1=auto)
		, intensity
		, p_chance
		, one_minus_p_chance
		, sqrt(exp(rows_dot_product( iZc_mat[iXc,1:nXc] , Xc )))[yXc] //threshold
		, sqrt(exp(rows_dot_product( iZc_mat[iXc,(nXc+1):nXc2] , Xc )))[yXc] //slope
		, inv_logit(iZc_mat[,nXc21])[iXc][yXc] //p_lapse
	) ;

}
