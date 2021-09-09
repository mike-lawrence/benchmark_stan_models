//aria: compile=1
//aria: syntax_ignore = c("has no priors","The variable partial_lupmf may not have been assigned a value before its use.","was declared but was not used in the density calculation.")

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
	matrix[nXc21,nXc21] iZc_cholcorr ;

	//next two parameters represent the same quantity, just with structures suitable for either
	//  centered or non-centered parameterization. The `centered` data variable will make one or the other zero-length.
	//  Hopefully Stan will soon allow matrix arguments to multi_normal_*(), which would obviate this hack.

	// iZc: by-individual coefficients (centered parameterization)
	array[nI] row_vector[nXc21] iZc ;


}

generated quantities{
	matrix[nI,nXc] threshold_coefs  ;
	matrix[nI,nXc] slope_coefs ;
	vector[rXc] threshold  ;
	vector[rXc] slope ;
	vector[rXc] p_lapse ;
	// Z_: copy of Z we can modify
	matrix[nXc21,nXg] Z_ = Z ; // would it be more efficient to not use to_vector(Z)?
	//replacing the entry in Z corresponding to the lapse intercept
	// with logit of p_lapse_intercept (which has an implicit uniform prior)
	Z_[nXc21,1] = logit(p_lapse_intercept) ;

	{ // local environment to avoid saving gZc & iZc_mat

		// using the group predictors and coefficients, compute condition coefficients for each group
		matrix[rXg,nXc21] gZc = to_matrix(
			rows_dot_product( Z_[i_arr_Z] , Xg_for_dot_Z )
			, rXg
			, nXc21
		) ;

		////
		// individual-level structure
		////
		matrix[nI,nXc21] iZc_mat ;
		for(this_nI in 1:nI){
			iZc_mat[this_nI] = iZc[this_nI] ;
		}

		threshold_coefs = iZc_mat[,1:nXc] ;
		slope_coefs = iZc_mat[,(nXc+1):nXc2] ;
		threshold = sqrt(exp(rows_dot_product( iZc_mat[iXc,1:nXc] , Xc ))) ;
		slope = sqrt(exp(rows_dot_product( iZc_mat[iXc,(nXc+1):nXc2] , Xc ))) ;
		p_lapse = inv_logit(iZc_mat[,nXc21])[iXc] ;
	}

}
