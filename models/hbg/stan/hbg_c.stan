//aria: compile=1
//aria: syntax_ignore = "Warning: The parameter iZ has no priors."
data{

	// nX: number of cols in the group-level predictor matrix
	int<lower=2> nX ;

	// rX: number of rows in the group-level predictor matrix
	int<lower=nX> rX ;

	// X: group-level predictor matrix
	matrix[rX,nX] X ;

	// nI: number of individuals
	int<lower=rX> nI ;

	// iX: which group each individual is associated with
	array[nI] int<lower=1,upper=rX> iX ;

	// nY: num entries in the observation vectors
	int<lower=nX> nY ;

	// Y: observations
	vector[nY] Y ;

	// iY: which individual is associated with each observation in Y
	array[nY] int<lower=1,upper=nI> iY ;

}
parameters{

	// Y_sd: magnitude of observation-level variability
	real<lower=0> Y_sd ;

	// Z: coefficients associated with predictors
	vector[nX] Z ;

	// iZ_sd: magnitude of variability among individuals within a group
	real<lower=0> iZ_sd ;

	// iZ: by-individual means
	vector[nI] iZ ;

}
model{

	////
	// group-level structure
	////

	// standard-normal priors on all group-level coefficients
	Z ~ std_normal() ;

	// using the group predictors and coefficients, compute condition coefficients for each group
	vector[rX] gZ = rows_dot_product(
		rep_matrix(to_row_vector(Z),rX)
		, X
	) ;

	////
	// individual-level structure
	////

	// positive-standard-normal priors on magnitude of variability among individuals within a group
	iZ_sd ~ std_normal() ;

	// iZ as normal (centered parameterization)
	iZ ~ normal(gZ[iX],iZ_sd) ;

	////
	// observation-level structure
	////

	// prior peaked around .8 for magnitude of observation-level variability
	Y_sd ~ weibull(2,1) ;

	Y ~ normal( iZ[iY] , Y_sd ) ;

}
