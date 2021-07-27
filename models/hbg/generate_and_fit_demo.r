# Preamble (options, installs, imports & custom functions) ----

options(warn=1) #really should be default in R
`%!in%` = Negate(`%in%`) #should be in base R!

# specify the packages used:
required_packages = c(
	'github.com/rmcelreath/rethinking' # for rlkjcorr & rmvrnom2
	, 'github.com/stan-dev/cmdstanr' #for Stan stuff
	, 'github.com/mike-lawrence/aria/aria' # for aria
	, 'tidyverse' #for all that is good and holy
)

# load the helper functions:
for(file in fs::dir_ls('r')){
	cat('Loading function: ',fs::path_ext_remove(fs::path_file(file)),'()\n',sep='')
	source(file)
}

#install any required packages not already present
install_if_missing(required_packages)

# load tidyverse & aria
library(tidyverse)
library(aria)

# Simulate data ----
set.seed(1) #change this to make different data

#setting the data simulation parameters

#parameters you can play with
nI_per_group = 3 #number of individuals, must be an integer >1
nG_vars = 3 #number of 2-level variables manipulated as crossed and across individuals, must be an integer >0
num_Y_per_iq = 3 #number of observations per individual/condition combo, must be an integer >1

#the rest of these you shouldn't touch
nX = 2^(nG_vars)
Z = rnorm(nX)
iZ_sd = rweibull(1,2,1)
Y_sd = rweibull(1,2,1)

#compute X
X = sim_contrasts(nG_vars,'G')

# compute gZq
gZq = rep(NA,nX)
for(this_nX in 1:nX){
	gZq[this_nX] = X[this_nX,] %*% Z
}

#compute iZ
(
	1:nX
	%>% map_dfr(
		.f = function(this_nX){(
			tibble(
				value = rnorm(
					n = nI_per_group
					, mean = gZq[this_nX]
					, sd = iZ_sd
				)
			)
			%>% mutate(
				g = this_nX
				, individual = paste(g,1:n(),sep='_') #important that individuals have unique identifiers across groups
			)
		)}
	)
) ->
	iZ


#compute Y
(
	iZ
	%>% group_by(individual,g)
	%>% summarise(
		value = rnorm(num_Y_per_iq,value,Y_sd)
		, .groups = 'keep'
	)
) ->
	Y

#join with group & condition columns
(
	Y
	%>% left_join((
		attr(X,'data')
		%>% mutate(g=1:n())
	))
	%>% ungroup()
	%>% select(everything(),-g,-value,value)
) ->
	dat

#dat is now what one would typically have as collected data
dat

# Compute inputs to Stan model ----

#if there's no missing data (as in this synthetic example), we could proceed with modelling,
#  but to demonstrate how to handle missing data, we'll do a bit of extra work:

#Uncomment next section to induce data-missingness (causes some Ss to be missing whole conditions)
(
	dat
	#first collapse to individual-by-condition stats
	%>% group_by(
		across(-value)
	)
	%>% summarise(.groups='drop')
	%>% slice_sample(prop=.9)
	%>% semi_join(
		x = dat
		, y = .
	)
) ->
	dat

#compute X from distinct combinations of groups
(
	dat
	%>% select(starts_with('G'))
	%>% distinct()
	#add the contrast matrix columns
	%>% mutate(
		contrasts = get_contrast_matrix_rows_as_list(
			data = .
			# the following complicated specification of the formula is a by-product of making this example
			# work for any number of variables =  normally you would do something like this
			# (for 2 variables for example):
			# formula = ~ G1*G2
			, formula = as.formula(paste0('~',paste0('G',1:nG_vars,collapse='*')))
			# half-sum contrasts are nice for 2-level variables bc they yield parameters whose value
			# is the difference between conditions
			, contrast_kind = halfsum_contrasts
		)
	)
) ->
	X_with_vars

#join X with dat to label individuals with corresponding row from X
(
	X_with_vars
	# add row identifier
	%>% mutate(X_row=1:n())
	%>% right_join((
		dat
		%>% select(individual,starts_with('G'))
		%>% distinct()
	))
	%>% arrange(individual)
	%>% pull(X_row)
) ->
	iX


# package for stan & sample ----

data_for_stan = lst( #lst permits later entries to refer to earlier entries

	####
	# Entries we need to specify ourselves
	####


	# X: group-level predictor matrix
	X = (
		X_with_vars
		%>% select(contrasts)
		%>% unnest(contrasts)
		%>% as.matrix()
	)

	# iX: which group each individual is associated with
	, iX = iX

	# Y: observations
	, Y = dat$value

	# iY: which individual is associated with each observation in Y
	, iY = as.numeric(factor(dat$individual))

	####
	# Entries computable from the above
	####

	# nX: number of cols in the group-level predictor matrix
	, nX = ncol(X)

	# rX: number of rows in the group-level predictor matrix
	, rX = nrow(X)

	# nI: number of individuals
	, nI = max(iY)

	# nY: num entries in the observation vectors
	, nY = length(Y)

)

# double-check:
glimpse(data_for_stan)

#set the model and posterior paths
mod_path = 'stan/hbg_c.stan'
post_path = 'nc/hbg_c.nc'

# ensure model is compiled
aria:::check_syntax_and_maybe_compile(mod_path)

# compose
aria::compose(
	data = data_for_stan
	, code_path = mod_path
	, out_path = post_path
	, overwrite = T
)

# check posterior diagnostics ----
post = aria::coda(post_path)

# Check treedepth, divergences, & rebfmi
(
	post$draws(group='sample_stats')
	%>% posterior::as_draws_df()
	%>% group_by(.chain)
	%>% summarise(
		max_treedepth = max(treedepth)
		, num_divergent = sum(divergent)
		, rebfmi = var(energy)/(sum(diff(energy)^2)/n()) #n.b. reciprocal of typical EBFMI, so bigger=bad, like rhat
	)
)

# gather summary for core parameters (inc. r̂ & ess)
(
	post$draws(group='parameters')
	%>% posterior::summarise_draws(.cores=parallel::detectCores())
) ->
	par_summary

# show the ranges of r̂/ess's
(
	par_summary
	%>% select(rhat,contains('ess'))
	%>% summary()
)

#View those with suspect r̂
(
	par_summary
	%>% filter(rhat>1.01)
	%>% View()
)

# Viz recovery of non-correlation parameters ----
(
	post$draws(variables=c('Y_sd','Z','iZ_sd'))
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = 'variable'
	)
	%>% left_join(
		bind_rows(
			tibble(
				true = Y_sd
				, variable = 'Y_sd'
			)
			, tibble(
				true = as.vector(Z)
				, variable = paste0('Z[',1:length(true),']')
			)
			, tibble(
				true = iZ_sd
				, variable = 'iZ_sd'
			)
		)
	)
	%>% group_by(variable)
	%>% arrange(variable,.chain,.iteration)
	%>% summarise(
		rhat = 1.01<posterior::rhat(matrix(value,ncol=length(unique(.chain))))
		, ess_bulk = 100>posterior::ess_bulk(matrix(value,ncol=length(unique(.chain))))
		, ess_tail = 100>posterior::ess_tail(matrix(value,ncol=length(unique(.chain))))
		, as_tibble(t(posterior::quantile2(value,c(.05,.25,.5,.75,.95))))
		, true = true[1]
	)
	%>% mutate(variable = factor_1d(variable))
	%>% ggplot()
	+ geom_hline(yintercept = 0)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q5
			, ymax = q95
			, colour = ess_tail
		)
		, alpha = .5
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q25
			, ymax = q75
			, colour = ess_bulk
		)
		, size = 3
		, alpha = .5
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = q50
			, fill = rhat
		)
		, shape = 21
		, size = 3
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = true
		)
		, size = 4
		, shape = 4
		, colour = 'blue'
	)
	+ coord_flip()
	+ scale_color_manual(
		values = lst(`TRUE`='red',`FALSE`='black')
		, labels = lst(`TRUE`='<100',`FALSE`='>=100')
	)
	+ scale_fill_manual(
		values = lst(`TRUE`='red',`FALSE`='white')
		, labels = lst(`TRUE`='>1.01',`FALSE`='<=1.01')
	)
	+ labs(
		y = 'True & Posterior Value'
		, x = 'Variable'
		, colour = 'ESS'
		, fill = 'Rhat'
	)
)


