# Glossary ----

# n.b. the following employs a mix of snake_case and camelCase that is sure to
#  vex many, but represents the author's best solution to the competing aims of
#  clarity & brevity.

# Y: observed outcome
# nY: number of observed outcomes
# X: predictor/contrast matrix
# nX: number of predictors (columns in the contrast matrix)
# rX: number of rows in the contrast matrix X
# (i)ndividual: a unit of observation within which correlated measurements may take place
# (c)ondition: a labelled set of observations within an individual that share some feature/predictor or conjunction of features/predictors
# Xc: condition-level contrast matrix
# nXc: number of predictors in the condition-level contrast matrix
# rXc: number of rows in the condition-level contrast matrix
# yXc: for each observation in y, an index indicating the associated row in Xc corresponding to that observation's individual/condition combo
# (g)roup: a collection of individuals that share some feature/predictor
# Xg: group-level contrast matrix
# nXg: number of predictors in the group-level contrast matrix
# rXg: number of rows in the group-level contrast matrix
# Z: matrix of coefficient row-vectors to be dot-product'd with a contrast matrix
# iZc: matrix of coefficient row-vectors associated with each individual

# Preamble (options, installs, imports & custom functions) ----

options(warn=1) # really should be default in R
`%!in%` = Negate(`%in%`) # should be in base R!
rstan::rstan_options(auto_write = FALSE)

# specify the packages used:
required_packages = c(
	'github.com/rmcelreath/rethinking' # for rlkjcorr & rmvrnom2
	, 'github.com/stan-dev/cmdstanr' # for Stan stuff
	, 'github.com/mike-lawrence/aria/aria' # for aria
	, 'tidyverse' # for all that is good and holy
)

# load the helper functions:
for(file in fs::dir_ls('r')){
	cat('Loading function: ',fs::path_ext_remove(fs::path_file(file)),'()\n',sep='')
	source(file)
}

# install any required packages not already present
install_if_missing(required_packages)

# load tidyverse & aria
library(tidyverse)
library(aria)

# Simulate data ----
set.seed(1) # change this to make different data

# setting the data simulation parameters

# parameters you can play with
nG_vars = 1 # number of 2-level variables manipulated as crossed and across individuals, must be an integer >0
nC_vars = 1 # number of 2-level variables manipulated as crossed and within each individual, must be an integer >0
nI_per_group = 100 # number of individuals, must be an integer >1
num_Y_per_ic = 100 # number of observations per individual/condition combo, must be an integer >1
# the latter two combine to determine whether centered or non-centered will sample better

# the rest of these you shouldn't touch
nXc = 2^(nC_vars)
nXc2 = nXc*2
nXc21 = nXc2+1
nXg = 2^(nG_vars)
Z = matrix(rnorm(nXg*nXc21),nrow=nXg,ncol=nXc21)
Z[1,nXc21] = qlogis(.1)
iZc_sd = rweibull(nXc21,2,1)
iZc_r_mat = rethinking::rlkjcorr(1,nXc21,eta=1)
iZc_r_vec = iZc_r_mat[lower.tri(iZc_r_mat)]


# compute Xg
Xg = sim_contrasts(nG_vars,'G')

# compute gZc
gZc = matrix(NA,nrow=nXg,ncol=nXc21)
for(this_nXc21 in 1:nXc21){
	for(this_nXg in 1:nXg){
		gZc[this_nXg,this_nXc21] = Xg[this_nXg,] %*% Z[,this_nXc21]
	}
}

# compute iZc
(
	1:nXg
	%>% map_dfr(
		.f = function(this_nXg){(
			rethinking::rmvnorm2(
				n = nI_per_group
				, Mu = gZc[this_nXg,]
				, sigma = iZc_sd
				, Rho = iZc_r_mat
			)
			%>% as_tibble(
				.name_repair = function(x) paste0('iZc[.,',1:length(x),']')
			)
			%>% rename()
			%>% mutate(
				g = this_nXg
				, individual = paste(g,1:n(),sep='_') # important that individuals have unique identifiers across groups
			)
		)}
	)
) ->
	iZc

# compute Xc
uXc = sim_contrasts(nC_vars,'C')

# compute iZc_dot_Xc
(
	iZc
	%>% group_by(g,individual)
	%>% summarise(
		{function(x){
			x = t(as.matrix(select(x,starts_with('iZc')),))
			threshold_this_dot = rep(NA,nXc)
			slope_this_dot = rep(NA,nXc)
			p_lapse_this_dot = rep(NA,nXc)
			for(this_nXc in 1:nXc){
				threshold_this_dot[this_nXc] = sqrt(exp(uXc[this_nXc,] %*% x[1:nXc,]))
				slope_this_dot[this_nXc] = sqrt(exp(uXc[this_nXc,] %*% x[(nXc+1):nXc2,]))
				p_lapse_this_dot[this_nXc] = plogis(x[nXc21,])
			}
			tibble(
				c = 1:nXc
				, threshold = threshold_this_dot
				, slope = slope_this_dot
				, p_lapse = p_lapse_this_dot
			)
		}}(cur_data())
		, .groups = 'keep'
	)
) ->
	iZc_dot_Xc

# compute Y
(
	iZc_dot_Xc
	%>% group_by(c,.add=T)
	%>% summarise(
		intensity = seq(0,max(iZc_dot_Xc$threshold)*2,length.out=num_Y_per_ic)
		# intensity = runif(num_Y_per_ic,0,1)
		, value = rbinom(
			n = num_Y_per_ic
			, size = 1
			, prob = pquick(
				x = intensity
				, chance = .5
				, lapse = p_lapse
				, threshold = threshold
				, slope = slope
			)
		)
		, .groups = 'keep'
	)
) ->
	Y

# join with group & condition columns
(
	Y
	%>% left_join((
		attr(Xg,'data')
		%>% mutate(g=1:n())
	))
	%>% left_join((
		attr(uXc,'data')
		%>% mutate(c=1:n())
	))
	%>% ungroup()
	%>% select(everything(),-g,-c,-value,value)
) ->
	dat

# dat is now what one would typically have as collected data
dat

(
	dat
	%>% ggplot()
	+ geom_point(
		aes(
			x = intensity
			, y = individual
			, colour = factor(value)
		)
		, alpha = .1
	)
	+ facet_grid(.~C1)
)

# Compute inputs to Stan model ----

# Some data prep
(
	dat
	# ungroup (necessary)
	%>% ungroup()
	# ensure individual is a sequential numeric
	%>% mutate(
		individual = as.numeric(factor(individual))
	)
	# arrange rows by individual
	%>% arrange(individual)
) ->
	dat

# compute Xg from distinct combinations of groups
(
	dat
	# select down to any G vars
	%>% select(starts_with('G'))
	# collapse to distinct set of rows (combinations of grouping variables)
	%>% distinct()
	# arrange (not really necessary, but why not)
	%>% arrange()
	# add the contrast matrix columns
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
	Xg_with_vars

# show contrasts
(
	Xg_with_vars
	%>% unnest(contrasts)
)

# join Xg with dat to label individuals with corresponding row from Xg
(
	Xg_with_vars
	# add row identifier
	%>% mutate(Xg_row=1:n())
	# join with dat, collapsed to 1-row per individual with their group info
	%>% right_join((# right-join to apply the row order from dat
		dat
		%>% select(individual,starts_with('G'))
		%>% distinct()
	))
	%>% arrange(individual)
	%>% pull(Xg_row)
) ->
	iXg

# compute complete_Xc from distinct combinations of individuals & conditions
(
	dat
	# select individual & any condition-defining columns
	%>% select(individual,starts_with('C'))
	# collapse down to distinct rows (1 per individual/conditions combo)
	%>% distinct()
	# expand (in case there's any missing data; doesn't hurt if not)
	# %>% exec(expand_grid,!!!.)
	%>% as.list()
	%>% map(unique)
	%>% cross_df()
	# arrange (not really necessary, but why not)
	%>% arrange()
	# add the contrast matrix columns
	%>% mutate(
		contrasts = get_contrast_matrix_rows_as_list(
			data = .
			# the following complicated specification of the formula is a by-product of making this example
			# work for any number of variables =  normally you would do something like this
			# (for 2 variables for example):
			# formula = ~ C1*C2
			, formula = as.formula(paste0('~',paste0('C',1:nC_vars,collapse='*')))
			# half-sum contrasts are nice for 2-level variables bc they yield parameters whose value
			# is the difference between conditions
			, contrast_kind = halfsum_contrasts
		)
	)
) ->
	complete_Xc_with_vars

# show the unique contrasts
(
	complete_Xc_with_vars
	%>% select(-individual)
	%>% distinct()
	%>% unnest(contrasts)
)

# subset down to just those individual-condition combos actually present in the data
#  it's ok if there's no missing data and nrow(complete_Xc_with_vars)==nrow(Xc_with_vars)
(
	complete_Xc_with_vars
	%>% semi_join(dat)
) ->
	Xc_with_vars


# join Xc with dat to label observations with corresponding row from Xc
(
	Xc_with_vars
	# add row identifier
	%>% mutate(Xc_row=1:n())
	# right-join with dat to preserve dat's row order
	%>% right_join(
		mutate(dat,dat_row=1:n())
	)
	%>% arrange(dat_row)
	# pull the Xc row identifier
	%>% pull(Xc_row)
) ->
	yXc


# package for stan & sample ----

data_for_stan = lst( # lst permits later entries to refer to earlier entries

	# # # #
	# Entries we need to specify ourselves
	# # # #


	# Xg: group-level predictor matrix
	Xg = (
		Xg_with_vars
		%>% select(contrasts)
		%>% unnest(contrasts)
		%>% as.matrix()
	)

	# iXg: which group each individual is associated with
	, iXg = iXg

	# Xc: condition-level predictor matrix
	, Xc = (
		Xc_with_vars
		%>% select(contrasts)
		%>% unnest(contrasts)
		%>% as.matrix()
	)

	# iXc: which individual is associated with each row in Xc
	, iXc = as.numeric(factor(Xc_with_vars$individual))

	# Y: observations
	, Y = dat$value

	# yXc: which row in Xc is associated with each observation in Y
	, yXc = yXc

	# intensity: intensity covariate for each observation
	, intensity = dat$intensity

	# p_chance: probability of success @ chance for each observation
	, p_chance = rep(.5,times=length(Y))

	# # # #
	# Entries computable from the above
	# # # #

	# nXg: number of cols in the group-level predictor matrix
	, nXg = ncol(Xg)

	# rXg: number of rows in the group-level predictor matrix
	, rXg = nrow(Xg)

	# nI: number of individuals
	, nI = max(iXc)

	# nXc: number of cols in the condition-level predictor matrix
	, nXc = ncol(Xc)

	# rXc: number of rows in the condition-level predictor matrix
	, rXc = nrow(Xc)

	# nY: num entries in the observation vectors
	, nY = length(Y)

)

# double-check:
glimpse(data_for_stan)

# set the model path (automated bc in this repo there's only one)
mod_path = fs::dir_ls(
	path = 'stan'
	, glob = '*.stan'
)

# set the model centered/non-centeredness
#  generally, if *either* nI_per_group *or* num_Y_per_ic is small, non-centered will sample better than centered
data_for_stan$centered = FALSE

# conversion to 1/0 for stan
data_for_stan$centered = as.numeric(data_for_stan$centered)

# set the posterior path (automated but you could do your own if you had multiple models)
(
	mod_path
	%>% fs::path_file()
	%>% fs::path_ext_remove()
	%>% paste0(
		ifelse(data_for_stan$centered,'_c','_nc')
	)
	%>% fs::path(
		'posteriors'
		, .
		, ext = 'netcdf4'
	)
) -> post_path

# ensure model is compiled
aria:::check_syntax_and_maybe_compile(mod_path)

# compose
aria::compose(
	data = data_for_stan
	, code_path = mod_path
	, out_path = post_path
	, overwrite = T
)

# how long did it take?
aria::marginalia()$time

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
		, rebfmi = var(energy)/(sum(diff(energy)^2)/n()) # n.b. reciprocal of typical EBFMI, so bigger=bad, like rhat
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

# View those with suspect r̂
(
	par_summary
	%>% filter(rhat>1.01)
	%>% (function(suspects){
		if(nrow(suspects)>1){
			View(suspects)
		}
		return(paste('# suspect parameters:',nrow(suspects)))
	})()
)

# Viz recovery of (some) non-correlation parameters ----
(
	post$draws(variables=c('Z','iZc_sd'))
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = 'variable'
	)
	%>% left_join(
		bind_rows(
			tibble(
				true = as.vector(t(Z))
				, variable = paste0(
					'Z['
					,rep(1:nXc2,times=nXg)
					,','
					,rep(1:nXg,each=nXc2)
					,']'
				)
			)
			, tibble(
				true = iZc_sd
				, variable = paste0('iZc_sd[',1:length(true),']')
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
	# %>% mutate(variable = factor_1d(variable))
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


# Viz recovery of correlations ----
(
	post$draws(variables='iZc_r_vec')
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = 'variable'
	)
	%>% left_join(
		tibble(
			true = iZc_r_vec
			, variable = case_when(
				length(true)==1 ~ 'r'
				, T ~ paste0('iZc_r_vec[',1:length(true),']')
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

