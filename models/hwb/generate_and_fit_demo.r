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
sim_pars = lst(

	#parameters you can play with
	num_subj = 30 #number of subjects, must be an integer >1
	, num_vars = 4 #number of 2-level variables manipulated as crossed and within each subject, must be an integer >0
	, num_trials = 10 #number of trials per subject/condition combo, must be an integer >1

	#the rest of these you shouldn't touch
	, num_coef = 2^(num_vars)
	, mu = rnorm(num_coef)
	, sigma = rweibull(num_coef,2,1)
	, r_mat = rethinking::rlkjcorr(1,num_coef,eta=1)

	, r = r_mat[lower.tri(r_mat)]
)

#compute the contrast matrix
contrast_matrix =
	(
		1:sim_pars$num_vars
		%>% map(.f=function(x){
			factor(c('lo','hi'))
		})
		%>% (function(x){
			names(x) = paste0('v',1:sim_pars$num_vars)
			return(x)
		})
		%>% cross_df()
		%>% (function(x){
			get_contrast_matrix(
				data = x
				, formula = as.formula(paste('~',paste0('v',1:sim_pars$num_vars,collapse='*')))
				, contrast_kind = halfsum_contrasts
			)
		})
	)

#get coefficients for each subject
subj_coef =
	(
		#subj coefs as mvn
		rethinking::rmvnorm2(
			n = sim_pars$num_subj
			, Mu = sim_pars$mu
			, sigma = sim_pars$sigma
			, Rho = sim_pars$r_mat
		)
		#add names to columns
		%>% (function(x){
			dimnames(x)=list(NULL,paste0('X',1:ncol(x)))
			return(x)
		})
		#make a tibble
		%>% as_tibble(.name_repair='unique')
		#add subject identifier column
		%>% mutate(
			subj = 1:sim_pars$num_subj
		)
	)

# get condition means implied by subject coefficients and contrast matrix
subj_cond =
	(
		subj_coef
		%>% group_by(subj)
		%>% summarise(
			(function(x){
				out = attr(contrast_matrix,'data')
				out$cond_mean = as.vector(contrast_matrix %*% t(x))
				return(out)
			})(cur_data())
			, .groups = 'drop'
		)
	)

# get noisy measurements in each condition for each subject
dat =
	(
		subj_cond
		%>% expand_grid(trial = 1:sim_pars$num_trials)
		%>% mutate(
			obs = rbinom(n(),1,plogis(cond_mean))
		)
		%>% select(-cond_mean)
	)


# Compute inputs to Stan model ----

(
	dat
	#first collapse to subject-by-condition stats
	%>% group_by(
		across(c(
			-obs
			, -trial
		))
	)
	%>% summarise(
		num_trials_total = n()
		, num_successes = sum(obs)
		, .groups = 'drop'
	)
	#add the contrast matrix columns
	%>% mutate(
		contrasts = get_contrast_matrix_rows_as_list(
			data = .
			# the following compilcated specification of the formula is a by-product of making this example
			# work for any value for sim_pars$num_vars; normally you would do something like this
			# (for 2 variables for example):
			# formula = ~ v1*v2
			, formula = as.formula(paste0('~',paste0('v',1:sim_pars$num_vars,collapse='*')))
			# half-sum contrasts are nice for 2-level variables bc they yield parameters whose value
			# is the difference between conditions
			, contrast_kind = halfsum_contrasts
		)
	)
) ->
	dat_summary_with_contrasts

#if there's no missing data (as in this synthetic example), we could proceed with modelling,
#  but to demonstrate how to handle missing data, we'll do a bit of extra work:

#Uncomment next line to induce data-missingness (causes some Ss to be missing whole conditions)
# dat_summary_with_contrasts %<>% slice_sample(prop=.9)

#We need a tbl containing the combinations of the unique contrasts and subjects

#first get the unique contrasts (will serve as input to Stan)
(
	dat_summary_with_contrasts
	%>% distinct(contrasts)
	%>% arrange(contrasts) #important
) ->
	unique_contrasts

#now cross with unique subjects
(
	unique_contrasts
	%>% expand_grid(
		tibble(subj = unique(dat$subj))
	)
	%>% arrange(contrasts, subj) #important
	# add a column containing the row index
	%>% mutate(
		row = 1:n()
	)
	#join with dat_summary
	%>% right_join(
		dat_summary_with_contrasts
		, by = c('contrasts','subj')
	)
) -> dat_summary_with_contrasts_and_row

#now get the unique contrast matrix
(
	unique_contrasts
	%>% unnest(contrasts)
	%>% as.matrix()
) -> uw


# package for stan & sample ----

data_for_stan = lst( #lst permits later entries to refer to earlier entries

	####
	# Entries we need to specify ourselves
	####

	# uw: unique within predictor matrix
	uw = uw

	# num_obs: number of observations in each subject/condition
	, num_total = dat_summary_with_contrasts_and_row$num_trials_total

	# sum_obs: number of 1's in each subject/condition
	, num_successes = dat_summary_with_contrasts_and_row$num_successes

	, row = dat_summary_with_contrasts_and_row$row

	####
	# Entries computable from the above
	####

	, k = ncol(uw)
	, n = length(unique(dat_summary_with_contrasts_and_row$subj))
	, m = length(row) # should be n*k, unless missing data

)

# double-check:
glimpse(data_for_stan)

# compose
aria::compose(
	data = data_for_stan
	, code_path = 'stan/hwb.stan'
	, out_path = 'nc/hwb.nc'
	, overwrite = T
)

# check posterior diagnostics ----

post = aria::coda('nc/hwb.nc')

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
	post$draws(variables=c('mu','sigma','r','z_z'))
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
	post$draws(variables=c('mu','sigma'))
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = 'variable'
	)
	%>% left_join(
		bind_rows(
			tibble(
				true = sim_pars$mu
				, variable = paste0('mu[',1:length(true),']')
			)
			, tibble(
				true = sim_pars$sigma
				, variable = paste0('sigma[',1:length(true),']')
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


# Viz recovery of correlations ----
(
	post$draws(variables='r')
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = 'variable'
	)
	%>% left_join(
		tibble(
			true = sim_pars$r_mat[lower.tri(sim_pars$r_mat)]
			, variable = case_when(
				length(true)==1 ~ 'r'
				, T ~ paste0('r[',1:length(true),']')
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
