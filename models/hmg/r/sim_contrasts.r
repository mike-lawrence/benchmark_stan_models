sim_contrasts = function(n,var_prefix){(
	1:n
	%>% map(.f=function(x){
		factor(c('lo','hi'))
	})
	%>% (function(x){
		names(x) = paste0(var_prefix,1:length(x))
		return(x)
	})
	%>% cross_df()
	%>% (function(x){
		get_contrast_matrix(
			data = x
			, formula = as.formula(paste('~',paste0(var_prefix,1:ncol(x),collapse='*')))
			, contrast_kind = halfsum_contrasts
		)
	})
)}
