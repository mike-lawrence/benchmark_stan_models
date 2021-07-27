factor_1d = function(x){(
	tibble(x=x)
	%>% mutate(
		y = str_remove(x,']')
		, y = str_replace(y,fixed('['),',')
	)
	%>% separate(
		y
		, into = c('variable','index')
		, sep = ','
		, fill = 'right'
		, convert = T
		, remove = F
	)
	%>% arrange(variable,index)
	%>% mutate(
		x = factor(x,levels=unique(x))
	)
	%>% pull(x)
)}
