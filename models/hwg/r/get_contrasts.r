get_contrast_matrix = function(
	data
	, formula
	, contrast_kind = NULL
){
	if (inherits(data, "tbl_df")) {
		data = as.data.frame(data)
	}
	vars = attr(terms(formula),'term.labels')
	vars = vars[!grepl(':',vars)]
	if(length(vars)==1){
		data = data.frame(data[,vars])
		names(data) = vars
	}else{
		data = data[,vars]
	}
	vars_to_rename = NULL
	for(i in vars){
		if(is.character(data[,i])){
			data[,i] = factor(data[,i])
		}
		if( is.factor(data[,i])){
			if(length(levels(data[,i]))==2){
				vars_to_rename = c(vars_to_rename,i)
			}
			if(!is.null(contrast_kind) ){
				contrasts(data[,i]) = contrast_kind
			}
		}
	}
	mm = model.matrix(data=data,object=formula)
	dimnames(mm)[[2]][dimnames(mm)[[2]]=='(Intercept)'] = '(I)'
	for(i in vars_to_rename){
		dimnames(mm)[[2]] = gsub(paste0(i,1),i,dimnames(mm)[[2]])
	}
	attr(mm,'formula') = formula
	attr(mm,'data') = data
	return(mm)
}
