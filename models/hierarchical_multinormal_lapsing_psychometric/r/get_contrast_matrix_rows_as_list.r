get_contrast_matrix_rows_as_list = function(
	data
	, formula
	, contrast_kind = NULL
){
	mm = get_contrast_matrix(data, formula, contrast_kind )
	mm_list = lapply(
		X = seq_len(nrow(mm))
		, FUN = function(i){
			array(
				mm[i,]
				, dim = c(1,ncol(mm))
				, dimnames = list(
					NULL
					, dimnames(mm)[[2]]
				)
			)
		}
	)
	return(mm_list)
}
