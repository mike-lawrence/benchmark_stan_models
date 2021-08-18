pow = function(x,y) x^y
pquick = function(x,chance,lapse,threshold,slope){(
	chance
	+ (
		(1-chance-lapse)
		* (
			1-pow(2,-pow(x/threshold,slope))
		)
	)
)}
