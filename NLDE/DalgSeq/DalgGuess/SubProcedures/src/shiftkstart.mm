#find the deltakshift order correspoding to a specific shift
shiftkstart := proc(n::nonnegint,k,$)
		return binomial(k+n,k)
	end proc: