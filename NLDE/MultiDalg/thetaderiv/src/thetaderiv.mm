
#The exported theta derivation
thetaderiv := proc(expr,V::list(name),k::nonnegint,$)::algebraic;
		return derivation(expr,V,k)
	end proc: