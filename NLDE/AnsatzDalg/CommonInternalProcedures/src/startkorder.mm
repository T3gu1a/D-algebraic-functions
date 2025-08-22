
#Interesting combinatorics: the deltakorder of the differential
#equation of degree k and maximum order n is binomial(k+n,k)

startkorder := proc(n::nonnegint,k,F,z,$)
		#local j:=k,df;
		#df:=deltakdiff(F,z,k,j);
		#while PDEtools:-difforder(df,z)<n do
		#	j:=j+1;
		#	df:=deltakdiff(F,z,k,binomial(j,k))
		#end do;
		#return binomial(j,k)
		return binomial(k+n,k)
	end proc: