
#Remark: minmaxorder ... ((n + 1)*(n^4 + 19*n^3 + 136*n^2 + 444*n + 600))/120

startkorder := proc(n::nonnegint,k,F,z,$)
		local j:=k+1,df;
		df:=deltakdiff(F,z,k,j);
		while PDEtools:-difforder(df,z)<n do
			j:=j+1;
			df:=deltakdiff(F,z,k,j)
		end do;
		return j
	end proc: