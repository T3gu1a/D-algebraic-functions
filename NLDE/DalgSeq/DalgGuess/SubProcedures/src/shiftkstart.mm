#find the order of deltakshift correspoding to a specific shift
shiftkstart := proc(n::nonnegint,k,F,z,$)
		local j:=0,df;
		df:=deltakshift(F,z,k,j);
		while REorders(df,F)[1]<n do
			j:=j+1;
			df:=deltakshift(F,z,k,j)
		end do;
		return j
	end proc: