
#comparing list of integers component wise
compwiseless := proc(a::list(nonnegint),b::list(nonnegint))::truefalse; 
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		local j::posint:=1,test::truefalse:=evalb(a[j]<=b[j]);
		while test and j<numelems(a) do:
			j:=j+1;
			test:=evalb(a[j]<=b[j]) 
		end do; 
		return test 
	end proc: