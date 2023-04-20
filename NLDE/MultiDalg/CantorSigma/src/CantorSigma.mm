
#Cantor (non algebraic) k-tuple function
#This may soon be changed with the exact algebraic formula
CantorSigma := proc(tuple::Or(list(nonnegint),NULL):=NULL)::nonnegint; option remember;
		local s::nonnegint, k::posint, N::posint;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		if tuple=NULL then
			return 0
		end if;
		k:=numelems(tuple);
		s:=add(tuple);
		CantorInvSigma(k,s);
		N:=numelems(thetakTuple[k])-1;
		if add(thetakTuple[k][N]) >= s then
			return ListTools:-Search(tuple,convert(thetakTuple[k],list))-1
		else
			CantorInvSigma(k,N+1);
			return CantorSigma(tuple)
		end if
	end proc: