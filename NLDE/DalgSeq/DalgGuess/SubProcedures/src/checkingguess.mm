#code to verify the guessed equation with the unused terms
checkingguess := proc(Sol::list,
				  REsol::algebraic,
				  Sinit::list,
				  M::posint,
				  nL::nonnegint,
				  A::anyfunc(name),
				  $)
		local S::list, RE::algebraic, checkL::list, checkset::set,i::nonnegint,
		      ord::nonnegint,maxind::nonnegint,n:=op(A);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		S:=map(simplify,Sol);
		RE:=subs(S,REsol);
		ord:=REorders(RE,A)[1];
		maxind:=nL-ord-1;
		checkset:={seq(simplify(subs(Sinit,eval(RE,n=i))),i=maxind..M,-1)};
		return RE, S, evalb(checkset = {0})
	end proc: