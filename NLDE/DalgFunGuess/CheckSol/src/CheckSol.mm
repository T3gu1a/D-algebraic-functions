
#checking the guess

checkSol:= proc(Sol::Or(list,set),
		      REsol::algebraic,
		     NegInd::list,
		      Sinit::list,
		          M::posint,
		         nL::nonnegint,
		          a::name,
			  n::name,
		          $)
		local S::list, RE::algebraic, checkL::list, checkset::set,i::nonnegint;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		S:=map(simplify,Sol);
		RE:=subs(S,REsol);
		checkL:=[op(NegInd),op(Sinit)];
		checkset:={seq(simplify(subs(checkL,eval(RE,[n=i,Sum=add]))),i=(nL-numelems(Sol)-1)..nL)};
		checkset:=remove(has,checkset,a);
		return RE, S, evalb(checkset in {{0},{}})
	end proc: