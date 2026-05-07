
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

polcheckSol:= proc(Sol::Or(list,set),
		ADEsol::algebraic,
		    Lf::algebraic,
		    nL::nonnegint,
		     y::name,
		     x::name,
		     $)
		local S::list, ADE::algebraic, i::nonnegint,
		      checkADE::algebraic, deg::extended_numeric;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		S:=map(simplify,Sol);
		ADE:=subs(S,ADEsol);
		checkADE:=expand(eval(ADE,y(x)=Lf));
		deg:= ldegree(checkADE,x);
		return ADE, S, evalb(checkADE=0 or deg>=nL-PDEtools:-difforder(ADE,x))
	end proc:	
	
