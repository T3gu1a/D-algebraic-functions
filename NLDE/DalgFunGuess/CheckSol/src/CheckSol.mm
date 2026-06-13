
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
		S:=map(normal,Sol);
		RE:=subs(S,REsol);
		checkL:=[op(NegInd),op(Sinit)];
		checkset:={seq(normal(subs(checkL,eval(RE,[n=i,Sum=add]))),i=(nL-numelems(Sol)-1)..nL)};
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
		S:=map(normal,Sol);
		ADE:=subs(S,ADEsol);
		checkADE:=expand(eval(ADE,y(x)=Lf));
		deg:= ldegree(checkADE,x);
		return ADE, S, evalb(checkADE=0 or deg>=nL-PDEtools:-difforder(ADE,x))
	end proc:

LetGenerateMatrix := proc(Eq::list,V::list,n::integer)
		local A, B,i,j;
		A := Matrix(n, n, [ seq([ seq( coeff(Eq[i], V[j]), j=1..n) ], i=1..n) ]);
		B := Vector(n, [ seq( (-subs(map(v -> v=0, V), Eq[i])), i=1..n) ]);
		return A,B
	end proc:	
	
