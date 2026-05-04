
#The D-algebraic function guesser (finite field)
modDalgFunGuess:= proc(L::list,
	         {degADE::posint:=2,
	         degPoly::nonnegint:=2,
          termsOfDegPoly::posint:=1,
	          devars::anyfunc(name):=NULL,
	    startfromord::nonnegint:=0,
	      allPolyDeg::truefalse:=false,
	        sparsity::Or(And(positive, fraction), identical(0)):=0,
	        approach::identical(recurrence,polynomialsubs):=polynomialsubs,
	    maxIteration::Or(posint,identical(infinity)):=infinity,
	         modulus::posint:=7,
          inputConstants::set(name):={}},
		      $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "Guessing D-algebraic functions (finding their differential equations)";
		local  Y::anyfunc(name),y::name,x::name,A::anyfunc(name),a::name,n::name,i::nonnegint,
		       c::nothing,N::posint,M::posint,V::list,j::nonnegint,Sinit::list(eq),k::name,Lf::algebraic,
		       nL::posint:=numelems(L),hasterm::truefalse:=true,ADE::algebraic,RE::algebraic,polEq::algebraic,	
		       Nmax:=ceil(nL/(degPoly+1)),K::list,Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),
		       NDE::algebraic,correct::boolean:=false,Arbconst::list,REcheck::algebraic,NRE::algebraic,NpolEq::algebraic,
		       Terms::list(integer),Meqs,beqs,Aindets,ADEcheck::algebraic,Param::set(name):=inputConstants;
		       
		Terms:= try map(term -> term mod modulus, L) catch : [] end try;
		if Terms = [] then
			return FAIL
		end if;
		#minimal deltak order for starting order startfromord
		N:=binomial(degADE+startfromord,degADE);
		#minimal number of unknown
		M:=(degPoly+1)*N;
		V:=[seq(c[i],i=0..M-1)];
		#ADE variables
		if devars=NULL then 
			`tools/genglobal`('y',{},'reset');
			`tools/genglobal`('x',{},'reset');
			y:=`tools/genglobal`('y');
			x:=`tools/genglobal`('x');
			Y:=y(x)
		else 
			y:=op(0,devars);
			x:=op(1,devars);
			#correctness of the input variables
			if y=x then
				error "invalid input: the variables of the differential equation must be distinct"
			end if;
			Y:=y(x)
		end if;
		#Recurrence variables
		if approach=recurrence then
			a:=`tools/genglobal`('a');
			n:=`tools/genglobal`('n');
			interface(warnlevel=0);
			for j from 0 to degADE-1 do
				cat(k,j):=`tools/genglobal`('k')
			end do;
			K:=[seq(cat(k,j),j=0..degADE-1)];
			interface(warnlevel=4);
			A:=a(n);
			Sinit:=[seq(a(i-1)=Terms[i],i=1..nL)]
		else
			Sinit:=[seq(a(i-1)=Terms[i],i=1..nL)];
			Lf:=add(a(i)*x^i,i=0..nL-1)
		end if;
		#underdetermined system
		if M > nL then
			if approach=recurrence then
				if nL-N<termsOfDegPoly*degPoly then
					N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
				end if;
				return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,maxIteration,maxmodulus,inputConstants),FAIL)
			else
				if sparsity<>0 then
					return ifelse(allPolyDeg,modFFixedOrdDegFunGuess(Lf,Sinit,a,degADE,degPoly,Y,N,y,x,maxIteration,modulus,inputConstants,sparsity),FAIL)
				else
					if nL-N<termsOfDegPoly*degPoly then
						N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
					end if;
					return ifelse(allPolyDeg,modFFixedOrdDegFunGuess2(Lf,Sinit,a,degADE,degPoly,Y,N,y,x,maxIteration,modulus,inputConstants),FAIL)
				end if
			end if
		end if;
		#initialization - ADE and RE of the first iteration
		ADE:=add(add(V[(degPoly+1)*(j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
		if approach=recurrence then
			RE:=ADEtoRE(ADE,Y,A,K);
			#write the RE for non-negative indices
			#build the linear system and solve it
			Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add]) mod modulus) mod modulus,i=0..M-1)];
			#Eq:=map(i -> subs(Sinit, eval(RE, [n = i, Sum = add]) mod modulus) mod modulus, [$0 .. M-1]);
			#NegInd: list for substituting terms with negative indices to zero
			NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
			Eq:=subs(NegInd,Eq);
			Aindets:=[op(indets(Eq,a('integer')))];
			Aindets:=[seq(Aindets[j]=cat(a,j),j=1..numelems(Aindets))];
			Eq:=subs(Aindets,Eq)
		else
			polEq:=expand(eval(ADE,Y=Lf) mod modulus) mod modulus;
			polEq:=subs(Sinit,polEq);
			Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M] #[seq(coeff(polEq,x,i),i=0..M-1)]
		end if;
		#solving the linear system
		Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
		S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
		S:= ifelse(type(S,list(algebraic)),S,NULL);
		#S:=try msolve({op(Eq)},modulus) catch : NULL  end try;
		if S<>NULL then
			if approach=recurrence then
				if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets)) then
					hasterm:=has(Eq,map(rhs,Aindets));
					S:=NULL;
					#too little initial values
					if hasterm then
						if nL-N<termsOfDegPoly*degPoly then
							N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
						end if;
						return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,modulus,inputConstants),FAIL)
					end if
				else
					#checking the solution
					S:=[seq(V[i]=S[i],i=1..M)];
					REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
				end if
			else
				if remove(v->v=0,S)=[] then
					S:=NULL
				else
					S:=[seq(V[i]=S[i],i=1..M)];
					ADEcheck, S, correct:=modpolcheckSol(S,ADE,Sinit,a,nL,y,x,modulus)
				end if
			end if
		end if;
		N:=N+1;
		while not(correct) and N<Nmax do
			V:=[op(V),seq(c[i],i=M..M+degPoly)];
			NDE:=add(V[M+i]*x^(i-1)*AnsatzDalg:-deltakdiff(Y,x,degADE,N),i=1..degPoly+1);
			ADE:=NDE+ADE;
			if approach=recurrence then
				NRE:=ADEtoRE(NDE,Y,A,K);
				RE:=NRE+RE;
				Eq:=[seq(Eq[i+1]+subs(Sinit,eval(NRE,[n=i,Sum=add]) mod modulus) mod modulus,i=0..M-1)];
				Eq:=[op(Eq),seq(subs(Sinit,eval(RE,[n=i,Sum=add]) mod modulus) mod modulus,i=M..M+degPoly+1)];
				NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
				Eq:=subs(NegInd,Eq);
				Aindets:=[op(indets(Eq,a('integer')))];
				Aindets:=[seq(Aindets[j]=cat(a,j),j=1..numelems(Aindets))];
				Eq:=subs(Aindets,Eq)
			else 
				NpolEq:=expand(eval(NDE,Y=Lf) mod modulus) mod modulus;
				polEq:=polEq+subs(Sinit,NpolEq) mod modulus;
				Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M+degPoly+1] #[seq(coeff(polEq,x,i),i=0..M+degPoly+1)]
			end if;
			M:=M+degPoly+1;
			#solving the linear system
			Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
			S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
			S:= ifelse(type(S,list(algebraic)),S,NULL);
			#S:=try msolve({op(Eq)},modulus) catch : NULL  end try;
			if S<>NULL then
				if approach=recurrence then 
					if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets)) then
						hasterm:=has(Eq,map(rhs,Aindets));
						S:=NULL;
						if hasterm then
							if nL-N<termsOfDegPoly*degPoly then
								N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
							end if;
							return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,maxIteration,modulus,inputConstants),FAIL)
						end if
					else
						#verification step
						S:=[seq(V[i]=S[i],i=1..M)];
						REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
					end if
				else
					if remove(v->v=0,S)=[] then
						S:=NULL
					else
						S:=[seq(V[i]=S[i],i=1..M)];
						ADEcheck, S, correct:=modpolcheckSol(S,ADE,Sinit,a,nL,y,x,modulus)
					end if
				end if
			end if;
			N:=N+1
		end do;
		if correct then 
			ADE:=subs(S,ADE);
			if approach=recurrence then
				Param:=Param union {n,op(K)};
				Arbconst:= sort([op(remove(has,indets(REcheck),a) minus Param)])
			else
				Param:=Param union {x,y,seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))};
				Arbconst:= sort([op(indets(ADEcheck) minus Param)])
			end if;
			if Arbconst <> [] then
				`tools/genglobal`('_C',{},'reset');
				Arbconst:=map(v->v=`tools/genglobal`('_C'),Arbconst);
				ADE:=subs(Arbconst,ADE)
			end if;
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))},'distributed');
			return ADE=0
		else
			if approach=recurrence then
				if nL-N<termsOfDegPoly*degPoly then
					N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
				end if;
				return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,maxIteration,modulus,inputConstants),FAIL)
			else
				if sparsity<>0 then
					return ifelse(allPolyDeg,modFFixedOrdDegFunGuess(Lf,Sinit,a,degADE,degPoly,Y,N,y,x,maxIteration,modulus,inputConstants,sparsity),FAIL)
				else
					if nL-N<termsOfDegPoly*degPoly then
						N:=ifelse(termsOfDegPoly<Nmax,nL-termsOfDegPoly*degPoly,nL-rand(1..floor(Nmax/2))()*(degPoly+1))
					end if;
					return ifelse(allPolyDeg,modFFixedOrdDegFunGuess2(Lf,Sinit,a,degADE,degPoly,Y,N,y,x,maxIteration,modulus,inputConstants),FAIL)
				end if
			end if
		end if
	end proc:
	
modcheckSol:= proc(Sol::Or(list,set),
	         REsol::algebraic,
	        NegInd::list,
	         Sinit::list,
		     M::posint,
		    nL::nonnegint,
		     a::name,
		     n::name,
		     m::posint,
		     $)
		local S::list, RE::algebraic, checkL::list, checkset::set,i::nonnegint;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		S:=map(simplify,Sol);
		RE:=subs(S,REsol);
		checkL:=[op(NegInd),op(Sinit)];
		checkset:={seq(simplify(subs(checkL,eval(RE,[n=i,Sum=add]) mod m) mod m),i=(nL-numelems(Sol)-1)..nL)};
		#print(checkset);
		checkset:=remove(has,checkset,a);
		return RE, S, evalb(checkset in {{0},{}})
	end proc:

modpolcheckSol:= proc(Sol::Or(list,set),
		   ADEsol::algebraic,
		    Sinit::list,
			a::name,
		       nL::nonnegint,
			y::name,
			x::name,
			m::posint,
			$)
	local S::list, ADE::algebraic, i::nonnegint,
	      solf::algebraic:=add(a(i)*x^i,i=0..nL-1),checkADE::algebraic, deg::extended_numeric;
	option `Copyright (c) 2022 Bertrand Teguia T.`;
	S:=map(simplify,Sol);
	ADE:=subs(S,ADEsol);
	checkADE:=expand(eval(ADE,y(x)=solf)) mod m;
	checkADE:=expand(subs(Sinit,checkADE)) mod m;
	#checkADE:=expand(checkADE) mod m;
	deg:= ldegree(checkADE,x);
	return ADE, S, evalb(checkADE=0 or deg>=nL-PDEtools:-difforder(ADE,x))
end proc:
	
$include <NLDE/DalgFunGuess/modDalgFunGuess/modFixedOrdDegFunGuess/src/modFixedOrdDegFunGuess.mm>