
#The D-algebraic function guesser (finite field)
modDalgFunGuess:= proc(L::list,
	         {degADE::posint:=2,
	         degPoly::nonnegint:=2,
	          devars::anyfunc(name):=NULL,
	    startfromord::nonnegint:=0,
	      allPolyDeg::truefalse:=false,
	         modulus::posint:=7},
		      $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "Guessing D-algebraic functions (finding their differential equations)";
		local  Y::anyfunc(name),y::name,x::name,A::anyfunc(name),a::name,n::name,i::nonnegint,
		       c::nothing,N::posint,M::posint,V::list,j::nonnegint,Sinit::list(eq),k::name,
		       nL::posint:=numelems(L),hasterm::truefalse:=true,ADE::algebraic,RE::algebraic,	
		       Nmax:=ceil(nL/(degPoly+1)),K::list,Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),
		       NDE::algebraic,correct::boolean:=false,Arbconst::list,REcheck::algebraic,NRE::algebraic,
		       Terms::list(integer),Meqs,beqs,Aindets;
		       
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
		a:=`tools/genglobal`('a');
		n:=`tools/genglobal`('n');
		interface(warnlevel=0);
		for j from 0 to degADE-1 do
			cat(k,j):=`tools/genglobal`('k')
		end do;
		K:=[seq(cat(k,j),j=0..degADE-1)];
		interface(warnlevel=4);
		A:=a(n);
		Sinit:=[seq(a(i-1)=Terms[i],i=1..nL)];
		#underdetermined system
		if M > nL then
			return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,modulus),FAIL)
		end if;
		#initialization - ADE and RE of the first iteration
		ADE:=add(add(V[(degPoly+1)*(j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
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
		Eq:=subs(Aindets,Eq);
		#solving the linear system
		Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
		S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
		S:= ifelse(type(S,list(algebraic)),S,NULL);
		#S:=try msolve({op(Eq)},modulus) catch : NULL  end try;
		if S<>NULL then
			if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets)) then
				hasterm:=has(Eq,map(rhs,Aindets));
				S:=NULL;
				#too little initial values
				if hasterm then
					return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,modulus),FAIL)
				end if
			else
				#checking the solution
				S:=[seq(V[i]=S[i],i=1..M)];
				REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
			end if
		end if;
		N:=N+1;
		while not(correct) and N<Nmax do
			V:=[op(V),seq(c[i],i=M..M+degPoly)];
			NDE:=add(V[M+i]*x^(i-1)*AnsatzDalg:-deltakdiff(Y,x,degADE,N),i=1..degPoly+1);
			ADE:=NDE+ADE;
			NRE:=ADEtoRE(NDE,Y,A,K);
			RE:=NRE+RE;
			Eq:=[seq(Eq[i+1]+subs(Sinit,eval(NRE,[n=i,Sum=add]) mod modulus) mod modulus,i=0..M-1)];
			Eq:=[op(Eq),seq(subs(Sinit,eval(RE,[n=i,Sum=add]) mod modulus) mod modulus,i=M..M+degPoly+1)];
			M:=M+degPoly+1;
			NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
			Eq:=subs(NegInd,Eq);
			Aindets:=[op(indets(Eq,a('integer')))];
			Aindets:=[seq(Aindets[j]=cat(a,j),j=1..numelems(Aindets))];
			Eq:=subs(Aindets,Eq);
			#solving the linear system
			Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
			S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
			S:= ifelse(type(S,list(algebraic)),S,NULL);
			#S:=try msolve({op(Eq)},modulus) catch : NULL  end try;
			if S<>NULL then
				if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets)) then
					hasterm:=has(Eq,map(rhs,Aindets));
					S:=NULL;
					if hasterm then
						return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,modulus),FAIL)
					end if
				else
					#verification step
					S:=[seq(V[i]=S[i],i=1..M)];
					REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
				end if
			end if;
			N:=N+1
		end do;
		if correct then 
			ADE:=subs(S,ADE);
			Arbconst:=sort([op(remove(has,indets(REcheck),a) minus {n,op(K)})]);
			if Arbconst <> [] then
				`tools/genglobal`('_C',{},'reset');
				Arbconst:=map(v->v=`tools/genglobal`('_C'),Arbconst);
				ADE:=subs(Arbconst,ADE)
			end if;
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))},'distributed');
			return ADE=0
		else
			return ifelse(allPolyDeg,modFixedOrdDegFunGuess(Sinit,degADE,degPoly,Y,A,N,y,x,a,n,K,modulus),FAIL)
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
	
$include <NLDE/DalgFunGuess/modDalgFunGuess/modFixedOrdDegFunGuess/src/modFixedOrdDegFunGuess.mm>