

modFixedOrdDegFunGuess:= proc(Sinit::list,
			     degADE::posint,
			    degPoly::nonnegint,
				  Y::anyfunc(name),
				  A::anyfunc(name),
				  N::nonnegint,
				  y::name,
				  x::name,
				  a::name,
				  n::name,
				  K::list,
			    modulus::posint,
				 $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
		       nL::posint:=numelems(Sinit),hasterm::truefalse:=false,ADE::algebraic,RE::algebraic,	
		       Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),
		       correct::truefalse:=false,Arbconst::list,REcheck::algebraic,Aindets,Meqs,beqs,
		       l::list(nonnegint),Ll::list(list),m::nonnegint,degCoeffs::list(nonnegint);
		#minimal number of unknown
		l:=[degPoly$N];
		l:=prevlistnumber(degPoly,l);
		while l<> FAIL and not(correct) do
			Ll:=AllListPermutations(l);
			j:=1;
			M:=add(Ll[j])+N;
			if M <= nL then
				while j<=numelems(Ll) and not(correct) and not(hasterm) do
					degCoeffs:=Ll[j];
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					RE:=ADEtoRE(ADE,Y,A,K);
					Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add]) mod modulus) mod modulus,i=0..M-1)];
					NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
					Eq:=subs(NegInd,Eq);
					Aindets:=[op(indets(Eq,a('integer')))];
					Aindets:=[seq(Aindets[j]=cat(a,j),j=1..numelems(Aindets))];
					Eq:=subs(Aindets,Eq);
					#solving the linear system
					Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
					S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
					S:= ifelse(type(S,list(algebraic)),S,NULL);
					#termfree:=evalb(indets(Eq,a('integer'))={});
					#S:=try msolve({op(Eq)},modulus) catch : NULL  end try;
					if S<>NULL then
						if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets))  then
							hasterm:=has(Eq,map(rhs,Aindets));
							S:=NULL
						else
							S:=[seq(V[i]=S[i],i=1..M)];
							REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
						end if
					end if;
					j:=j+1
				end do
			end if;
			l:=prevlistnumber(degPoly,l);
			hasterm:=false
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
			return FAIL
		end if          
	end proc: