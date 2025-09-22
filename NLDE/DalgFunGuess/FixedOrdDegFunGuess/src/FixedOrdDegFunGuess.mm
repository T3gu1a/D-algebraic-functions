
FixedOrdDegFunGuess:= proc(Sinit::list,
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
		       linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense):=AlgebraicFunction,
			      $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
		       nL::posint:=numelems(Sinit),termfree::truefalse:=true,ADE::algebraic,RE::algebraic,	
		       Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),
		       correct::truefalse:=false,Arbconst::list,REcheck::algebraic,
		       l::list(nonnegint),Ll::list(list),m::nonnegint,degCoeffs::list(nonnegint);
		#minimal number of unknown
		l:=[degPoly$N];
		l:=prevlistnumber(degPoly,l);
		while l<> FAIL and not(correct) do
			Ll:=AllListPermutations(l);
			j:=1;
			M:=add(Ll[j])+N;
			while j<=numelems(Ll) and not(correct) and termfree do
				degCoeffs:=Ll[j];
				V:=[seq(c[i],i=0..M-1)];
				ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j)
								  ,i=0..degCoeffs[j]),j=1..N);
				RE:=ADEtoRE(ADE,Y,A,K);
				Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=0..M-1)];
				NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
				Eq:=subs(NegInd,Eq);
				termfree:=evalb(indets(Eq,a('integer'))={});
				if termfree then
					S:=SolveTools:-Linear(Eq,V,method=linsolver);
					if S<>NULL then
						if remove(v->rhs(v)=0,S)={} or has(map(rhs,S),a) then
							S:=NULL
						else
							REcheck, S, correct:=checkSol(S,RE,NegInd,Sinit,M,nL,a,n)
						end if
					end if
				end if;
				j:=j+1
			end do;
			l:=prevlistnumber(degPoly,l);
			termfree:=true
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
	
prevlistnumber := proc(maxn::nonnegint, L::list(nonnegint))
		local L1::list(nonnegint), N::posint; 
		if L = [0 $ (numelems(L) - 1), maxn] then 
			return FAIL 
		end if; 
		L1 := convert(L, base, maxn + 1, 10); 
		N := parse(StringTools:-Join(map(n -> convert(n, string), ListTools:-Reverse(L1)), "")); 
		convert(N - 1, base, maxn + 1)
	end proc;
	
AllListPermutations := proc(L::list(nonnegint))
		local p, S; 
		S := [seq(StringTools:-Permute(StringTools:-Join(map(n -> convert(n, string), L), ""), p),
			p = combinat['permute'](numelems(L)))]; 
		map(n -> map(parse, convert(n, list)), S) 
	end proc:
	
