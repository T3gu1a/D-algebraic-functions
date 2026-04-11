
FixedOrdDegFunGuess:= proc(        Sinit::list,
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
			       linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense),
			  inputConstants::set(name),
			              $)::Or(identical(FAIL),`=`);
			option `Copyright (c) 2025 Bertrand Teguia T.`;
			description "Looking for an equation among all possible equations of the given maximum polynomial degree";
			local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
			       nL::posint:=numelems(Sinit),hasterm::truefalse:=false,ADE::algebraic,RE::algebraic,	
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
				if M <= nL then
					while j<=numelems(Ll) and not(correct) and not(hasterm) do
						degCoeffs:=Ll[j];
						V:=[seq(c[i],i=0..M-1)];
						ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j)
										  ,i=0..degCoeffs[j]),j=1..N);
						RE:=ADEtoRE(ADE,Y,A,K);
						Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=0..M-1)];
						NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
						Eq:=subs(NegInd,Eq);
						S:=SolveTools:-Linear(Eq,V,method=linsolver);
						if S<>NULL then
							if remove(v->rhs(v)=0,S)={} or has(map(rhs,S),a) then
								hasterm:=has(Eq,a);
								S:=NULL
							else
								REcheck, S, correct:=checkSol(S,RE,NegInd,Sinit,M,nL,a,n)
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
				Arbconst:=sort([op(remove(has,indets(REcheck),a) minus ({n,op(K)} union inputConstants))]);
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
	
FFixedOrdDegFunGuess:= proc(           Lf::algebraic,
			           degADE::posint,
			          degPoly::nonnegint,
			                Y::anyfunc(name),
			                N::nonnegint,
			                y::name,
				        x::name,
			        linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense),
		           inputConstants::set(name),
				 sparsity::fraction,
			               $)::Or(identical(FAIL),`=`);
			option `Copyright (c) 2026 Bertrand Teguia T.`;
			description "Looking for an equation among all possible equations of the given maximum polynomial degree";
			local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
			       nL::posint:=degree(Lf,x)+1,ADE::algebraic,polEq::algebraic,	
			       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),
			       correct::truefalse:=false,Arbconst::list,ADEcheck::algebraic,
			       MnL::posint,ZerosV,zV,zzV::list,unkV::list;
			#minimal number of unknown
			M:=(degPoly+1)*N;
			V:=[seq(c[i],i=0..M-1)];
			if M>nL then
				MnL:=ceil((1+sparsity)*M-nL);
				ZerosV:=Iterator:-Combination(M, MnL);
				for zV in ZerosV do
					zzV := [seq(c[i], i in zV)];
					zzV:=map(t->t=0,zzV);
					unkV:=subs(zzV,V);
					ADE:=add(add(unkV[(degPoly+1)*(j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
					polEq:=expand(eval(ADE,Y=Lf));
					unkV:=remove(t->t=0,unkV);
					Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..numelems(unkV)]; #[seq(coeff(polEq,x,i),i=0..M-1)];
					S:=SolveTools:-Linear(Eq,unkV,method=linsolver);
					if S<>NULL then
						if remove(v->rhs(v)=0,S)={} then
							S:=NULL
						else
							ADEcheck, S, correct:=polcheckSol(S,ADE,Lf,nL,y,x)
						end if
					end if;
					if correct then
						break
					end if
				end do
			else
				ADE:=add(add(V[(degPoly+1)*(j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
				polEq:=expand(eval(ADE,Y=Lf));
				Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; #[seq(coeff(polEq,x,i),i=0..M-1)];
				S:=SolveTools:-Linear(Eq,V,method=linsolver);
				if S<>NULL then
					if remove(v->rhs(v)=0,S)={} then
						S:=NULL
					else
						ADEcheck, S, correct:=polcheckSol(S,ADE,Lf,nL,y,x)
					end if
				end if
			end if;
			if correct then
				ADE:=subs(S,ADE);
				Arbconst:=sort([op(indets(ADEcheck) 
					minus (inputConstants union {x,y,seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))}))]);
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
	
FFixedOrdDegFunGuess2:= proc(          Lf::algebraic,
			           degADE::posint,
			          degPoly::nonnegint,
			                Y::anyfunc(name),
			                N::nonnegint,
			                y::name,
				        x::name,
			        linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense),
		           inputConstants::set(name),
			               $)::Or(identical(FAIL),`=`);
			option `Copyright (c) 2026 Bertrand Teguia T.`;
			description "Looking for an equation among all possible equations of the given maximum polynomial degree";
			local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
			       nL::posint:=degree(Lf,x)+1,ADE::algebraic,polEq::algebraic,	
			       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),
			       correct::truefalse:=false,Arbconst::list,ADEcheck::algebraic,
			       l::list(nonnegint),Ll::list(list),m::nonnegint,degCoeffs::list(nonnegint);
			#minimal number of unknown
			l:=[degPoly$N];
			l:=prevlistnumber(degPoly,l);
			while l<> FAIL and not(correct) do
				Ll:=AllListPermutations(l);
				j:=1;
				M:=add(Ll[j])+N;
				if M <= nL then
					while j<=numelems(Ll) and not(correct) do
						degCoeffs:=Ll[j];
						V:=[seq(c[i],i=0..M-1)];
						ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j)
										  ,i=0..degCoeffs[j]),j=1..N);
						polEq:=expand(eval(ADE,Y=Lf));
						Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; #[seq(coeff(polEq,x,i),i=0..M-1)];
						S:=SolveTools:-Linear(Eq,V,method=linsolver);
						if S<>NULL then
							if remove(v->rhs(v)=0,S)={} then
								S:=NULL
							else
								ADEcheck, S, correct:=polcheckSol(S,ADE,Lf,nL,y,x)
							end if
						end if;
						j:=j+1
					end do
				end if;
				l:=prevlistnumber(degPoly,l)
			end do;
			if correct then
				ADE:=subs(S,ADE);
				Arbconst:=sort([op(indets(ADEcheck) 
					minus (inputConstants union {x,y,seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))}))]);
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
	
