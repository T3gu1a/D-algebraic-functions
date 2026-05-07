
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
		       maxIteration::Or(posint,identical(infinity)),
			    modulus::posint,
		     inputConstants::set(name),
				 $)::Or(identical(FAIL),`=`);
			option `Copyright (c) 2025 Bertrand Teguia T.`;
			description "Looking for an equation among all possible equations of the given maximum polynomial degree";
			local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
			       nL::posint:=numelems(Sinit),hasterm::truefalse:=false,ADE::algebraic,RE::algebraic,total_perms,	
			       Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),ul,tl,freqs,val,
			       correct::truefalse:=false,Arbconst::list,REcheck::algebraic,Aindets,Meqs,beqs,randpick,
			       l::list(nonnegint),Ll::list(list),m::nonnegint,degCoeffs::list(nonnegint);
			
			l:=GenMaxlistnumber(N,degPoly,nL-N);
			ul:=sort([op({op(l)})]);
			tl:=Statistics:-Tally(l);
			freqs:=Array([seq(eval(val, tl), val = ul)], datatype = integer);
			total_perms:=MultinomialCount(freqs);
			while l<> FAIL and not(correct) do
				if total_perms < maxIteration then
					Ll:=Iterator:-Permute(l,N);
					M:=add(l)+N;
					for degCoeffs in Ll do:
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
						Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
						S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
						S:= ifelse(type(S,list(algebraic)),S,NULL);
						if S<>NULL then
							if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets))  then
								hasterm:=has(Eq,map(rhs,Aindets));
								S:=NULL
							else
								S:=[seq(V[i]=S[i],i=1..M)];
								REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
							end if
						end if;
						if correct or hasterm then
							break
						end if
					end do
				else
					randpick:=rand(1..total_perms);
					M:=add(l)+N;
					to maxIteration do
						degCoeffs:=UnrankMultiset(randpick(), Array(ul), freqs, N);
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
						Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
						S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
						S:= ifelse(type(S,list(algebraic)),S,NULL);
						if S<>NULL then
							if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets))  then
								hasterm:=has(Eq,map(rhs,Aindets));
								S:=NULL
							else
								S:=[seq(V[i]=S[i],i=1..M)];
								REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
							end if
						end if;
						if correct or hasterm then
							break
						end if
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


modFFixedOrdDegFunGuess2:= proc(Lf::algebraic,
			    degADE::posint,
			   degPoly::nonnegint,
				 Y::anyfunc(name),
				 N::nonnegint,
				 y::name,
				 x::name,
		      maxIteration::Or(posint,identical(infinity)),
			   modulus::posint,
		    inputConstants::set(name),
				$)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2026 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,total_perms,
		       nL::posint:=degree(Lf,x)+1,ADE::algebraic,polEq::algebraic,val,randpick,	
		       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),ul,tl,freqs,
		       correct::truefalse:=false,Arbconst::list,ADEcheck::algebraic,Meqs,beqs,
		       l::list(nonnegint),Ll::list(list),m::nonnegint,degCoeffs::list(nonnegint);
		       
		l:=GenMaxlistnumber(N,degPoly,nL-N);
		while l<> FAIL and not(correct) do
			ul:=sort([op({op(l)})]);
			tl:=Statistics:-Tally(l);
			freqs:=Array([seq(eval(val, tl), val = ul)], datatype = integer);
			total_perms:=MultinomialCount(freqs);
			if total_perms < maxIteration then
				Ll:=Iterator:-Permute(l,N);
				M:=add(l)+N;
				for degCoeffs in Ll do
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					polEq:=expand(eval(ADE,Y=Lf) mod modulus) mod modulus;
					Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; 
					Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
					S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
					S:= ifelse(type(S,list(algebraic)),S,NULL);
					if S<>NULL then
						if remove(v->v=0,S)=[] then
							S:=NULL
						else
							S:=[seq(V[i]=S[i],i=1..M)];
							ADEcheck, S, correct:=modpolcheckSol(S,ADE,Lf,nL,y,x,modulus)
						end if
					end if;
					if correct then
						break
					end if
				end do
			else
				randpick:=rand(1..total_perms);
				M:=add(l)+N;
				to maxIteration do
					degCoeffs:=UnrankMultiset(randpick(), Array(ul), freqs, N);
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					polEq:=expand(eval(ADE,Y=Lf) mod modulus) mod modulus;
					Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; #[seq(coeff(polEq,x,i),i=0..M-1)];
					Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
					S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
					S:= ifelse(type(S,list(algebraic)),S,NULL);
					if S<>NULL then
						if remove(v->v=0,S)=[] then
							S:=NULL
						else
							S:=[seq(V[i]=S[i],i=1..M)];
							ADEcheck, S, correct:=modpolcheckSol(S,ADE,Lf,nL,y,x,modulus)
						end if
					end if;
					if correct then
						break
					end if
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
		       maxIteration::Or(posint,identical(infinity)),
			    modulus::posint,
	             inputConstants::set(name),
				 $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
		       nL::posint:=numelems(Sinit),hasterm::truefalse:=false,ADE::algebraic,RE::algebraic,total_perms,	
		       Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),ul,tl,freqs,val,
		       correct::truefalse:=false,Arbconst::list,REcheck::algebraic,Aindets,Meqs,beqs,randpick,
		       l::list(nonnegint),Ll::list(list),m::nonnegint,degCoeffs::list(nonnegint);
		
		l:=GenMaxlistnumber(N,degPoly,max(nL-N,0));
		ul:=sort([op({op(l)})]);
		tl:=Statistics:-Tally(l);
		freqs:=Array([seq(eval(val, tl), val = ul)], datatype = integer);
		total_perms:=MultinomialCount(freqs);
		while l<> FAIL and not(correct) do
			if total_perms < maxIteration then
				Ll:=Iterator:-Permute(l,N);
				M:=add(l)+N;
				for degCoeffs in Ll do:
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
					Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
					S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
					S:= ifelse(type(S,list(algebraic)),S,NULL);
					if S<>NULL then
						if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets))  then
							hasterm:=has(Eq,map(rhs,Aindets));
							S:=NULL
						else
							S:=[seq(V[i]=S[i],i=1..M)];
							REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
						end if
					end if;
					if correct or hasterm then
						break
					end if
				end do
			else
				randpick:=rand(1..total_perms);
				degCoeffs:=copy(l);
				to maxIteration do
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
					Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
					S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
					S:= ifelse(type(S,list(algebraic)),S,NULL);
					if S<>NULL then
						if remove(v->v=0,S)=[] or has(S,map(rhs,Aindets))  then
							hasterm:=has(Eq,map(rhs,Aindets));
							S:=NULL
						else
							S:=[seq(V[i]=S[i],i=1..M)];
							REcheck, S, correct:=modcheckSol(S,RE,NegInd,Sinit,M,nL,a,n,modulus)
						end if
					end if;
					if correct or hasterm then
						break
					end if;
					degCoeffs:=UnrankMultiset(randpick(), Array(ul), freqs, N)
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