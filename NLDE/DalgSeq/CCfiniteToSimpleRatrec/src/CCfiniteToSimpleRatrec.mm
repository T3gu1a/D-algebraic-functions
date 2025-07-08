#efficient algorithm for converting CC-finite equations into Dalgebraic rational recursions

CCfiniteToSimpleRatrec:=proc(DE::`=`,S::anyfunc(name),Leq::list(`=`),Lvars::list(anyfunc(name)),$)::`=`;
			local p::`=`, C::list,r::nonnegint,j::nonnegint,n::name,s::name,k::nonnegint,
			      l::nonnegint,R::list(nonnegint),Elim::list(`=`):=[], eqelim::`=`;
			option `Copyright (c) 2025 Bertrand Teguia T.`;
			description "efficient algorithm (compared to the classical Groebner bases method) for converting"
				     "CC-finite equations into D-algebraic equations, which moreover are rationalizing";
			n:=op(S);
			s:=op(0,S);
			r:=REorders(DE,s(n))[1];
			p:=DE;
			l:=numelems(Lvars);
			R:=[seq(REorders(Leq[j],Lvars[j])[1],j=1..l)];
			C:=[seq(subs(n=n+R[j],Lvars[j])=solve(Leq[j],subs(n=n+R[j],Lvars[j])),j=1..l)];
			#decrement the order by 1 for the last c-finite equation
			R[-1]:=R[-1]-1;
			for j to l do
				for k from 0 to R[j]-1 do
					p:=subs(C,p);
					if has(p,subs(n=n+k,Lvars[j])) then
						eqelim:=subs(n=n+k,Lvars[j])=solve(p,subs(n=n+k,Lvars[j]));
						Elim:=subs(eqelim,Elim);
						Elim:= [eqelim,op(Elim)]
					end if;
					C:=subs(Elim,C);
					p:=subs(Elim,subs(n=n+1,p))
				end do
			end do;
			r:=r+add(R);
			return s(n+r)=solve(subs(C,p),s(n+r))
		end proc: