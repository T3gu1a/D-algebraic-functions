

ComputDegkDE := proc(f::algebraic,z::name,k::posint,degkNmax::posint,Subdiff::list(`=`),startfromord::posint,$)
		local a::nothing, A, N:=startfromord, Eq, n, S, Sumds, nisol, j, tmp, rmS, Eqs, s, factSumds, polfact, Coef:=[], i,
		      updateSub::list;		
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Computational part of DegreekDE";
		while REorders(deltakshift(a(z),z,k,N),a(z))[1]<=degkNmax and Coef=[] do
			A:=[seq(a[i],i=0..N-1)];
			Eq:=deltakshift(f,z,k,N)+add(A[i+1]*deltakshift(f,z,k,i),i=0..N-1);
			n:=REorders(deltakshift(a(z),z,k,N),a(z))[1];
			updateSub:=[op(Subdiff),seq(op(LREtools:-shift(Subdiff, z, j)), j = 1..n)];
			to n do
				Eq:=evala(eval(Eq,updateSub))
			end do;
			Eq:=expand(numer(normal(Eq)));
			if Eq=0 then
				return [1], 1
			end if;
			S:=[op(Eq)];
			Sumds:=[];
			nisol:=true;
			while S<>[] and nisol do
				tmp:=S[1];
				rmS:=[tmp];
				if numelems(S)>1 then
					for j from 2 to numelems(S) do
						if type(normal(S[j]/S[1]), ratpoly(anything,z)) then
							tmp:=normal(tmp+S[j]);
							rmS:=[op(rmS),S[j]]
						end if
					end do
				end if;
				Sumds:=[op(Sumds),tmp];
				nisol:=has(tmp,A);
				S:=remove(member,S,rmS)
			end do;
			if nisol then 
				Sumds:=map(r->numer(factor(r)),Sumds);
				Eqs:=[];
				for s in Sumds do
					if type(s,polynom(anything,z)) then
						if has(s,A) then
							Eqs:=[op(Eqs),s]
						end if
					else
						factSumds:=map(t->exp(t),[op(simplify(ln(s),ln,'symbolic'))]);
						polfact:=select(type,factSumds,polynom(anything,z));
						polfact:=select(has,polfact,A);
						Eqs:=[op(Eqs),op(polfact)]
					end if
				end do;
				#Coef:= solve(Eqs,A);
				Coef:=[SolveTools:-Linear(Eqs,A)];
				if Coef<>[] then
					Coef:=map(rhs,convert(Coef[1],list));
					A:=map(r->r=1,A);
					Coef:=factor(subs(A,Coef))
				else
					N:=N+1
				end if
			else 
				N:=N+1
			end if
		end do;
		Coef, N
	end proc: