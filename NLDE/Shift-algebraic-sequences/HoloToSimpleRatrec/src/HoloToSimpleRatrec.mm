#Conversion of holonomic equation to simple ratrec equations

HoloToSimpleRatrec := proc(heq::Or(`=`,algebraic),v::anyfunc(name),{userbound::nonnegint:=0},$)::`=`;
		local P::algebraic, n::name, a::name, l::nonnegint, d::nonnegint, L::list,
		x::nothing, J, par::set, E::list:=[], j:=nonnegint, i::nonnegint:=0, A::list, R::`=`;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		P:=ifelse(type(heq,`=`),lhs(heq)-rhs(heq),heq);
		par:=indets(heq) minus {indets(v)};
		n:=op(v);
		a:=op(0,v);
		l:=HyperTypeSeq:-AlgebraHolonomicSeq:-REorder(P=0, a(n));
		L:=[subs([seq(a(n + j) = a[n + j], j = 0..l), n = x[n]], P)];
		d:=max(degree(L[1],x[n]),userbound);
		while E = [] and i<=d do:
			i:=i+1;
			L:=[subs([seq(a(n + j) = a[n + j], j = i .. l+i), n = x[n]], normal(LREtools:-shift(P,n,i))),op(L)];
			J:=PolynomialIdeals:-PolynomialIdeal(L,'parameters'=par);
			J:=PolynomialIdeals:-EliminationIdeal(J, {seq(a[n + j], j = 0 .. l + i)});
			J:=select(type, convert(J, list), polynom);
			E:=select(r -> degree(r, a[n + l + i]) = 1, J);
		end do;
		if E<>[] then
			A:=[seq(a[n + j], j = 0 .. l + i)];
			E:=sort(E,(s,t)->degree(s,A)<=degree(t,A));
			R:=a[n+l+i]=solve(E[1],a[n+l+i]);
			return subs([seq(a[n+j]=a(n+j), j=0..l+i)],R)
		end if;
		return heq
	end proc:
	
	



