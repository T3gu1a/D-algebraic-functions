#Conversion of holonomic equation to simple ratrec equations

HoloToSimpleRatrec := proc(heq::Or(`=`,algebraic),
			   v::anyfunc(name),
			   {userbound::nonnegint:=0,
			   method::identical(GB,LA):=LA},
			   $)::`=`;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		if method=LA then
			return HoloToSimpleRatrecLA(heq,v)
		elif method = GB then
			return HoloToSimpleRatrecGB(heq,v,':-userbound'=userbound)
		end if
	end proc:
	
HoloToSimpleRatrecGB := proc(heq::Or(`=`,algebraic),v::anyfunc(name),{userbound::nonnegint:=0},$)::`=`;
		local P::algebraic, n::name, a::name, l::nonnegint, d::nonnegint, L::list,
		x::nothing, J, par::set, E::list:=[], j:=nonnegint, i::nonnegint:=0, A::list, R::`=`;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		P:=ifelse(type(heq,`=`),lhs(heq)-rhs(heq),heq);
		par:=indets(heq) minus {indets(v)};
		n:=op(v);
		a:=op(0,v);
		l:=HyperTypeSeq:-AlgebraHolonomicSeq:-REorder(P=0, a(n));
		L:=[subs([seq(a(n + j) = a[n + j], j = 0..l), n = x[n]], P)];
		d:=min(degree(L[1],x[n]),userbound);
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

HoloToSimpleRatrecLA := proc(heq::Or(`=`,algebraic),v::anyfunc(name),$)::`=`;
		local P::algebraic, n::name, a::name, l::nonnegint, d::nonnegint, L::list,
		x::nothing, E::list:=[], j:=nonnegint, i::nonnegint:=0, R::algebraic;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		P:=ifelse(type(heq,`=`),lhs(heq)-rhs(heq),heq);
		n:=op(v);
		a:=op(0,v);
		l:=HyperTypeSeq:-AlgebraHolonomicSeq:-REorder(P=0, a(n));
		R:=subs([seq(a(n+j) = a[n+j], j=0..l), n=x[n]], collect(P,n,'distributed'));
		d:=degree(R,x[n]);
		L:=[R,op(subs([seq(a(n+j) = a[n+j], j=1..l+d-1), n=x[n]],[seq(collect(LREtools:-shift(P,n,i),n),i=1..d-1)]))];
		L:=subs([seq(x[n]^i=x[n+i-1],i=2..d)],L);
		E:=SolveTools:-Linear(L,[seq(x[n+j],j=0..d)]);
		R:=subs([seq(a(n+j) = a[n+j], j=d..d+l),n=x[n]],collect(LREtools:-shift(P,n,d),n));
		R:=subs([seq(x[n]^i=x[n+i-1],i=1..d)],R);
		R:=subs(E,R);
		R:=a[n+d+l]=solve(R,a[n+l+d]);
		return subs([seq(a[n+j]=a(n+j), j=0..l+d)],R)
	end proc:
	
	



