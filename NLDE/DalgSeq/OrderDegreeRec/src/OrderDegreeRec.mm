
OrderDegreeRec:= proc(RE::`=`,
		       s::anyfunc(name),
		      $)::list(nonnegint,nonnegint);
		local n::name:=op(s),X::name:=op(0,s),r::integer,r0::integer,d::nonnegint,j::posint,
		      P::polynom;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		description "Returns [order, degree] of the input recurrence equation";
		P:=numer(normal(lhs(RE) - rhs(RE)));
		(r,r0):=DalgSeq:-REorders(P=0,X(n));
		r:=r;
		P:=subs([seq(X(n+j)=X[j],j=r0..r)],P);
		d:=degree(P,[seq(X[j],j=r0..r)]);
		return [r-r0,d]
	end proc: