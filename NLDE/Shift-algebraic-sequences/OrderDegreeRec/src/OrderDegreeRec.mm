
OrderDegreeRec:= proc(RE::`=`,
		       s::anyfunc(name),
		       $)::list(nonnegint,nonnegint);
		local n::name:=op(s),X::name:=op(0,s),r::nonnegint,d::nonnegint,j::posint,
		      P::polynom;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		description "Returns [order, degree] of the input recurrence equation";
		P:=numer(normal(lhs(RE) - rhs(RE)));
		r:=HyperTypeSeq:-AlgebraHolonomicSeq:-REorder(P=0,X(n));
		P:=subs([seq(X(n+j)=X[j],j=0..r)],P);
		d:=degree(P,[seq(X[j],j=0..r)]);
		return [r,d]
	end proc: