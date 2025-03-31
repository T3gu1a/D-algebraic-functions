
OrderDegreeADE:= proc(ADE::`=`,
		        y::anyfunc(name),
		       $)::list(nonnegint,nonnegint);
		local x::name:=op(y),X::name:=op(0,y),r::nonnegint,d::nonnegint,j::posint,
		      P::polynom;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Returns [order, degree] of the input ADE";
		P:=lhs(ADE)-rhs(ADE);
		r:=PDEtools:-difforder(P,x);
		P:=subs([seq(diff(y,[x$j])=X[j],j=0..r)],P);
		d:=degree(P,{seq(X[j],j=0..r)});
		return [r,d]
	end proc: