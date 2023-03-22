
OrderDegreeADE:= proc(ADE::`=`,
		        y::anyfunc(name),
		       $)::list(nonnegint,nonnegint);
		local x::name:=op(y),X::name:=op(0,y),r::nonnegint,d::nonnegint,j::posint,
		      P::polynom, L::list;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Returns [order, degree] of the input ADE";
		P:=lhs(ADE)-rhs(ADE);
		r:=PDEtools:-difforder(P,x);
		P:=subs([seq(diff(y,[x$j])=X[j],j=0..r)],P);
		#convert P to the list of its monomials
		#to avoid cancellation after the substitution (x=1)
		L:=subs(x=1,convert(P,list));
		d:=max(map(degree,L));
		return [r,d]
	end proc: