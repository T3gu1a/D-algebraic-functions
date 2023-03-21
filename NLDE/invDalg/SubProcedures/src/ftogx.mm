
ftogx:= proc(f::name,g::name,x::name,N::posint,$)::set(`=`); option remember;
	   local j::nonnegint,Lderiv::list(algebraic);
	   option `Copyright (c) 2023 Bertrand Teguia T.`;
	   description  "subprocedure of invDalg for expressing the derivatives of f    "
			"in terms of those of g and y                                   "
			"INPUT:  - the name for g                                       "
			"        - the name for t (always use with x for remembrance)   "
			"        - the integer N (representing n from invDalg)          "
			"OUPUT: the list with the f[j] in terms of g[j] and x           "
			"       j=0..N                                                  ";
		#this procedure uses the method which solve the triangular linear system
		#of dimention N. It turns out to be more efficient in general compare to
		#the recursive approach, even though the latter is more effective 
		#for remembrance.
		Lderiv:= [seq(diff(x,[x$j])=diff(f(g(x)),[x$j]),j=0..N)];
		Lderiv:= map(r->subs([seq(diff(g(x),[x$j])=g[j],j=0..N)],r),Lderiv);
		Lderiv:= map(r->subs(g[0]=x,r),Lderiv);
		Lderiv:= map(r->convert(r,diff),Lderiv);
		Lderiv:= map(r->subs([seq(diff(f(x),[x$j])=f[j],j=0..N)],r),Lderiv);
		return SolveTools:-Linear(Lderiv,[seq(f[j],j=0..N)])
	end proc:

