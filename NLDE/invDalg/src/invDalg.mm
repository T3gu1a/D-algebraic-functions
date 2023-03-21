
invDalg:= proc(DE::`=`,
		y::anyfunc(name),
	        z::anyfunc(name),
	       $)::`=`;
		local t::name:=op(y),n::posint,fgx::set(`=`),R::algebraic,
		      f::nothing,x::nothing,g::nothing,j::nonnegint,DE1::`=`;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description "inverse a D-algebraic function from its differential        "
			    "equation in the dependent variable y(t)                     "
			    "INPUT: - a differential equation DE(y(t))                   "
			    "       - the dependent variable y(t)                        "
			    "	    - the dependent variable of the output z(t)          "
			    "OUPUT: a differential equations for the inverse f^{-1}(t)   ";
		DE1:=lhs(DE)-rhs(DE);
		n:=PDEtools:-difforder(DE1,t);
		#main computation: linear algebra (see ftogx)
		fgx:=subs(x=t,ftogx(f,g,x,n));
		#substitution to remove derivatives and t from DE1
		R:=subs(t=g[0],subs([seq(diff(y,[t$j])=f[j],j=0..n)],DE1));
		#the elimination of y
		R:=numer(normal(subs(fgx,R)));
		R:=collect(R,[seq(g[j],j=0..n)],'distributed');
		#Return the numerator of R after substitution
		R:=subs([seq(g[j]=diff(z,[t$j]),j=0..n)],R);
		return R=0
	end proc:
	
$include <NLDE/invDalg/SubProcedures/src/ftogx.mm>