
mergesystem:= proc(L::list(`=`),
		   V::list(anyfunc(name)),
		  $)::`=`;
		local l::posint:=numelems(L), j::posint, Sys::list, vars::list, deriv::list, 
		      n::posint, x::nothing, X::list, i::posint, Ind::list;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description     "Merge the dynamical systems of a list of differential equations. "
				"INPUT: -a list of differential equations,                        "
				"       -their dependent variable like y(t)                       "
				"OUPUT: A list of three lists:                                    "
				"       - the list of derivatives with the variables of           " 
				"         the new system                                          "
				"	- the list of variables representing the solutions        "
				"	  of the input equations                                  "
				"	- the variables of the system                             ";
		Sys:=[seq(buildsystem(L[j],V[j],cat(x,j)),j=1..l)];
		vars:=map(r->op(r[2]),Sys);
		deriv:=map(r->op(r[1]),Sys);
		
		#deriv:=[1,op(deriv)]; #constant coefficients
		
		n:=numelems(vars);
		X:=[seq(vars[j]=x[j],j=1..n)];
		
		#X:=[t=x[0],seq(vars[j]=x[j],j=1..n)]; #constant coefficients
		
		#indices of the variables representing the solutions of the input DEs
		Ind:=[seq(1+add(numelems(Sys[i][2]),i=1..(j-1)),j=1..l)];
		return [subs(X,deriv),map(r->x[r],Ind),map(rhs,X)]
	end proc:

$include <NLDE/CommonInternalProcedures/src/buildsystem.mm>