NLDE:= module()

option `Copyright (c) 2022 Bertrand Teguia Tabuguia, Max Planck Institute for MiS, Leipzig`, package;

export unaryDalg, SystoMinDiffPoly, composeDalg, arithmeticDalg;

local buildsystem, mergesystem, ftogh, subsgfurther;


buildsystem:= proc(DE::`=`,
		    y::anyfunc(name),
	            x::name,
		   $)::list(`=`);
		local  t::name, d::posint, SubL::list, PolDE::polynom, j::nonnegint;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description     "Build a dynamical system (or model) from a differential equation. "
				"If the differential equation is not LEF, its derivatives is used. "
				"INPUT: -A differential equation DE,                               "
				"       -its dependent variable like y(t)                          "
				"	-a name x for the variable of the system                   "
				"OUPUT: A list of two lists:                                       "
				"       - the list of derivatives of the variables of the system   "
				"         in in terms of these variables                           " 
				"	- the variables of the system                              ";
		t:=op(y);
		d:=PDEtools:-difforder(DE,t);
		#variables of substitution for the model, the input x with indices
		SubL:=[seq(diff(y,[t$j])=x[j],j=0..d)];
		PolDE:=subs(SubL,lhs(DE));
		#the differential equation is not LEF
		if degree(PolDE,x[d])>1 then
			d:=d+1;
			SubL:=[op(SubL),diff(y,t$d)=x[d]];
			PolDE:=subs(SubL,lhs(diff(DE,t)));
			return [[seq(x[j],j=1..(d-1)),solve(PolDE,x[d])],[seq(x[j],j=0..(d-1))]]
		else
			return [[seq(x[j],j=1..(d-1)),solve(PolDE,x[d])],[seq(x[j],j=0..(d-1))]]
		end if		
	end proc:

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


SystoMinDiffPoly:= proc(f::list(algebraic),
                        g::algebraic,
                        X::Or(list,set),
                        z::anyfunc(name),$)::algebraic;
		local F,G,q1,q2,Q,Svars,J1,J2,J,n,Xt,t,DE,Sub:=[],j,k,y,alpha;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description     "Compute the minimal order non-linear DE of y=g(p,x)             "
				"from the system {x'=f(p,x),y=g(p,x)}, for any parametric        "
				"vector p, and a variable x=(x_1,...x_n), where                  "
				"f and g are rational functions in x_1,...,x_n                   "
				"INPUT:  - the list of derivatives of the variables of the system"
				"          in in terms of them.                                  "
				"        - the rational expession representing g                 "
				"	 - the list of variables of the system                   "
				"        - the dependent variable (like y(t)) for the            "
				"          output differential equation                          "
				"OUPUT: a differential equations for g                           ";
		t:=op(1,z);
		y:=op(0,z);
		alpha:=indets([f,g]) minus {op(X)};
		n:=numelems(X);
		F:=normal(f);
		q1:=mul(map(denom,F));
		G:=normal(g);
		q2:=denom(G);
		#least common multiple of the denominators of the system
		Q:=lcm(q1,q2);
		#to differentiate, the variables should be functions of the 
		#independent variable t
		Xt:=map(x->x=x(t),X);
		Q:=subs(Xt,Q);
		F:=subs(Xt,F);
		G:=subs(Xt,G);
		Xt:=map(rhs,Xt);
		J1:=[seq(Q*diff(Xt[j],t)-normal(Q*F[j]),j=1..n)];
		#differentiating n-1 times the polynomials Q*x'-Q*f
		for j to n do:
			J1:=[op(J1),seq(diff(J1[j],t$k),k=1..(n-1))]
		end do;
		J2:=[seq(diff(Q*y(t)-normal(Q*G),[t$j]),j=0..n)];
		J:=[op(J1),op(J2)];
		#differentiating n times the polynomials Q*y - Q*g
		for j to n do:
			Sub:=[op(Sub),seq(diff(Xt[j],[t$k])=x[j,k],k=0..n)]
		end do;
		Sub:={op(Sub),seq(diff(y(t),[t$j])=y[j],j=0..n)};
		#elimination and saturation with Groebner bases
		J:=PolynomialIdeals:-PolynomialIdeal(subs(Sub,J),parameters=alpha);
		J:=PolynomialIdeals:-Saturate(J,subs(Sub,Q));
		J:=PolynomialIdeals:-EliminationIdeal(J,select(has,map(rhs,Sub),y));
		DE:=op(1,J);
		DE:=collect(DE,[seq(y[j],j=0..n)],'distributed');
		return subs(select(has,map(e->rhs(e)=lhs(e),Sub),y),DE)=0
	end proc:

subsgfurther :=proc(gm1::algebraic,g::name,t::name,m::posint,n::posint,$)::list;
		local k::posint,j::nonnegint,Subdiff::list,rSubdiff::list,eqg,Sub::list;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description     "subprocedure of composeDalg for substituting higher derivatives"
		                "of g in the rational relation R                                "
				"INPUT:  - an algebraic relation representing the link from     "
				"          the differential equation                            "
				"        - the name for g                                       "
				"        - the name for the independent variable t              "
				"        - the order of the equation of g, m                    "
				"        - the order n >= m of the other equation               "
				"OUPUT: the list required to make the substitution              ";
		k:=0;
		eqg:=gm1;
		#the first substition
		Sub:=[g[m]=eqg];
		if k<n-k then
			Subdiff:=[seq(g[j]=diff(g(t),[t$j]),j=0..m-1)];
			rSubdiff:=[op(map(r->rhs(r)=lhs(r),Subdiff)),diff(g(t),t$m)=g[m]];
			#the other substitutions are obtained by
			#differentiating and substituting the first substitution
			while k<n-m do
				eqg:=subs(Subdiff,eqg);
				eqg:=diff(eqg,t);
				eqg:=subs(rSubdiff,eqg);
				eqg:=normal(subs(Sub[1],eqg));
				k:=k+1;
				Sub:=[op(Sub),g[m+k]=eqg]
			end do
		end if;
		return Sub
	end proc:

ftogh:= proc(f::name,g::name,h::name,x::name,N::posint,$)::set(`=`); option remember;
	   local j::nonnegint,Lderiv::list(algebraic);
	   option `Copyright (c) 2022 Bertrand Teguia T.`;
	   description  "subprocedure of composeDalg for expressing the derivatives of f"
		        "in terms of those of h and g                                   "
			"INPUT:  - the name for f                                       "
			"        - the name for g                                       "
			"        - the name for h                                       "
			"        - the name for t (always use with x for remembrance)   "
			"        - the integer N (representing n from composeDalg)      "
			"OUPUT: the list with the f[j] in terms of g[j] and h[j]        "
			"       j=0..N                                                  ";
	   #this procedure uses the method which solve the triangular linear system
	   #of dimention N. It turns out to be more efficient in general compare to
	   #the recursive approach, even though the latter is more effective 
	   #for remembrance.
	   Lderiv:= [seq(h[j]=diff(f(g(x)),[x$j]),j=0..N)];
	   Lderiv:= map(r->subs([seq(diff(g(x),[x$j])=g[j],j=0..N)],r),Lderiv);
	   Lderiv:= map(r->subs(g[0]=x,r),Lderiv);
	   Lderiv:= map(r->convert(r,diff),Lderiv);
	   Lderiv:= map(r->subs([seq(diff(f(x),[x$j])=f[j],j=0..N)],r),Lderiv);
	   return SolveTools:-Linear(Lderiv,[seq(f[j],j=0..N)])
	end proc:

composeDalg:= proc(L::[`=`,`=`],
	           V::[anyfunc(name),anyfunc(name)],
	           z::anyfunc(name),
	          $)::`=`;
		local t::name:=op(z),n::posint,fgh::set(`=`),R,f::nothing,Sys::list,x::nothing,
		      Sysh::list,g::nothing,h::nothing,j::nonnegint,m::posint,Subg::list,
		      DE1::`=`,DE2::`=`;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Compose two D-algebraic functions from their differential   "
			    "equations given in the same order with the second argument, "
			    "representing the list of dependent variables.               "
			    "INPUT: - a list of two differential equations [DE(f),DE(g)] "
			    "       - a list of their dependent variables  [f(t),g(t)]   "
			    "	    - the dependent variable of the output h(t)          "
			    "OUPUT: a differential equations for the composition f(g)    "
			    "       of the second g by the first f.                      ";
		DE1:=lhs(L[1])-rhs(L[1]);
		DE2:=lhs(L[2])-rhs(L[2])=0;
		n:=PDEtools:-difforder(DE1,t);
		#substitution to remove derivatives and t from the first differential equation (DE)
		R:=subs(t=g[0],subs([seq(diff(V[1],[t$j])=f[j],j=0..n)],DE1));
		#in case on non-LEF, one derivation is needed
		if degree(R,f[n])>1 then
			n:=n+1;
			R:=subs(t=g[0],subs([seq(diff(V[1],[t$j])=f[j],j=0..n)],diff(DE1,t)))
		end if;
		#build the system with the second DE
		Sys:=buildsystem(DE2,V[2],g);
		m:=numelems(Sys[1]);
		#expressing derivatives of f in terms of those of h and g
		fgh:=ftogh(f,g,h,x,n);
		if m<=n then
			#build a list of further substitutions in R if m<=n
			Subg:=subsgfurther(Sys[1][m],g,t,m,n);
			fgh:=subs(Subg,fgh)
		end if;
		#the full substitution (f is already eliminated here)
		R:=normal(subs(fgh,R));
		#the system for h is easily obtained from R
		Sysh:=[[seq(h[j],j=1..(n-1)),solve(R,h[n])],[seq(h[j],j=0..(n-1))]];
		#use SystoMinDiffPoly to return the desired output
		return SystoMinDiffPoly([op(Sys[1]),op(Sysh[1])],h[0],[op(Sys[2]),op(Sysh[2])],z)
	end proc:			
			
unaryDalg:= proc(DE::`=`,
		  y::anyfunc(name),
		  z::name=ratpoly,
		 $)::`=`;
		local t::name:=op(y),var::name:=op(0,y),dvar::name=lhs(z),
		      r::ratpoly:=rhs(z),eq::algebraic,j::nonnegint,Sys::list,x::nothing;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Compute a differential equation for                                   "
		            "a rational expression of a D-algebraic function from a                "
			    "differential equation that it satisfies.                              "
			    "INPUT: - a differential equation                                      "
			    "       - its dependent variable, say f(t)                             "
			    "	    - an equation h=r(f) (a rational expression in f)              "
			    "	      h is the name for the dependent variable in the output.      "
			    "OUPUT: a differential equation satisfied by r(f)                      ";
		var:=op(0,y);
		dvar:=lhs(z);
		r:=normal(rhs(z));
		#Simple case: y appears linearly in r(y)
		if degree(numer(r),var)<=1 and degree(denom(r),var)<=1 then
			eq:=subs(dvar=dvar(t),solve(dvar-r,var));
			eq:=eval(lhs(DE) - rhs(DE), y=eq);
			eq:=numer(normal(eq));
			return collect(eq,[seq(diff(dvar(t),[t$j]),
			       j=0..PDEtools:-difforder(DE,t))],'distributed')=0
		end if;
		#General case
		#build the system using buildsystem
		Sys:=buildsystem(lhs(DE) - rhs(DE)=0,y,x);
		#use SystoMinDiffPoly to return the desired output
		return SystoMinDiffPoly(Sys[1],subs(var=Sys[2][1],r),Sys[2],dvar(t))
	end proc:
	
arithmeticDalg:=proc(L::list(`=`),
		     V::list(anyfunc(name)),
		     z::name=ratpoly,
		    $)::`=`;
		local t:=op(1,V[1]),DEs::list(`=`),Sys::list,j::posint,subV::list;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Compute a differential equaton for a rational expression of        "
		            "D-algebraic functions from differential equations given in the     "
			    "same order with the second argument, representing the list of      "
			    "dependent variables.                                               "
			    "One may supply more than two differential equations.               "
			    "INPUT: - a list of differential equations [DE(f1(t)),...,DE(fn(t))]"
			    "       - a list of their dependent variables [f1(t),...,fn(t)]     "
			    "	    - an equation z=r(f1,...,fn) where z is the name of the     "
			    "         dependent variable for the output, and r(f1,...,fn) is    "
			    "         is a rational expression in f1,...,fn.                    "
			    "OUPUT: a differential equation satisfied by the rational expresion " 
	                    "       r(f1,f2,...,fn)                           		        ";
		if numelems(L)=1 then
			return L
		end if;
		DEs:=map(r->lhs(r) - rhs(r)=0,L);
		#build the systems and merge them using mergesystem
		Sys:=mergesystem(DEs,V);
		#prepare the list for the change of variables 
		#in r according to Sys
		subV:=[seq(op(0,V[j])=Sys[2][j],j=1..numelems(V))];
		#use SystoMinDiffPoly to return the desired output
		return SystoMinDiffPoly(Sys[1],subs(subV,rhs(z)),Sys[3],lhs(z)(t))
	end proc:

end module:

savelib('NLDE',"C:/Users/bertr/maple/toolbox/personal/lib/NLDE.mla"):





