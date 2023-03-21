NLDE:= module()

option `Copyright (c) 2022 Bertrand Teguia Tabuguia, Max Planck Institute for MiS, Leipzig`, package;

export unaryDalg, diffDalg, invDalg, SystoMinDiffPoly, composeDalg, arithmeticDalg, AnsatzDalg, DDfiniteToDalg;

local buildsystem, mergesystem, ftogh, subsgfurther, ftogx, NLDE_nlho;

buildsystem:= proc(DE::`=`,
		    y::anyfunc(name),
		    x::name,
		   $)::list(`=`);
		local  t::name, d::posint, SubL::list, PolDE::polynom, j::nonnegint;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description     "Build a dynamical system (or model) from a differential equation.    "
				"If the differential equation is not l.h.o, its derivatives is used.  "
				"INPUT: -A differential equation DE,                                  "
				"       -its dependent variable like y(t)                             "
				"	-a name x for the variable of the system                      "
				"OUPUT: A list of two lists:                                          "
				"       - the list of derivatives of the variables of the system      "
				"         in in terms of these variables                              " 
				"	- the variables of the system                                 ";
		t:=op(y);
		d:=PDEtools:-difforder(DE,t);
		#variables of substitution for the model, the input x with indices
		SubL:=[seq(diff(y,[t$j])=x[j],j=0..d)];
		PolDE:=subs(SubL,lhs(DE));
		#the differential equation is not l.h.o
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
                        z::anyfunc(name),
			{ordering::identical(plex,lexdeg):=plex},
			$)::algebraic;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description     "Compute the minimal order non-linear DE of y=g(p,x)             "
				"from the system {x'=f(p,x),y=g(p,x)}, for any parametric        "
				"vector p, and a variable x=(x_1,...x_n), where                  "
				"f and g are rational functions in x_1,...,x_n                   "
				"INPUT:  - the list of derivatives of the variables of the system"
				"          in in terms of them.                                  "
				"        - the rational expession representing g                 "
				"	 - the list of variables of the system                	 "
				"        - the dependent variable (like y(t)) for the            "
				"          output differential equation                          "
				"OUPUT: a differential equations for g                           ";
		local F,G,q1,q2,Q,Svars,J1,J2,J,n,Xt,t,DE,Sub:=[],allvars,yvars,ord,j,k,y,alpha;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
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
		#differentiating n times the polynomials Q*y - Q*g
		J2:=[seq(diff(Q*y(t)-normal(Q*G),[t$j]),j=0..n)];
		J:=[op(J1),op(J2)];
		#build the list of substitution to see derivatives as variables
		for j to n do:
			Sub:=[op(Sub),seq(diff(Xt[j],[t$k])=x[j,k],k=0..n)]
		end do;
		if ordering=plex then
			#elimination and saturation with Groebner bases
			#w.r.t. pure lex monomial ordering
			Sub:=[op(Sub),seq(diff(y(t),[t$j])=y[j],j=0..n)];
			allvars:=ListTools:-Reverse(map(rhs,Sub));
			yvars:=select(has,allvars,y);
			allvars:=allvars[numelems(yvars)+1..-1];
			J:=PolynomialIdeals:-PolynomialIdeal(subs(Sub,J),parameters=alpha);
			J:=PolynomialIdeals:-Saturate(J,subs(Sub,Q));
			J:=Groebner:-Basis(J,plex(op(allvars),op(yvars)));
			J:=remove(has,J,allvars)
		else
			#elimination and saturation with Groebner bases
			#w.r.t. lexdeg elimination ordering
			Sub:={op(Sub),seq(diff(y(t),[t$j])=y[j],j=0..n)};
			J:=PolynomialIdeals:-PolynomialIdeal(subs(Sub,J),parameters=alpha);
			J:=PolynomialIdeals:-Saturate(J,subs(Sub,Q));
			yvars:=select(has,map(rhs,Sub),y);
			J:=PolynomialIdeals:-EliminationIdeal(J,yvars);
			J:=select(type,convert(J,list),polynom)
		end if;
		#Taking a diff polynomial of minimal total degree
		# among those of the minimal order
		J:=map(de->collect(de,[seq(y[j],j=0..n)],'distributed'),J);
		Sub:=select(has,map(e->rhs(e)=lhs(e),Sub),y);
		J:=map(de->subs(Sub,de),J);
		#order
		ord:=min(map(de->PDEtools:-difforder(de,t),J));
		DE:=select(de->PDEtools:-difforder(de,t)=ord,J);
		Sub:=map(e->rhs(e)=lhs(e),Sub);
		DE:=map(de->subs(Sub,de),DE);
		#degree
		DE:=sort(DE,(a,b)->degree(a,yvars)<=degree(b,yvars));
		DE:=DE[1];
		Sub:=map(e->rhs(e)=lhs(e),Sub);
		return subs(Sub,DE)=0
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
		   {ordering::identical(plex,lexdeg):=plex},
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
		#in case of non-l.h.o, one derivation is needed
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
		return SystoMinDiffPoly([op(Sys[1]),op(Sysh[1])],h[0],[op(Sys[2]),op(Sysh[2])],z,':-ordering'=ordering)
	end proc:			
			
unaryDalg:= proc(DE::`=`,
		y::anyfunc(name),
		z::name=ratpoly,
		{ordering::identical(plex,lexdeg):=plex},
		$)::`=`;
		local t::name:=op(y),var::name:=op(0,y),dvar::name:=lhs(z),
		      r::ratpoly:=rhs(z),Sys::list,x::nothing;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Compute a differential equation for                                   "
		            "a rational expression of a D-algebraic function from a                "
			    "differential equation that it satisfies.                              "
			    "INPUT: - a differential equation                                      "
			    "       - its dependent variable, say f(t)                             "
			    "	    - an equation h=r(f) (a rational expression in f)              "
			    "	      h is the name for the dependent variable in the output.      "
			    "OUPUT: a differential equation satisfied by r(f)                      ";
		#build the system using buildsystem
		Sys:=NLDE_nlho:-buildsystem(lhs(DE) - rhs(DE)=0,y,x);
		#use SystoMinDiffPoly to return the desired output
		return NLDE_nlho:-SystoMinDiffPoly(Sys[1],subs(var=Sys[2][1][1],r),
		                          Sys[2],dvar(t),':-ordering'=ordering)
	end proc:
	
arithmeticDalg:=proc(L::list(`=`),
		     V::list(anyfunc(name)),
		     z::name=ratpoly,
		    {ordering::identical(plex,lexdeg):=plex,
		    lho::truefalse:=true,
		    lhoplex::truefalse:=false},
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
		Sys:=ifelse(lho,mergesystem(DEs,V),NLDE_nlho:-mergesystem(DEs,V));
		#prepare the list for the change of variables 
		#in r according to Sys
		subV:=[seq(op(0,V[j])=Sys[2][j],j=1..numelems(V))];
		#use SystoMinDiffPoly to return the desired output
		if lho then
			return SystoMinDiffPoly(Sys[1],subs(subV,rhs(z)),Sys[3],lhs(z)(t),':-ordering'=ordering)
		else
			if lhoplex then 
				return NLDE_nlho:-SystoMinDiffPoly(Sys[1],subs(subV,rhs(z)),Sys[3],lhs(z)(t),':-ordering'=plex)
			else
				return NLDE_nlho:-SystoMinDiffPoly(Sys[1],subs(subV,rhs(z)),Sys[3],lhs(z)(t),':-ordering'=lexdeg)
			end if
		end if
	end proc:
	
diffDalg :=proc(DE::`=`,
		 y::anyfunc(name),
		 n::posint:=1,
		{ordering::identical(plex,lexdeg):=plex},
		$)::`=`;
		local t::name:=op(y),var::name:=op(0,y),
		      d::nonnegint,j::nonnegint,p,q,R,V;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description "Compute a differential equation for                                   "
		            "the derivative of a D-algebraic function from a                       "
			    "differential equation that it satisfies.                              "
			    "INPUT: - a differential equation                                      "
			    "       - its dependent variable, say f(t)                             "
			    "OUPUT: a differential equation satisfied by diff(f(t),t)              ";
		if n=1 then	    
			p:=lhs(DE)-rhs(DE);
			d:=PDEtools:-difforder(p,t);
			q:=diff(p,t);
			p:=subs([seq(diff(var(t),[t$j])=var[j],j=0..d)],p);
			q:=subs([seq(diff(var(t),[t$j])=var[j],j=0..d+1)],q);
			R:=resultant(p,q,var[0]);
			if type(R,`+`) then
				return subs([seq(var[j]=diff(var(t),[t$(j-1)]),j=1..d+1)],R)=0
			elif type(R,`*`) then
				R:=[select(has,R,var[d+1])];
				V:=[seq(var[j],j=0..d+1)];
				R:=sort(R,(a,b)->degree(a,V)<=degree(b,V));
				return subs([seq(var[j]=diff(var(t),[t$(j-1)]),j=1..d+1)],R[1])=0
			else
				return subs([seq(var[j]=diff(var(t),[t$(j-1)]),j=1..d+1)],R)=0
			end if
		else
			return diffDalg(diffDalg(DE,y,n-1),y)
		end if
	end proc:

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
	
	DDfiniteToDalg:= proc(DE::`=`,y::anyfunc(name),
			     Leq::list(`=`),
			   Lvars::list(anyfunc(name)),
		       {ordering::identical(plex,lexdeg):=plex},
			      $)::`=`;
		local DEs, V, Sys, j, subV;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description  "Convert a DD-finite equation into an ADE    "
			     "INPUT: - The DD-finite ODE DE               "
			     "       - Its dependent variable y(x)        "
			     "       - The list of holonomic DEs          "
			     "       - Their dependent variables in Lvars "
			     "OUTPUT: An ADE fulfills by the solutions of "
			     "        the given DD-finite ODE             ";
		DEs:=map(eq->lhs(eq) - rhs(eq)=0,[DE,op(Leq)]);
		V:=[y,op(Lvars)];
		Sys:=mergesystem(DEs,V);
		subV:=[seq(op(0,Lvars[j])=Sys[2][j+1],j=1..numelems(Lvars))];
		return SystoMinDiffPoly(subs(subV,Sys[1]),Sys[2][1],Sys[3],y,':-ordering'=ordering)
	end proc:
	
	#submodule used to overcome the non-l.h.o situation in the unary case
	#This generalizes to the general computation (proof to be done)
	#However, generally the minimal ADE in the non-l.h.o case is not "good-looking"
	#Hence the reason for not using it in general.
	NLDE_nlho:= module()

	option `Copyright (c) 2023 Bertrand Teguia Tabuguia, Max Planck Institute for MiS, Leipzig`, package;

	export buildsystem, mergesystem, SystoMinDiffPoly;


	buildsystem:= proc(DE::`=`,
			   y::anyfunc(name),
			   x::name,
			   $)::list(`=`);
			local  t::name, r::posint, SubL::list, PolDE::polynom, d::posint, j::nonnegint;
			option `Copyright (c) 2023 Bertrand Teguia T.`;
			description     "The non-lho anologue of buildsystem";
			t:=op(y);
			r:=PDEtools:-difforder(DE,t);
			#variables of substitution for the model, the input x with indices
			SubL:=[seq(diff(y,[t$j])=x[j],j=0..r)];
			PolDE:=subs(SubL,lhs(DE));
			d:=degree(PolDE,x[r]);
			#the differential equation is not l.h.o
			if d>1 then
				PolDE:=subs(x[r]^d=x[r],PolDE);
				return [[seq(x[j],j=1..(r-1)),solve(PolDE,x[r])],[seq([x[j],1],j=0..(r-2)),[x[r-1],d]]]
			else
				return [[seq(x[j],j=1..(r-1)),solve(PolDE,x[r])],[seq([x[j],1],j=0..(r-1))]]
			end if	
		end proc:

	mergesystem:= proc(L::list(`=`),
			   V::list(anyfunc(name)),
			   $)::`=`;
			local l::posint:=numelems(L), j::posint, Sys::list, vars::list, deriv::list, 
			      n::posint, x::nothing, X::list, i::posint, Ind::list;
			option `Copyright (c) 2023 Bertrand Teguia T.`;
			description     "The non-lho analogue of mergesystem";
			Sys:=[seq(buildsystem(L[j],V[j],cat(x,j)),j=1..l)];
			vars:=map(r->op(r[2]),Sys);
			deriv:=map(r->op(r[1]),Sys);
			n:=numelems(vars);
			X:=[seq(vars[j][1]=x[j],j=1..n)];
			#indices of the variables representing the solutions of the input DEs
			Ind:=[seq(1+add(numelems(Sys[i][2]),i=1..(j-1)),j=1..l)];
			return [subs(X,deriv),map(r->x[r],Ind),subs(X,vars)]
		end proc:


	SystoMinDiffPoly:= proc(f::list(algebraic),
				g::algebraic,
				X::Or(list,set),
				z::anyfunc(name),
				{ordering::identical(plex,lexdeg):=plex},
				$)::algebraic;
			option `Copyright (c) 2023 Bertrand Teguia T.`;
			description     "The non-lho analogue of NLDE:-SystoMinDiffPoly";
			local F,G,q1,q2,Q,Svars,J1,J2,J,n,Xt,nlho_pow,t,DE,Sub:=[],allvars,yvars,ord,j,k,y,alpha;
			t:=op(1,z);
			y:=op(0,z);
			alpha:=indets([f,g]) minus {op(map(x-x[1],X))};
			n:=numelems(X);
			F:=normal(f);
			q1:=mul(map(denom,F));
			G:=normal(g);
			q2:=denom(G);
			#least common multiple of the denominators of the system
			Q:=lcm(q1,q2);
			#to differentiate, the variables should be functions of the 
			#independent variable t
			Xt:=map(x->x[1]=x[1](t),X);
			Q:=subs(Xt,Q);
			F:=subs(Xt,F);
			G:=subs(Xt,G);
			Xt:=map(rhs,Xt);
			nlho_pow:=map(x->x[2],X);
			J1:=[seq(Q*diff(Xt[j],t)^nlho_pow[j]-normal(Q*F[j]),j=1..n)];
			#differentiating n-1 times the polynomials Q*x'-Q*f
			for j to n do:
				J1:=[op(J1),seq(diff(J1[j],t$k),k=1..(n-1))]
			end do;
			#differentiating n times the polynomials Q*y - Q*g
			J2:=[seq(diff(Q*y(t)-normal(Q*G),[t$j]),j=0..n)];
			J:=[op(J1),op(J2)];
			#build the list of substitution to see derivatives as variables
			for j to n do:
				Sub:=[op(Sub),seq(diff(Xt[j],[t$k])=x[j,k],k=0..n)]
			end do;
			if ordering=plex then
				#elimination and saturation with Groebner bases
				#w.r.t. pure lex monomial ordering
				Sub:=[op(Sub),seq(diff(y(t),[t$j])=y[j],j=0..n)];
				allvars:=ListTools:-Reverse(map(rhs,Sub));
				yvars:=select(has,allvars,y);
				allvars:=allvars[numelems(yvars)+1..-1];
				J:=PolynomialIdeals:-PolynomialIdeal(subs(Sub,J),parameters=alpha);
				J:=PolynomialIdeals:-Saturate(J,subs(Sub,Q));
				J:=Groebner:-Basis(J,plex(op(allvars),op(yvars)));
				J:=remove(has,J,allvars)
			else
				#elimination and saturation with Groebner bases
				#w.r.t. lexdeg elimination ordering
				Sub:={op(Sub),seq(diff(y(t),[t$j])=y[j],j=0..n)};
				J:=PolynomialIdeals:-PolynomialIdeal(subs(Sub,J),parameters=alpha);
				J:=PolynomialIdeals:-Saturate(J,subs(Sub,Q));
				yvars:=select(has,map(rhs,Sub),y);
				J:=PolynomialIdeals:-EliminationIdeal(J,yvars);
				J:=select(type,convert(J,list),polynom)
			end if;
			#Taking a diff polynomial of minimal total degree
			# among those of the minimal order
			J:=map(de->collect(de,[seq(y[j],j=0..n)],'distributed'),J);
			Sub:=select(has,map(e->rhs(e)=lhs(e),Sub),y);
			J:=map(de->subs(Sub,de),J);
			#order
			ord:=min(map(de->PDEtools:-difforder(de,t),J));
			DE:=select(de->PDEtools:-difforder(de,t)=ord,J);
			Sub:=map(e->rhs(e)=lhs(e),Sub);
			DE:=map(de->subs(Sub,de),DE);
			#degree
			DE:=sort(DE,(a,b)->degree(a,yvars)<=degree(b,yvars));
			DE:=DE[1];
			Sub:=map(e->rhs(e)=lhs(e),Sub);
			return subs(Sub,DE)=0
		end proc:	
	
	end module: #end NLDE_nlho
	
	AnsatzDalg:= module()

	option `Copyright (c) 2022 Bertrand Teguia Tabuguia, Max Planck Institute for MiS, Leipzig`, package;

	export deltakdiff, unaryDeltak, arithmeticDeltak;

	local buildsystem, mergesystem, ComputDegkDE, DegreekDE, startkorder, ordertoktuple;


	unaryDeltak:= proc(DE::`=`,
			    y::anyfunc(name),
			    z::name=algebraic,
			    {degreeDE::posint:=2,
			    deorder::posint:=10},
			    $)::`=`;
			local t::name:=op(y),start::posint,Sys::list,x::nothing,var::name,subvars::list,SubL::list,j::posint;
			option `Copyright (c) 2022 Bertrand Teguia T.`;
			start:=PDEtools:-difforder(DE,t);
			Sys:=buildsystem(lhs(DE) - rhs(DE)=0,y,x);
			var:=op(0,y);
			subvars:=map(r->r=r(t),Sys[2]);
			Sys:=subs(subvars,Sys);
			SubL:=[seq(diff(Sys[2][j],t)=Sys[1][j],j=1..numelems(Sys[1]))];
			DegreekDE(subs(var=Sys[2][1],normal(rhs(z))),lhs(z)(t),SubL,maxdeorder=deorder,
				  ':-degreeDE'=degreeDE,startfromord=start)
		end proc:

	arithmeticDeltak:=proc(L::list(`=`),
			       V::list(anyfunc(name)),
			       z::name=ratpoly,
			       {degreeDE::posint:=2,
			       deorder::posint:=10},
			       $)::`=`;
			local t::name:=op(1,V[1]), start::posint, DEs::list(`=`), 
			      j::posint, Sys::list, subvars::list, SubL::list, subV::list;
			option `Copyright (c) 2022 Bertrand Teguia T.`;
			description "Ansatz method for the arithmetic of D-algebraic functions ";
			if numelems(L)=1 then
				return L
			end if;
			start:=min(map(r->PDEtools:-difforder(r,t),L));
			DEs:=map(r->lhs(r) - rhs(r)=0,L);
			Sys:=mergesystem(DEs,V);
			subvars:=map(r->r=r(t),Sys[3]);
			Sys:=subs(subvars,Sys);
			SubL:=[seq(diff(Sys[3][j],t)=Sys[1][j],j=1..numelems(Sys[1]))];
			subV:=[seq(op(0,V[j])=Sys[2][j],j=1..numelems(V))];
			return DegreekDE(subs(subV,rhs(z)),lhs(z)(t),SubL,
				maxdeorder=deorder,':-degreeDE'=degreeDE,startfromord=start)
		end proc:

	buildsystem:= proc(DE::`=`,
			    y::anyfunc(name),
			    x::name,
			   $)::list(`=`);
			local d::posint, t::name, SubL::list, PolDE::polynom, Xd, j::posint;
			option `Copyright (c) 2022 Bertrand Teguia T.`;
			t:=op(y);
			d:=PDEtools:-difforder(DE,t);
			SubL:=[seq(diff(y,[t$j])=x[j],j=0..d)];
			PolDE:=subs(SubL,lhs(DE));
			if degree(PolDE,x[d])>1 then
				Xd:=SolveTools:-AbstractRootOfSolution([PolDE],[x[d]]);
				return [[seq(x[j],j=1..(d-1)),rhs(op(Xd))],[seq(x[j],j=0..(d-1))]]
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
			Sys:=[seq(buildsystem(L[j],V[j],cat(x,j)),j=1..l)];
			vars:=map(r->op(r[2]),Sys);
			deriv:=map(r->op(r[1]),Sys);
			n:=numelems(vars);
			X:=[seq(vars[j]=x[j],j=1..n)];
			Ind:=[seq(1+add(numelems(Sys[i][2]),i=1..(j-1)),j=1..l)];
			[subs(X,deriv),map(r->x[r],Ind),map(rhs,X)]
		end proc:
		
	ComputDegkDE := proc(f::algebraic,z::name,k::posint,degkNmax::posint,Subdiff::list(`=`),startfromord::posint,$)
			local a::nothing, A, N:=startfromord, Eq, n, S, Sumds, nisol, j, tmp, rmS, Eqs, s, factSumds, polfact, Coef:=[], i;		
			option `Copyright (c) 2022 Bertrand Teguia T.`;
			description "Computational part of DegreekDE";
			while PDEtools:-difforder(deltakdiff(a(z),z,k,N))<=degkNmax and Coef=[] do
				A:=[seq(a[i],i=0..N-1)];
				Eq:=deltakdiff(f,z,k,N)+add(A[i+1]*deltakdiff(f,z,k,i),i=0..N-1);
				n:=PDEtools:-difforder(Eq,z);
				to n do
					Eq:=evala(eval(Eq,Subdiff))
				end do;
				Eq:=expand(numer(normal(Eq)));
				if Eq=0 then
					return [1], 1
				end if;
				S:=[op(Eq)];
				Sumds:=[];
				nisol:=true;
				while S<>[] and nisol do
					tmp:=S[1];
					rmS:=[tmp];
					if numelems(S)>1 then
						for j from 2 to numelems(S) do
							if type(normal(S[j]/S[1]), ratpoly(anything,z)) then
								tmp:=normal(tmp+S[j]);
								rmS:=[op(rmS),S[j]]
							end if
						end do
					end if;
					Sumds:=[op(Sumds),tmp];
					nisol:=has(tmp,A);
					S:=remove(member,S,rmS)
				end do;
				if nisol then 
					Sumds:=map(r->numer(factor(r)),Sumds);
					Eqs:=[];
					for s in Sumds do
						if type(s,polynom(anything,z)) then
							if has(s,A) then
								Eqs:=[op(Eqs),s]
							end if
						else
							factSumds:=map(t->exp(t),[op(simplify(ln(s),ln,'symbolic'))]);
							polfact:=select(type,factSumds,polynom(anything,z));
							polfact:=select(has,polfact,A);
							Eqs:=[op(Eqs),op(polfact)]
						end if
					end do;
					Coef:= solve(Eqs,A);
					if Coef<>[] then
						Coef:=map(rhs,Coef[1]);
						A:=map(r->r=1,A);
						Coef:=factor(subs(A,Coef))
					else
						N:=N+1
					end if
				else 
					N:=N+1
				end if
			end do;
			Coef, N
		end proc:

	DegreekDE := proc(expr::algebraic,
			F::anyfunc(name),
			sublistdiff::list(`=`),
			{maxdeorder::posint:=4,
			degreeDE::posint:=2,
			startfromord::posint:=1},
			$)::Or(`=`,identical(FAIL));
			local z, f, N, Coef, Dde, i;
			option `Copyright (c) 2022 Bertrand Teguia T.`;
			description "Compute a homogeneous degree k differential equation";
			z:=op(1,F);
			#some normalization on the input
			f:=expand(evala(expr));
			if f=0 then 
				return F=0
			end if;
			#compute the coeficients of the DE sought
			Coef,N:=ComputDegkDE(f,z,degreeDE,maxdeorder,sublistdiff,startkorder(startfromord,degreeDE,F,z));
			#if there is a solution (Coef is not empty)
			if Coef<>[] then
				#clearing denominators
				Dde:=lcm(op(denom(Coef)));
				Coef:=map(r-> factor(normal(Dde*r)), Coef);
				Coef:=[op(Coef),Dde];
				return add(Coef[i+1]*deltakdiff(F,z,degreeDE,i),i=0..N)=0
			else
				userinfo(2,DegreekDE,printf("No homogeneous degree %d DE of order at most %d found\n", k, maxdeorder));
				return FAIL
			end if
		end proc:

	#minmaxorder ... ((n + 1)*(n^4 + 19*n^3 + 136*n^2 + 444*n + 600))/120

	startkorder := proc(n::nonnegint,k,F,z,$)
			local j:=k+1,df;
			df:=deltakdiff(F,z,k,j);
			while PDEtools:-difforder(df,z)<n do
				j:=j+1;
				df:=deltakdiff(F,z,k,j)
			end do;
			return j
		end proc:

	#Derivative operator used to compute product of k derivatives
	#with respect to a certain 'natural ordering'
	deltakdiff := proc(expr,z::name,k::posint:=2,n::nonnegint:=1,$)
			local tuple;
			option `Copyright (c) 2022 Bertrand Teguia T.`;
			tuple:= select(type,ordertoktuple(k,n)-~1,nonnegint);
			return mul(map(d-> diff(expr,[seq(z,1..d)]), tuple))
		end proc:
		
		
	#A bijection between a subset of N^k and N is used 
	#to define the ordering for differential monomials
	ordertoktuple := proc(k::posint,n::nonnegint) option remember;
			local m, im, Tkn, j;
			option `Copyright (c) 2022 Bertrand Teguia T.`;
			if n=0 then 
			    return [seq(0,j=1..k)]
			else   
			    Tkn:=ordertoktuple(k,n-1);
			    m:=min(Tkn);
			    im:=ListTools:-Search(m,Tkn);
			    Tkn[im]:=m+1;
			    if im=k then
				return Tkn
			    else
				return [op(Tkn[1..im]),seq(0,j=im+1..k)]
			    end if     
			end if    
		end proc:

	end module: #end AnsatzDalg

end module:

savelib('NLDE',"path_to_your_folder_for_softwarepackages/NLDE.mla"):
