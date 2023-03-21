
SysToMinDiffPoly:= proc(f::list(algebraic),
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

