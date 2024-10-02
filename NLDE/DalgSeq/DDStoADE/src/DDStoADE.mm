#Discrete dynamical system to algebraic difference equation

DDStoADE:= proc(f::list(algebraic),
                        g::algebraic,
                        X::Or(list,set),
                        z::anyfunc(name),
			$)::algebraic;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		local F,G,q1,q2,Q,x::nothing,J1,J2,J,n,Xt,t,DE,Sub:=[],yvars,ord,j,k,y,alpha;
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
		J1:=[seq(Q*LREtools:-shift(Xt[j],t)-normal(Q*F[j]),j=1..n)];
		#differentiating n-1 times the polynomials Q*x'-Q*f
		for j to n do:
			J1:=[op(J1),seq(LREtools:-shift(J1[j],t,k),k=1..(n-1))]
		end do;
		#differentiating n times the polynomials Q*y - Q*g
		J2:=[seq(LREtools:-shift(Q*y(t)-normal(Q*G),t,j),j=0..n)];
		J:=[op(J1),op(J2)];
		#build the list of substitution to see derivatives as variables
		for j to n do:
			Sub:=[op(Sub),seq(LREtools:-shift(Xt[j],t,k)=x[j,k],k=0..n)]
		end do;
		#elimination and saturation with Groebner bases
		#w.r.t. lexdeg elimination ordering
		Sub:={op(Sub),seq(LREtools:-shift(y(t),t,j)=y[j],j=0..n)};
		J:=PolynomialIdeals:-PolynomialIdeal(subs(Sub,J),'parameters'=alpha);
		J:=PolynomialIdeals:-Saturate(J,subs(Sub,[seq(LREtools:-shift(Q,t,j),j=0..n)]));
		yvars:=select(has,map(rhs,Sub),y);
		J:=PolynomialIdeals:-EliminationIdeal(J,yvars);
		J:=select(type,convert(J,list),polynom);
		#Taking a diff polynomial of minimal total degree
		# among those of the minimal order
		J:=map(de->collect(de,[seq(y[j],j=0..n)],'distributed'),J);
		Sub:=select(has,map(e->rhs(e)=lhs(e),Sub),y);
		J:=map(de->subs(Sub,de),J);
		#order
		ord:=min(map(de->REorders(de,y(t))[1],J));
		DE:=select(de->REorders(de,y(t))[1]=ord,J);
		Sub:=map(e->rhs(e)=lhs(e),Sub);
		DE:=map(de->subs(Sub,de),DE);
		#degree
		DE:=sort(DE,(a,b)->degree(a,yvars)<=degree(b,yvars));
		DE:=DE[1];
		Sub:=map(e->rhs(e)=lhs(e),Sub);
		return subs(Sub,DE)=0
	end proc: