

SystoDE:= proc(f::list(algebraic),
	       g::algebraic,
	       X::Or(list,set),
	       z::anyfunc(name),
	      $)::algebraic;
			option `Copyright (c) 2024 Bertrand Teguia T.`;
			description  "The non-lho analogue of NLDE:-SysToMinDiffPoly for sequences";
			local F,G,q1,q2,Q,Svars,J1,J2,J,n,Xt,nlho_pow,nlho_lead,
			      t,DE,Sub:=[],allvars,yvars,ord,j,k,y,alpha;
			t:=op(1,z);
			y:=op(0,z);
			alpha:=indets([f,g]) minus {op(map(x->x[1],X))};
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
			J1:=[seq(Q*LREtools:-shift(Xt[j],t)^nlho_pow[j]-normal(Q*F[j]),j=1..n)];
			#non-lho case with many different degrees on the leader
			nlho_lead:=select(v->numelems(v)=3,X);
			if nlho_lead <> [] then
				nlho_lead:=map(v->v[3]=LREtools:-shift(v[1](t),t),nlho_lead);
				J1:=subs(nlho_lead,J1)
				#J1[n]:=subs(X[n][3]=diff(Xt[n],t),J1[n])
			end if;
			#shifting n-1 times the polynomials Q*x'-Q*f
			for j to n do:
				J1:=[op(J1),seq(LREtools:-shift(J1[j],t,k),k=1..(n-1))]
			end do;
			#differentiating n times the polynomials Q*y - Q*g
			J2:=[seq(LREtools:-shift(Q*y(t)-normal(Q*G),t,j),j=0..n)];
			J:=[op(J1),op(J2)];
			#build the list of substitution to see shifts as variables
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
			ord:=min(map(de->HyperTypeSeq:-AlgebraHolonomicSeq:-REorder(de,y(t)),J));
			DE:=select(de->HyperTypeSeq:-AlgebraHolonomicSeq:-REorder(de,y(t))=ord,J);
			Sub:=map(e->rhs(e)=lhs(e),Sub);
			DE:=map(de->subs(Sub,de),DE);
			#degree
			DE:=sort(DE,(a,b)->degree(a,yvars)<=degree(b,yvars));
			DE:=DE[1];
			Sub:=map(e->rhs(e)=lhs(e),Sub);
			return subs(Sub,DE)=0
		end proc: