
#submodule used to overcome the non-l.h.o situation in the unary case
#This generalizes to the general computation (proof to be done)
#However, generally the minimal ADE in the non-l.h.o case is not "good-looking"
#Hence the reason for not using it in general.
NLDE_nlho:= module()

option `Copyright (c) 2023 Bertrand Teguia Tabuguia, Max Planck Institute for MiS, Leipzig`, package;

export buildsystem, mergesystem, SysToMinDiffPoly;


buildsystem:= proc(DE::`=`,
		   y::anyfunc(name),
		   x::name,
		   $)::list(`=`);
		local  t::name, r::posint, SubL::list, PolDE::polynom, d::posint, j::nonnegint, sep;
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
			sep:=diff(PolDE,x[r]);
			PolDE:=subs(x[r]^d=x[r+1],collect(PolDE,x[r],'distributed'));
			return [[seq(x[j],j=1..(r-1)),solve(PolDE,x[r+1])],[seq([x[j],1],j=0..(r-2)),[x[r-1],d,x[r]]],sep]
		else
			return [[seq(x[j],j=1..(r-1)),solve(PolDE,x[r])],[seq([x[j],1],j=0..(r-1))],1]
		end if	
	end proc:

mergesystem:= proc(L::list(`=`),
		   V::list(anyfunc(name)),
		   $)::`=`;
		local l::posint:=numelems(L), j::posint, Sys::list, vars::list, deriv::list, 
		      n::posint, Sep::algebraic, x::nothing, X::list, i::posint, Ind::list;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description     "The non-lho analogue of mergesystem";
		Sys:=[seq(buildsystem(L[j],V[j],cat(x,j)),j=1..l)];
		vars:=map(r->op(r[2]),Sys);
		deriv:=map(r->op(r[1]),Sys);
		Sep:=lcm(op(map(r->r[3], Sys)));
		n:=numelems(vars);
		X:=[seq(vars[j][1]=x[j],j=1..n)];
		#indices of the variables representing the solutions of the input DEs
		Ind:=[seq(1+add(numelems(Sys[i][2]),i=1..(j-1)),j=1..l)];
		return [subs(X,deriv),map(r->x[r],Ind),subs(X,vars),subs(X,Sep)]
	end proc:


SysToMinDiffPoly:= proc(f::list(algebraic),
			g::algebraic,
			X::Or(list,set),
			z::anyfunc(name),
			Sep::algebraic,
			{sepsols::truefalse:=false,
			ordering::identical(plex,lexdeg):=plex},
			$)::algebraic;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description     "The non-lho analogue of NLDE:-SysToMinDiffPoly";
		local F,G,q1,q2,Q,sep,Svars,J1,J2,J,n,Xt,nlho_pow,nlho_lead,
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
		sep:=subs(Xt,Sep);
		Xt:=map(rhs,Xt);
		nlho_pow:=map(x->x[2],X);
		J1:=[seq(Q*diff(Xt[j],t)^nlho_pow[j]-normal(Q*F[j]),j=1..n)];
		#non-lho case with many different degrees on the leader
		nlho_lead:=select(v->numelems(v)=3,X);
		if nlho_lead <> [] then
			nlho_lead:=map(v->v[3]=diff(v[1](t),t),nlho_lead);
			sep:=subs(nlho_lead,sep);
			J1:=subs(nlho_lead,J1)
			#J1[n]:=subs(X[n][3]=diff(Xt[n],t),J1[n])
		end if;
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
			J:=ifelse(sepsols,PolynomialIdeals:-Saturate(J,subs(Sub,Q)),
					  PolynomialIdeals:-Saturate(J,subs(Sub,lcm(Q,sep))));
			J:=Groebner:-Basis(J,plex(op(allvars),op(yvars)));
			J:=remove(has,J,allvars)
		else
			#elimination and saturation with Groebner bases
			#w.r.t. lexdeg elimination ordering
			Sub:={op(Sub),seq(diff(y(t),[t$j])=y[j],j=0..n)};
			J:=PolynomialIdeals:-PolynomialIdeal(subs(Sub,J),parameters=alpha);
			J:=ifelse(sepsols,PolynomialIdeals:-Saturate(J,subs(Sub,Q)),
					  PolynomialIdeals:-Saturate(J,subs(Sub,lcm(Q,sep))));
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