
#arithmetic of multivariate D-algebraic functions
arithmeticMDalg := proc(L::list(`=`),
			V::list(function(name)),
			rat::name=ratpoly,
			IV::list(name):=[],
			{ordering::identical(plex,lexdeg):=plex,
			maxord::list(nonnegint):=NULL,
			elimvars::set(name):={}},
			$)::`=`;
		local DEs, indvars, n, orders, N, i, j, Nu, nu, z, R, E, sigmaorders,
		      d, Subz, dm, SubV, Sub, P, J, backSubz, mthord;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		#remove right-hand sides
		DEs:=map(ade->lhs(ade)-rhs(ade),L);
		#the independent variables
		if IV=[] then
			indvars:=[op(V[1])]
		else
			indvars:=IV
		end if;
		n:=numelems(indvars);
		#the maximum order w.r.t. to each independent variables
		orders:=map(ade->[seq(PDEtools:-difforder(ade,indvars[j]),j=1..n)],DEs);
		#the sum of these orders component wise
		N:=numelems(DEs);
		Nu:=[seq(add(orders[i][j],i=1..N),j=1..n)];
		#nu is the max-theta 
		nu:=max(CantorSigma(Nu),CantorSigma(maxord));
		Nu:=CantorInvSigma(n,nu);
		z:=lhs(rat)(op(indvars));
		#writing the dependent variables as functions of the independent ones
		SubV:=map(v->op(0,v)=v,V);
		R:=numer(z-subs(SubV,rhs(rat)));
		#compute the first theta derivatives of R and select those of the desired
		# order in z
		E:=[seq(derivation(R,indvars,j),j=0..nu)];
		E:=select(ade->compwiseless([seq(PDEtools:-difforder(ade,indvars[j]),j=1..n)],Nu),E);
		#compute d as the theta-order at which derivatives of the input
		#differential polynomials are of order at least nu
		sigmaorders:=map(r->CantorSigma(r),map(r->Nu-r,orders));
		d:=max(sigmaorders);
		#add the first d theta-derivatives of the input diff. poly
		E:=[op(E),seq(seq(derivation(DEs[i],indvars,j),j=0..d),i=1..N)];
		#build a list of substitution to have everything as polynomials
		Subz:=[seq(derivation(z,indvars,i)=op(0,z)[seq(PDEtools:-difforder(derivation(z,indvars,i),indvars[j]),j=1..n)],i=0..nu)];
		dm:=CantorSigma(CantorInvSigma(n,d-min(sigmaorders))+Nu);
		SubV:=map(v->seq(derivation(v,indvars,i)=op(0,v)[seq(PDEtools:-difforder(derivation(v,indvars,i),indvars[j]),j=1..n)],i=0..dm),V);
		SubV:=remove(v->lhs(v)=0,SubV);
		Sub:=[op(SubV),op(Subz)];
		P:=subs(Sub,E);
		#elimination with Groebner basis
		Sub:=[op(elimvars),op(ListTools:-MakeUnique(map(rhs,Sub)))];
		J:=PolynomialIdeals:-PolynomialIdeal(P,parameters={op(indvars)} minus elimvars);
		if ordering = plex then
			J:=Groebner:-Basis(J,plex(op(Sub)));
			J:=remove(has,J,[op(elimvars),op(map(rhs,SubV))])
		else
			J:=PolynomialIdeals:-EliminationIdeal(J,{op(map(rhs,Subz))});
			J:=select(type,convert(J,list),polynom);
			J:=map(de->collect(de,map(rhs,Subz),'distributed'),J)
		end if;
		d:=d+1;
		#if the elimination ideal is empty, repeat the computation with d+1,..,nu
		while J=[] and d<=nu do:
			E:=[op(E),seq(derivation(DEs[i],indvars,d),i=1..N)];
			dm:=CantorSigma(CantorInvSigma(n,d-min(sigmaorders))+Nu);
			SubV:=map(v->seq(derivation(v,indvars,i)=op(0,v)[seq(PDEtools:-difforder(derivation(v,indvars,i),indvars[j]),j=1..n)],i=0..dm),V);
			SubV:=remove(v->lhs(v)=0,SubV);
			Sub:=[op(SubV),op(Subz)];
			P:=subs(Sub,E);
			Sub:=[op(elimvars),op(ListTools:-MakeUnique(map(rhs,Sub)))];
			J:=PolynomialIdeals:-PolynomialIdeal(P,parameters={op(indvars)} minus elimvars);
			if ordering = plex then
				J:=Groebner:-Basis(J,plex(op(Sub)));
				J:=remove(has,J,[op(elimvars),op(map(rhs,SubV))])
			else
				J:=PolynomialIdeals:-EliminationIdeal(J,{op(map(rhs,Subz))});
				J:=select(type,convert(J,list),polynom);
				J:=map(de->collect(de,map(rhs,Subz),'distributed'),J)
			end if;
			d:=d+1
		end do;
		#No ADE of order at most Nu or maxord (component wise) found if J=[]
		if J=[] then
			return 0
		end if;
		#selecting the differential polynomial of lowest degree among those
		#among those of the lowest order possible
		backSubz:=map(v->rhs(v)=lhs(v),Subz);
		J:=subs(backSubz,J);
		if numelems(J)=1 then
			return(J[1])=0
		end if;
		#the least theta-order differential polynomial
		orders:=map(ade->[seq(PDEtools:-difforder(ade,indvars[j]),j=1..n)],J);
		orders:=map(CantorSigma,orders);
		mthord:=min(orders);
		J:=map(pos->J[pos],[ListTools:-SearchAll(mthord,orders)]);
		#the least degree differential polynomial
		J:=subs(Subz,J);
		J:=sort(J,(p1,p2)->degree(p1,map(rhs,Subz))<=degree(p2,map(rhs,Subz)));
		J:=subs(backSubz,J);
		return J[1]=0
	end proc:	