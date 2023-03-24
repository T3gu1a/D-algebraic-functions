
arithmeticDeltak:=proc(L::list(`=`),
		       V::list(anyfunc(name)),
		       z::name=ratpoly,
		       {degreeDE::posint:=2,
		       startingorder::posint:=1,
		       maxdeorder::posint:=2},
		       $)::`=`;
		local t::name:=op(1,V[1]), start::posint, Ords::list(posint), DEs::list(`=`), 
		      j::posint, Sys::list, subvars::list, SubL::list, subV::list;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Ansatz method for the arithmetic of D-algebraic functions ";
		if numelems(L)=1 then
			return L
		end if;
		Ords:=map(r->PDEtools:-difforder(r,t),L);
		start:=max(min(Ords),startingorder);
		DEs:=map(r->lhs(r) - rhs(r)=0,L);
		Sys:=mergesystem(DEs,V);
		subvars:=map(r->r=r(t),Sys[3]);
		Sys:=subs(subvars,Sys);
		SubL:=[seq(diff(Sys[3][j],t)=Sys[1][j],j=1..numelems(Sys[1]))];
		subV:=[seq(op(0,V[j])=Sys[2][j],j=1..numelems(V))];
		return DegreekDE(subs(subV,rhs(z)),lhs(z)(t),SubL,
			':-maxdeorder'=max(maxdeorder,start),
			':-degreeDE'=degreeDE,startfromord=start)
	end proc: