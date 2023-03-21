
unaryDeltak:= proc(DE::`=`,
		    y::anyfunc(name),
		    z::name=algebraic,
		    {degreeDE::posint:=2,
		    maxdeorder::posint:=10},
		    $)::`=`;
		local t::name:=op(y),start::posint,Sys::list,Ords::list(posint),
		      x::nothing,var::name,subvars::list,SubL::list,j::posint;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		Ords:=map(r->PDEtools:-difforder(r,t),L);
		start:=max(min(Ords),startingorder);
		Sys:=buildsystem(lhs(DE) - rhs(DE)=0,y,x);
		var:=op(0,y);
		subvars:=map(r->r=r(t),Sys[2]);
		Sys:=subs(subvars,Sys);
		SubL:=[seq(diff(Sys[2][j],t)=Sys[1][j],j=1..numelems(Sys[1]))];
		DegreekDE(subs(var=Sys[2][1],normal(rhs(z))),lhs(z)(t),SubL,
			':-maxdeorder'=max(add(Ords),maxdeorder),
			':-degreeDE'=degreeDE,startfromord=start)
	end proc:

