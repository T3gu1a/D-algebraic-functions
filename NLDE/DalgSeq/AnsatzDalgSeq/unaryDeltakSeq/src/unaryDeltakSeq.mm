
unaryDeltakSeq:= proc(DE::`=`,
		    y::anyfunc(name),
		    z::name=algebraic,
		    {degreeDE::posint:=2,
		    startingorder::posint:=1,
		    maxdeorder::posint:=2},
		    $)::`=`;
		local t::name:=op(y),Sys::list,ord::posint,
		      x::nothing,var::name,subvars::list,SubL::list,j::posint;
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		Sys:=buildsystem(lhs(DE) - rhs(DE)=0,y,x);
		var:=op(0,y);
		ord:=REorders(DE,var(t))[1];
		subvars:=map(r->r=r(t),Sys[2]);
		Sys:=subs(subvars,Sys);
		SubL:=[seq(diff(Sys[2][j],t)=Sys[1][j],j=1..numelems(Sys[1]))];
		DegreekDE(subs(var=Sys[2][1],normal(rhs(z))),lhs(z)(t),SubL,
			':-maxdeorder'=max(maxdeorder,ord,startingorder),
			':-degreeDE'=degreeDE,startfromord=max(ord,startingorder))
	end proc: