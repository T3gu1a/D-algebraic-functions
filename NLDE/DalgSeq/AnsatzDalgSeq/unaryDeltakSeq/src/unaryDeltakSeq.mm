
unaryDeltakSeq:= proc(DE::`=`,
		       y::anyfunc(name),
		       z::name=algebraic,
		 {degADE::posint:=2,
	    startfromord::posint:=1,
	      maxdeorder::posint:=2},
		      $)::Or(`=`,identical(FAIL));
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
			':-maxdeorder'=max(maxdeorder,ord,startfromord),
			':-degreeDE'=degADE,':-startfromord'=max(ord,startfromord))
	end proc: