
unaryDalg:= proc(DE::`=`,
		y::anyfunc(name),
		z::name=ratpoly,
		{ordering::identical(plex,lexdeg):=lexdeg},
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