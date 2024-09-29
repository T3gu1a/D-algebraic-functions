
#arithmetic of a single D-algebraic sequence

unaryDalgSeq:= proc(DE::`=`,
	             y::anyfunc(name),
		     z::name=ratpoly,
		    $)::`=`;
			local t::name:=op(y),var::name:=op(0,y),dvar::name:=lhs(z),
			      r::ratpoly:=rhs(z),Sys::list,x::nothing;
			option `Copyright (c) 2024 Bertrand Teguia T.`;
			description " Analog of unaryDalg for sequences";
			#build the system using buildsystem
			Sys:=buildsystemseq(lhs(DE) - rhs(DE)=0,y,x);
			#use SysToMinDiffPoly to return the desired output
			return SystoDE(Sys[1],subs(var=Sys[2][1][1],r),Sys[2],dvar(t))
		end proc: