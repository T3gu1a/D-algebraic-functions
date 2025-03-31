
#arithmetic of a single D-algebraic sequence

unaryDalgSeq:= proc(DE::`=`,
		y::anyfunc(name),
		z::name=algebraic,
		$)::`=`;
		local t::name:=op(y),var::name:=op(0,y),dvar::name:=lhs(z),
		      r::ratpoly:=rhs(z),Sys::list,x::nothing,shiftL::list,j::nonnegint;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		description " Analog of unaryDalg for sequences";
		#build the system using buildsystem
		Sys:=buildsystemseq(numer(normal(lhs(DE) - rhs(DE)))=0,y,x);
		#when the ratpoly contains shifts of y(t)
		if has(r,t) then
			shiftL:=[seq(var(t+j)=Sys[2][j+1][1],j=0..numelems(Sys[1])-1)];
			r:=subs(shiftL,r)
		end if; 
		#use SysToMinDiffPoly to return the desired output
		return SystoDE(Sys[1],subs(var=Sys[2][1][1],r),Sys[2],dvar(t))
	end proc: