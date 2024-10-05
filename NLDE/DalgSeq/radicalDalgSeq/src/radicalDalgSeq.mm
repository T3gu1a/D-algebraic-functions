
#radicals of D-algebraic sequences

radicalDalgSeq:= proc(DE::`=`,
			y::anyfunc(name),
			z::name,
			N::nonnegint,
			$)::`=`;
			local t::name:=op(y),Sys::list,x::nothing;
			option `Copyright (c) 2024 Bertrand Teguia T.`;
			description " Analog of unaryDalg for sequences";
			#build the system using buildsystem
			Sys:=buildsystemseq(numer(normal(lhs(DE) - rhs(DE)))=0,y,x);
			#use SysToMinDiffPoly to return the desired output
			return SystoDE(Sys[1],Sys[2][1][1],Sys[2],z(t),N)
		end proc: