
#Finite series of a D-algebraic sequence

FiniteSeriesDalgSeq:=proc(DE::`=`,
			   y::anyfunc(name),
			   z::anyfunc(name),
		          $)::`=`;
			local t::name:=op(y),var::name:=op(0,y),Sys::list,x::nothing;
			option `Copyright (c) 2024 Bertrand Teguia T.`;
			description "Analog of unaryDalg for sequences";
			#build the system using buildsystemseq
			Sys:=buildsystemseq(lhs(DE) - rhs(DE)=0,y,x);
			#use SystoDE to return the desired output
			return SystoDE([x[-1]+Sys[2][1][1],op(Sys[1])],x[-1],[[x[-1],1],op(Sys[2])],z)
		end proc: