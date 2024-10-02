
#Converting CC-finite to D-algebraic

CCfiniteToDalg := proc(DE::`=`,
			y::anyfunc(name),
		      Leq::list(`=`),
		    Lvars::list(anyfunc(name)),
		       $)::`=`;
		local DEs, V, Sys, j, subV;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		description  "Analog of DDfiniteToDalg for sequences";
		DEs:=map(eq->lhs(eq) - rhs(eq)=0,[DE,op(Leq)]);
		V:=[y,op(Lvars)];
		Sys:=mergesystemseq(DEs,V);
		subV:=[seq(op(0,Lvars[j])=Sys[2][j+1],j=1..numelems(Lvars))];
		return SystoDE(subs(subV,Sys[1]),Sys[2][1],Sys[3],y)
	end proc:
