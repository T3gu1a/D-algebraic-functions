
DDfiniteToDalg:= proc(DE::`=`,y::anyfunc(name),
		     Leq::list(`=`),
		   Lvars::list(anyfunc(name)),
	       {ordering::identical(plex,lexdeg):=plex},
		      $)::`=`;
		local DEs, V, Sys, j, subV;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description  "Convert a DD-finite equation into an ADE    "
			     "INPUT: - The DD-finite ODE DE               "
			     "       - Its dependent variable y(x)        "
			     "       - The list of holonomic DEs          "
			     "       - Their dependent variables in Lvars "
			     "OUTPUT: An ADE fulfills by the solutions of "
			     "        the given DD-finite ODE             ";
		DEs:=map(eq->lhs(eq) - rhs(eq)=0,[DE,op(Leq)]);
		V:=[y,op(Lvars)];
		Sys:=mergesystem(DEs,V);
		subV:=[seq(op(0,Lvars[j])=Sys[2][j+1],j=1..numelems(Lvars))];
		return SysToMinDiffPoly(subs(subV,Sys[1]),Sys[2][1],Sys[3],y,':-ordering'=ordering)
	end proc:
	