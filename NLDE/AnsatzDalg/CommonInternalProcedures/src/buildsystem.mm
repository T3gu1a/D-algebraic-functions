
buildsystem:= proc(DE::`=`,
		    y::anyfunc(name),
		    x::name,
		   $)::list(`=`);
		local d::posint, t::name, SubL::list, PolDE::polynom, Xd, j::posint;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		t:=op(y);
		d:=PDEtools:-difforder(DE,t);
		SubL:=[seq(diff(y,[t$j])=x[j],j=0..d)];
		PolDE:=subs(SubL,lhs(DE));
		if degree(PolDE,x[d])>1 then
			Xd:=SolveTools:-AbstractRootOfSolution([PolDE],[x[d]]);
			return [[seq(x[j],j=1..(d-1)),rhs(op(Xd))],[seq(x[j],j=0..(d-1))]]
		else
			return [[seq(x[j],j=1..(d-1)),solve(PolDE,x[d])],[seq(x[j],j=0..(d-1))]]
		end if	
	end proc:
	