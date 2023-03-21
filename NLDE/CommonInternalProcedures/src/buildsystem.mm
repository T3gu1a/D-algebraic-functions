

buildsystem:= proc(DE::`=`,
		    y::anyfunc(name),
		    x::name,
		   $)::list(`=`);
		local  t::name, d::posint, SubL::list, PolDE::polynom, j::nonnegint;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description     "Build a dynamical system (or model) from a differential equation.    "
				"If the differential equation is not l.h.o, its derivatives is used.  "
				"INPUT: -A differential equation DE,                                  "
				"       -its dependent variable like y(t)                             "
				"	-a name x for the variable of the system                      "
				"OUPUT: A list of two lists:                                          "
				"       - the list of derivatives of the variables of the system      "
				"         in in terms of these variables                              " 
				"	- the variables of the system                                 ";
		t:=op(y);
		d:=PDEtools:-difforder(DE,t);
		#variables of substitution for the model, the input x with indices
		SubL:=[seq(diff(y,[t$j])=x[j],j=0..d)];
		PolDE:=subs(SubL,lhs(DE));
		#the differential equation is not l.h.o
		if degree(PolDE,x[d])>1 then
			d:=d+1;
			SubL:=[op(SubL),diff(y,t$d)=x[d]];
			PolDE:=subs(SubL,lhs(diff(DE,t)));
			return [[seq(x[j],j=1..(d-1)),solve(PolDE,x[d])],[seq(x[j],j=0..(d-1))]]
		else
			return [[seq(x[j],j=1..(d-1)),solve(PolDE,x[d])],[seq(x[j],j=0..(d-1))]]
		end if		
	end proc:
