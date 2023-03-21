
arithmeticDalg:=proc(L::list(`=`),
		     V::list(anyfunc(name)),
		     z::name=ratpoly,
		    {ordering::identical(plex,lexdeg):=plex,
		    lho::truefalse:=true,
		    lhoplex::truefalse:=false},
		    $)::`=`;
		local t:=op(1,V[1]),DEs::list(`=`),Sys::list,j::posint,subV::list;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Compute a differential equaton for a rational expression of        "
		            "D-algebraic functions from differential equations given in the     "
			    "same order with the second argument, representing the list of      "
			    "dependent variables.                                               "
			    "One may supply more than two differential equations.               "
			    "INPUT: - a list of differential equations [DE(f1(t)),...,DE(fn(t))]"
			    "       - a list of their dependent variables [f1(t),...,fn(t)]     "
			    "	    - an equation z=r(f1,...,fn) where z is the name of the     "
			    "         dependent variable for the output, and r(f1,...,fn) is    "
			    "         is a rational expression in f1,...,fn.                    "
			    "OUPUT: a differential equation satisfied by the rational expresion " 
	                    "       r(f1,f2,...,fn)                           		        ";
		if numelems(L)=1 then
			return L
		end if;
		DEs:=map(r->lhs(r) - rhs(r)=0,L);
		#build the systems and merge them using mergesystem
		Sys:=ifelse(lho,mergesystem(DEs,V),NLDE_nlho:-mergesystem(DEs,V));
		#prepare the list for the change of variables 
		#in r according to Sys
		subV:=[seq(op(0,V[j])=Sys[2][j],j=1..numelems(V))];
		#use SystoMinDiffPoly to return the desired output
		if lho then
			return SystoMinDiffPoly(Sys[1],subs(subV,rhs(z)),Sys[3],lhs(z)(t),':-ordering'=ordering)
		else
			if lhoplex then 
				return NLDE_nlho:-SystoMinDiffPoly(Sys[1],subs(subV,rhs(z)),Sys[3],lhs(z)(t),':-ordering'=plex)
			else
				return NLDE_nlho:-SystoMinDiffPoly(Sys[1],subs(subV,rhs(z)),Sys[3],lhs(z)(t),':-ordering'=lexdeg)
			end if
		end if
	end proc:
	