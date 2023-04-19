
diffDalg :=proc(DE::`=`,
		 y::anyfunc(name),
		 n::posint:=1,
		{ordering::identical(plex,lexdeg):=plex},
		$)::`=`;
		local t::name:=op(y),var::name:=op(0,y),C::nothing,
		      d::nonnegint,j::nonnegint,p,q,R,V;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		description "Compute a differential equation for                                   "
		            "the derivative of a D-algebraic function from a                       "
			    "differential equation that it satisfies.                              "
			    "INPUT: - a differential equation                                      "
			    "       - its dependent variable, say f(t)                             "
			    "OUPUT: a differential equation satisfied by diff(f(t),t)              ";
		if n=1 then	
			p:=lhs(DE)-rhs(DE);
			d:=PDEtools:-difforder(p,t);
			q:=diff(p,t);
			p:=subs([seq(diff(var(t),[t$j])=var[j],j=0..d)],p);
			#no need for elimination
			if not has(p, var[0]) then
				return subs([seq(var[j]=diff(var(t),[t$(j-1)]),j=1..d)],p)=0
			end if;
			q:=subs([seq(diff(var(t),[t$j])=var[j],j=0..d+1)],q);
			R:=resultant(p,q,var[0]);
			if type(R,`+`) then
				return subs([seq(var[j]=diff(var(t),[t$(j-1)]),j=1..d+1)],R)=0
			elif type(R,`*`) then
				R:=[select(has,R,var[d+1])];
				V:=[seq(var[j],j=0..d+1)];
				R:=sort(R,(a,b)->degree(a,V)<=degree(b,V));
				return subs([seq(var[j]=diff(var(t),[t$(j-1)]),j=1..d+1)],R[1])=0
			else
				return subs([seq(var[j]=diff(var(t),[t$(j-1)]),j=1..d+1)],R)=0
			end if
		else
			return diffDalg(diffDalg(DE,y,n-1),y)
		end if
	end proc:
	