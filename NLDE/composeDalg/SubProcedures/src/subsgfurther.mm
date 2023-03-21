
subsgfurther :=proc(gm1::algebraic,g::name,t::name,m::posint,n::posint,$)::list;
		local k::posint,j::nonnegint,Subdiff::list,rSubdiff::list,eqg,Sub::list;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description     "subprocedure of composeDalg for substituting higher derivatives"
		                "of g in the rational relation R                                "
				"INPUT:  - an algebraic relation representing the link from     "
				"          the differential equation                            "
				"        - the name for g                                       "
				"        - the name for the independent variable t              "
				"        - the order of the equation of g, m                    "
				"        - the order n >= m of the other equation               "
				"OUPUT: the list required to make the substitution              ";
		k:=0;
		eqg:=gm1;
		#the first substition
		Sub:=[g[m]=eqg];
		if k<n-k then
			Subdiff:=[seq(g[j]=diff(g(t),[t$j]),j=0..m-1)];
			rSubdiff:=[op(map(r->rhs(r)=lhs(r),Subdiff)),diff(g(t),t$m)=g[m]];
			#the other substitutions are obtained by
			#differentiating and substituting the first substitution
			while k<n-m do
				eqg:=subs(Subdiff,eqg);
				eqg:=diff(eqg,t);
				eqg:=subs(rSubdiff,eqg);
				eqg:=normal(subs(Sub[1],eqg));
				k:=k+1;
				Sub:=[op(Sub),g[m+k]=eqg]
			end do
		end if;
		return Sub
	end proc:
