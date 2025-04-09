

DegreekDE := proc(expr::algebraic,
		F::anyfunc(name),
		sublistdiff::list(`=`),
		{maxdeorder::posint:=4,
		degreeDE::posint:=2,
		startfromord::posint:=1},
		$)::Or(`=`,identical(FAIL));
		local z, f, N, Coef, Dde, i;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Compute a homogeneous degree k differential equation";
		z:=op(1,F);
		#some normalization on the input
		f:=expand(evala(expr));
		if f=0 then 
			return F=0
		end if;
		#compute the coeficients of the DE sought
		Coef,N:=ComputDegkDE(f,z,degreeDE,maxdeorder,sublistdiff,shiftkstart(startfromord,degreeDE,F,z));
		#if there is a solution (Coef is not empty)
		if Coef<>[] then
			#clearing denominators
			Dde:=lcm(op(denom(Coef)));
			Coef:=map(r-> factor(normal(Dde*r)), Coef);
			Coef:=[op(Coef),Dde];
			return add(Coef[i+1]*deltakshift(F,z,degreeDE,i),i=0..N)=0
		else
			userinfo(2,DegreekDE,printf("No homogeneous degree %d DE of order at most %d found\n", k, maxdeorder));
			return FAIL
		end if
	end proc:
	
$include <NLDE/DalgSeq/AnsatzDalgSeq/CommonInternalProcedures/src/ComputDegkDE.mm>