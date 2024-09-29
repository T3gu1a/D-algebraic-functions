

buildsystemseq:= proc(DE::`=`,
		       y::anyfunc(name),
		       x::name,
		      $)::list(`=`);
			local  t::name, r::integer, r0::integer, SubL::list, PolDE::polynom, 
			       d::posint, j::nonnegint;
			option `Copyright (c) 2024 Bertrand Teguia T.`;
			description     "The non-lho anologue of buildsystem";
			t:=op(y);
			(r,r0):=REorders(DE, y);
			PolDE:=LREtools:-shift(lhs(DE),t,-r0);
			r:=r-r0;
			#variables of substitution for the model, the input x with indices
			SubL:=[seq(LREtools:-shift(y,t,j)=x[j],j=0..r)];
			PolDE:=subs(SubL,PolDE);
			d:=degree(PolDE,x[r]);
			#the differential equation is not l.h.o
			if d>1 then
				PolDE:=subs(x[r]^d=x[r+1],collect(PolDE,x[r],'distributed'));
				return [[seq(x[j],j=1..(r-1)),solve(PolDE,x[r+1])],[seq([x[j],1],j=0..(r-2)),[x[r-1],d,x[r]]]]
			else
				return [[seq(x[j],j=1..(r-1)),solve(PolDE,x[r])],[seq([x[j],1],j=0..(r-1))]]
			end if	
		end proc: