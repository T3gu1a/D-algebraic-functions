
#convert differential equations to recurrences
ADEtoRE := proc(ADE::algebraic,Y::anyfunc(name),A::anyfunc(name),K::list(name),$)::algebraic;
		local terms, RE;
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "convert algebraic differential equations to recurrence equations";
		#Let terms be the terms in the computed differential equation ADE
		if type(ADE,`+`) then
			terms:=[op(ADE)]
		else
			terms:=[ADE]
		end if;
		#Transform each terms of the differential equation to its correspondent for the recurrence equation
		RE:=map(T->ADEtermToREterm(T,Y,A,K),terms);
		#Return the sum recurrence terms
		return add(RE)
	end proc:
	
ADEtermToREterm:= proc(term::algebraic,Y::anyfunc(name),A::anyfunc(name),K::list(name),$)::algebraic;
		local x,n,j,q,mterm,Ldiff,Lrec,xpow,c,Cauchyterm,i;
		option `Copyright (c) 2020 Bertrand Teguia T.`;
		description "Conversion of a differential term to a recurrence term.";
		x:=op(Y);
		n:=op(A);
		j:=PDEtools['difforder'](term,x);
		q:=degree(term,diff(Y,[x$j]));
		mterm:=term;
		Ldiff:=[];
		#collect the derivative data iteratively until it remains the monomial: C*x^l
		do
			Ldiff:=[op(Ldiff),`$`(j,q)];
			mterm:=subs(diff(Y,[x$j])^q=1,mterm);
			j:=PDEtools['difforder'](mterm,x);
			q:=degree(mterm,diff(Y,[x$j]))
		until type(mterm,polynom(anything,x));
		xpow:=degree(mterm,x);
		c:=coeff(mterm,x,xpow);
		#Apply the conversion formula
		Lrec:=map(r->poch(n+1,r)*subs(n=n+r,A), Ldiff);
		Cauchyterm:=Lrec[1];
		#Apply the Cauchy product formula for the derivative part
		for i from 2 to numelems(Lrec) do
			Cauchyterm:=ADECauchyprod(Cauchyterm,Lrec[i],n,K[i-1])
		end do;
		return c*subs(n=n-xpow,Cauchyterm)
	end proc:

ADECauchyprod := proc(t1::algebraic,t2::algebraic,n::name,k::name,$)::algebraic;
			return Sum(subs(n=k,t1)*subs(n=n-k,t2),k=0..n)
		end proc:
		
poch := proc(p::algebraic,k::nonnegint,$)::algebraic; 
		local j;
		return mul(p+j,j=0..k-1) 
	end proc: