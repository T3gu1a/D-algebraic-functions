
FixedOrdDegFunGuess:= proc(Sinit::list,
			  degADE::posint,
			 degPoly::nonnegint,
			       Y::anyfunc(name),
			       A::anyfunc(name),
			       N::nonnegint,
			       y::name,
			       x::name,
			       a::name,
			       n::name,
			       K::list,
		       linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense),
		    maxIteration::posint,
		  inputConstants::set(name),
			      $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,
		       nL::posint:=numelems(Sinit),hasterm::truefalse:=false,ADE::algebraic,RE::algebraic,	
		       Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),randpick,val,
		       correct::truefalse:=false,Arbconst::list,REcheck::algebraic,Meqs,beqs,ul::list,tl,
		       l::list(nonnegint),Ll,m::nonnegint,degCoeffs::list(nonnegint),freqs,total_perms;
		
		l:=GenMaxlistnumber(N,degPoly,max(nL-N,0));
		ul:=sort([op({op(l)})]);
		tl:=Statistics:-Tally(l);
		freqs:=Array([seq(eval(val, tl), val = ul)], datatype = integer);
		total_perms:=MultinomialCount(freqs);	
		while l<> FAIL and not(correct) do
			if total_perms < maxIteration then
				Ll:=Iterator:-Permute(l,N);
				M:=add(l)+N;
				for degCoeffs in Ll do:
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					RE:=ADEtoRE(ADE,Y,A,K);
					Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=0..M-1)];
					NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
					Eq:=subs(NegInd,Eq);
					#Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq, V);
					#S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
					#S:=ifelse(S<>NULL,[seq(V[i] = S[i], i = 1 .. numelems(V))],NULL);
					S:=SolveTools:-Linear(Eq,V,method=linsolver);
					if S<>NULL then
						if remove(v->rhs(v)=0,S)={} or has(map(rhs,S),a) then
							hasterm:=has(Eq,a);
							S:=NULL
						else
							REcheck, S, correct:=checkSol(S,RE,NegInd,Sinit,M,nL,a,n)
						end if
					end if;
					if correct or hasterm then
						break
					end if
				end do
			else
				randpick:=rand(1..total_perms);
				degCoeffs:=copy(l);
				M:=add(l)+N;
				to maxIteration do
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					RE:=ADEtoRE(ADE,Y,A,K);
					Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=0..M-1)];
					NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
					Eq:=subs(NegInd,Eq);
					Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq, V);
					S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
					S:=ifelse(S<>NULL,[seq(V[i] = S[i], i = 1 .. numelems(V))],NULL);
					#S:=SolveTools:-Linear(Eq,V,method=linsolver);
					if S<>NULL then
						if remove(v->rhs(v)=0,S)=[] or has(map(rhs,S),a) then
							hasterm:=has(Eq,a);
							S:=NULL
						else
							REcheck, S, correct:=checkSol(S,RE,NegInd,Sinit,M,nL,a,n)
						end if
					end if;
					if correct or hasterm then
						break
					end if;
					degCoeffs:=UnrankMultiset(randpick(), Array(ul), freqs, N)
				end do
			end if;
			l:=prevlistnumber(degPoly,l);
			hasterm:=false
		end do;
		if correct then
			ADE:=subs(S,ADE);
			Arbconst:=sort([op(remove(has,indets(REcheck),a) minus ({n,op(K)} union inputConstants))]);
			if Arbconst <> [] then
				`tools/genglobal`('_C',{},'reset');
				Arbconst:=map(v->v=`tools/genglobal`('_C'),Arbconst);
				ADE:=subs(Arbconst,ADE)
			end if;
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))},'distributed');
			return ADE=0
		else
			return FAIL
		end if     
			       
	end proc:
	
FFixedOrdDegFunGuess:= proc(   Lf::algebraic,
			    Sinit::list,
				a::name,
			   degADE::posint,
			  degPoly::nonnegint,
				Y::anyfunc(name),
				N::nonnegint,
				y::name,
				x::name,
			linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense),
		     maxIteration::posint,
		   inputConstants::set(name),
			 sparsity::fraction,
			       $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2026 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,randpick,
		       nL::posint:=degree(Lf,x)+1,ADE::algebraic,polEq::algebraic,total_combs::posint,	
		       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),pcentge,
		       correct::truefalse:=false,Arbconst::list,ADEcheck::algebraic,nleft,
		       MnL::posint,ZerosV,zV,zzV::list,unkV::list,idx,Meqs,beqs,pool;
		M:=(degPoly+1)*N;
		pcentge:=max(sparsity,1-nL/M);
		MnL:=ceil(pcentge*M);
		total_combs:=binomial(M,MnL);
		V:=[seq(c[i],i=0..M-1)];
		if total_combs < maxIteration then
			ZerosV:=Iterator:-Combination(M, MnL);
			for zV in ZerosV do
				zzV := [seq(c[idx], idx in zV)];
				zzV:=map(t->t=0,zzV);
				unkV:=subs(zzV,V);
				ADE:=add(add(unkV[(degPoly+1)*(j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
				polEq:=expand(eval(ADE,Y=Lf));
				polEq:=subs(Sinit,polEq);
				unkV:=remove(t->t=0,unkV);
				Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..numelems(unkV)]; 
				#Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq, unkV);
				#S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
				#S:=ifelse(S<>NULL,[seq(unkV[i] = S[i], i = 1 .. numelems(unkV))],NULL);
				S:=SolveTools:-Linear(Eq,unkV,method=linsolver);
				if S<>NULL then
					if remove(v->rhs(v)=0,S)={} then
						S:=NULL
					else
						ADEcheck, S, correct:=polcheckSol(S,ADE,Sinit,a,nL,y,x)
					end if
				end if;
				if correct then
					break
				end if
			end do
		else
			randpick:=rand(1..M);
			to maxIteration do
				pool:=table([seq(i=i,i=1..M)]);
				nleft:=M;
				if 2*MnL<M then
					zV:=Array(1..MnL);
					for j from 1 to MnL do
						idx:=rand(1..nleft)();
						zV[j]:=pool[idx];
						pool[idx]:=pool[nleft];
						unassign('pool[nleft]');
						nleft:=nleft-1
					end do
				else
					zV:=Array(1..M-MnL);
					for j from 1 to M-MnL do
						idx:=rand(1..nleft)();
						zV[j]:=pool[idx];
						pool[idx]:=pool[nleft];
						unassign('pool[nleft]');
						nleft:=nleft-1
					end do;
					zV:=convert({seq(1..M)} minus convert(zV,set),list)
				end if;
				zzV:=[seq(c[idx], idx in zV)];
				zzV:=map(t->t=0,zzV);
				unkV:=subs(zzV,V);
				ADE:=add(add(unkV[(degPoly+1)*(j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j),i=0..degPoly),j=1..N);
				polEq:=expand(eval(ADE,Y=Lf));
				polEq:=subs(Sinit,polEq);
				unkV:=remove(t->t=0,unkV);
				Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..numelems(unkV)]; 
				Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq, unkV);
				S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
				S:=ifelse(S<>NULL,[seq(unkV[i] = S[i], i = 1 .. numelems(unkV))],NULL);
				if S<>NULL then
					if remove(v->rhs(v)=0,S)=[] then
						S:=NULL
					else
						ADEcheck, S, correct:=polcheckSol(S,ADE,Sinit,a,nL,y,x)
					end if
				end if;
				if correct then
					break
				end if		
			end do
			
		end if;
		if correct then
			ADE:=subs(S,ADE);
			Arbconst:=sort([op(indets(ADEcheck) 
				minus (inputConstants union {x,y,seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))}))]);
			if Arbconst <> [] then
				`tools/genglobal`('_C',{},'reset');
				Arbconst:=map(v->v=`tools/genglobal`('_C'),Arbconst);
				ADE:=subs(Arbconst,ADE)
			end if;
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))},'distributed');
			return ADE=0
		else
			return FAIL
		end if     
			       
	end proc:
	
FFixedOrdDegFunGuess2:= proc(  Lf::algebraic,
			    Sinit::list,
				a::name,
			   degADE::posint,
			  degPoly::nonnegint,
				Y::anyfunc(name),
				N::nonnegint,
				y::name,
				x::name,
			linsolver::identical(AlgebraicFunction,Rational,AlgebraicNumber,RadicalFunction,RationalDense),
		     maxIteration::posint,
			constants::set(name),
			       $)::Or(identical(FAIL),`=`);
		option `Copyright (c) 2026 Bertrand Teguia T.`;
		description "Looking for an equation among all possible equations of the given maximum polynomial degree";
		local  i::nonnegint,c::nothing,M::posint,V::list,j::nonnegint,total_perms::nonnegint,
		       nL::posint:=degree(Lf,x)+1,ADE::algebraic,polEq::algebraic,randpick,	
		       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),tl,freqs,val,
		       correct::truefalse:=false,Arbconst::list,ADEcheck::algebraic,ul::list,Meqs,beqs,
		       l::list(nonnegint),Ll::list(list),m::nonnegint,degCoeffs::list(nonnegint);
		l:=GenMaxlistnumber(N,degPoly,max(nL-N,0));
		while l<> FAIL and not(correct) do
			ul:=sort([op({op(l)})]);
			tl:=Statistics:-Tally(l);
			freqs:=Array([seq(eval(val, tl), val = ul)], datatype = integer);
			total_perms:=MultinomialCount(freqs);
			if total_perms < maxIteration then
				Ll:=Iterator:-Permute(l,N);
				M:=add(l)+N;
				for degCoeffs in Ll do
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					polEq:=expand(eval(ADE,Y=Lf));
					polEq:=subs(Sinit,polEq);
					Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; 
					#Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq, V);
					#S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
					#S:=ifelse(S<>NULL,[seq(V[i] = S[i], i = 1 .. numelems(V))],NULL);
					S:=SolveTools:-Linear(Eq,V,method=linsolver);
					if S<>NULL then
						if remove(v->rhs(v)=0,S)={} then
							S:=NULL
						else
							ADEcheck, S, correct:=polcheckSol(S,ADE,Sinit,a,nL,y,x)
						end if
					end if;
					if correct then
						break
					end if
				end do
			else
				randpick:=rand(1..total_perms);
				degCoeffs:=copy(l);
				M:=add(l)+N;
				to maxIteration do
					V:=[seq(c[i],i=0..M-1)];
					ADE:=add(add(V[add(degCoeffs[m]+1,m=1..j-1)+i+1]*x^i*AnsatzDalg:-deltakdiff(Y,x,degADE,j)
									  ,i=0..degCoeffs[j]),j=1..N);
					polEq:=expand(eval(ADE,Y=Lf));
					polEq:=subs(Sinit,polEq);
					Eq:=PolynomialTools:-CoefficientList(polEq,x)[1..M]; 
					Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq, V);
					S:=try LinearAlgebra:-LinearSolve(Meqs, beqs) catch: NULL end try;
					S:=ifelse(S<>NULL,[seq(V[i] = S[i], i = 1 .. numelems(V))],NULL);
					if S<>NULL then
						if remove(v->rhs(v)=0,S)=[] then
							S:=NULL
						else
							ADEcheck, S, correct:=polcheckSol(S,ADE,Sinit,a,nL,y,x)
						end if
					end if;
					if correct then
						break
					end if;
					degCoeffs:=UnrankMultiset(randpick(), Array(ul), freqs, N)
				end do
			end if;
			l:=prevlistnumber(degPoly,l)
		end do;
		if correct then
			ADE:=subs(S,ADE);
			Arbconst:=sort([op(indets(ADEcheck) 
				minus (constants union {x,y,seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))}))]);
			if Arbconst <> [] then
				`tools/genglobal`('_C',{},'reset');
				Arbconst:=map(v->v=`tools/genglobal`('_C'),Arbconst);
				ADE:=subs(Arbconst,ADE)
			end if;
			ADE:=collect(ADE,{seq(diff(Y,[x$i]),i=0..PDEtools:-difforder(ADE,x))},'distributed');
			return ADE=0
		else
			return FAIL
		end if     
	end proc:
	
prevlistnumber := proc(k::integer, L::list)
	local A, temp_L, n, m, i, j, current_sum, val;

	A := Array(ListTools:-Reverse(L),datatype=integer);
	m:=ArrayNumElems(A);
	n:=0;
	for i from 1 to m do n:=n+A[i] od;

	# Step 1: Find the rightmost index i that we can decrement
	# To keep it a 'partition' (descending), L[i] must remain >= L[i+1]
	for i from m-1 by -1 to 1 do
		if A[i] > 0 and A[i] > A[i+1] then
		    # Decrease this pivot
		    temp_L:=copy(A);
		    temp_L[i] := temp_L[i] - 1;
		    
		    # Step 2: Recalculate the sum used so far
		    current_sum := 0;
		    for j from 1 to i do current_sum:=current_sum+temp_L[j] od;
		    
		    # Step 3: Fill the rest greedily from left to right
		    # We must respect: digit <= k, digit <= previous_digit, and total sum <= n
		    for j from i+1 to m do
			# The digit cannot exceed k, the remaining budget, or the digit to its left
			val := min(k, n - current_sum, temp_L[j-1]);
			temp_L[j] := val;
			current_sum := current_sum + val
		    end do;
		    
		    if current_sum = n and has(temp_L,k) then
			return ListTools:-Reverse(convert(temp_L,list))
		    end if
		end if
	end do;

	# If we exit the loop, we've exhausted all partitions for this sum N.
	# Now we must drop the sum to N-1 and start from the Max again.
	if n > k then
		return GenMaxlistnumber(m, k, n-1);
	else
		return FAIL 
	end if
end proc:

MultinomialCount := proc(counts::Array)
	local total, den, c;
	total := 0;
	for c in counts do total := total + c; od;
	den := 1;
	for c in counts do den := den * factorial(c); od;
	factorial(total) / den;
end proc:

UnrankMultiset := proc(R_in, elements::Array, initial_counts::Array, N::integer)
    local R, P, i, j, counts, sub_total;
    R := R_in;
    P := Array(1..N, datatype=integer);
    counts := copy(initial_counts);
    
    for i from 1 to N do
	for j from 1 to ArrayNumElems(elements) do
	    if counts[j] > 0 then
		# If we put elements[j] at position i, 
		# how many unique ways to finish the rest?
		counts[j] := counts[j] - 1;
		sub_total := MultinomialCount(counts);
		
		if R <= sub_total then
		    # This is our element! Move to the next position i.
		    P[i] := elements[j];
		    break
		else
		    # R is further down the list; skip this group
		    R := R - sub_total;
		    # Restore count to try the next unique element for this position
		    counts[j] := counts[j] + 1
		end if
	    end if
	end do
    end do;
    P
end proc:

GenMaxlistnumber := proc(m::posint, k::integer, n::nonnegint)
	local digits, remaining_sum, i, current_digit;

	remaining_sum := n;
	digits := Array(1..m, datatype=integer);

	for i from 1 to m do
		# Greedy choice: take as much as possible, up to k
		current_digit := min(k, remaining_sum);
		digits[i] := current_digit;
		remaining_sum := remaining_sum - current_digit
	end do;

	# Convert list of digits to a single integer
	ListTools:-Reverse(convert(digits,list))
end proc:
	
