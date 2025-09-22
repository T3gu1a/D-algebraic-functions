
#The D-algebraic sequence guesser (finite field)
modDalgGuess := proc(L::list,
		    {degADE::posint:=2,
		     revars::anyfunc(name):=NULL,
	       startfromord::nonnegint:=0,
		    modulus::posint:=7},
		 $)::Or(identical(FAIL),list(equation));
		local  k:=degADE,A,a,n,N::posint:=k+1,M::nonnegint,c::nothing,Arbconst::list,
		       i::nonnegint,V::list:=[seq(c[i],i=0..N-1)],RE::algebraic,REcheck::algebraic,
		       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),Sinit::list(algebraic),
		       correct::boolean:=false,termfree::truefalse:=true,nL::nonnegint:=numelems(L),
		       Terms::list(integer), Meqs, beqs, rank::name;
		option `Copyright (c) 2025 Bertrand Teguia T.`;
			
		Terms:= try map(term -> term mod modulus, L) catch : [] end try;
		if Terms = [] then
			return FAIL
		end if;
		#ADE variables
		if revars=NULL then 
			`tools/genglobal`('a',{},'reset');
			`tools/genglobal`('n',{},'reset');
			a:=`tools/genglobal`('a');
			n:=`tools/genglobal`('n');
			A:=a(n)
		else 
			a:=op(0,revars);
			n:=op(1,revars);
			#correctness of the input variables
			if a=n then
				error "invalid input: the variables of the differential equation must be distinct"
			end if;
			A:=a(n)
		end if;

		M:=shiftkstart(startfromord,k);
		if M>N then
			V:=[op(V),seq(c[i],i=N..M-1)];
			N:=M
		end if;
		
		#initial values
		Sinit:=[seq(a(i-1)=Terms[i],i=1..nL)];			
		#initialization - DE and RE of the first iteration
		RE:=add(V[i+1]*deltakshift(A,n,k,i),i=0..N-1);
		#build the linear system
		Eq:=[seq(subs(Sinit,eval(RE,n=i) mod modulus) mod modulus,i=0..N-1)];
		#too few initial values
		termfree:=evalb(indets(Eq,a('integer'))={});
		if not(termfree) then
			return FAIL
		end if;
		#solve the linear system
		Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
		S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
		S:= ifelse(type(S,list(algebraic)),S,NULL);
		#S:=try msolve(Eq,modulus) catch : NULL  end try;
		if S<>NULL then
			if remove(v->v=0,S)=[] then
				S:=NULL
			else
				#checking the solution
				S:=[seq(V[i]=S[i],i=1..N)];
				REcheck, S, correct:=checkingguessmod(S,RE,Sinit,N,nL,A,modulus)
			end if
		end if;
		N:=N+1;
		while termfree and not(correct) do
			V:=[op(V),c[N-1]];
			RE:=RE+V[N]*deltakshift(A,n,k,N-1);
			Eq:=[seq(subs(Sinit,eval(RE,n=i) mod modulus) mod modulus,i=0..N-1)];
			termfree:=evalb(indets(Eq,a('integer'))={});
			#S:=try msolve(Eq,modulus) catch : NULL  end try;
			Meqs, beqs := LinearAlgebra:-GenerateMatrix(Eq,V);
			S:= try convert(Linsolve(Meqs,beqs) mod modulus, list) catch : NULL end try;
			S:= ifelse(type(S,list(algebraic)),S,NULL);
			if S<>NULL then
				if remove(v->v=0,S)=[] then
					S:=NULL
				else
					#verification step
					S:=[seq(V[i]=S[i],i=1..N)];
					REcheck, S, correct:=checkingguessmod(S,RE,Sinit,N,nL,A,modulus)
				end if
			end if;
			N:=N+1
		end do;
		if correct then
			RE:=REcheck;
			Arbconst:=sort([op(remove(has,indets(RE),a) minus {n})]);
			RE:=collect(RE,Arbconst,'distributed');
			if Arbconst <> [] then
				`tools/genglobal`('_C',{},'reset');
				Arbconst:=map(v->v=`tools/genglobal`('_C'),Arbconst);
				RE:=subs(Arbconst,RE)
			end if;
			return RE=0
		else
			return FAIL
		end if
	end proc:
	
#code to verify the guessed equation with the unused terms
checkingguessmod := proc(Sol::Or(set,list),
		       REsol::algebraic,
		       Sinit::list,
		           M::posint,
		          nL::nonnegint,
		           A::anyfunc(name),
		           m::posint,
		           $)
		local S::list, RE::algebraic, checkL::list, checkset::set,i::nonnegint,
		      ord::nonnegint,maxind::nonnegint,n:=op(A);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		S:=map(simplify,Sol);
		RE:=subs(S,REsol);
		ord:=REorders(RE,A)[1];
		maxind:=nL-ord-1;
		checkset:={seq(simplify(subs(Sinit,eval(RE,n=i) mod m)) mod m,i=maxind..M,-1)};
		return RE, S, evalb(checkset = {0})
	end proc: