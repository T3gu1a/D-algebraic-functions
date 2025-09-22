#Guessing D-alg sequences

#The D-algebraic guesser
DalgGuess := proc(L::list,
	    {degADE::posint:=2,
	     revars::anyfunc(name):=NULL,
       startfromord::nonnegint:=0},
		 $)::Or(identical(FAIL),list(equation));
		local  k:=degADE,A,a,n,N::posint:=k+1,M::nonnegint,c::nothing,Arbconst::list,
		       i::nonnegint,V::list:=[seq(c[i],i=0..N-1)],RE::algebraic,REcheck::algebraic,
		       Eq::list(algebraic),S::Or(identical(NULL),list(algebraic)),Sinit::list(algebraic),
		       correct::boolean:=false,termfree::truefalse:=true,nL::nonnegint:=numelems(L);
		option `Copyright (c) 2025 Bertrand Teguia T.`;
			
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
		Sinit:=[seq(a(i-1)=L[i],i=1..nL)];			
		#initialization - DE and RE of the first iteration
		RE:=add(V[i+1]*deltakshift(A,n,k,i),i=0..N-1);
		#build the linear system
		Eq:=[seq(subs(Sinit,eval(RE,n=i)),i=0..N-1)];
		#too few initial values
		termfree:=evalb(indets(Eq,a('integer'))={});
		if not(termfree) then
			return FAIL
		end if;
		#solve the linear system
		S:=try op(solve(Eq,V)) catch : NULL  end try;
		if S<>NULL then
			if remove(v->rhs(v)=0,S)=[] then
				S:=NULL
			else
				#checking the solution
				REcheck, S, correct:=checkingguess(S,RE,Sinit,N,nL,A)
			end if
		end if;
		N:=N+1;
		while termfree and not(correct) do
			V:=[op(V),c[N-1]];
			RE:=RE+V[N]*deltakshift(A,n,k,N-1);
			Eq:=[seq(subs(Sinit,eval(RE,n=i)),i=0..N-1)];
			termfree:=evalb(indets(Eq,a('integer'))={});
			S:=try op(solve(Eq,V)) catch : NULL  end try;
			if S<>NULL then
				if remove(v->rhs(v)=0,S)=[] then
					S:=NULL
				else
					#verification step
					REcheck, S, correct:=checkingguess(S,RE,Sinit,N,nL,A)
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
	
$include <NLDE/DalgSeq/DalgGuess/SubProcedures/src/deltakshift.mm>
$include <NLDE/DalgSeq/DalgGuess/SubProcedures/src/checkingguess.mm>
$include <NLDE/DalgSeq/DalgGuess/SubProcedures/src/shiftkstart.mm>
$include <NLDE/DalgSeq/DalgGuess/modDalgGuess/src/modDalgGuess.mm>