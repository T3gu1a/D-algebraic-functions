
#Guessing D-algebraic functions

DalgFunGuess:= proc(L::list,
		      {degADE::posint:=2,
		      degPoly::nonnegint:=2,
		       devars::anyfunc(name):=NULL,
	         startfromord::nonnegint:=0},
			   $)::Or(identical(FAIL),`=`);
			option `Copyright (c) 2025 Bertrand Teguia T.`;
			description "Guessing D-algebraic functions (finding their differential equations)";
			local  Y::anyfunc(name),y::name,x::name,A::anyfunc(name),a::name,n::name,i::nonnegint,
			       c::nothing,N::posint,M::posint,V::list,j::nonnegint,Sinit::list(eq),k::name,
			       nL::posint:=numelems(L),Nmax:=ceil(nL/(degPoly+1)),ADE::algebraic,RE::algebraic,	
			       K::list,Eq::list(algebraic),NegInd::list,S::Or(identical(NULL),list(algebraic)),
			       NDE::algebraic,correct::boolean:=false,Arbconst::list,REcheck::algebraic,NRE::algebraic;
			#minimal deltak order for starting order startfromord
			N:=binomial(degADE+startfromord,degADE);
			#minimal number of unknown
			M:=(degPoly+1)*N;
			#too few initial values
			if M>nL then
				return FAIL
			end if;
			V:=[seq(c[i],i=0..M-1)];
			#ADE variables
			if devars=NULL then 
				`tools/genglobal`('y',{},'reset');
				`tools/genglobal`('x',{},'reset');
				y:=`tools/genglobal`('y');
				x:=`tools/genglobal`('x');
				Y:=y(x)
			else 
				y:=op(0,devars);
				x:=op(1,devars);
				#correctness of the input variables
				if y=x then
					error "invalid input: the variables of the differential equation must be distinct"
				end if;
				Y:=y(x)
			end if;
			#Recurrence variables
			a:=`tools/genglobal`('a');
			n:=`tools/genglobal`('n');
			interface(warnlevel=0);
			for j from 0 to degADE-1 do
				cat(k,j):=`tools/genglobal`('k')
			end do;
			K:=[seq(cat(k,j),j=0..degADE-1)];
			interface(warnlevel=4);
			A:=a(n);
			Sinit:=[seq(a(i-1)=L[i],i=1..nL)];
			#initialization - ADE and RE of the first iteration
			ADE:=expand(add(add(V[(degPoly+1)*(j-1)+i+1]*x^i,i=0..degPoly)*AnsatzDalg:-deltakdiff(Y,x,degADE,j),j=1..N));
			RE:=ADEtoRE(ADE,Y,A,K);
			#write the RE for non-negative indices
			#build the linear system and solve it
			Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=0..M-1)];
			#NegInd: list for substituting terms with negative indices to zero
			NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
			Eq:=subs(NegInd,Eq);
			#S:=try op(solve(Eq,V)) catch : NULL  end try;
			S:=SolveTools:-Linear(Eq,V,method='AlgebraicFunction');
			if S<>NULL then
				if remove(v->rhs(v)=0,S)={} or has(map(rhs,S),a) then
					S:=NULL
				else
					#checking the solution
					REcheck, S, correct:=checkSol(S,RE,NegInd,Sinit,M,nL,a,n)
				end if
			end if;
			while N<Nmax and not(correct) do
				N:=N+1;
				#print(N,Nmax);
				V:=[op(V),seq(c[i],i=M..M+degPoly)];
				NDE:=add(V[M+i]*x^(i-1)*AnsatzDalg:-deltakdiff(Y,x,degADE,N),i=1..degPoly+1);
				ADE:=NDE+ADE;
				NRE:=ADEtoRE(NDE,Y,A,K);
				RE:=NRE+RE;
				Eq:=[seq(Eq[i+1]+subs(Sinit,eval(NRE,[n=i,Sum=add])),i=0..M-1)];
				Eq:=[op(Eq),seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=M..M+degPoly+1)];
				M:=M+degPoly+1;
				#Eq:=[seq(subs(Sinit,eval(RE,[n=i,Sum=add])),i=0..M-1)];
				NegInd:=map(v->v=0,[op(indets(Eq,a(negint)))]);
				Eq:=subs(NegInd,Eq);
				#print("ok1");
				S:=SolveTools:-Linear(Eq,V,method='AlgebraicFunction');
				#print("ok2");
				if S<>NULL then
					if remove(v->rhs(v)=0,S)={} or has(map(rhs,S),a) then
						S:=NULL
					else
						#verification step
						REcheck, S, correct:=checkSol(S,RE,NegInd,Sinit,M,nL,a,n)
					end if
				end if
			end do;
			if correct then
				ADE:=subs(S,ADE);
				Arbconst:=sort([op(remove(has,indets(REcheck),a) minus {n,op(K)})]);
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

$include <NLDE/DalgFunGuess/ADEtoRE/src/ADEtoRE.mm>
$include <NLDE/DalgFunGuess/CheckSol/src/CheckSol.mm>