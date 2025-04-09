

mergesystem:= proc(L::list(`=`),
		   V::list(anyfunc(name)),
		  $)::`=`;
		local l::posint:=numelems(L), j::posint, Sys::list, vars::list, deriv::list, 
		      n::posint, x::nothing, X::list, i::posint, Ind::list;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		Sys:=[seq(buildsystem(L[j],V[j],cat(x,j)),j=1..l)];
		vars:=map(r->op(r[2]),Sys);
		deriv:=map(r->op(r[1]),Sys);
		n:=numelems(vars);
		X:=[seq(vars[j]=x[j],j=1..n)];
		Ind:=[seq(1+add(numelems(Sys[i][2]),i=1..(j-1)),j=1..l)];
		[subs(X,deriv),map(r->x[r],Ind),map(rhs,X)]
	end proc:
	
$include <NLDE/DalgSeq/AnsatzDalgSeq/CommonInternalProcedures/src/buildsystem.mm>