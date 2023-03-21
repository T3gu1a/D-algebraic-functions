
composeDalg:= proc(L::[`=`,`=`],
		   V::[anyfunc(name),anyfunc(name)],
	           z::anyfunc(name),
		   {ordering::identical(plex,lexdeg):=plex},
	           $)::`=`;
		local t::name:=op(z),n::posint,fgh::set(`=`),R,f::nothing,Sys::list,x::nothing,
		      Sysh::list,g::nothing,h::nothing,j::nonnegint,m::posint,Subg::list,
		      DE1::`=`,DE2::`=`;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		description "Compose two D-algebraic functions from their differential   "
			    "equations given in the same order with the second argument, "
			    "representing the list of dependent variables.               "
			    "INPUT: - a list of two differential equations [DE(f),DE(g)] "
			    "       - a list of their dependent variables  [f(t),g(t)]   "
			    "	    - the dependent variable of the output h(t)          "
			    "OUPUT: a differential equations for the composition f(g)    "
			    "       of the second g by the first f.                      ";
		DE1:=lhs(L[1])-rhs(L[1]);
		DE2:=lhs(L[2])-rhs(L[2])=0;
		n:=PDEtools:-difforder(DE1,t);
		#substitution to remove derivatives and t from the first differential equation (DE)
		R:=subs(t=g[0],subs([seq(diff(V[1],[t$j])=f[j],j=0..n)],DE1));
		#in case of non-l.h.o, one derivation is needed
		if degree(R,f[n])>1 then
			n:=n+1;
			R:=subs(t=g[0],subs([seq(diff(V[1],[t$j])=f[j],j=0..n)],diff(DE1,t)))
		end if;
		#build the system with the second DE
		Sys:=buildsystem(DE2,V[2],g);
		m:=numelems(Sys[1]);
		#expressing derivatives of f in terms of those of h and g
		fgh:=ftogh(f,g,h,x,n);
		if m<=n then
			#build a list of further substitutions in R if m<=n
			Subg:=subsgfurther(Sys[1][m],g,t,m,n);
			fgh:=subs(Subg,fgh)
		end if;
		#the full substitution (f is already eliminated here)
		R:=normal(subs(fgh,R));
		#the system for h is easily obtained from R
		Sysh:=[[seq(h[j],j=1..(n-1)),solve(R,h[n])],[seq(h[j],j=0..(n-1))]];
		#use SystoMinDiffPoly to return the desired output
		return SystoMinDiffPoly([op(Sys[1]),op(Sysh[1])],h[0],[op(Sys[2]),op(Sysh[2])],z,':-ordering'=ordering)
	end proc:

$include <NLDE/composeDalg/SubProcedures/src/subsgfurther.mm>
$include <NLDE/composeDalg/SubProcedures/src/ftogh.mm>
