
#arithmetic of D-algebraic sequences: addition, product, ratio; more generaly, rational functions of
#D-algebraic sequences

arithmeticDalgSeq:=proc(L::list(`=`),
			V::list(anyfunc(name)),
			z::name=ratpoly,
		       $)::`=`;
			local t:=op(1,V[1]),DEs::list(`=`),Sys::list,j::posint,subV::list;
			option `Copyright (c) 2024 Bertrand Teguia T.`;
			description "Analog of arithmeticDalg for sequences";
			if numelems(L)=1 then
				return L
			end if;
			DEs:=map(r->lhs(r) - rhs(r)=0,L);
			#build the systems and merge them using mergesystemseq
			Sys:=mergesystemseq(DEs,V);
			#prepare the list for the change of variables 
			#in r according to Sys
			subV:=[seq(op(0,V[j])=Sys[2][j],j=1..numelems(V))];
			#use SystoDE to return the desired output
			return SystoDE(Sys[1],subs(subV,rhs(z)),Sys[3],lhs(z)(t))
		end proc: