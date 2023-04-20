
#deglex-sort of an integer partition
sortedkpartition := proc(k::posint,s::nonnegint,$)::list(list(nonnegint));
		local P::list(nonnegint);
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		P:=select(l->numelems(l)<=k,combinat:-partition(s));
		P:=map(l->op(combinat:-permute([op(l),$(0,k-numelems(l))],k)),P);
		return sort(P,(a,b)->deglexisless(a,b))
	end proc: