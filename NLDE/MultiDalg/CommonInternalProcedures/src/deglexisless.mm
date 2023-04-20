
#degree lexicographic comparisons of fixed-length-k lists
deglexisless:=proc(l1::list(nonnegint),l2::list(nonnegint),$)::truefalse;
		local n1::nonnegint,n2::nonnegint;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		if add(l1) = add(l2) then
			n1 := parse(cat(op(ListTools:-Reverse(l1))));
			n2 := parse(cat(op(ListTools:-Reverse(l2))));
			return ifelse(n1<n2,true,false)
		else
			return ifelse(add(l1)<add(l2),true,false)
		end if
	end proc: