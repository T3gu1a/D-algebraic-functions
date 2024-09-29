

REorders := proc(RE::Or(algebraic,`=`),a::anyfunc(name),$)::(integer,integer);
		local n::name,aterms::list,A;
		option `Copyright (c) 2024 Bertrand Teguia T.`;
		n:=op(1,a);
		A:=op(0,a);
		aterms:=indets(RE,A(`+`)) union indets(RE,A(name));
		aterms:=subs(n=0,map(r-> op(r), aterms));
		return (max(aterms),min(aterms))
	end proc: