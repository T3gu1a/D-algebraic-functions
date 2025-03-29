#shift operator used to compute product of shifts
#with respect to the ordering from AnsatzDalg:-ordertoktuple
deltakshift := proc(expr,z::name,k::posint:=2,n::nonnegint:=1,$)
		local tuple;
		option `Copyright (c) 2025 Bertrand Teguia T.`;
		tuple:= select(type,AnsatzDalg:-ordertoktuple(k,n)-~1,nonnegint);
		return mul(map(d-> LREtools:-shift(expr,z,d), tuple))
	end proc: