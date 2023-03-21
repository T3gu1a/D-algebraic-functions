
#Derivative operator used to compute product of k derivatives
#with respect to a certain 'natural ordering'
deltakdiff := proc(expr,z::name,k::posint:=2,n::nonnegint:=1,$)
		local tuple;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		tuple:= select(type,ordertoktuple(k,n)-~1,nonnegint);
		return mul(map(d-> diff(expr,[seq(z,1..d)]), tuple))
	end proc:

$include <NLDE/AnsatzDalg/CommonInternalProcedures/src/ordertoktuple.mm>
	