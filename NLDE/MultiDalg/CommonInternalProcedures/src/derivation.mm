
#the derivation
derivation := proc(expr,Z::list(name),n::nonnegint,$)
		local tuple, k::posint:=numelems(Z),j;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		tuple:=CantorInvSigma(k,n);
		return diff(expr,[seq(Z[j]$tuple[j],j=1..k)])
	end proc: