
#A bijection between a subset of N^k and N is used 
#to define the ordering for differential monomials
ordertoktuple := proc(k::posint,n::nonnegint) option remember;
		local m, im, Tkn, j;
		option `Copyright (c) 2022 Bertrand Teguia T.`;
		if n=0 then 
		    return [seq(0,j=1..k)]
		else   
		    Tkn:=ordertoktuple(k,n-1);
		    m:=min(Tkn);
		    im:=ListTools:-Search(m,Tkn);
		    Tkn[im]:=m+1;
		    if im=k then
			return Tkn
		    else
			return [op(Tkn[1..im]),seq(0,j=im+1..k)]
		    end if     
		end if    
	end proc:
	