
#inverse of the Cantor k-tuple  function
CantorInvSigma := proc(k::posint,n::nonnegint,$) option remember;
		local N::posint,  s::posint, L::list(list(nonnegint)), j::integer;
		option `Copyright (c) 2023 Bertrand Teguia T.`;
		if assigned(thetakTuple[k]) then
			N:=numelems(thetakTuple[k])-1;
			if n <= N then
				return thetakTuple[k][n]
			else
				s:=add(thetakTuple[k][N])+1;
				L:=sortedkpartition(k,s);
				assign(seq(thetakTuple[k][N+j]=L[j],j=1..numelems(L)));
				while numelems(thetakTuple[k]) < n+1 do:
					N:=numelems(thetakTuple[k])-1;
					s:=s+1;
					L:=sortedkpartition(k,s);
					assign(seq(thetakTuple[k][N+j]=L[j],j=1..numelems(L)))
				end do;
				return thetakTuple[k][n]
			end if
		else
			thetakTuple[k]:=table([]);
			thetakTuple[k][0]:=[$(0,k)];
			s:=1;
			while numelems(thetakTuple[k])<n+1 do:
				L:=sortedkpartition(k,s);
				N:=numelems(thetakTuple[k])-1;
				assign(seq(thetakTuple[k][N+j]=L[j],j=1..numelems(L)));
				s:=s+1
			end do;
			return thetakTuple[k][n]
		end if
	end proc: