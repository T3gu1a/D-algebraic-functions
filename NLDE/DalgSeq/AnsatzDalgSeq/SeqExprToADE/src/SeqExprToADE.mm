
SeqExprToADE := proc(expr::algebraic,
		        F::anyfunc(name),
	      {maxdeorder::posint:=10,
	     startfromord::posint:=0,
		   degADE::posint:=1},
		       $)::Or(`=`,identical(FAIL));
		return DegreekDE(expr,F,[],':-maxdeorder'=maxdeorder,degreeDE=degADE,':-startfromord'=startfromord)
	end proc: