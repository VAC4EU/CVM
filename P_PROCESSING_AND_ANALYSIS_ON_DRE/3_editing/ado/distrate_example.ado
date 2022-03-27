*! version 1.0.0  10mar2006
program distrate_example
	if (_caller() < 8) {
		di as err "This example requires version 8"
		exit 198
	}
	if (_caller() < 8.2)  version 8
	else		      version 8.2
	gettoken dsn 0 : 0, parse(" :")
	gettoken null 0 : 0, parse(" :")
	di as txt
	di as txt "-> " as res "preserve"
	preserve
	di as txt
	cap findfile SuffolkCounty.dta
	if _rc {
		di as err "file SuffolkCounty.dta not found"
		exit 601 
	}
	local fileful `"`r(fn)'"'
	cap use `"`fileful'"'
	if _rc>900 { 
		window stopbox stop ///   
		"Dataset used in this example" ///
		"too large for Small Stata"
		exit _rc 
	}
	di 
	di as txt "-> " as res "Example Direct Standardization" 
	di as txt "-> " as res "use `dsn', clear"
	de
	di
	collapse (sum)  deaths pop,by(cindpov agegr)
	di as txt "-> " as res "collapse (sum)  deaths pop,by(cindpov agegr)"
	di
	di as txt "-> " as res "distrate deaths pop using year2000st, stand(agegr) by(cindpov) mult(100000)"
	di
	distrate deaths pop using year2000st, stand(agegr) by(cindpov) mult(100000) 
	di
	di as txt "Further options"
	di
	di as txt "-> " as res "distrate deaths pop using year2000st, stand(agegr) by(cindpov) dobson" ///
		" saving(DirectSuffolk,replace) format(%8.2f) mult(100000) level(90) list(rateadj lb_f ub_f lb_d ub_d)"
	di
	distrate deaths pop using year2000st, stand(agegr) by(cindpov) dobson ///
		saving(DirectSuffolk,replace) format(%8.2f) mult(100000) level(90) list(rateadj lb_f ub_f lb_d ub_d)    
	di
	di as txt "-> " as res `"restore"'
end
