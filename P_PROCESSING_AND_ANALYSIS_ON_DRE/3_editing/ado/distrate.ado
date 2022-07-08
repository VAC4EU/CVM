*! version 1.0.1 10mar2006
*! Direct standardized rates with improved confidence intervals
program define distrate, rclass
	version 8.0
	syntax varlist(min=2 max=2) [if] [in] using/ , STANDstrata(varlist) ///
                  [ BY(varlist) POPSTAND(name) DOBSON MULT(numlist integer max=1) Format(string) /// 
                    LIst(namelist) SAving(string asis) Level(integer $S_level) ]

/* Relevant sample */
	marksample touse
	markout `touse' `by' `standstrata', strok
	tokenize `varlist'
	local death `1'
	local studypop `2'
	cap assert `death' <= `studypop'
	if _rc {
		di as err "Number events in `death' cannot be greater than population specified in `studypop'."
		exit 459 
	}
	confirm file `"`using'.dta"'

/* Check that format is valid */
	if "`format'" != "" {
		if index("`format'",",") local format = subinstr("`format'", "," , "." , 1) /* european numeric format */
		local fmt = substr("`format'",index("`format'",".")-1,3) 
		capture {
			assert substr("`format'",1,1)=="%" & substr("`format'",2,1)!="d" ///
				& substr("`format'",2,1)!="t" & index("`format'","s")==0
			confirm number `fmt'
		}
		if _rc {
			di as err "invalid format. format has been set to default %6.4f"
		}
	}

	if "`saving'" != "" {
		gettoken stfile saving : saving, parse(",")
		gettoken comma saving  : saving, parse(",") 
		if `"`comma'"' == "," { 
			gettoken outrgro saving : saving, parse(" ,")
			gettoken comma saving : saving, parse(" ,")
			if `"`outrgro'"' != "replace" | `"`comma'"'!="" { 
				di as err "option saving() invalid"
			exit 198
			}
		}
		else if `"`comma'"' != "" {
			di as err "option saving() invalid"
			exit 198
		}
		else 	confirm new file `"`stfile'.dta"'
	}
	preserve
	qui keep if `touse'
	qui keep `death' `studypop' `standstrata' `by'
	if "`by'"== "" {
		tempvar by
		g byte `by' = 1
		local nby y
	}

	if "`mult'"!="" local mult "*`mult'"
	tempvar ckstr rate sumpop 
	g double `rate' = `death'/`studypop'

	if "`popstand'"=="" {
		local popstand `studypop'
		tempvar stpop
		g `stpop' = `studypop'
		drop `studypop'
		local studypop `stpop'
	}

/* Check that all standstrata strata are present in each by level */
	bysort `by' : g `ckstr' = _N
	cap assert `ckstr' = `ckstr'[_n-1] if _n>1
	if _rc {
		fillin `by' `standstrata'
		qui replace `rate'=0 if `rate'==.
		drop _fillin
	}

	sort `standstrata'
	qui merge `standstrata' using "`using'"
	
/* Print a warning message if any records do not match with standard population file and exit */
	qui count if _merge==1
	if r(N) {
		di in red "`r(N)' records fail to match with standard population file " ///
			"(records who do not match are saved to _merge_error_.dta)."
		qui keep if _merge==1
		qui save _merge_error_.dta, replace
		exit 459 
	}
	
	tempvar w var wfay vardob yl yu 
	qui {
		replace `rate' = 0 if _m==2
		drop _m
		egen `sumpop' = sum(`popstand'), by(`by') 
		g double `w' = `popstand'/ `sumpop'
		g double rateadj = `rate'*`w' 
		g double `wfay' = `w'/`studypop'
		g double `var' = `death' * (`w'/`studypop')^2
*		g double `vardob'=`popstand'^2*`rate'/`studypop'
		collapse (sum) `death' `studypop' `popstand' rateadj `var' (max) `wfay', by(`by')
		g double crude = `death'/`studypop' `mult' 
		g double lb_fay = (`var'/(2*rateadj))*invchi2(2*rateadj^2/`var',1-(.5+`level'/200)) `mult'
		replace lb_fay = 0 if `death'== 0 
		g double ub_fay = ((`var' + `wfay'^2)/(2*(rateadj+`wfay'))) * ///
			invchi2((2*(rateadj+`wfay')^2)/(`var'+`wfay'^2),.5+`level'/200) `mult'
*		replace `vardob' = `vardob'/`popstand'^2
		g double `yl' = .
		g double `yu' = .
		count
		local n = r(N)
		forval i = 1/`r(N)' {
			_crccip `death'[`i'] `level'
			replace `yl' = `r(lower)' if _n == `i'
			replace `yu' = `r(upper)' if _n == `i'
		}
		g double ub_dob= (rateadj+(`yu'-`death')*sqrt(`var'/`death')) `mult' 
		g double lb_dob= (rateadj+(`yl'-`death')*sqrt(`var'/`death')) `mult'
		replace rateadj = rateadj `mult'
	}
	if "`nby'"!="" local by
	g double se_F = sqrt(`var')`mult'

	if "`format'" != "" format crude rateadj lb_f ub_f lb_dob ub_dob se_F `format'
	
	tempname NDeath Nobs Crude Adj Se_F Ub_F Lb_F Ub_D Lb_D 
	mat `NDeath'=J(1,`n',0)
	mat `Nobs'=J(1,`n',0)
	mat `Crude'=J(1,`n',0)
	mat `Adj'=J(1,`n',0)
	mat `Se_F'=J(1,`n',0)
	mat `Ub_F'=J(1,`n',0)
	mat `Lb_F'=J(1,`n',0)
	mat `Ub_D'=J(1,`n',0)
	mat `Lb_D'=J(1,`n',0)
	forval i = 1 / `n' {
			mat `NDeath'[1,`i']=`death'[`i']
			mat `Nobs'[1,`i']=`studypop'[`i']
			mat `Crude'[1,`i']=  cond(crude[`i']==.,9,crude[`i']) 
			mat `Adj'[1,`i']= cond(rateadj[`i']==.,9,rateadj[`i'])
			mat `Se_F'[1,`i']= cond(`var'[`i']==.,9,sqrt(`var'[`i'])) 
			mat `Ub_F'[1,`i']= cond(ub_fay[`i']==.,9,ub_fay[`i']) 
			mat `Lb_F'[1,`i']= cond(lb_fay[`i']==.,9,lb_fay[`i']) 
			mat `Ub_D'[1,`i']= cond(ub_dob[`i']==.,9,ub_dob[`i']) 
			mat `Lb_D'[1,`i']= cond(lb_dob[`i']==.,9,lb_dob[`i']) 
	}
	return mat NDeath `NDeath'
	return mat Nobs `Nobs'
	return mat crude `Crude'
	return mat adj `Adj'
	return mat se_F `Se_F'
	return mat ub_F `Ub_F'
	return mat lb_F `Lb_F'
	return mat ub_D `Ub_D'
	return mat lb_D `Lb_D'
	rename `studypop' N
	if "`dobson'"!=""{
		local dd "and as proposed by Dobson et al."
		local dlist "lb_dob ub_dob"
	}
	di as txt "Directly standardized rates"
	di as txt "CI based on the gamma distibution (Fay and Feuer) `dd'"
	if "`list'"=="" list `by' `death' N crude rateadj lb_f ub_f se_F `dlist', table noobs 
	else {
		foreach name of local list {
			cap confirm var `name'
			if _rc    di as err "WARNING: `name' invalid or ambiguous in list option" 
			else	  unab ilist: `name'
			local tolist "`tolist' `ilist'"
		}
		if "`by'" != "" local tolist : list tolist - by
		list `by' `tolist', table noobs
	}
	if "`stfile'" != "" {
		di
		keep `by' `death' N crude rateadj lb_f ub_f lb_dob ub_dob se_F
		label var `death'	"Events"
		label var N		"Study Population"
		label var crude		"Crude Rate"
		label var rateadj	"Direct Adjusted Rate"
		label var lb_f		"Lower  Fay and Feuer CI" 
		label var ub_f		"Higher Fay and Feuer CI" 
		label var lb_dob	"Lower  DKES (Dobson et al) CI" 
		label var ub_dob	"Higher DKES (Dobson et al) CI" 
		label var se_F		"Standard Error"
		label data		"Directly Standardized Rates `mult'"
		save `stfile', `outrgro'
	}
end
