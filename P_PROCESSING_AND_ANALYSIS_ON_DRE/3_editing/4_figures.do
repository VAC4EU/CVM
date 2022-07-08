/* this step creates figures, including their labels, titles, and legends */ 

set more off

qui do par.do
qui do colors.do

// {}
// `

/*

// distribution of PERSONTIME

use ${dirdtain}PY_week_str,clear
bysort year month: gen tobelabelled = 1 if _n ==1
//set trace on
local xlab ""
forvalues mon = 1/12{
	qui levelsof week_labelled if month == `mon' & tobelabelled == 1,local(firstweek) clean
	foreach week of local firstweek{
		local monname = word("$months",`mon')
		di `"label define week_lab `week' "`monname'",modify"'
		label define week_lab `week' "`monname'",modify
		local xlab `"`xlab'  `week' "`monname'""'
		label list week_lab
		}
	}

label list week_lab
label values week_labelled week_lab
qui levelsof week_labelled if tobelabelled == 1,local(xtobelab) clean
gen ord = .
local j = 1
foreach datasource of global datasources{
	replace  ord = `j' if  datasource == `"`datasource'"'
	label define ord_lab `j' `"${name`datasource'}"',modify
	local j = `j'  +1
	}
label values ord ord_lab
di  `"`xlab'"'

//local graph `"graph bar (asis) PY,over(ord,descending) over(week_labelled,relabel(`xlab') label(labsize(vsmall)) )  asy stack  graphregion(${graphregion}`marginregion') ylabel(#6,angle(hor) labsize(small) format("%20.0gc")) legend(cols(1) pos(3) $legendstyle) xsize(10) ytitle(Person years)"'
//local graph `"graph bar (asis) PY,over(ord,descending) over(week_labelled,label(nolabels))  asy stack  graphregion(${graphregion}`marginregion') ylabel(#6,angle(hor) labsize(small) format("%20.0gc")) legend(cols(1) pos(3) $legendstyle) xsize(10) ytitle(Person years)"'
//di `"`graph'"'
//`graph'

gen minusord = 6-ord
sort week_lab minusord
bysort week_lab : gen cumPY = PY if _n == 1 
bysort week_lab : replace cumPY = cumPY[_n-1] + PY if _n > 1
local bars ""
local legend ""
local j = 1
levelsof ord, local(ords) clean
foreach ord of local ords{
	local bars  `"`bars' (bar cumPY week_labelled if ord== `ord', bcolor("${color`j'}") barwidth(.8) )"'  
	qui levelsof datasource if ord == `ord',local(datasource) clean
	local legend  `"`legend' `j' `"${name`datasource'}"' "'  
	local j = `j' + 1
	}
local graph `"twoway `bars' , graphregion(${graphregion}`marginregion') ylabel(#6,angle(hor) labsize(small)  format("%20.0gc")) xlabel(`xlab', labsize(small)) legend(order(`legend') cols(1) pos(3) $legendstyle) xsize(10) ytitle(Person years)"'
di `"`graph'"'
`graph'
graph export "${dirpdf}PYs_overall.pdf",replace



// CUMULATIVE COVERAGE

set more off
use  ${dirdtain}cumulative_coverage_by_datasource_and_ageband_and_vax,clear
sort week_labelled
bysort year month: gen tobelabelled = 1 if _n ==1
//set trace on
forvalues mon = 1/12{
	local monname = word("$months",`mon')
	qui levelsof week_labelled if month == `mon' & tobelabelled == 1,local(firstweek) clean
	di `"label define week_lab `firstweek' "`monname'""'
	capture label define week_lab `firstweek' "`monname'",modify
	}
label list week_lab
label values week_labelled week_lab
qui levelsof week_labelled if tobelabelled == 1,local(xtobelab) clean
label define dose_lab 1 "Dose 1" 2 "Dose 2"
label values dose dose_lab

foreach datasource of global datasources{
//	local datasource "ARS"
	local graphs ""
	foreach ageband of global agebands{
		if real(word(subinstr("`ageband'","-"," ",.),1)) >= 12{
			//local ageband  "60-69"
			di "`ageband'"
			preserve
			keep if datasource == "`datasource'" & ageband == "`ageband'"
			encode vx_manufacturer,gen(ordvax)
			gen minusordvax = 6-ordvax	
			
			// cum_vx: cumulative number of doses of all vaccines upt to minus alphabetic order vax, by week x (this is just to make the graph work)
			sort dose week_str minusordvax
			/*
			bysort  dose week_str: gen cum = cum_vx if _n == 1 
			bysort  dose week_str: replace cum = cum[_n-1] + cum_vx if _n > 1
			gen cum_perc = 100*cum/pop
			*/
			bysort  dose week_str: gen cum_perc = percentage if _n == 1 
			bysort  dose week_str: replace cum_perc = cum_perc[_n-1] + percentage if _n > 1
			local bars ""
			local legend ""
			local j = 1

			levelsof ordvax , local(ordvaxes) clean
			foreach ordvax of local ordvaxes{
				local bars  `"`bars' (bar cum_perc week_labelled if ordvax == `ordvax', bcolor("${color`j'}") )"'  
				qui levelsof vx_manufacturer if ordvax == `ordvax',local(vax) clean
				local vaxfixed = cond("`vax'"!="J&J","`vax'","Janssen")
				local legend  `"`legend' `j' `"${name`vaxfixed'}"' "'  
				local j = `j' + 1
				}
			local ysize "2"
			local legendstyle `" region(lcolor("${bianco}") )  keygap(1) symxsize(4) size(vsmall)"'
			local plotregion `""'
			// //local graph `"twoway `bars',by(dose , title(`ageband',size(small)) rows(1) colfirst graphregion(${graphregion}`marginregion') `plotregion'   legend(position(3) ) note("")) $subtitle legend(order(`legend')  cols(1) `legendstyle' title("Vaccine",size(vsmall)) )  xlabel(`xtobelab',valuelabel angle(hor) labsize(vsmall))  ytitle("") ylabel(0(10)100,angle(hor) labsize(vsmall)) ysize(`ysize')  xtitle("")"'  //xtitle(Weeks,size(vsmall)) title(${name`datasource'}) ytitle("Cumulative" "coverage",size(small)) 
			local graph `"twoway `bars',by(dose , title(`ageband',size(small)) rows(1) colfirst graphregion(${graphregion}`marginregion') `plotregion'   legend(off) note("")) $subtitle legend(off)  xlabel(`xtobelab',valuelabel angle(hor) labsize(vsmall))  ytitle("") ylabel(0(10)100,angle(hor) labsize(vsmall)) ysize(`ysize')  xtitle("")"'  //xtitle(Weeks,size(vsmall)) title(${name`datasource'}) ytitle("Cumulative" "coverage",size(small)) 
			di `"`graph'"'
			`graph'
			graph save "${dirfig}cum_coverage_`ageband'_`datasource'.gph",replace
			if "`ageband'" == word("$agebands",1){
				local bars ""
				local legend ""
				local j = 1
				foreach ordvax of local ordvaxes{
					local bars  `"`bars' (bar cum_perc week_labelled if ordvax == `ordvax' & _n == 0, bcolor("${color`j'}") )"'  
					qui levelsof vx_manufacturer if ordvax == `ordvax',local(vax) clean
					local vaxfixed = cond("`vax'"!="J&J","`vax'","Janssen")
					local legend  `"`legend' `j' `"${name`vaxfixed'}"' "'  
					local j = `j' + 1
					}
				local legendstyle `" region(lcolor("${bianco}") )  keygap(1) symxsize(8) size(large) "'
				local plotregion `"plotregion(style(none) istyle(none)) yscale(off) xscale(off) ylabel("") xlabel("") xtitle("") ytitle("")"'
				local legendin `"legend(order(`legend')  cols(1) `legendstyle' title("Vaccine",size(large)) bplacement(c) )"'
				local graph `"twoway `bars', `plotregion' graphregion(${graphregion}`marginregion') `legendin' "'
				di `"`graph'"'
				`graph'
				// pause on
				// pause
				graph save "${dirfig}cum_coverage_legend_`datasource'.gph",replace
				}
			restore
			local graphs `"`graphs' "${dirfig}cum_coverage_`ageband'_`datasource'.gph""'
			}
		}
	local graphs `"`graphs' "${dirfig}cum_coverage_legend_`datasource'.gph" "'
	di `"`graphs' "'
	graph combine `graphs', cols(2) ysize(5) graphregion(${graphregion}`marginregion') imargin(tiny)
	graph export "${dirpdf}cum_coverage_`datasource'.pdf",replace
	}



/*levels of covid per cohort and calendartime and coarse ageband */
use ${dirdtain}IR_covid_calendartime_vax_cohort_labelled,clear
local bandlist "young old"
set more off
sort week_labelled
bysort year month: gen tobelabelled = 1 if _n ==1
//set trace on
forvalues mon = 1/12{
	local monname = word("$months",`mon')
	qui levelsof week_labelled if month == `mon' & tobelabelled == 1,local(firstweek) clean
	foreach ord of local firstweek{
		di `"label define week_lab `ord' "`monname'""'
		capture label define week_lab `ord' "`monname'",modify
		}
	}
label list week_lab
label values week_labelled week_lab
qui levelsof week_labelled if tobelabelled == 1,local(xtobelab) clean
destring dose, replace force
label define dose_lab 1 "Dose 1" 2 "Dose 2"
label values dose dose_lab
local doses1 = "1 2"
local doses2 = "1 2"
local doses3 = "1 2"
local doses4 = "1"
local legendstyle `" region(lcolor("${bianco}") )  keygap(1) symxsize(8) size(medium) "'
local plotregion `"plotregion(style(none) istyle(none)) yscale(off) xscale(off) ylabel("") xlabel("") xtitle("") ytitle("")"'
local j = 1
local lines "(line IR week_labelled if ordvac == 0 & _n <0)" 
local legend `"`j' "No vaccine" "'
forvalues ordvac= 1/4{
	foreach dose of local doses`ordvac'{
		local j = 	`j' + 1
		local lpattern = cond(`dose' == 2 | `ordvac' == 4,"","lpattern(dash)")
		local correct_color = cond(`dose' == 2 | `ordvac' == 4,"","*.2")
		qui levelsof vx_manufacturer if ordvac == `ordvac',local(vax) clean
		local vaxfixed = cond("`vax'"!="J&J","`vax'","Janssen")
		local lab = cond(`ordvac' == 4,"${nameJanssen}","${name`vax'} dose `dose'")
		local lines `"`lines'(line IR week_labelled if ordvac == `ordvac' & dose == `dose' & _n < 0,lcolor("${color`ordvac'}`correct_color'") `lpattern' )"'
		local legend `"`legend' `j' "`lab'" "'
		}
	}
local legendin `"legend(order(`legend')  cols(1) `legendstyle' title("Vaccine",size(large)) bplacement(c) )"'
local xsize "10"
local graph `"twoway  `lines' , xtitle("")  xlabel(`xtobelab', labsize(small) valuelabel)  graphregion(${graphregion}`marginregion') `plotregion' `legendin'  ylabel(,angle(hor) labsize(small)) xsize(`xsize') ytitle("")"' // title(`"${name`AESI'}"')
di `"`graph'"'
`graph'
graph save "${dirfig}covidlevels_legend.gph",replace
foreach band of local bandlist{
	local Age = cond("`band'" == "young","0-59","60+")
	foreach datasource of global datasources{
		local graphs `""'
		foreach level of global levels`datasource'{
			set more off
			preserve
			keep if datasource == "`datasource'" & COVID == "L`level'" & ageband == "`Age'"
			gen IRfixed = IR
			qui sum IR if ordvac == 0,det
			local threshold = `r(max)'*1.5
			replace IRfixed = `threshold' if IRfixed > `threshold'
			gen extra = IRfixed <  IR & ordvac > 0
			drop if ordvac == 4 & dose == 2
			local lines "(line IR week_labelled if ordvac == 0)" 
			forvalues ordvac= 1/4{
				foreach dose of local doses`ordvac'{
					qui count if ordvac == `ordvac' & dose == `dose' 
					if `r(N)' > 0 {
						local lpattern = cond(`dose' == 2 | `ordvac' == 4,"","lpattern(dash)")
						local correct_color = cond(`dose' == 2 | `ordvac' == 4,"","*.2")
						local lines `"`lines'(line IRfixed week_labelled if ordvac == `ordvac' & dose == `dose',lcolor("${color`ordvac'}`correct_color'") `lpattern' )"'
						}
					}
				}
			local lines `"`lines'(line IRfixed week_labelled if extra == 1,lcolor(white) ) (scatter IRfixed week_labelled if extra == 1, msymbol(X) mcolor(red))"'
			local xsize "10"
			// local title = cond(`level' == 1, "Any level of severity", cond(`level' == 5,"Death with COVID-19", "Level of severity `level' or worse"))
			local title = "${namelevelseverity`level'}"
			local graph `"twoway  `lines' , xtitle("")  xlabel(`xtobelab', labsize(small) valuelabel)  graphregion(${graphregion}`marginregion') legend(off) ylabel(,angle(hor) labsize(small)) xsize(`xsize') ytitle("") title("`title'")"' // title(`"${name`AESI'}"')
			di `"`graph'"'
			`graph'
			graph save "${dirfig}covidlevel_`level'_`datasource'_`band'.gph",replace
			local graphs `"`graphs' "${dirfig}covidlevel_`level'_`datasource'_`band'.gph" "'
			restore
			}
		local graphs `"`graphs' "${dirfig}covidlevels_legend.gph" "'
		graph combine `graphs', cols(1) ysize(8) graphregion(${graphregion}`marginregion') imargin(tiny)
		graph export "${dirpdf}covidlevels_`datasource'_`band'.pdf",replace
		}
	}


	
/*distribution of severity levels of covid per cohort and time since vaccination by coarse ageband */
use ${dirdtain}IR_covid_timesincevaccination_vax_cohort_labelled,clear
local bandlist "young old"
set more off
destring dose, replace force
label define dose_lab 1 "Dose 1" 2 "Dose 2"
label values dose dose_lab
local doses1 = "1 2"
local doses2 = "1 2"
local doses3 = "1 2"
local doses4 = "1"
foreach datasource of global datasources{
	local legendstyle `" region(lcolor("${bianco}") )  keygap(1) symxsize(8) size(medium) "'
	local plotregion `"plotregion(style(none) istyle(none)) yscale(off) xscale(off) ylabel("") xlabel("") xtitle("") ytitle("")"'
	local lines ""
	local graph ""
	local legend ""
	local j = 0
	foreach level of global levels`datasource'{
		local j = 	`j' + 1
		local lines `"`lines' (bar IR week_since_vaccination if  _n < 0,color("${color`level'}") )"'
		local legend `"`legend' `j' "${namelevelseverity`level'}" "'
		}
	local legendin `"legend(order(`legend')  cols(1) `legendstyle' title("Level of COVID severity",size(large)) bplacement(c) )"'
	local xsize "10"
	local graph `"twoway  `lines' , xtitle("")  xlabel(`xtobelab', labsize(small) valuelabel)  graphregion(${graphregion}`marginregion') `plotregion' `legendin'  ylabel(,angle(hor) labsize(small)) xsize(`xsize') ytitle("")"' // title(`"${name`AESI'}"')
	di `"`graph'"'
	`graph'
	graph save "${dirfig}covidlevels_timesincevax_`datasource'_legend.gph",replace
	foreach band of local bandlist{
		local Age = cond("`band'" == "young","0-59","60+")
		forvalues ordvac= 1/4{			
			local graphs `""'
			foreach dose of local doses`ordvac'{
				set more off
				local graph ""
				local lines ""
				preserve
				keep if datasource == "`datasource'" & ordvac == `ordvac' & dose == `dose' & ageband == "`Age'"
				qui count 
				if `r(N)' > 0 {
					gen IRfixed = IR
					// qui sum IR if ordvac == 0,det
					// local threshold = `r(max)'*1.5
					// replace IRfixed = `threshold' if IRfixed > `threshold'
					// gen extra = IRfixed <  IR & ordvac > 0
					foreach level of global levels`datasource'{
						local lines `"`lines'(bar IRfixed week_since if COVID == "L`level'",color("${color`level'}") )"'
						}
	//				local lines `"`lines'(line IRfixed week_labelled if extra == 1,lcolor(white) ) (scatter IRfixed week_labelled if extra == 1, msymbol(X) mcolor(red))"'
					local xsize "10"
					local title = cond(`ordvac' == 4, "", "Dose `dose'")
					local graph `"twoway  `lines' , xtitle("")  xlabel(1(1)34, labsize(small) valuelabel)  graphregion(${graphregion}`marginregion') legend(off) ylabel(,angle(hor) labsize(small)) xsize(`xsize') ytitle("") title("`title'")"' // title(`"${name`AESI'}"')
					di `"`graph'"'
					`graph'
					graph save "${dirfig}covidlevel_timesincevax_`ordvac'_`dose'_`datasource'_`band'.gph",replace
					local graphs `"`graphs' "${dirfig}covidlevel_timesincevax_`ordvac'_`dose'_`datasource'_`band'.gph" "'
					}
				restore
				}

			if `"`graphs'"'! != "" {
				local graphs `"`graphs' "${dirfig}covidlevels_timesincevax_`datasource'_legend.gph" "'
				graph combine `graphs', cols(1) ysize(8) graphregion(${graphregion}`marginregion') imargin(tiny)
				qui levelsof vx_manufacturer if ordvac == `ordvac',local(vax) clean
				local vaxfixed = cond("`vax'"!="J&J","`vax'","Janssen")
				graph export "${dirpdf}covidlevels_timesincevax_`vaxfixed'_`datasource'_`band'.pdf",replace
				}
			}
		
		}
	}
	
	
	
// {}
// `


*/
	
/*************/
/*figure AESI BGR*/
/*************/

set more off
use ${dirdtain}BGR_IR_AESI_ageband, clear
gen ordageband = .
local j = 1
foreach ageband of global agebands{
	replace  ordageband = `j' if ageband == "`ageband'"
	label define ordageband_lab `j' `"`ageband'"',modify
	local j = `j' + 1
	}
label values ordageband ordageband_lab
gen DAPname = "" 
local j = 1
foreach datasource of global datasources{
	replace DAPname = `"${name`datasource'}"' if datasource == `"`datasource'"'
	}
local j = 1
gen ord = .
foreach datasource of global datasources{
	foreach ageband of global agebands{
		replace  ord = `j' if ageband == "`ageband'" & datasource == `"`datasource'"'
		label define ord_lab `j' `"${name`datasource'} - `ageband'"',modify
		local j = `j' + 1
		}
	}
label values ord ord_lab
drop if mi(ord)
//replace IR = .0000001 if IR == 0
export excel using "${dirtab}for_inspection_age_specific_BGR_2020",replace first(var)
foreach var of varlist IR* ub* lb*{
	tostring `var',force gen(`var'_BGR) format("%9.1f")
	}
gen IRCI_BGR = IR_BGR + " (" + lb_BGR + ";" + ub_BGR + ")" 

gen spot_label = .
gen spot_aux = .
foreach AESI of global AESIs {
	if strmatch("${events_narrow_broad}","*`AESI'*"){
		local condAESI = "`AESI'_narrow"
		}
	else{
		local condAESI = "`AESI'"
		}
	di "`condAESI'"
//	local AESI "CAD"
	qui count if AESI == "`AESI'_narrow" | AESI == "`AESI'"
	local howmuch = `r(N)'
	if `howmuch' > 0 { 
		qui sum ub if AESI == "`AESI'_narrow" | AESI == "`AESI'"
		local spot = cond(`r(max)' > 0,`r(max)',0)
		}
	else{
		local spot = 0
		}
	replace spot_label = `spot' if AESI == "`AESI'_narrow" | AESI == "`AESI'"

	replace spot_aux = `spot'*1.5 if AESI == "`AESI'_narrow" | AESI == "`AESI'"
	local scatterlabel `"(scatter ord spot_label if AESI == "`AESI'_narrow" | AESI == "`AESI'", msize(0) mlabel(IRCI_BGR)  mlabposition(3) mlabsize(vsmall) mlabcolor(black))"'
	local scatteraux `"(scatter ord spot_aux if AESI == "`AESI'_narrow" | AESI == "`AESI'", msize(0))"'
	local rcap `"(rcap ub lb ord  if AESI == "`AESI'_narrow" | AESI == "`AESI'", lcolor(gray*.5) color(gray*.5) hor)"'
	local scatter `"(scatter  ord IR if AESI == "`AESI'_narrow" | AESI == "`AESI'",mcolor(black))"'
	local plotregion " plotregion(margin( r = 20))"
	// local marginregion " margin( r = 5) "
	local ysize "10"
	local log "xscale(log)"
	local log ""
//	local graph `"twoway (rcap ub lb ord  if AESI == "`AESI'_narrow" | AESI == "`AESI'", lcolor(gray*.5) color(gray*.5) hor)  (scatter  ord IR if AESI == "`AESI'_narrow" | AESI == "`AESI'",mcolor(black)) `scatterlabel' , title(`"${name`AESI'}"')  xtitle("Background rate in 2020" "per 100,000 person-years")  xlabel(, labsize(small) )  graphregion(${graphregion}`marginregion') legend(off) ylabel(1(1)42,angle(hor) labsize(small) valuelabel) ysize(10)"'
	local graph `"twoway `rcap'  `scatter' `scatterlabel' ,   `log'  xtitle("Background rate in 2020" "per 100,000 person-years")  xlabel(, labsize(small) )  graphregion(${graphregion}`marginregion') `plotregion' legend(off) ylabel(1(1)55,angle(hor) labsize(small) valuelabel) ysize(`ysize') ytitle("")"' // title(`"${name`AESI'}"')

	di `"`graph'"'
	`graph'
	graph export "${dirpdf}labelled_forest_plot_`AESI'.pdf",replace
	}


	

/*************/
/*figure AESI vax cohorts*/
/*************/

set more off
//use ${dirdtain}IR_std_and_observed_and_expected_28_labelled, clear
use ${dirdtain}StandardizedRates_labelled, clear
/*
drop if type_vax == "J&J" & Dose == 2
drop if type_vax == "UKN"
*/
drop if ordvac == 4 & Dose ==2
drop if ordvac == 5
keep datasource AESI ordvac Dose IR_* lb_* ub_* IRCI_*
//reshape long IR_ lb_ ub_ IRCI_,i(datasource AESI ordvac Dose) j(type_rate) string

qui levelsof ordvac,local(ordvaxes) clean
sort AESI ordvac Dose datasource
local j = 35
gen ord = .
foreach ordvac of local ordvaxes{
	local firstvac = 1
	foreach Dose in 1 2{
		local firstdose = 1
		foreach datasource of global datasources{
			qui count if ordvac == `ordvac' & datasource == `"`datasource'"' & Dose == `Dose' 
			{
					replace  ord = `j' if ordvac == `ordvac' & datasource == `"`datasource'"' & Dose == `Dose' // & type_rate == "`type_rate'"
					local namedose = cond(`firstdose' == 1,"dose `Dose' -  ","")
					local nameordvac: label (ordvac) `ordvac'
					local namevac = cond(`firstvac' == 1 | `firstdose' == 1,"`nameordvac' - ","")
					label define ord_lab `j' `"`namevac'`namedose'${name`datasource'}"',modify
					local j = `j' - 1
					local firstdose = 2
					local firstvac = 2
				}
			}
		}
	}
label values ord ord_lab
drop if mi(ord)
foreach var of varlist IR_* lb_* ub_* {
	destring `var',replace
	}
bysort AESI: egen max  = max( IR_std)
bysort AESI: egen max2  = max( IR_BGR_std)
replace max = max(max,max2)
gen ub_std_fix = 3*max 
//replace IR = .0000001 if IR == 0
//export excel using "${dirtab}for_inspection_",replace first(var)
gen spot_label = .
gen spot_aux = .
foreach AESI of global AESIs {
	qui count if AESI == "`AESI'"
	local howmuch = `r(N)'
	if `howmuch' > 0 { 
		qui sum ub_std_fix if AESI == "`AESI'_narrow" | AESI == "`AESI'"
		local spot = cond(`r(max)' > 0,`r(max)',0)
		}
	else{
		local spot = 0
		}
	replace spot_label = `spot' if AESI == "`AESI'_narrow" | AESI == "`AESI'"

	replace spot_aux = `spot'*1.5 if AESI == "`AESI'_narrow" | AESI == "`AESI'"
	local scatterlabel `"(scatter ord spot_label if  AESI == "`AESI'", msize(0) mlabel(IRCI_std)  mlabposition(3) mlabsize(vsmall) mlabcolor(black))"'
	local scatteraux `"(scatter ord spot_aux if AESI == "`AESI'_narrow" | AESI == "`AESI'", msize(0))"'
	local scatterBGR `"(scatter ord IR_BGR_std if AESI == "`AESI'_narrow" | AESI == "`AESI'", mcolor(red) msymbol(X))"'
	local rcap `"(rcap ub_std_fix lb_std ord  if AESI == "`AESI'_narrow" | AESI == "`AESI'", lcolor(gray*.5) color(gray*.5) hor)"'
	local scatter `"(scatter  ord IR_std if AESI == "`AESI'_narrow" | AESI == "`AESI'",mcolor(black))"'
	local plotregion " plotregion(margin( r = 20))"
	// local marginregion " margin( r = 5) "
	local ysize "10"
	local log "xscale(log)"
	local log ""
//	local graph `"twoway (rcap ub lb ord  if AESI == "`AESI'_narrow" | AESI == "`AESI'", lcolor(gray*.5) color(gray*.5) hor)  (scatter  ord IR if AESI == "`AESI'_narrow" | AESI == "`AESI'",mcolor(black)) `scatterlabel' , title(`"${name`AESI'}"')  xtitle("Background rate in 2020" "per 100,000 person-years")  xlabel(, labsize(small) )  graphregion(${graphregion}`marginregion') legend(off) ylabel(1(1)42,angle(hor) labsize(small) valuelabel) ysize(10)"'
	local graph `"twoway `rcap'  `scatter' `scatterBGR' `scatterlabel' , `log'  xtitle("Standardised incidence rate" "per 100,000 person-years")  xlabel(, labsize(small) )  graphregion(${graphregion}`marginregion') `plotregion' legend(off) ylabel(1(1)35,angle(hor) labsize(small) valuelabel) ysize(`ysize') ytitle("")"' // title(`"${name`AESI'}"')

	di `"`graph'"'
	`graph'
	graph export "${dirpdf}AESI_std_`AESI'.pdf",replace
	}

	
	
	
	
	
/*figure cumulative IR*/
	
use ${dirdta}IR_fup_weeks_cumulative_and_BGR,clear
local manufacturerlist "AstraZeneca Janssen Moderna Pfizer UKN"
local bandlist "young old"
local doselistAstraZeneca "1"
local doselistModerna "1 2"
local doselistPfizer "1 2"
local doselistJanssen "1"
gen orddays = cumulative_week_fup
label define dose_lab 1 "Dose 1" 2 "Dose 2"
label values dose dose_lab
drop if Vaccine == "Janssen" & dose== 2
qui levelsof AESI,local(AESIs)
foreach AESI of local AESIs {
	foreach manufacturer in `manufacturerlist'{
		foreach band in `bandlist'{
			di `"`AESI' ${name`AESI'} - "`band'" `manufacturer'"'
			local Age = cond("`band'" == "young","0-59","60+")
			levelsof ordvac if Vaccine == "`manufacturer'",local( ordvac) clean
			qui count if AESI == "`AESI'" & ordvac == `ordvac' & Age =="`Age'"
			if `r(N)' > 0{	
				preserve
				keep if AESI == "`AESI'" & ordvac == `ordvac' & Age =="`Age'"
				egen max = max(IR)
				replace UL = min(UL,1.2*max) if max > 0
 	
			
				local lineBGR `"(line  IR_BGR orddays if ordvac == `ordvac' & Age =="`Age'",lcolor(red)  lwidth(thin)  sort)"'
				local ysize "ysize(7)"
				local numcol = cond("`manufacturer'" == "Janssen",1,2)
				local graph `"twoway (rarea UL LL orddays if ordvac == `ordvac' & Age =="`Age'", lcolor(gray*.2) fcolor(gray*.2)  color(gray*.2) sort) (line IR orddays if ordvac == `ordvac' & Age =="`Age'",lcolor(black) sort) `lineBGR' , by( ordDAP dose, cols(`numcol')  title("${name`manufacturer'}, age `Age'",size(small)) graphregion(${graphregion}`marginregion') note("") legend(off))  xline(4) xtitle("Week after vaccination") $subtitle xlabel(1(4)32) xticks(1(1)35)   ylabel(,angle(hor) labsize(vsmall)) xtitle("") ytitle("") `ysize' `xsize'"' // `"ytitle("Incidence per 100,000 person-years and 95% CI") xtitle(Weeks since vaccination)  xlabel(`orddays', valuelabel labsize(vsmall) angle(vert))"'
				di `"`graph'"'
				`graph'
					//graph save "${dirfig}IR_cumulative_per_week_after_vaccination_`manufacturer' `band'.gph",replace
				graph export "${dirpdf}IR_cumulative_per_week_after_vaccination_`AESI'_`manufacturer' `band'.pdf",replace
				restore
				}
 			}
 		}
	

	}
