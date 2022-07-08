/* this step creates tables in both LaTeX and xls format */ 
/* moreover it masks small counts*/

set more off
qui do par.do


/*
// ATTRITION

use ${dirdtain}attrition_table_edited, clear
foreach datasource of global datasources {
	gen `datasource'_r = real(`datasource')
	format `datasource'_r %20.0gc
	tostring `datasource'_r, replace use force
	replace `datasource' = `datasource'_r if !mi(`datasource'_r) & real(`datasource') != . & strmatch(Condition,"Recomm*") == 0
	replace `datasource' = "<5" if real(`datasource') > 0 & real(`datasource') < 5
	drop `datasource'_r
	}
export excel using "${dirtab}attrition_table",replace first(var)

use ${dirdtain}attrition_table_edited, clear
gen ord = _n
foreach datasource of global datasources {
	gen `datasource'_r = real(`datasource')
	format `datasource'_r %20.0gc
	tostring `datasource'_r, replace use force
	replace `datasource' = `datasource'_r if !mi(`datasource'_r) & real(`datasource') != . & strmatch(Condition,"Recomm*") == 0
	replace `datasource' = "$<5$" if real(`datasource') > 0 & real(`datasource') < 5
	replace `datasource' = "\multicolumn{1}{p{2cm}}{"+`datasource'+"}" if Condition == "Type of datasource" & "`datasource'" != "BIFAP_PC_HOSP"
	replace `datasource' = "\multicolumn{1}{p{3.5cm}}{"+`datasource'+"}" if Condition == "Type of datasource" & "`datasource'" == "BIFAP_PC_HOSP"
	}
local columns "p{4cm}"
local intest1 "\hline\multicolumn{1}{l}{Condition}"
local intest2 ""
local listvarlisttex "Condition  "
local n = 1
foreach datasource of global datasourcesdescr {
	local columns "`columns'rr"
	local intest1 "`intest1' & \multicolumn{1}{c}{${name`datasource'}}"
	local intest2 ""
	local listvarlisttex "`listvarlisttex' `datasource' "
	local n = `n' + 1
	}
replace Condition = "\hline\multicolumn{`n'}{c}{\bf Attrition} \\ \hline" if Condition == "Attrition"
sort ord
listtex `listvarlisttex' using ${dirtabtex}attrition.tex,rstyle(tabular) headlines("\begin{tabular}{`columns'}" "  `intest1' \\  " " `intest2' \\  \hline") footlines("\hline\end{tabular}") replace
 

// COHORT CHARACTERISTICS

use "${dirdtain}Table 6, Cohort characteristics at start of study (1-1-2020)",replace
// Variable Values ARS_N ARS_perc BIFAP_PC_N BIFAP_PC_perc BIFAP_PC_HOSP_N BIFAP_PC_HOSP_perc PHARMO_N PHARMO_perc CPRD_N CPRD_perc
// `"Age in categories"' `"Age in years"' `"At risk population at January 1-2020"' `"Person years across age categories"' `"Person years across sex"' `"Person years of follow-up"' `"Study population"'


foreach datasource of global datasources {
	local var `datasource'_N
	gen `var'_r = real(`var')
	format `var'_r %20.0gc
	tostring `var'_r, replace use force
	replace `var' = `var'_r if !mi(`var'_r) & real(`var') != . 
	drop `var'_r
	replace `var' = "<5" if real(`var') > 0 & real(`var') < 5 & (Variable == `"At risk population at January 1-2020"' | Variable == `"Person years across age categories"' |Variable == `"Person years across sex"'  |Variable == `"Person years of follow-up"' |Variable == `"Study population"')
	}
export excel using "${dirtab}Table 6, Cohort characteristics at start of study (1-1-2020)",replace first(var)

use "${dirdtain}Table 6, Cohort characteristics at start of study (1-1-2020)",replace
gen ord = _n
foreach datasource of global datasources {
	local var `datasource'_N
	gen `var'_r = real(`var')
	format `var'_r %20.0gc
	tostring `var'_r, replace use force
	replace `var' = `var'_r if !mi(`var'_r) & real(`var') != . 
	drop `var'_r
	replace `var' = "$<5$" if real(`var') > 0 & real(`var') < 5 & (Variable == `"At risk population at January 1-2020"' | Variable == `"Person years across age categories"' |Variable == `"Person years across sex"'  |Variable == `"Person years of follow-up"' |Variable == `"Study population"')
	}
foreach var of varlist *_perc{
	gen `var'_tex = subinstr(`var', "%","\%",.)
	}
sort Variable ord
bysort Variable: gen num = _N
bysort Variable: gen varnametex = "\hline \multirow{ "+string(num)+"}{2cm}{\begin{minipage}{2cm}" + Variable + "\end{minipage}}" if _n == 1
replace varnametex = "\hline " + Variable  if num == 1
local columns "p{2cm}p{2cm}"
local intest1 "\hline\multicolumn{1}{l}{Variables}&\multicolumn{1}{l}{Values}"
local intest2 "& "
local listvarlisttex "varnametex Values  "
foreach datasource of global datasourcesdescr {
	local columns "`columns'rr"
	local intest1 "`intest1' & \multicolumn{2}{c}{${name`datasource'}}"
	local intest2 "`intest2' & N & \%"
	local listvarlisttex "`listvarlisttex' `datasource'_N `datasource'_perc_tex "
	}
sort ord
listtex `listvarlisttex' using ${dirtabtex}cohort_characteristics.tex,rstyle(tabular) headlines("\begin{tabular}{`columns'}" "  `intest1' \\  " " `intest2' \\  \hline") footlines("\hline\end{tabular}") replace
 

 
// COHORT CHARACTERISTICS AT FIRST AND SECOND DOSE

foreach datasource of global datasourcesdescr {
	foreach dose in first second{
		local b = cond("`dose'" == "first","","_b")
		use "${dirdtain}Table ${numtab10`datasource'}`b', Cohort characteristics at `dose' dose of COVID-19 vaccine `datasource'",replace
		capture rename AZ AZ_N
		capture rename UKN UKN_N
		replace Variable = "Persons with a `dose' dose" if Variable == "Study population"
		//`"2020 December"' `"2021 April"' `"2021 February"' `"2021 January"' `"2021 June"' `"2021 March"' `"2021 May"' `"Age in categories"' `"Age in years"' `"At risk population at date of vaccination"' `"Month of first vaccination"' `"Person years across age categories"' `"Person years across sex"' `"Person-years of follow-up between first and second dose"' `"Study population"'
		preserve
		foreach vax in Pfizer Moderna AZ Janssen UKN {
			local var `vax'_N
			capture gen `var'_r = real(`var')
			capture format `var'_r %20.0gc
			capture tostring `var'_r, replace use force
			capture replace `var' = `var'_r if !mi(`var'_r) & real(`var') != . 
			capture drop `var'_r
			capture replace `var' = "<5" if real(`var') > 0 & real(``var') < 5 & (Variable != `"Age in years"' & Variable !=  `"Month of first vaccination"')
			}

		export excel using "${dirtab}Table ${numtab10`datasource'}, Cohort characteristics at `dose' dose of COVID-19 vaccine `datasource'",replace first(var)
		restore
		
		// Variable Values Pfizer_N Pfizer_perc Moderna_N Moderna_perc AZ AZ_perc Janssen_N Janssen_perc
		gen ord = _n
		foreach vax in Pfizer Moderna AZ Janssen UKN {
			local var `vax'_N
			capture gen `var'_r = real(`var')
			capture format `var'_r %20.0gc
			capture tostring `var'_r, replace use force
			capture replace `var' = `var'_r if !mi(`var'_r) & real(`var') != . 
			capture drop `var'_r
			capture replace `var' = "$<5$" if real(`var') > 0 & real(`var') < 5 & (Variable != `"Age in years"' & Variable !=  `"Month of first vaccination"' & Variable !=  `"Month of second vaccination"')
			}
		foreach var of varlist *_perc{
			gen `var'_tex = subinstr(`var', "%","\%",.)
			}
		sort Variable ord
		bysort Variable: gen num = _N
		bysort Variable: gen varnametex = "\hline \multirow{ "+string(num)+"}{2cm}{\begin{minipage}{2cm}" + Variable + "\end{minipage}}" if _n == 1
		replace varnametex = "\hline " + Variable  if num == 1
		local columns "p{2cm}p{2cm}"
		local intest1 "\hline\multicolumn{1}{l}{Variables}&\multicolumn{1}{l}{Values}"
		local intest2 "& "
		local listvarlisttex "varnametex Values  "
		foreach vax in Pfizer Moderna AZ Janssen UKN {
			capture confirm var `vax'_N
			if _rc == 0 {
				local columns "`columns'rr"
				local intest1 "`intest1' & \multicolumn{2}{c}{${name`vax'}}"
				local intest2 "`intest2' & N & \%"
				local listvarlisttex "`listvarlisttex' `vax'_N `vax'_perc_tex "
				}
			}
		sort ord
		listtex `listvarlisttex' using ${dirtabtex}`dose'dose`datasource'.tex,rstyle(tabular) headlines("\begin{tabular}{`columns'}" "  `intest1' \\  " " `intest2' \\  \hline") footlines("\hline\end{tabular}") replace
		}
	}


// DOSES AND GAPS

use "${dirdtain}Table 11, COVID-19 vaccination by dose and time period between first and second dose (days)",replace
gen ord = _n
append using ${dirdtain}N_study_population
replace ord = 0 if mi(ord)
replace Dose = "Study population" if mi(Dose)
replace Measure = "N" if mi(Measure)
foreach num of varlist *_N {
	local var `num'_N
	capture gen `num'_r = real(`num')
	capture 
	format `num'_r %20.0gc
	capture tostring `num'_r, replace use force
	capture replace `num' = `num'_r if !mi(`num'_r) & real(`num') != . 
	capture drop `num'_r
	capture replace `num' = "<5" if real(`num') > 0 & real(`num') < 5 & (strmatch(Dose,"*distance*") == 0)
	}
drop Variable Value
sort ord 
drop ord
// // Dose Measure ARS_N  ARS_N ARS_perc BIFAP_PC_N BIFAP_PC_perc BIFAP_PC_HOSP_N BIFAP_PC_HOSP_perc PHARMO_N PHARMO_perc CPRD_N CPRD_perc
export excel using "${dirtab}Table 11, COVID-19 vaccination by dose and time period between first and second dose (days)",replace first(var)

use "${dirdtain}Table 11, COVID-19 vaccination by dose and time period between first and second dose (days)",replace
gen ord = _n
append using ${dirdtain}N_study_population
replace ord = 0 if mi(ord)
replace Dose = "Study population" if mi(Dose)
replace Measure = "N" if mi(Measure)
drop Variable Value
sort ord 
foreach num of varlist *_N {
	local var `num'_N
	capture gen `num'_r = real(`num')
	capture 
	format `num'_r %20.0gc
	capture tostring `num'_r, replace use force
	capture replace `num' = `num'_r if !mi(`num'_r) & real(`num') != . 
	capture drop `num'_r
	capture replace `num' = "$<5$" if real(`num') > 0 & real(`num') < 5 & (strmatch(Dose,"*distance*") == 0)
	}
gen group = word(Dose,1) 
local manufacturers = "Pfizer Moderna AstraZeneca Janssen"
replace group = "" if subinstr("`manufacturers'",group,"",1) == "`manufacturers'"
replace group = group[_n-1] if mi(group)
foreach var of varlist *_perc{
	gen `var'_tex = subinstr(`var', "%","\%",.)
	}
egen ggroup = group(group Dose)
// // bysort group: gen num = _N
// // bysort group: gen varnametex = "\hline \multirow{ "+string(num)+"}{2cm}{\begin{minipage}{2cm}" + Dose + "\end{minipage}}" if _n == 1
// // replace varnametex = "\hline " + Variable  if num == 1
gen varnametex = Dose
sort ggroup ord
bysort ggroup: replace  varnametex = "" if _n > 1
// replace varnametex = "" if varnametex == varnametex[_n-1]
sort group ord
replace varnametex = "\multirow{ "+string(5)+"}{2cm}{\begin{minipage}{2cm}" + varnametex + "\end{minipage}}" if varnametex[_n+1] == "" & varnametex != ""
bysort group: replace varnametex = "\hline "+varnametex if _n == 1
local columns "p{2cm}p{2cm}"
local intest1 "\hline\multicolumn{1}{l}{Doses}&\multicolumn{1}{l}{Measures}"
local intest2 "& "
local listvarlisttex "varnametex Measure  "
foreach datasource of global datasourcesdescr {
	local columns "`columns'rr"
	local intest1 "`intest1' & \multicolumn{2}{c}{${name`datasource'}}"
	local intest2 "`intest2' & N & \%"
	local listvarlisttex "`listvarlisttex' `datasource'_N `datasource'_perc_tex "
	}
sort ord
listtex `listvarlisttex' using ${dirtabtex}dosesgaps.tex,rstyle(tabular) headlines("\begin{tabular}{`columns'}" "  `intest1' \\  " " `intest2' \\  \hline") footlines("\hline\end{tabular}") replace
 

/*COVID 19*/

use ${dirdtain}IR_covid_calendartime_vax_cohort_labelled,clear
sort  datasource week_str ordvac dose ageband COVID
order  datasource week_str ordvac dose ageband COVID
keep datasource week_str ordvac dose ageband COVID IR ub lb
export excel using ${dirtab}IR_covid_calendartime_vax_cohort,replace first(var)


use ${dirdtain}IR_covid_timesincevaccination_vax_cohort_labelled,clear
sort datasource week_since ordvac dose ageband COVID
order datasource week_since ordvac dose ageband COVID Numerator IR ub lb
keep datasource week_since ordvac dose ageband COVID Numerator IR ub lb
foreach num of varlist Numerator* {
	capture tostring `num', replace
	capture gen `num'_r = real(`num')
	capture format `num'_r %20.0gc
	capture tostring `num'_r, replace use force
	capture replace `num' = `num'_r if !mi(`num'_r) & real(`num') != . 
	capture drop `num'_r
	capture 
	replace `num' = "<5" if real(`num') > 0 & real(`num') < 5 // & (strmatch(Dose,"*distance*") == 0)
	}
export excel using ${dirtab}IR_covid_timesincevaccination_vax_cohort,replace first(var)


// INCIDENT CASES

use ${dirdtain}incident_cases_AESI_labelled,clear
foreach num of varlist count* {
	local var `num'_N
	capture gen `num'_r = real(`num')
	capture 
	format `num'_r %20.0gc
	capture tostring `num'_r, replace use force
	capture replace `num' = `num'_r if !mi(`num'_r) & real(`num') != . 
	capture drop `num'_r
	capture replace `num' = "<5" if real(`num') > 0 & real(`num') < 5 // & (strmatch(Dose,"*distance*") == 0)
	}
keep system_num AESI_num year count*
keep if !mi(AESI_num)
export excel using ${dirtab}incident_cases_AESI_exported,replace first(var)


use ${dirdtain}incident_cases_AESI_labelled,clear
keep system_num AESI_num year count*
// system_num AESI_num year countARS countBIFAP_PC countBIFAP_PC_HOSP countPHARMO countCPRD
gen ord = _n
// // foreach var of varlist *_perc{
// // 	gen `var'_tex = subinstr(`var', "%","\%",.)
// // 	}
sort system_num ord
foreach num of varlist count* {
	local var `num'_N
	capture gen `num'_r = real(`num')
	capture 
	format `num'_r %20.0gc
	capture tostring `num'_r, replace use force
	capture replace `num' = `num'_r if !mi(`num'_r) & real(`num') != . 
	capture drop `num'_r
	capture replace `num' = "$<5$" if real(`num') > 0 & real(`num') < 5 //B& (strmatch(Dose,"*distance*") == 0)
	}
bysort system_num: gen num = _N
decode system_num,gen(system)
decode AESI_num,gen(AESI)
drop if mi(AESI)
local len = 3 + wordcount("${datasourcedescr}") -1 
bysort system_num: gen varnametex = "\hline \multirow{ "+string(num)+"}{2cm}{\begin{minipage}{2cm}" + system + "\end{minipage}}" if _n == 1
bysort system_num AESI_num: gen AESInametex = "\multirow{ "+string(2)+"}{4cm}{\begin{minipage}{4cm}" +  AESI + "\end{minipage}}" if _n == 1
// replace AESInametex = " " if mi(AESInametex)
local columns "p{2cm}p{4cm}l"
local intest1 "\hline\multicolumn{1}{l}{System}&\multicolumn{1}{l}{AESI} & Year"
local intest2 "& "
local listvarlisttex "varnametex AESInametex year  "
foreach datasource of global datasourcesdescr {
	local columns "`columns'r"
	local intest1 "`intest1' & ${name`datasource'}"
	local intest2 ""
	local listvarlisttex "`listvarlisttex' count`datasource' "
	replace count`datasource' = "$<5$" if count`datasource' == "<5"
	}
sort ord
listtex `listvarlisttex' using ${dirtabtex}incident_cases.tex,rstyle(tabular) headlines("\begin{tabular}{`columns'}" "  `intest1' \\  " " `intest2' \\  \hline") footlines("\hline\end{tabular}") replace

*/

// code counts
set more off
use  ${dirdtain}AESI_code_counts_labelled,clear
sort system_num AESI_num DAPname meaning coding code year
order system_num AESI_num DAPname meaning coding code year 
foreach num of varlist count* {
	local var `num'_N
	capture gen `num'_r = real(`num')
	capture 
	format `num'_r %20.0gc
	capture tostring `num'_r, replace use force
	capture replace `num' = `num'_r if !mi(`num'_r) & real(`num') != . 
	capture drop `num'_r
	capture replace `num' = "<5" if real(`num') > 0 & real(`num') < 5 // & (strmatch(Dose,"*distance*") == 0)
	}
keep system_num AESI_num DAPname meaning coding code year count AESI
save ${dirdtain}AESI_code_counts_labelled_FOR_APPENDIX,replace

/*

/*


// INCIDENCE RATES AFTER VACCINATION


local AESIlist " $eventsnoCOVID"
foreach AESI in `AESIlist'{
	use ${dirdtain}table_IR_std_and_observed_and_expected_28_labelled if AESI == "`AESI'",clear
	di "`AESI'"
	qui count
	if `r(N)'>0{
		gen namelong = "Background crude incidence rate" if strmatch(name,"1_*") 
		replace namelong = "Background age-standardised incidence rate" if strmatch(name,"2_*") 
		replace namelong = "Observed person-years after vaccination" if strmatch(name,"3_*") 
		//replace namelong = "Crude number of cases expected" if strmatch(name,"4_*") 
		replace namelong = "Expected cases" if strmatch(name,"4_*") 
		replace namelong = "Observed cases" if strmatch(name,"5_*") 
		replace namelong = "Age-standardised incidence rate" if strmatch(name,"6_*") 
		replace namelong = "Age-standardised rate difference" if strmatch(name,"7_*") 
		preserve
	
		foreach num of varlist row_* {
			local var `num'_N
			capture gen `num'_r = real(`num')
			format `num'_r %20.0gc
			capture tostring `num'_r, replace use force
			capture replace `num' = `num'_r if !mi(`num'_r) & real(`num') != . & strmatch(name,"*_PY")
			capture drop `num'_r
			capture replace `num' = "<5" if real(`num') > 0 & real(`num') < 5 & strmatch(name, "*_obs")
			}
		order AESI_num ordvac Dose namelong row*
		keep AESI_num ordvac Dose namelong row*
		export excel using "${dirtab}Table XX Crude and standardised incidence rates for `AESI'",replace first(var)
		restore
		foreach num of varlist row_* {
			local var `num'_N
			capture gen `num'_r = real(`num')
			format `num'_r %20.0gc
			capture tostring `num'_r, replace use force
			capture replace `num' = `num'_r if !mi(`num'_r) & real(`num') != . & strmatch(name,"*_PY")
			capture drop `num'_r
			capture replace `num' = "$<5$" if real(`num') > 0 & real(`num') < 5 &  strmatch(name, "*_obs")
			}
		gen ddose = Dose
		replace ddose = 0 if mi(ddose)
		egen group = group(ordvac ddose)

		decode ordvac,gen(vax)
		gen varnametex = vax
		sort group name
		bysort group: replace  varnametex = varnametex + cond(!mi(Dose)," dose " +string(Dose),"") if _n == 1
		bysort group: replace  varnametex = "" if _n > 1
		// replace varnametex = "" if varnametex == varnametex[_n-1]
		sort group name
		bysort group: gen numrows = _N
		replace varnametex = "\multirow{ "+string(numrows)+"}{1.5cm}{\begin{minipage}{2cm}" + varnametex + "\end{minipage}}" if varnametex[_n+1] == "" & varnametex != ""
		bysort group: replace varnametex = "\hline "+varnametex if _n == 1
		local columns "p{1.5cm}p{4cm}"
		local intest1 "\hline\multicolumn{1}{l}{Vaccine}&\multicolumn{1}{l}{Measure}"
		// local intest2 "& "
		local listvarlisttex "varnametex namelong  "
		foreach datasource of global datasourcesdescr {
			local columns "`columns'r"
			local intest1 "`intest1' & \multicolumn{1}{c}{${name`datasource'}}"
			// local intest2 "`intest2' & Value"
			local listvarlisttex "`listvarlisttex' row_`datasource' "
			}
		sort ord
		listtex `listvarlisttex' using "${dirtabtex}Table XX Crude and standardised incidence rates for `AESI'.tex",rstyle(tabular) headlines("\begin{tabular}{`columns'}" "  `intest1' \\  " " `intest2' \\  \hline") footlines("\hline\end{tabular}") replace
 
		}
	}

