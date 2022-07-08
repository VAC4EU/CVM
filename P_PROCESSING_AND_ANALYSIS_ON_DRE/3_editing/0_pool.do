/*this step gathers data from DAPs' input and from the T4 step executed on DRE (folder 2)*/
/*data is pooled across DAPs and stored in Stata format*/

set more off

qui do par.do
qui do folders.do
qui do colors.do

// {}
// `

/********/
/*import*/
/********/



/*input from DAPs: CDM_SOURCE*/
/*output: CDM_SOURCE  */
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}CDM_SOURCE.csv",clear case(preserve) delim(",")
	qui count
	if `r(N)' > 1{
		import delimited using "${dirinput`datasource'}CDM_SOURCE.csv",clear case(preserve) delim(";")
		}
	gen datasource = "`datasource'"
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
save ${dirdtain}CDM_SOURCE,replace


/*input from DAPs: RISK_BIRTHCOHORTS_TIMESINCEVACCINATION*/
/*output: BGR_IR_AESI_ageband including stratfied BGR IRs */
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}dashboard tables/RISK_BIRTHCOHORTS_TIMESINCEVACCINATION${suffix`datasource'}.csv",clear case(preserve)
	keep if dose == "no_dose"
	replace datasource = "`datasource'"
	di "`datasource'"
	foreach AESI of global AESIexcluded`datasource'{
		di "`AESI'"
		drop if AESI == "`AESI'"
		drop if AESI == "`AESI'_narrow"
		drop if AESI == "`AESI'_possible"
		}
	isid datasource ageband AESI
	keep  datasource  ageband AESI Numerator IR lb ub
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
save ${dirdtain}BGR_IR_AESI_ageband,replace



/*input from DAPs: RISK_BIRTHCOHORTS_TIMESINCEVACCINATION*/
/*output: IR_AESI_after_vx_dose_week_ageban */
set more off
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}dashboard tables/RISK_BIRTHCOHORTS_TIMESINCEVACCINATION${suffix`datasource'}.csv",clear case(preserve)
	drop if dose == "no_dose"
	replace datasource = "`datasource'"
	di "`datasource'"
	foreach AESI of global AESIexcluded`datasource'{
		di "`AESI'"
		drop if AESI == "`AESI'"
		drop if AESI == "`AESI'_narrow"
		drop if AESI == "`AESI'_possible"
		}
	drop if ageband == "all_birth_cohorts" | ageband == "60+" | ageband == "0-59" 
	isid datasource vx_manufacturer dose ageband AESI week

	keep  datasource vx_manufacturer dose ageband AESI week IR lb ub
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
save ${dirdtain}IR_AESI_after_vx_dose_week_ageband,replace

/*input from DAPs: BENEFIT_BIRTHCOHORTS_CALENDARTIME*/
/*output: IR_covid_calendartime_vax_cohort */
set more off
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}dashboard tables/BENEFIT_BIRTHCOHORTS_CALENDARTIME${suffix`datasource'}.csv",clear case(preserve)
	keep if ageband  == "60+" | ageband == "0-59"
	replace datasource = "`datasource'"
	keep  datasource vx_manufacturer dose COVID week IR lb ub ageband
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
isid datasource vx_manufacturer dose COVID week ageband
save ${dirdtain}IR_covid_calendartime_vax_cohort,replace


/*input from DAPs: BENEFIT_BIRTHCOHORTS_CALENDARTIME*/
/*output: IR_covid_calendartime_vax_cohort */
set more off
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}dashboard tables/BENEFIT_BIRTHCOHORTS_TIMESINCEVACCINATION${suffix`datasource'}.csv",clear case(preserve)
	keep if ageband  == "60+" | ageband == "0-59"
	drop if vx == "none"
	replace datasource = "`datasource'"
	keep  datasource vx_manufacturer dose COVID week_since_vaccination IR lb ub ageband Numerator
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
isid datasource vx_manufacturer dose COVID week_since_vaccination ageband
save ${dirdtain}IR_covid_timesincevaccination_vax_cohort,replace



/*input from DAPs: COVERAGE_BIRTHCOHORTS*/
/*output: coverage_weeks_manufacturer_ageband including age-stratfied coverage per week */
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}dashboard tables/COVERAGE_BIRTHCOHORTS${suffix`datasource'}.csv",clear case(preserve)
	replace datasource = "`datasource'"
	di "`datasource'"
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
save ${dirdtain}coverage_weeks_manufacturer_ageband,replace


/*input from DAPs: DOSES_BIRTHCOHORTS*/
/*output: doses_weeks_manufacturer_ageband including age-stratfied doses per week */
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}dashboard tables/DOSES_BIRTHCOHORTS${suffix`datasource'}.csv",clear case(preserve)
	replace datasource = "`datasource'"
	di "`datasource'"
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
save ${dirdtain}doses_weeks_manufacturer_ageband,replace


/*input from DAPs: D4_persontime_risk_week.csv*/
/*output: Persontime_ageband_week_sex */
set more off
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}D4_persontime_risk_week.csv",clear case(preserve)
	keep sex ageband_at_study_entry week Persontime
	gen datasource = "`datasource'"
	collapse (sum) Persontime,by(datasource ageband_at_study_entry week sex)
	di "`datasource'"
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
save ${dirdtain}Persontime_ageband_week_sex,replace

/*input from DAPs: D4_persontime_risk_year.csv*/
/*output: counts_AESI_and_persontime_dose_type_vax_week_fup_sex_ageband */
set more off
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}D4_persontime_risk_year.csv",clear case(preserve)
	drop if Dose == 0
	drop *broad*
	collapse (sum) Persontime* *_b,by(sex ageband_at_study_entry Dose type_vax week_fup)
	isid sex ageband_at_study_entry Dose type_vax week_fup
	gen datasource = "`datasource'"
	di "`datasource'"
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
save ${dirdtain}counts_AESI_and_persontime_dose_type_vax_week_fup_sex_ageband,replace


/*input from DAPs: D4_persontime_risk_year.csv*/
/*output: counts_AESI_and_persontime_dose_type_vax_week_fup_sex_ageband_28 */
set more off
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}D4_persontime_risk_year.csv",clear case(preserve)
	drop if Dose == 0
	drop *broad*
	keep if week_fup <=4
	collapse (sum) Persontime* *_b,by(sex ageband_at_study_entry Dose type_vax week_fup)
	isid sex ageband_at_study_entry Dose type_vax week_fup
	gen datasource = "`datasource'"
	di "`datasource'"
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
save ${dirdtain}counts_AESI_and_persontime_dose_type_vax_week_fup_sex_ageband_28,replace


/*input from DAPs: D4_persontime_risk_year.csv*/
/*output: counts_AESI_and_persontime_BGR_sex_ageband */
set more off
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}D4_persontime_risk_year.csv",clear case(preserve)
	keep if Dose == 0 & year == 2020
	drop *broad*
	collapse (sum) Persontime* *_b,by(sex ageband_at_study_entry )
	isid sex ageband_at_study_entry 
	gen datasource = "`datasource'"
	di "`datasource'"
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
save ${dirdtain}counts_AESI_and_persontime_BGR_sex_ageband,replace


/*input from DAPs: RISK_BIRTHCOHORTS_CALENDARTIME.csv*/
/*output: incident_cases_AESI_DEATH, */
set more off
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}dashboard tables/RISK_BIRTHCOHORTS_CALENDARTIME${suffix`datasource'}.csv",clear case(preserve)
	keep if AESI == "DEATH"
	keep if ageband == "all_birth_cohorts"
	tostring week,replace
	gen year = substr(week, 1,4)
	tab year
	collapse (sum) Numerator ,by(year AESI)
	rename Numerator count
	tostring count,replace
	isid year 
	gen datasource = "`datasource'"
	di "`datasource'"
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
reshape wide count, i(year AESI) j(datasource) string
save ${dirdtain}incident_cases_AESI_DEATH,replace




/*input from DAPs: Table 5, Attrition diagram 1${suffix`datasource'}.cs*/
/*output: attrition_table */
local files ""
foreach datasource of global datasourcesdescr {
	import delimited using "${dirinput`datasource'}Dummy tables October/Table 5, Attrition diagram 1${suffix`datasource'}.csv",clear case(preserve)  varnames(nonames)    
	if "`datasource'" == word("$datasources",1){
		gen n = _n
		}
	rename v1 Condition
	rename v2 `datasource'
	replace `datasource' = "<5" if real(`datasource') > 0 & real(`datasource') < 5 
	tempfile xxx`datasource'
	save `xxx`datasource'',replace
	local files "`files' `xxx`datasource''"
	}
foreach datasource of global datasourcesdescr {
	if "`datasource'" == word("$datasources",1){
		use `xxx`datasource'',clear
		}
	else{
		merge 1:1  Condition using `xxx`datasource'',keep(master match)
		drop _merge
		}
	}
sort n
drop n
save ${dirdtain}attrition_table,replace

/*input from DAPs: Table 6, Cohort characteristics at start of study (1-1-2020)${suffix`datasource'}.csv*/
/*output: Table 6, Cohort characteristics at start of study (1-1-2020) */
local files ""
foreach datasource of global datasourcesdescr {
	import delimited using "${dirinput`datasource'}Dummy tables October/Table 6, Cohort characteristics at start of study (1-1-2020)${suffix`datasource'}.csv",clear case(preserve)  varnames(nonames)    
	if "`datasource'" == word("$datasources",1){
		gen n = _n
		}
	rename v1 Variable
	rename v2 Values
	rename v3 `datasource'_N
	rename v4 `datasource'_perc
	drop in 1/2
	local PY = `datasource'_N[2]
	levelsof `datasource'_N if Variable == `"Person years across sex"'  & Values == `"Male"', local(PYm) clean
	levelsof `datasource'_N if Variable == `"Person years across sex"'  & Values == `"Female"', local(PYf) clean
	replace `datasource'_perc = string(100*`PYm'/`PY',"%3.1f") if Variable == `"Person years across sex"'  & Values == `"Male"'
	replace `datasource'_perc = string(100*`PYf'/`PY',"%3.1f") if Variable == `"Person years across sex"'  & Values == `"Female"'
//	replace `datasource'_N = "<5" if real(`datasource'_N) > 0 & real(`datasource'_N) < 5 
	tempfile xxx`datasource'
	save `xxx`datasource'',replace
	local files "`files' `xxx`datasource''"
	}
foreach datasource of global datasourcesdescr {
	if "`datasource'" == word("$datasources",1){
		use `xxx`datasource'',clear
		}
	else{
		merge 1:1  Variable Values using `xxx`datasource'',keep(master match)
		drop _merge
		}
	}
sort n 
drop n
save "${dirdtain}Table 6, Cohort characteristics at start of study (1-1-2020)",replace


/*input from DAPs: from each datasource Cohort characteristics at first dose of COVID-19 vaccine*/
/*output: Cohort characteristics at first dose of COVID-19 vaccine `datasource'*/
local files ""
foreach datasource of global datasourcesdescr {
	import delimited using "${dirinput`datasource'}Dummy tables October/Table ${numtab10`datasource'}, Cohort characteristics at first dose of COVID-19 vaccine ${namesuff`datasource'}${suffix`datasource'}.csv",clear case(preserve)  varnames(nonames)    
	
		gen n = _n
		bysort v1: egen ordvar = min(n)
		sort ordvar n
		bysort ordvar: gen ordvalues = _n
		replace ordvalues = real(word(subinstr(v2,"-"," ",1),1)) if (strmatch(v1,"Age*") | strmatch(v1,"* age*")) & v1 != "Age in years"
		replace ordvalues = 88 if v2 == "80+" 
		replace ordvalues = 99 if v2 == "60+" 
		egen newn = group(ordvar ordvalues)
		replace n = newn
		sort n
		drop ordvar ordvalues newn
	rename v1 Variable
	rename v2 Values
	rename v3 Pfizer_N
	rename v4 Pfizer_perc
	rename v5 Moderna_N
	rename v6 Moderna_perc
	rename v7 AZ
	rename v8 AZ_perc
	capture rename v9 Janssen_N
	capture rename v10 Janssen_perc
	capture rename v11 UKN_N
	capture rename v12 UKN_perc
	
	drop in 1/2
	// replace `datasource'_N = "<5" if real(`datasource'_N) > 0 & real(`datasource'_N) < 5 
	tempfile xxx`datasource'
	save `xxx`datasource'',replace
	local files "`files' `xxx`datasource''"
	capture sort n 
	capture drop n
	save "${dirdtain}Table ${numtab10`datasource'}, Cohort characteristics at first dose of COVID-19 vaccine `datasource'",replace
	

	}

/*input from DAPs: Cohort characteristics at second dose of COVID-19 vaccine*/
/*output for each datasource: Cohort characteristics at second dose of COVID-19 vaccine `datasource'*/
local files ""
foreach datasource of global datasourcesdescr {
	import delimited using "${dirinput`datasource'}Dummy tables October/Table ${numtab10`datasource'}, Cohort characteristics at second dose of COVID-19 vaccine ${namesuff`datasource'}${suffix`datasource'}.csv",clear case(preserve)  varnames(nonames)    
		gen n = _n
		bysort v1: egen ordvar = min(n)
		sort ordvar n
		bysort ordvar: gen ordvalues = _n
		replace ordvalues = real(word(subinstr(v2,"-"," ",1),1)) if (strmatch(v1,"Age*") | strmatch(v1,"* age*")) & v1 != "Age in years"
		replace ordvalues = 88 if v2 == "80+" 
		replace ordvalues = 99 if v2 == "60+" 
		egen newn = group(ordvar ordvalues)
		replace n = newn
		sort n
		drop ordvar ordvalues newn
	rename v1 Variable
	rename v2 Values
	foreach j in 3 5 7 9 11{
		capture local name = v`j'[1]
		capture rename v`j' `name'_N
		local jp = `j' + 1
		capture rename v`jp' `name'_perc
		}
	
	drop in 1/2
	
	tempfile xxx`datasource'
	save `xxx`datasource'',replace
	local files "`files' `xxx`datasource''"
	capture sort n 
	capture drop n
	save "${dirdtain}Table ${numtab10`datasource'}_b, Cohort characteristics at second dose of COVID-19 vaccine `datasource'",replace
	}	
	//Table 10, Cohort characteristics at first dose of COVID-19 vaccine ES_BIFAP_PC_HOSP

/*input from DAPs: Table 11, COVID-19 vaccination by dose and time period between first and second dose (days)${suffix`datasource'}.csv*/
/*output: Table 11, COVID-19 vaccination by dose and time period between first and second dose (days) */
local files ""
foreach datasource of global datasourcesdescr {
	import delimited using "${dirinput`datasource'}Dummy tables October/Table 11, COVID-19 vaccination by dose and time period between first and second dose (days)${suffix`datasource'}.csv",clear case(preserve)  varnames(nonames)    
	if "`datasource'" == word("$datasources",1){
		gen n = _n
		}
	rename v1 Dose
	rename v2 Measure
	rename v3 `datasource'_N
	rename v4 `datasource'_perc
	drop in 1/2
	gen check = ""
	foreach vaccine in AstraZeneca Janssen Moderna Pfizer UKN{
		replace check = "`vaccine'" if strmatch(Dose,"*`vaccine'*")
		}
	replace check = check[_n-1] if mi(check)
	isid check Dose Measure
	replace `datasource'_N = "<5" if real(`datasource'_N) > 0 & real(`datasource'_N) < 5 
	tempfile xxx`datasource'
	save `xxx`datasource'',replace
	local files "`files' `xxx`datasource''"
	}
foreach datasource of global datasources {
	if "`datasource'" == word("$datasources",1){
		use `xxx`datasource'',clear
		}
	else{
		merge 1:1  check Dose Measure using `xxx`datasource'',keep(master match)
		drop _merge
		}
	}
sort n 
drop n check
save "${dirdtain}Table 11, COVID-19 vaccination by dose and time period between first and second dose (days)",replace


/*input from DAPs: QC_code_counts_in_study_population_`AESI'_narrow_`year'.csv*/
/*output: AESI_coe_counts */
set more off
local files ""
foreach datasource of global datasources {
	foreach AESI in $eventsnoCOVID MISCC{
		foreach year in 2020 2021{
			capture confirm file "${dirinput`datasource'}QC_code_counts_in_study_population_`AESI'_narrow_`year'.csv"
			if _rc == 0{
				import delimited using "${dirinput`datasource'}QC_code_counts_in_study_population_`AESI'_narrow_`year'.csv",clear case(preserve)  stringcols(_all)   
				gen datasource = "`datasource'"
				gen AESI = "`AESI'"
				gen year = `year'
				tempfile xxx
				save `xxx',replace
				local files "`files' `xxx'"
				}
			}
		}
	}
dsconcat `files'
order  AESI datasource meaning coding code year
sort AESI datasource meaning coding code year
save ${dirdtain}AESI_code_counts,replace




/*input from DAPs: Table XX, Incidence rates of AESI by vaccine and datasource${suffix`datasource'}.csv*/
/*output: IR_AESI*/
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}Dummy tables October/Table XX, Incidence rates of AESI by vaccine and datasource${suffix`datasource'}.csv",clear case(preserve)
	gen datasource = "`datasource'"
	rename IR IR 
	rename lb lb
	rename ub ub
	rename Event AESI
	foreach AESI of global AESIexcluded`datasource'{
		drop if AESI == "`AESI'"
		}
	tempfile xxx
	save `xxx',replace
	local files "`files' `xxx'"
	}
dsconcat `files'
save ${dirdtain}IR_AESI,replace


/*input from DAPs: Table 12, Number of incident cases entire study period${suffix`datasource'}.csv*/
/*output: incident_cases_AESI*/
local files ""
foreach datasource of global datasources {
	import delimited using "${dirinput`datasource'}Dummy tables October/Table 12, Number of incident cases entire study period${suffix`datasource'}.csv",clear case(preserve)
	rename v1 AESIyear
	rename v2 count`datasource'
	drop in 1/2 
	gen AESI = AESIyear[_n-1] if AESIyear == "2020"
	replace AESI = AESIyear[_n-2] if AESIyear == "2021"
	drop if count`datasource' == ""
	foreach AESI of global AESIexcluded`datasource'{
		drop if AESI == "`AESI'"
		}
	rename AESIyear year
	tempfile xxx`datasource'
	save `xxx`datasource'',replace
	local files "`files' `xxx`datasource''"
	}
foreach datasource of global datasources {
	if "`datasource'" == word("$datasources",1){
		use `xxx`datasource'',clear
		}
	else{
		merge 1:1  AESI year using `xxx`datasource'',keep(master match)
		drop _merge
		}
	}
save ${dirdtain}incident_cases_AESI,replace



/*standardised rates from folder 2*/


local files ""
foreach datasource in ARS BIFAP_PC BIFAP_PC_HOSP CPRD PHARMO{
	local namedat = cond(strmatch("`datasource'","BIFAP*") == 0,"`datasource'", subinstr("`datasource'","BIFAP_PC","BIFAPPC",1))
	import delimited using ${dirinputpoisson}RCodes\ObjectUsed/DAP_`namedat'_StandardizedRates.csv,clear case(preserve) varnames(1)
	gen datasource = "`datasource'"
	foreach AESI of global AESIexcluded`datasource'{
		drop if strmatch(AESI, "`AESI'*")
		}
	tempfile xxx
	save `xxx', replace
	local files "`files' `xxx'"
	}
dsconcat `files'
capture drop v1
save ${dirdtain}StandardizedRates,replace

/*rate differences from folder 2*/

local files ""
foreach datasource in ARS BIFAP_PC BIFAP_PC_HOSP CPRD PHARMO{
	local namedat = cond(strmatch("`datasource'","BIFAP*") == 0,"`datasource'", subinstr("`datasource'","BIFAP_PC","BIFAPPC",1))
	import delimited using ${dirinputpoisson}RCodes\ObjectUsed/DAP_`namedat'_RateDiff.csv,clear case(preserve) varnames(1)
	gen datasource = "`datasource'"
		foreach AESI of global AESIexcluded`datasource'{
		drop if strmatch(AESI, "`AESI'*")
		
	tempfile xxx
	save `xxx', replace
	local files "`files' `xxx'"
	}
dsconcat `files'
capture drop v1
save ${dirdtain}RateDiff,replace
