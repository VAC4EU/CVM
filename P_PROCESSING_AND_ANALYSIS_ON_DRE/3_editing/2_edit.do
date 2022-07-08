/* this step edits datasets by adding labels to vaccines, datasource, AESIs, and other categorical variables; */ 
/* moreover it splits and reshapes datasets to form tables or datasets for figures*/

set more off

qui do par.do
qui do folders.do
qui do colors.do

// {}
// `

/*

/*input: attrition_table, CDM_SOURCE*/
/*output: attrition_table_edited*/
use ${dirdtain}attrition_table, clear
gen ord  = 100 + _n
drop in 1
tempfile attrition
save `attrition'
replace ord = _n
keep if ord < 5
drop if ord == 2
replace Condition = "Type of datasource" if ord == 1 
replace ARS = "Record linkage of multiple data banks" if ord == 1 
replace BIFAP_PC = "Primary care medical records" if ord == 1 
replace BIFAP_PC_HOSP = "Primary care medical records linked to hospital administrative records" if ord == 1 
replace CPRD = "Primary care medical records" if ord == 1 
replace PHARMO = "Primary care medical records" if ord == 1 
replace Condition = "Coding systems for diagnoses" if ord == 3 
replace ARS = "ICD9CM" if ord == 3 
replace BIFAP_PC = "ICPC, SNOMED" if ord == 3 
replace BIFAP_PC_HOSP = "ICPC, SNOMED, ICD10CM" if ord == 3 
replace CPRD = "SNOMED, RCD2" if ord == 3 
replace PHARMO = "ICPC" if ord == 3 
replace Condition = "Attrition" if ord == 4
replace ARS = "" if ord == 4
replace BIFAP_PC = "" if ord == 4 
replace BIFAP_PC_HOSP = "" if ord == 4 
replace CPRD = "" if ord == 4
replace PHARMO = "" if ord == 4
tempfile descr
save `descr'
use ${dirdtain}CDM_SOURCE, clear
keep recomme datasource
rename reco _
gen aux = 1
reshape wide _,i(aux) j(datasource) string
drop aux
foreach datasource of global datasources {
//	rename _`datasource' `datasource'
	tostring _`datasource',replace force
	gen `datasource' = mdy(real(substr(_`datasource',5,2)),real(substr(_`datasource',7,2)),real(substr(_`datasource',1,4)))
	format `datasource' %td
	tostring `datasource', replace use force
	drop _`datasource' 
	}
gen ord = 2
gen Condition = "Recommended end date"
append using `attrition'
append using `descr'
order Condition
sort ord 
drop ord
save ${dirdtain}attrition_table_edited, replace


/*input: Persontime_ageband_week_sex*/
/*output: PY_ageband_week_str_sex PY_ageband_week_str PYweek*/
use	${dirdtain}Persontime_ageband_week_sex,replace
replace week = subinstr(week,"-","",.)
destring week,replace
tostring week,gen(week_str)
rename ageband ageband
gen popapprox = Persontime/7
save ${dirdtain}Persontime_ageband_week_str_sex,replace
collapse (sum) popapprox Persontime,by(datasource ageband week*)
save ${dirdtain}Persontime_ageband_week_str,replace
use ${dirdtain}Persontime_ageband_week_str,clear
collapse (sum) Persontime,by(week* datasource)
gen PY = Persontime/365.25
format PY %20.0gc
encode week_str,gen(week_labelled)
gen year = substr( week_str,1,4)
gen month = real(substr( week_str,5,2))
sort week_labelled
save  ${dirdtain}PY_week_str,replace

/*input: Table 6, Cohort characteristics at start of study (1-1-2020)*/
/*output: N_study_population*/
use "${dirdtain}Table 6, Cohort characteristics at start of study (1-1-2020)",clear
keep if Variable == "Study population"
drop *_perc
save  ${dirdtain}N_study_population,replace


/*input: Table 6, Cohort characteristics at start of study (1-1-2020, doses_weeks_manufacturer_ageband, Persontime_ageband_week_str, coverage_weeks_manufacturer_ageband*/
/*output: cumulative_coverage_by_datasource_and_ageband_and_vax*/
use "${dirdtain}Table 6, Cohort characteristics at start of study (1-1-2020)",clear
keep if Variable == `"Age in categories"'
drop Variable
drop *_perc
rename Values ageband
foreach datasource of global datasources {
	rename `datasource'_N pop`datasource'
	}
reshape long pop, i(ageband) j(datasource) string
destring pop,replace
tempfile pop
save `pop'
use ${dirdtain}doses_weeks_manufacturer_ageband,clear
rename N doses
tostring week, gen(week_str)
keep if week_str > "20201130"
capture drop _merge
merge m:1 datasource ageband using `pop',keep(master match)
drop _merge
merge m:1 datasource ageband week using ${dirdtain}Persontime_ageband_week_str,keep(match)
drop _merge
// cum_vx: cumulative number of doses of a vaccine by week x
sort  datasource ageband dose vx_manufacturer  week_str
bysort datasource ageband  dose vx_manufacturer: gen cum_vx = doses if _n == 1 
bysort datasource ageband  dose vx_manufacturer: replace cum_vx = cum_vx[_n-1] + doses if _n > 1 
gen cum_vx_perc = 100*cum_vx/pop
label variable cum_vx_perc "Cumulative coverage, with fixed population" 
gen cum_vx_perc_timedep = 100*cum_vx/popapprox
label variable cum_vx_perc "Cumulative coverage, with time-dependent population approximated using weekly persontime" 
preserve
use ${dirdtain}coverage_weeks_manufacturer_ageband,clear
tostring week, gen(week_str)
keep if week_str > "20201130"
label variable percentage "Cumulative coverage, with time-dependent population as computed by the ECVM script" 
tempfile coverage_from_ECVM_script
save `coverage_from_ECVM_script'
restore
merge 1:1 datasource ageband week_str vx_manufacturer dose using `coverage_from_ECVM_script',keep(match)
drop _merge
encode week_str,gen(week_labelled)
gen year = substr( week_str,1,4)
gen month = real(substr( week_str,5,2))
sort week_labelled
save ${dirdtain}cumulative_coverage_by_datasource_and_ageband_and_vax, replace

*/

use ${dirdtain}IR_covid_calendartime_vax_cohort,clear
gen ordvac = 0
label define ordvac_lab 0 `"None"',modify
local j = 1
foreach vaccine in Pfizer AstraZeneca Moderna J&J UKN{
	replace ordvac = `j' if strmatch(vx,"*`vaccine'*")
	label define ordvac_lab `j' `"${name`vaccine'}"',modify
	local j = `j' +1
	}
label value ordvac ordvac_lab
tostring week,gen(week_str)
encode week_str,gen(week_labelled)
gen year = substr( week_str,1,4)
gen month = real(substr( week_str,5,2))
sort week_labelled
gen ordDAP = .
local j = 1
foreach datasource of global datasources{
	replace ordDAP = `j' if datasource == `"`datasource'"'
	label define ordDAP_lab `j' `"${name`datasource'}"' ,modify
	local j  = `j' + 1
	}
label values ordDAP ordDAP_lab
save ${dirdtain}IR_covid_calendartime_vax_cohort_labelled,replace


/*
use ${dirdtain}IR_covid_timesincevaccination_vax_cohort,clear
gen ordvac = 0
label define ordvac_lab 0 `"None"',modify
local j = 1
foreach vaccine in Pfizer AstraZeneca Moderna J&J UKN{
	replace ordvac = `j' if strmatch(vx,"*`vaccine'*")
	label define ordvac_lab `j' `"${name`vaccine'}"',modify
	local j = `j' +1
	}
save ${dirdtain}IR_covid_timesincevaccination_vax_cohort_labelled,replace





use  ${dirdtain}AESI_code_counts,clear
gen system_num = .
gen AESI_num = .
gen AESI_name = ""
gen tokeepAESI = 0
local j = 1
local s = 1
foreach system of global systems{
	foreach AESI of  global events`system'{
		replace system_num = `s' if AESI == `"`AESI'"'
		replace AESI_num = `j' if AESI == `"`AESI'"'
		replace AESI_name = `"${name`AESI'}"' if AESI == `"`AESI'"'
		replace tokeepAESI = 1 if AESI == `"`AESI'"'
		label define AESI_num_lab `j' `"${name`AESI'}"',modify
		local j = `j' +1
		}
	label define system_num_lab `s' `"${sysname`system'}"',modify
	local s = `s' +1
	}
label values system_num system_num_lab
label values AESI_num AESI_num_lab
gen DAPname = "" 
foreach datasource of global datasources{
	replace DAPname = `"${name`datasource'}"' if datasource == `"`datasource'"'
	}
save  ${dirdtain}AESI_code_counts_labelled,replace



/*input incident_cases_AESI*/
/*output: incident_cases_AESI_labelled*/
use ${dirdtain}incident_cases_AESI,replace
append using ${dirdtain}incident_cases_AESI_DEATH
gen system_num = .
gen AESI_num = .
gen AESI_name = ""
gen tokeepAESI = 0
local j = 1
local s = 1
foreach system of global systems{
	foreach AESI of  global events`system'{
		replace system_num = `s' if AESI == `"`AESI'"'
		replace AESI_num = `j' if AESI == `"`AESI'"'
		replace AESI_name = `"${name`AESI'}"' if AESI == `"`AESI'"'
		replace tokeepAESI = 1 if AESI == `"`AESI'"'
		label define AESI_num_lab `j' `"${name`AESI'}"',modify
		local j = `j' +1
		}
	label define system_num_lab `s' `"${sysname`system'}"',modify
	local s = `s' +1
	}
label values system_num system_num_lab
label values AESI_num AESI_num_lab
sort system_num AESI_num year
order system_num AESI_num year count*
foreach datasource of global datasources{
	foreach AESI of global AESIexcluded`datasource'{
		di "`AESI'"
		replace count`datasource' = "" if AESI == "`AESI'"
		}
	}

save ${dirdtain}incident_cases_AESI_labelled,replace



// INCIDENCE RATES AFTER VACCINATION




/*for figure: cumulative IR per week*/
use ${dirdtain}BGR_IR_AESI_ageband,clear
drop if strmatch(AESI, "*broad*")
replace AESI = subinstr(AESI,"_narrow","",.)
rename ageband Age
rename IR IR_BGR
rename lb lb_BGR
rename ub ub_BGR
merge 1:m  datasource Age AESI using ${dirdtain}IR_fup_weeks_cumulative,keep(using match)
drop _merge
gen ordDAP = .
local j = 1
foreach datasource of global datasources{
	replace ordDAP = `j' if datasource == `"`datasource'"'
	label define ordDAP_lab `j' `"${name`datasource'}"' ,modify
	local j  = `j' + 1
	}
label values ordDAP ordDAP_lab
replace manufacturer = "Janssen" if manufacturer == "J&J"
rename manufacturer Vaccine
gen ordvac = .
local j = 1
foreach vaccine in Pfizer AstraZeneca Moderna Janssen   UKN{
	replace ordvac = `j' if strmatch(Vaccine,"*`vaccine'*")
	label define ordvac_lab `j' `"${name`vaccine'}"',modify
	local j = `j' +1
	}
label values ordvac ordvac_lab
local j = 1
local s = 1
gen system_num = .
gen AESI_num = .
gen AESI_name = ""
gen tokeepAESI = 0
foreach system of global systems{
	foreach AESI of  global events`system'{
		replace system_num = `s' if AESI == `"`AESI'"'
		replace AESI_num = `j' if AESI == `"`AESI'"'
		replace AESI_name = `"${name`AESI'}"' if AESI == `"`AESI'"'
		replace tokeepAESI = 1 if AESI == `"`AESI'"'
		label define AESI_num_lab `j' `"${name`AESI'}"',modify
		local j = `j' +1
		}
	label define system_num_lab `s' `"${sysname`system'}"',modify
	local s = `s' +1
	}
label values system_num system_num_lab
label values AESI_num AESI_num_lab
sort system_num AESI_num 
order system_num AESI_num
foreach datasource of global datasources {
	foreach AESI of global AESIexcluded`datasource'{
		di "`AESI'"
		drop if AESI == "`AESI'" & datasource == "`datasource'"
		}
	}
keep if !mi(AESI_num)
isid system_num AESI_num ordvac dose datasource Age cumulat
save ${dirdta}IR_fup_weeks_cumulative_and_BGR,replace



/*during 28 days after vaccination, IR per age and sex; then per age only; then all lumped together*/
set more off
local collapse_sex_ageband "sex ageband"
local collapse_ageband "ageband"
local collapse_all ""
foreach comb in sex_ageband ageband all{
	use ${dirdtain}IRs_AESI_and_persontime_dose_type_vax_`comb'_28,clear
	foreach datasource of global datasources {
	foreach AESI of global AESIexcluded`datasource'{
		di "`AESI'"
		drop if AESI == "`AESI'" & datasource == "`datasource'"
		}
	}
	gen system_num = .
	gen AESI_num = .
	gen AESI_name = ""
	gen tokeepAESI = 0
	local j = 1
	local s = 1
	foreach system of global systems{
		foreach AESI of  global events`system'{
			replace system_num = `s' if AESI == `"`AESI'"'
			replace AESI_num = `j' if AESI == `"`AESI'"'
			replace AESI_name = `"${name`AESI'}"' if AESI == `"`AESI'"'
			replace tokeepAESI = 1 if AESI == `"`AESI'"'
			label define AESI_num_lab `j' `"${name`AESI'}"',modify
			local j = `j' +1
			}
		label define system_num_lab `s' `"${sysname`system'}"',modify
		local s = `s' +1
		}
	label values system_num system_num_lab
	label values AESI_num AESI_num_lab
	keep if tokeepAESI
	drop tokeepAESI
	gen DAPname = ""
	foreach datasource of global datasources{
		replace DAPname = `"${name`datasource'}"' if datasource == `"`datasource'"'
		}
	gen ordvac = .
	local j = 1
	rename type_vax Vaccine
	foreach vaccine in Pfizer AstraZeneca Moderna J&J   UKN{
		replace ordvac = `j' if strmatch(Vaccine,"*`vaccine'*")
		label define ordvac_lab `j' `"${name`vaccine'}"',modify
		if "`vaccine'" == "J&J"{
			label define ordvac_lab `j' `"${nameJanssen}"',modify
			}
		local j = `j' +1
		}
	label values ordvac ordvac_lab
	sort DAPname AESI ordvac
	order system_num AESI_num DAPname ordvac Dose `collapse_`comb''
	sort  system_num AESI_num DAPname ordvac Dose `collapse_`comb''
	save ${dirdtain}IRs_AESI_and_persontime_dose_type_vax_`comb'_28_labelled,replace
	}



/*input observed_and_expected_28*/
/*output: observed_and_expected_28_labelled*/
set more off
use ${dirdtain}observed_and_expected_28,clear
foreach datasource of global datasources{
	foreach AESI of global AESIexcluded`datasource'{
		di "drop `AESI' in `datasource'"
		drop if AESI == "`AESI'" & datasource == "`datasource'"
		}
	}
gen system_num = .
gen AESI_num = .
gen AESI_name = ""
gen tokeepAESI = 0
local j = 1
local s = 1
foreach system of global systems{
	foreach AESI of  global events`system'{
		replace system_num = `s' if AESI == `"`AESI'"'
		replace AESI_num = `j' if AESI == `"`AESI'"'
		replace AESI_name = `"${name`AESI'}"' if AESI == `"`AESI'"'
		replace tokeepAESI = 1 if AESI == `"`AESI'"'
		label define AESI_num_lab `j' `"${name`AESI'}"',modify
		local j = `j' +1
		}
	label define system_num_lab `s' `"${sysname`system'}"',modify
	local s = `s' +1
	}
label values system_num system_num_lab
label values AESI_num AESI_num_lab
keep if tokeepAESI
drop tokeepAESI
gen DAPname = ""
foreach datasource of global datasources{
	replace DAPname = `"${name`datasource'}"' if datasource == `"`datasource'"'
	}
gen ordvac = .
local j = 1
rename type_vax Vaccine
foreach vaccine in Pfizer AstraZeneca Moderna J&J   UKN{
	replace ordvac = `j' if strmatch(Vaccine,"*`vaccine'*")
	label define ordvac_lab `j' `"${name`vaccine'}"',modify
	if "`vaccine'" == "J&J"{
		label define ordvac_lab `j' `"${nameJanssen}"',modify
		}
	local j = `j' +1
	}
label values ordvac ordvac_lab
rename Vaccine type_vax
gen expected_crude_rounded = round(expected_crude)
gen expected_perageband_rounded =round(expected_perageband)
order system_num AESI_num DAPname ordvac Dose expected_crude_rounded expected_perageband_rounded  observed
sort  system_num AESI_num DAPname ordvac Dose 
save ${dirdtain}observed_and_expected_28_labelled,replace


/*merge with std IRs calculated with stata*/
use ${dirdtain}BGR_IR_AESI_ageband,clear
replace AESI = subinstr(AESI,"_narrow","",.)
drop if strmatch(AESI,"*broad*")
keep if ageband == "all_birth_cohorts"
keep datasource  AESI IR lb ub
isid datasource  AESI
tempfile BGR
save `BGR'
use ${dirdtain}observed_and_expected_28_labelled,clear
isid datasource AESI type_vax Dose
merge 1:1 datasource AESI type_vax Dose using ${dirdtain}ageband_standardised_IR_type_vax_Dose_28,keepus(rateadj lb_fay ub_fay) keep(master match)
rename rateadj IR_std
rename lb_fay lb_std
rename ub_fay ub_std
drop _merge
merge m:1 datasource AESI using ${dirdtain}ageband_standardised_BGR,keepus(rateadj lb_fay ub_fay) keep(master match)
rename rateadj IR_BGR_std
rename lb_fay lb_BGR_std
rename ub_fay ub_BGR_std
drop _merge
merge m:1 datasource AESI using `BGR', keep(master match)
drop _merge
rename IR IR_BGR
rename lb lb_BGR
rename ub ub_BGR
gen diff_from_direct_stand =  IR_std - IR_BGR_std
foreach var of varlist IR* ub* lb*{
	tostring `var',force replace format("%9.1f")
	}
gen IRCI_BGR = IR_BGR + " (" + lb_BGR + ";" + ub_BGR + ")" 
gen IRCI_BGR_std = IR_BGR_std + " (" + lb_BGR_std + ";" + ub_BGR_std + ")" 
gen IRCI_std = IR_std + " (" + lb_std + ";" + ub_std + ")" 
drop if type_vax == "J&J" & Dose == 2
save ${dirdtain}IR_std_and_observed_and_expected_28_labelled,replace

/*RATE DIFFERENCES*/

// {}
// `


use ${dirdtain}StandardizedRates,clear
replace AESI = subinstr(AESI,"_b","",.)
replace AESI = subinstr(AESI,"_narrow","",.)
replace Dose = word(subinstr(Dose,"Dose","",.),1)
destring Dose,replace	
gen ordvac = .
local j = 1
foreach vaccine in Pfizer AZ Moderna JJ   UKN{
	replace ordvac = `j' if strmatch(Subgroup,"*`vaccine'*")
	label define ordvac_lab `j' `"${name`vaccine'}"',modify
	if "`vaccine'" == "JJ"{
		label define ordvac_lab `j' `"${nameJanssen}"',modify
		}
	local j = `j' +1
	}
label values ordvac ordvac_lab
replace ordvac = 0 if mi(ordvac)
drop if Dose == 1 & Subgroup == "Baseline"
drop if Dose == 2 & Subgroup == "Baseline"
drop if Dose == 12
preserve
keep if Subgroup == "Baseline"
rename StdRateper1e05 IR_BGR_std_R
rename LCLStd lb_BGR_std_R
rename UCLStd ub_BGR_std_R
keep AESI datasource *BGR*
contract *
drop _freq
isid AESI datasource
tempfile BGR
save `BGR'
restore
drop if Subgroup == "Baseline"
drop if mi(Dose)
rename StdRateper1e05 IR_std_R
rename LCLStd lb_std_R
rename UCLStd ub_std_R
drop if IR_std == "NA"
keep AESI Dose ordvac datasource *std*
isid AESI Dose ordvac datasource
merge m:1 AESI datasource using `BGR',keep(master match using)
drop _merge
foreach var of varlist IR* lb* ub*{
	destring `var', replace force
	tostring `var', replace force format(%9.1f)
	replace `var' = "" if `var' == "."
	}
replace IR_BGR = "" if mi(lb_BGR)
replace IR_std = "" if mi(lb_std) 
isid AESI datasource Dose ordvac
gen IRCI_BGR_std_R = IR_BGR_std_R + " (" + lb_BGR_std_R + ";" + ub_BGR_std_R + ")" if !mi(lb_BGR_std_R)
gen IRCI_std_R = IR_std_R + " (" + lb_std_R + ";" + ub_std_R + ")" if !mi(lb_std_R)
save ${dirdtain}StandardizedRates_labelled,replace


use ${dirdtain}RateDiff,clear
replace AESI = subinstr(AESI,"_b","",.)
replace AESI = subinstr(AESI,"_narrow","",.)
gen Dose = real(word(subinstr(DoseBackground,"Dose","",.),1))
// replace Dose = word(subinstr(Dose,"Dose","",.),1)
// destring Dose,replace
gen ordvac = .
local j = 1
foreach vaccine in Pfizer AZ Moderna JJ   UKN{
	replace ordvac = `j' if strmatch(Vaccine,"*`vaccine'*")
	label define ordvac_lab `j' `"${name`vaccine'}"',modify
	if "`vaccine'" == "JJ"{
		label define ordvac_lab `j' `"${nameJanssen}"',modify
		}
	local j = `j' +1
	}
label values ordvac ordvac_lab
rename lower lb_RD
rename upper ub_RD
merge 1:1 AESI datasource Dose ordvac using ${dirdtain}StandardizedRates_labelled,keep(master match using)
drop _merge
/*set to missing RD if IR_std is missing*/
foreach var in RD lb_RD ub_RD{
	replace `var' = "" if mi(IR_std_R)
	destring `var', replace force
	tostring `var', replace force format(%9.1f)
	replace `var' = "" if `var' == "."
	}
drop DoseBac
gen RDCI = RD + " (" + lb_RD + ";" + ub_RD + ")" if !mi(lb_RD)
drop if Dose == 12
isid AESI datasource Dose ordvac

save ${dirdtain}RateDiff_labelled,replace



set more off
use ${dirdtain}BGR_IR_AESI_ageband,clear
replace AESI = subinstr(AESI,"_narrow","",.)
drop if strmatch(AESI,"*broad*")
keep if ageband == "all_birth_cohorts"
keep datasource  AESI IR lb ub
isid datasource  AESI
rename IR IR_BGR
rename lb lb_BGR
rename ub ub_BGR
foreach var of varlist IR* ub* lb*{
	tostring `var',force replace format("%9.1f")
	}
gen IRCI_BGR = IR_BGR + " (" + lb_BGR + ";" + ub_BGR + ")" 
tempfile BGR
save `BGR'
use ${dirdtain}observed_and_expected_28_labelled,clear
isid datasource AESI type_vax Dose
merge 1:1 AESI datasource ordvac  Dose using ${dirdtain}RateDiff_labelled,keep(master match using)
drop _merge
merge m:1 datasource AESI using `BGR', keep(master match )
drop _merge
//keep AESI ordvac Dose IRCI_BGR IRCI_BGR_std IRCI_std   Persontime expected_perageband_rounded expected_crude_rounded observed datasource
gen PY = round(Persontime/365.35)
drop Persontime
drop if ordvac == 4 & Dose == 2
rename IRCI_BGR row_1_IRCI_BGR
rename IRCI_BGR_std_R row_2_IRCI_BGR_std 
bysort AESI datasource: replace row_1 = row_1[1] if mi(row_1_IRCI_BGR)
bysort AESI datasource: replace row_2 = row_2[1] if mi(row_2)
rename PY row_3_PY 
tostring row_3_PY,format(%20.0gc) replace force
//rename expected_crude_rounded row_4_exp_crude
//tostring row_4,format(%9.0gc) replace force
drop  expected_crude_rounded
rename expected_perageband_rounded row_4_exp_std
tostring row_4,format(%9.0gc) replace force
rename observed row_5_obs
tostring row_5,format(%9.0gc) replace force
rename IRCI_std_R row_6_IRCI_std
rename RDCI row_7_RDCI
keep row* AESI ordvac Dose datasource
reshape long row_,i(AESI ordvac Dose datasource) j(name) string
preserve
keep if strmatch(name,"1*") | strmatch(name,"2*")	
drop ordvac Dose
contract *
drop _freq
drop if mi(row)
reshape wide row_,i(AESI name) j(datasource) string
gen ordvac = 0
tempfile header
save `header'
restore
label define ordvac_lab 0 "",modify
label values ordvac ordvac_lab
drop if strmatch(name,"1*") | strmatch(name,"2*")	
reshape wide row_,i(AESI name ordvac Dose) j(datasource) string
append using `header'
gen system_num = .
gen AESI_num = .
gen AESI_name = ""
gen tokeepAESI = 0
local j = 1
local s = 1
foreach system of global systems{
	foreach AESI of  global events`system'{
		replace system_num = `s' if AESI == `"`AESI'"'
		replace AESI_num = `j' if AESI == `"`AESI'"'
		replace AESI_name = `"${name`AESI'}"' if AESI == `"`AESI'"'
		replace tokeepAESI = 1 if AESI == `"`AESI'"'
		label define AESI_num_lab `j' `"${name`AESI'}"',modify
		local j = `j' +1
		}
	label define system_num_lab `s' `"${sysname`system'}"',modify
	local s = `s' +1
	}
label values system_num system_num_lab
label values AESI_num AESI_num_lab
sort AESI_num ordvac Dose name
order AESI_num ordvac Dose name
save ${dirdtain}table_IR_std_and_observed_and_expected_28_labelled,replace






	
//{}

//` 
