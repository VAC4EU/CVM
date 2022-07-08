/* this steps executes final steps of statistical analysis, namely*/
/* (1) computation of background persontime per ageband */
/* (2) cumulation of person time per week of follow up after vaccination and estimation of IRs on such cumulative persontime, including both wholeobservetion epriod and 28 dasy*/
/* (3) observed and expected events at 28 days, using age-specific IRs*/

set more off

qui do par.do
qui do folders.do
qui do colors.do

// {}
// `


// BACKGROUND PERSONTIME

set more off
local collapse_sex_ageband "sex ageband"
local collapse_ageband "ageband"
local collapse_all ""
foreach comb in sex_ageband ageband all{
	use ${dirdtain}counts_AESI_and_persontime_BGR_sex_ageband,clear
	rename ageband ageband
	collapse (sum) Persontime* *_b,by(datasource  `collapse_`comb'')
	foreach var of varlist Persontime_*{
		local core: subinstr local var "Persontime_" ""
		rename `var' den`core'
		}
	capture rename denACUASEARTHRITIS_narro denACUASEARTHRITIS_narrow
	foreach var of varlist *_b{
		local core: subinstr local var "_b" ""
		rename `var' num`core'
		}
	reshape long den num, i(datasource Persontime  `collapse_`comb'') j(AESI) string
	replace AESI = subinstr(AESI,"_narrow","",.)
	foreach AESI of global AESIexcluded`datasource'{
		di "`AESI'"
		drop if AESI == "`AESI'"
		drop if AESI == "`AESI'_narrow"
		drop if AESI == "`AESI'_possible"
		}
	save ${dirdtain}count_AESI_and_persontime_BGR_`comb',replace
	}

// INCIDENCE RATES AFTER VACCINATION

/*both during all persontime after vaccination, and truncating at 28 days: per age and sex; then per age only; then all lumped together; also compute persontime per vacination cohort*/
set more off
local collapse_sex_ageband "sex ageband"
local collapse_ageband "ageband"
local collapse_all ""
local filter_entire ""
local filter_28 "keep if week <= 4"
local suffix_entire ""
local suffix_28 "_28"
foreach comb in sex_ageband ageband all{
	foreach filter in entire 28{
		use ${dirdtain}counts_AESI_and_persontime_dose_type_vax_week_fup_sex_ageband,clear
		rename ageband ageband
		`filter_`filter''
		collapse (sum) Persontime* *_b,by(datasource  Dose type_vax `collapse_`comb'')
		preserve
		keep datasource  Dose type_vax `collapse_`comb'' Persontime
		save ${dirdtain}Persontime_type_vax_dose_`comb'`suffix_`filter'',replace
		restore
		foreach var of varlist Persontime_*{
			local core: subinstr local var "Persontime_" ""
			rename `var' den`core'
			}
		capture rename denACUASEARTHRITIS_narro denACUASEARTHRITIS_narrow
		foreach var of varlist *_b{
			local core: subinstr local var "_b" ""
			rename `var' num`core'
			}
		reshape long den num, i(datasource Persontime Dose type_vax `collapse_`comb'') j(AESI) string
		global K = 100000*365.25
		gen IR = ${K}*(num/den)
		gen invl =invchi2(2*num,.025)/2
		gen invu =invchi2(2*(num+1),.975)/2
		gen LL =(${K}/den)*invl
		gen UL = (${K}/den)*invu
		replace LL = 0 if mi(LL) & !mi(UL)
		order AESI datasource  IR LL UL
		gen Personyears = round(den/365.25)
		replace Personyears = . if Personyears == 0
		// tostring num,gen(Cases)
		foreach var in IR LL UL{
			rename `var' `var'
			format `var' %9.2f
			tostring `var', force replace use
			replace `var' = "" if `var' == "."
			}
		drop inv*
		replace AESI = subinstr(AESI,"_narrow","",.)
		foreach AESI of global AESIexcluded`datasource'{
			di "`AESI'"
			drop if AESI == "`AESI'"
			drop if AESI == "`AESI'_narrow"
			drop if AESI == "`AESI'_possible"
			}
		save ${dirdtain}IRs_AESI_and_persontime_dose_type_vax_`comb'`suffix_`filter'',replace
		}
	}


/*CREATE EUROPEAN STANDARD POPULATION*/

clear
set obs 96
gen age = _n -1
local agebandsESP "0 1 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 150"
cut_ageband age, at(`agebandsESP') gen(agebandESP)
gen popESP = 1000 if agebandESP == 0
replace popESP = 4000 if agebandESP == 1
replace popESP = 5500 if agebandESP == 5
replace popESP = 5500 if agebandESP == 10
replace popESP = 5500 if agebandESP == 15
replace popESP = 6000 if agebandESP == 20
replace popESP = 6000 if agebandESP == 25
replace popESP = 6500 if agebandESP == 30
replace popESP = 7000 if agebandESP == 35
replace popESP = 7000 if agebandESP == 40
replace popESP = 7000 if agebandESP == 45
replace popESP = 7000 if agebandESP == 50
replace popESP = 6500 if agebandESP == 55
replace popESP = 6000 if agebandESP == 60
replace popESP = 5500 if agebandESP == 65
replace popESP = 5000 if agebandESP == 70
replace popESP = 4000 if agebandESP == 75
replace popESP = 2500 if agebandESP == 80
replace popESP = 1500 if agebandESP == 85
replace popESP = 800 if agebandESP == 90
replace popESP = 200 if agebandESP == 95

bysort agebandESP: gen numyears = _N
gen popESPyear = popESP/numyears

//global agebands `"80+ 70-79 60-69 50-59 40-49 30-39 25-29 18-24 12-17 5-11 0-4"'

cut_ageband age, at(0 5 12 18 25 30 40 50 60 70 80 150) gen(ageband_num)
decode ageband_num,gen(ageband) 
save ${dirdtain}ESP_age,replace
use ${dirdtain}ESP_age,clear
collapse (sum) popESPyear,by(ageband)
rename popESPyear popESP
save ${dirdtain}ESP_ageband,replace
use ${dirdtain}ESP_ageband,clear
gen sex = 1
append using ${dirdtain}ESP_ageband
replace sex = 0 if mi(sex)
replace pop = pop/2
save ${dirdtain}ESP_ageband_sex,replace



use ${dirdtain}ESP_age,clear
outsheet using ${dircsv}ESP_age.csv,replace

use ${dirdtain}ESP_ageband,clear
outsheet using ${dircsv}ESP_ageband.csv,replace

use ${dirdtain}ESP_ageband_sex,clear
outsheet using ${dircsv}ESP_ageband_sex.csv,replace



/*calculate expected cases of AESI within 28 days, both crude and per age*/

use ${dirdtain}BGR_IR_AESI_ageband,clear
replace AESI = subinstr(AESI,"_narrow","",.)
drop if strmatch(AESI,"*broad*")
preserve
keep if ageband == "all_birth_cohorts"
keep datasource  AESI IR
isid datasource AESI
joinby datasource using ${dirdtain}Persontime_type_vax_dose_all_28

gen expected_crude = (Persontime/365.25)*(IR/100000)
rename Persontime Persontime_tot
tempfile crude
save `crude'
restore
joinby datasource ageband using ${dirdtain}Persontime_type_vax_dose_ageband_28
gen expectedstratified = (Persontime/365.25)*(IR/100000)
collapse (sum) expected Persontime,by(datasource AESI Dose type_vax )
rename expected expected_perageband
merge 1:1 datasource AESI Dose type_vax using `crude'
drop _merge IR
merge 1:1 datasource AESI Dose type_vax using ${dirdtain}IRs_AESI_and_persontime_dose_type_vax_all_28,keepus(num)
drop _merge
rename num observed
format expecte* observed %9.0gc
drop Persontime_tot
order datasource AESI  type_vax Dose Persontime
sort datasource AESI  type_vax Dose
save ${dirdtain}observed_and_expected_28,replace

/*compute incidence cumuulative per week of fup after vaccination*/
set more off
use ${dirdtain}counts_AESI_and_persontime_dose_type_vax_week_fup_sex_ageband,clear
drop Persontime sex 
gen ageband_coarse = age == "80+" | age == "70-79" | age == "60-69"
collapse (sum) Persontime* *_b,by(ageband_coarse type_vax week_fup Dose  datasource)
gen Age = "0-59" if ageband_coarse == 0
replace Age = "60+" if ageband_coarse == 1
foreach var of varlist Persontime*{
	local core: subinstr local var "Persontime_" ""
	local core = subinstr("`core'","_narrow","",.)
	rename `var' den`core'
	}
foreach var of varlist *_b{
	local core: subinstr local var "_b" ""
	local core = subinstr("`core'","_narrow","",.)
	rename `var' num`core'
	}
reshape long den num, i(Age type_vax week_fup Dose   datasource) j(AESI) string
	foreach AESI of global AESIexcluded`datasource'{
		di "`AESI'"
		drop if AESI == "`AESI'"
		drop if AESI == "`AESI'_narrow"
		drop if AESI == "`AESI'_possible"
		}

/*discard weeks where <5 people are observed*/
drop if den < 35
bysort AESI Age  type_vax Dose   datasource:egen maxweek = max(week)
local listfiles  ""
qui sum week
local N = `r(max)'
forvalues week = 1/`N'{
	preserve
	keep if week <=`week'
	gen extraweek = `week' > maxweek  
	bysort AESI Age  type_vax Dose   datasource: egen todrop = max(extraweek) 
	drop if todrop
	collapse (sum) den* num* ,by(AESI Age  type_vax Dose datasource)
	gen cumulative_week_fup = `week'
	tempfile week
	save `week'
	local listfiles "`listfiles' `week'"
	restore
	}
dsconcat `listfiles'
rename type_vax manufacturer
rename Dose dose
global K = 100000*365.25
gen IR = ${K}*(num/den)
gen invl =invchi2(2*num,.025)/2
gen invu =invchi2(2*(num+1),.975)/2
gen LL =(${K}/den)*invl
gen UL = (${K}/den)*invu
replace LL = 0 if mi(LL) & !mi(UL)
order AESI   Age IR LL UL
gen Personyears = round(den/365.25)
foreach datasource of global datasources {
	foreach AESI of global AESIexcluded`datasource'{
		di "`AESI'"
		drop if AESI == "`AESI'" & datasource == "`datasource'"
		}
	}
drop inv* 
save ${dirdtain}IR_fup_weeks_cumulative,replace
 	


/*age-standardised BGR*/

use ${dirdtain}BGR_IR_AESI_ageband,clear
replace AESI = subinstr(AESI,"_narrow","",.)
drop if strmatch(AESI,"*broad*")
drop if ageband == "all_birth_cohorts" | ageband == "60+" | ageband == "0-59" 
merge 1:1 datasource AESI ageband using ${dirdtain}count_AESI_and_persontime_BGR_ageband,keepus(num den) keep(match)
foreach datasource of global datasources {
	foreach AESI of global AESIexcluded`datasource'{
		di "`AESI'"
		drop if AESI == "`AESI'" & datasource == "`datasource'"
		}
	}
drop _merge
rename den persontime
gen PY = persontime/365.25
distrate num PY using ${dirdtain}ESP_ageband,stand(ageband) popstand(popESP) by(datasource AESI) mult(100000) format(%9.2f) saving(${dirdtain}ageband_standardised_BGR,replace)

/*age-standardised rates within 28 days from vaccine*/

use ${dirdtain}IRs_AESI_and_persontime_dose_type_vax_ageband_28,clear
rename Personyears PY
drop if strmatch(AESI, "CONTR*")
foreach datasource of global datasources {
	foreach AESI of global AESIexcluded`datasource'{
		di "`AESI'"
		drop if AESI == "`AESI'" & datasource == "`datasource'"
		}
	}
levelsof AESI,local(AESIs)
local listfiles ""
foreach AESI of local AESIs{
	preserve
	keep if AESI == "`AESI'"
	tempfile xxx
	distrate num PY using ${dirdtain}ESP_ageband,stand(ageband) popstand(popESP) by(datasource AESI type_vax Dose) mult(100000) format(%9.2f) saving(`xxx',replace)
	restore
	local listfiles "`listfiles' `xxx'"

	}
dsconcat `listfiles' 
save ${dirdtain}ageband_standardised_IR_type_vax_Dose_28,replace


// {}
// `






// {}
// `

