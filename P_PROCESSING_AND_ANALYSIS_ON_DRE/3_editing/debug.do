set more off
qui do par.do
qui do colors.do

use ${dirdtain}Persontime_ageband_week_str,clear
keep if datasource == "BIFAP_PC_HOSP"
sort week ageband
bysort week:egen totpop = sum( popapprox )
gen distr = 100*popapprox/ totpop
format distr %3.1f
gen ordage = word(subinstr(ageband,"-"," ",.),1)
replace ordage=substr(ordage,1,2)
destring ordage, replace
sort ageband week
sort week ordage
graph bar distr,over( ordage) by(week) graphregion(${graphregion}`marginregion') $subtitle
graph bar popapprox ,over( ordage) by(week) graphregion(${graphregion}`marginregion') $subtitle

