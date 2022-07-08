*! Author: Rosa Gini
*! Date: 6 August 2008
program define cut_ageband
version 9.0
syntax varlist(max=1) [if] [in] , AT(numlist integer ascending) GENerate(string) [CLOSEDend LABVAR(string) LABelname(string)] 
local age "`varlist'"
tempvar touse
mark `touse' `if' `in'

confirm new variable `generate'
capture confirm numeric variable `age'
if _rc!=0{
	di as err "`age' must be a numeric variable"
	}

if "`labelname'"==""{
	local labelname "`generate'_label"
	}
	
numlist "`at'"
local cutatnew="`r(numlist)'"
qui egen byte `generate'=cut(`age') if `touse',at(`at') 

if "`labvar'"!=""{
	label variable `generate' "`labvar'"	
	}

local j=1
foreach cut of local cutatnew{
	local `j'=`cut'
	local last=`cut'
	local j=`j'+1
	}

local j=1
foreach cut of local cutatnew{
	if ``j''!=`last'{
		local eta=`cut'
		local jp=`j'+1
		local etasec=``jp''-1
		local label=cond(``jp''!=`last'|"`closedend'"!="","`eta'-`etasec'","`eta'+")
		label define `labelname' `eta' "`label'",modify
		}
	local j=`j'+1
	}

label values `generate' `generate'_label

end;
