*this file is to make the latent to active analysis, basically, using Stata

clear all
set more off
set maxvar 10000


forvalues runNum = 1/1 {

if `runNum' == 1 {
	clear
	cd C:\Users\Sze\Desktop\localTBproj\Mar2\inf0p0020_lat2p55_aveTreat_fullCatIV_empUptake_noTB_noSmokingStrats\2012-03-06_17-21-32
	disp "aveTreat_fullCatIV_empUptake"
}
if `runNum' == 2 {
	clear
	cd C:\Users\Sze\Desktop\localTBproj\Mar2\inf0p0020_lat2p55_aveTreat_fullCatIV_empUptake_noTB_noSmokingStrats\2012-03-06_15-51-12
	disp "C:\Users\Sze\Desktop\localTBproj\Feb27\inf0p0020_lat2p55_aveTreat_fullCatIV_empUptake_noTB"
	
}
if `runNum' == 3 {
	clear
	cd C:\Users\Sze\Desktop\localTBproj\Mar2\inf0p0020_lat2p55_aveTreat_fullCatIV_empUptake_noTB_noSmokingStrats\2012-03-06_15-58-45
	disp "C:\Users\Sze\Desktop\localTBproj\Feb27\inf0p0020_lat2p55_aveTreat_fullCatIV_empUptake_noSmoke"
}
if `runNum' == 4 {
	clear
	cd 
	disp "C:\Users\Sze\Desktop\localTBproj\Feb27\inf0p0020_lat2p55_aveTreat_fullCatIV_empUptake_noTB_noSmoke"
}


log using lifeExpStata_noSmokingStrat.log, replace

*cd C:\Users\Sze\Dropbox\TBproject\code\outputs\Feb27\inf0p0020_lat2p55_aveTreat_fullCatIV_empUptake_noTB\2012-02-27_12-02-02

insheet using "recNo.csv", comma
count
local totNumVars = 1199


*name the variables
rename varnamesgoeshere t_1
forvalues varNum = 2/`totNumVars' {
	rename v`varNum' t_`varNum'
}


*generate health Matrix
/*
forvalues varNum = 1/`totNumVars'{
	gen health_`varNum' = floor( (  mod(t_`varNum',48) )/6 )
}
*/

forvalues varNum = 1/`totNumVars'{
	gen health_`varNum' = floor( (  mod(t_`varNum',48) )/6 )
	gen smoking_`varNum' = floor(mod(t_`varNum',43632)/14544)
	gen age_`varNum' =   floor( mod(t_`varNum',14544)/144 ) 
	gen trtmtVal_`varNum' = floor( mod(t_`varNum',6) )
	gen sexMat_`varNum' = (1/48)*( floor((mod(t_`varNum',144))) - 6*health_`varNum' - trtmtVal_`varNum')
}


*mkdir graphs

*forvalues varNum = 1(60)`totNumVars'{
*	histogram health_`varNum', width(0.5) start(0) xscale( range(0 7))
*	gr export graph_`varNum'.png, replace
*}

*create a cohort of people who were born in the first three years
gen cohort1 = 0
forvalues cohortNum = 1/35 {
	local cohortNumPlus = `cohortNum'+1
	replace cohort1 = 1 if health_`cohortNum' == 5 & health_`cohortNumPlus' == 0
}

keep if cohort1 == 1

*make some proportion graphs
*preserve
*keep if cohort1 == 1
*forvalues varNum = 1(60)`totNumVars'{
*	histogram health_`varNum', width(0.5) start(0) xscale( range(0 7)) yscale(range( 0 2))
*	gr export graph_`varNum'.png, replace
*}
*restore

*num people with TB
gen everGetActTB = 0
gen everGetLatTB = 0
foreach varName of varlist health*  {
	replace everGetActTB = 1 if `varName' == 3 | `varName' == 3
	replace everGetLatTB = 1 if `varName' == 1 | `varName' == 2
}

count if everGetActTB
count if everGetLatTB
count
count if everGetActTB & cohort1
count if everGetLatTB & cohort1


*when did they die, and what of
gen bornTime = 0
gen deathTime = 0
gen diedOf = 0
gen deathAge = 0
gen sex1 = 0
local totVarsOneLess = `totNumVars' - 1
forvalues varNum = 2/`totVarsOneLess' {
	local varBefore = `varNum' - 1
	local varAfter = `varNum' + 1
	replace bornTime = `varNum' if health_`varBefore' == 5 & inlist(health_`varNum',0,1,2,3,4)
	replace sex1 = sexMat_`varNum' if bornTime == `varNum'
	replace deathTime = `varNum' if inlist(health_`varNum',0,1,2,3,4) & health_`varAfter' == 6 
	replace diedOf = health_`varNum' if deathTime == `varNum'
	replace deathAge = age_`varNum' if deathTime == `varNum'
}


*browse bornTime deathTime diedOf health* if cohort1 == 1

*everyone who had actTB dies with it, which makes sense since there is no cure
count if everGetActTB == 1 & cohort1 == 1
count if diedOf == 3 & everGetActTB == 1 & cohort1 == 1

gen lifeLength = deathTime-bornTime
mean(lifeLength) if cohort1 == 1
*668.4285 is about 56 years, which is a bit low

*generate a mortality table by age
tab deathAge if cohort1 == 1
tab deathAge if cohort1 == 1 & sex1 == 1
tab deathAge if cohort1 == 1 & sex1 == 2


save recNoStata.dta, replace


log close


}
