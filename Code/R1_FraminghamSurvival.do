//global db3275 "db3275"
global db3275 "danielbelsky"
//***************************************************************************//
//Survival Analysis 
//***************************************************************************//
use "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham.dta", clear


foreach Y in cvd stroke_tia death { 
stset T_`Y', failure(`Y'=1)
}

/*
foreach Y in [list test outcome] {
	local y zpace 
	stset T_`Y', failure(`Y'=1)
	stcox `y' age8 sex B1 B2, cluster(familyid) robust 
	stcox `y' zar_grimage age8 sex B1 B2, cluster(familyid) robust 
	}
*/

//Clock Var macro
global clocks "zpace zpoam zar_horvath zar_hannum zar_phenoage zar_grimage" 
//Cell control macro
global cells "cd4t cd8t nk bcell mono gran cd8pcd28ncd45ran cd8naive cd4naive plasmablast"

foreach Y in cvd chf chd afix dem death cvddeath chddeath stroke stroke_tia { 
	
stset T_`Y', failure(`Y'=1)

matrix M1 = J(1,5,999) 
matrix M2 = J(1,5,999) 
matrix M3 = J(1,5,999) 
foreach y in $clocks{ 
	//M1
	stcox `y' age8 sex B1 B2, cluster(familyid) robust 
		matrix A = exp(_b[`y']) , exp(_b[`y'] - invnormal(0.975)*_se[`y']), exp(_b[`y'] + invnormal(0.975)*_se[`y']), 2*normal(-abs(_b[`y']/_se[`y'])), e(N) 
		matrix rownames A = `y'
		matrix M1 = M1 \ A
	//M2
	stcox `y' age8 sex B1 B2 $cells, cluster(familyid) robust
		matrix A = exp(_b[`y']) , exp(_b[`y'] - invnormal(0.975)*_se[`y']), exp(_b[`y'] + invnormal(0.975)*_se[`y']), 2*normal(-abs(_b[`y']/_se[`y'])), e(N)  
		matrix rownames A = `y'
		matrix M2 = M2 \ A	
	//M3
	stcox `y' age8 sex B1 B2 smk cpd, cluster(familyid) robust
		matrix A = exp(_b[`y']) , exp(_b[`y'] - invnormal(0.975)*_se[`y']), exp(_b[`y'] + invnormal(0.975)*_se[`y']), 2*normal(-abs(_b[`y']/_se[`y'])), e(N) 
		matrix rownames A = `y'
		matrix M3 = M3 \ A	
	}
foreach x in 1 2 3{
	matrix `Y'M`x'=M`x'[2...,1...]
	matrix colnames `Y'M`x'=HR lb ub p N
	matrix list `Y'M`x'
	}

//Clock-adjusted PACE Effects 	
matrix X = J(1,5,999) 	
foreach x in zpoam zar_horvath zar_hannum zar_phenoage zar_grimage { 
	stcox zpace `x' age8 sex B1 B2, cluster(familyid)
		matrix A = (exp(_b[zpace]) , exp(_b[zpace] - invnormal(0.975)*_se[zpace]), exp(_b[zpace] + invnormal(0.975)*_se[zpace]), 2*normal(-abs(_b[zpace]/_se[zpace])), e(N) )  \ (  exp(_b[`x']) , exp(_b[`x'] - invnormal(0.975)*_se[`x']), exp(_b[`x'] + invnormal(0.975)*_se[`x']), 2*normal(-abs(_b[`x']/_se[`x'])), e(N) )
		matrix rownames A = pace `x'
		matrix X = X \ A
		}
	matrix `Y'X=X[2...,1...]
	matrix colnames `Y'X=HR lb ub p N
	matrix list `Y'X
		
putexcel set "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/R1_Framingham.xlsx", sheet(`Y') modify

count if `Y'==1 & zpace<. & T_`Y'>0
local N = r(N)
capture drop temp 
gen temp = age8+T_`Y'
sum temp if `Y'==1 & T_`Y'>0 & zpace<.
local M = round(r(mean), .01)
local SD = round(r(sd), .01)
sum T_`Y'
local X = round(r(max),.01)
drop temp 

putexcel A1 = `"`Y': N=`N' cases at mean age = `M' (SD=`SD') over up to `X' years of follow-up"'
putexcel B2 = matrix(`Y'M1), names 
putexcel B11 = matrix(`Y'M2), names 
putexcel B21 = matrix(`Y'M3), names 
putexcel B31 = matrix(`Y'X), names

}

//***************************************************************************//

//Figure showing main effects of all clocks 
matrix B = (1,2,3,4,5,6)'
matrix colnames B = clock
foreach x in 1 2 3 4 { 
	matrix A`x' = J(6,1,`x'), B 
	}
matrix X = (cvdM1,A2) \ (stroke_tiaM1,A3) \ (demM1, A4)
clear 
svmat2 X, names(col)
list 
rename c1 Y 
capture label drop clock
#delimit ;
label define clock
	1 "DunedinPACE"
	2 "DunedinPoAm"
	3 "Horvath Clock"
	4 "Hannum Clock"
	5 "PhenoAge Clock"
	6 "GrimAge Clock" ; #delimit cr 
label values clock clock
capture label drop survival 
#delimit ;
label define survival 
	1 "Mortality"
	2 "CVD"
	3 "Stroke/TIA"
	4 "Dementia"
	; #delimit cr 
label values Y survival 

#delimit ; 
twoway rcap lb ub clock if clock==1, lcolor(orange_red) by(Y, cols(3) legend(off) note("")) scheme(plotplain)
	|| rcap lb ub clock if clock==2, lcolor(orange)
	|| rcap lb ub clock if clock==3, lcolor(maroon)
	|| rcap lb ub clock if clock==4, lcolor(purple)
	|| rcap lb ub clock if clock==5, lcolor(navy)
	|| rcap lb ub clock if clock==6, lcolor(dknavy)
	|| scatter HR clock if clock==1, msymbol(O) mfcolor(gold) mlcolor(orange_red) msize(large)
	|| scatter HR clock if clock==2, msymbol(O) mcolor(orange) msize(large)
	|| scatter HR clock if clock==3, msymbol(O) mcolor(maroon) msize(large)
	|| scatter HR clock if clock==4, msymbol(O) mcolor(purple) msize(large)
	|| scatter HR clock if clock==5, msymbol(O) mcolor(navy) msize(large)
	|| scatter HR clock if clock==6, msymbol(O) mcolor(dknavy) msize(large)
	xscale(range(.75 6.25)) 
	xlabel(1(1)6,valuelabels angle(50) labsize(medlarge))
	xtitle("")
	ytitle(Time-to-Event Effect-size (HR), margin(medium) size(medlarge))
	ylabel(.75(.25)2,labsize(medlarge) format(%9.2f))
	yline(1)
name(survival, replace)
; #delimit cr
 
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/Figures/FraminghamSurvival.pdf" ,replace

drop if Y==1 | Y==4
#delimit ; 
twoway rcap lb ub clock if clock==1, lcolor(orange_red) by(Y, cols(3) legend(off) note("")) scheme(plotplain)
	|| rcap lb ub clock if clock==2, lcolor(orange)
	|| rcap lb ub clock if clock==3, lcolor(maroon)
	|| rcap lb ub clock if clock==4, lcolor(purple)
	|| rcap lb ub clock if clock==5, lcolor(navy)
	|| rcap lb ub clock if clock==6, lcolor(dknavy)
	|| scatter HR clock if clock==1, msymbol(O) mfcolor(gold) mlcolor(orange_red) msize(large)
	|| scatter HR clock if clock==2, msymbol(O) mcolor(orange) msize(large)
	|| scatter HR clock if clock==3, msymbol(O) mcolor(maroon) msize(large)
	|| scatter HR clock if clock==4, msymbol(O) mcolor(purple) msize(large)
	|| scatter HR clock if clock==5, msymbol(O) mcolor(navy) msize(large)
	|| scatter HR clock if clock==6, msymbol(O) mcolor(dknavy) msize(large)
	xscale(range(.75 6.25)) 
	xlabel(1(1)6,valuelabels angle(50) labsize(medlarge))
	xtitle("")
	ytitle(Time-to-Event Effect-size (HR), margin(medium) size(medlarge))
	ylabel(.75(.25)2,labsize(medlarge) format(%9.2f))
	yline(1)
name(survival_l, replace)
; #delimit cr
 
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/Figures/FraminghamSurvival_limited.pdf" ,replace

//Figure showing multivariate-adjusted effect-sizes for PACE 
matrix fx = J(1,5,999)
foreach x in dem stroke_tia cvd death  { 
	matrix fx = `x'M1[1,1...] \ fx 
	}
matrix fx = fx[1..4,1...]
matrix fx1 = (1,0,1\2,0,1\3,0,1\4,0,1)
matrix colnames fx1 = c1 model clock
matrix fx = fx , fx1 
matrix list fx 

matrix C = (1,1,2,2,3,3,4,4,5,5)'
matrix colnames C = model
matrix D = (1,2,1,3,1,4,1,5,1,6)'
matrix colnames D = clock 
foreach x in 1 2 3 4{ 
	matrix A`x' = J(10,1,`x'), C ,D 
	}
	
matrix XX = fx \ (deathX,A1) \ (cvdX,A2) \ (stroke_tiaX,A3) \ (demX,A4) 
matrix list XX
clear 
svmat2 XX, names(col)
rename c1 Y 
capture label drop model
#delimit ;
label define model
	0 "DunedinPACE"
	1 " "
	2 "adj Horvath Clock"
	3 "adj Hannum Clock"
	4 "adj PhenoAge Clock"
	5 "adj GrimAge Clock" ; #delimit cr 
label values model model
capture label drop survival 
#delimit ;
label define survival 
	1 "Mortality - Framingham"
	2 "CVD"
	3 "Stroke/TIA"
	4 "Dementia"
	5 "Mortality - NAS"
	6 "Incidence"
	7 "Prevalence"
	; #delimit cr 
label values Y survival 
list 
sort Y clock model
save temp, replace 

use temp, clear 
capture label drop Z
#delimit ;
label define Z 
	1 "Mortality"
	2 "CVD"
	3 "Stroke/TIA"
	; #delimit cr 
label values Y Z 
drop if Y == 4
#delimit ; 
twoway rcap lb ub model if model==0 & clock==1, lcolor(orange_red) by(Y, cols(3) legend(off) note("")) scheme(plotplain)
	|| rcap lb ub model if model==2 & clock==1, lcolor(maroon)
	|| rcap lb ub model if model==3 & clock==1, lcolor(purple)
	|| rcap lb ub model if model==4 & clock==1, lcolor(navy)
	|| rcap lb ub model if model==5 & clock==1, lcolor(dknavy)
	|| scatter HR model if model==0 & clock==1, msymbol(O) mfcolor(gold) mlcolor(orange_red) msize(large)
	|| scatter HR model if model==2 & clock==1, msymbol(O) mfcolor(gold) mlcolor(maroon) msize(large)
	|| scatter HR model if model==3 & clock==1, msymbol(O) mcolor(gold) mlcolor(purple) msize(large)
	|| scatter HR model if model==4 & clock==1, msymbol(O) mcolor(gold) mlcolor(navy) msize(large)
	|| scatter HR model if model==5 & clock==1, msymbol(O) mcolor(gold) mlcolor(dknavy) msize(large)
	xscale(range(.75 6.25)) 
	xlabel(0(1)5,valuelabels angle(50) labsize(medlarge))
	xtitle("")
	ytitle(Time-to-Event Effect-size (HR), margin(medium) size(medlarge))
	ylabel(.75(.25)2,labsize(medlarge) format(%9.2f))
	yline(1)
name(adjSurvival, replace)
; #delimit cr
 
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/Figures/FraminghamSurvival_adj.pdf" ,replace

use temp, clear 
drop if Y==1
#delimit ; 
twoway rcap lb ub model if model==0 & clock==1, lcolor(orange_red) by(Y, cols(3) legend(off) note("")) scheme(plotplain)
	|| rcap lb ub model if model==2 & clock==1, lcolor(maroon)
	|| rcap lb ub model if model==3 & clock==1, lcolor(purple)
	|| rcap lb ub model if model==4 & clock==1, lcolor(navy)
	|| rcap lb ub model if model==5 & clock==1, lcolor(dknavy)
	|| scatter HR model if model==0 & clock==1, msymbol(O) mfcolor(gold) mlcolor(orange_red) msize(large)
	|| scatter HR model if model==2 & clock==1, msymbol(O) mfcolor(gold) mlcolor(maroon) msize(large)
	|| scatter HR model if model==3 & clock==1, msymbol(O) mcolor(gold) mlcolor(purple) msize(large)
	|| scatter HR model if model==4 & clock==1, msymbol(O) mcolor(gold) mlcolor(navy) msize(large)
	|| scatter HR model if model==5 & clock==1, msymbol(O) mcolor(gold) mlcolor(dknavy) msize(large)

	xscale(range(.75 6.25)) 
	xlabel(0(1)5,valuelabels angle(50) labsize(medlarge))
	xtitle("")
	ytitle(Time-to-Event Effect-size (HR), margin(medium) size(medlarge))
	ylabel(.75(.25)2,labsize(medlarge) format(%9.2f))
	yline(1)
name(adjSurvival_l, replace)
; #delimit cr

 
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/Figures/FraminghamSurvival_adj_limited.pdf" ,replace


use temp, clear 
preserve 
clear 
input HR	lb	ub	Y	model	clock
1.26	1.14	1.4	5	0	1
1.257	1.129	1.398	5	1	1
1.011	0.833	1.158	5	1	3
1.271	1.141	1.416	5	3	1
0.96	0.82	1.124	5	3	4
1.23	1.102	1.373	5	4	1
1.094	0.96	1.247	5	4	5
1.066	0.922	1.231	5	5	1
1.261	1.098	1.447	5	5	6
1.15	1.01	1.309	5	2	1
1.182	1.028	1.359	5	2	2
1.16	1.12	1.2	7	0	1
1.153	1.113	1.194	7	1	1
1.042	0.997	1.09	7	1	3
1.158	1.118	1.199	7	3	1
1.008	0.963	1.055	7	3	4
1.136	1.095	1.179	7	4	1
1.069	1.022	1.119	7	4	5
1.073	1.025	1.122	7	5	1
1.128	1.079	1.18	7	5	6
1.121	1.075	1.168	7	2	1
1.072	1.023	1.124	7	2	2
1.23	1.07	1.42	6	0	1
1.234	1.069	1.424	6	1	1
0.973	0.822	1.152	6	1	3
1.203	1.04	1.39	6	3	1
1.118	0.942	1.328	6	3	4
1.221	1.053	1.415	6	4	1
1.028	0.859	1.23	6	4	5
1.123	0.943	1.336	6	5	1
1.172	0.981	1.4	6	5	6
1.176	0.998	1.387	6	2	1
1.093	0.92	1.299	6	2	2
end 
save temp1, replace 
restore 
append using temp1 
sort Y clock model

//Mortality Figure 
preserve 
keep if Y==1 | Y==5
#delimit ; 
twoway rcap lb ub model if model==0 & clock==1, lcolor(orange_red) by(Y, cols(3) legend(off) note("")) scheme(plotplain)
	|| rcap lb ub model if model==2 & clock==1, lcolor(maroon)
	|| rcap lb ub model if model==3 & clock==1, lcolor(purple)
	|| rcap lb ub model if model==4 & clock==1, lcolor(navy)
	|| rcap lb ub model if model==5 & clock==1, lcolor(dknavy)
	|| scatter HR model if model==0 & clock==1, msymbol(O) mfcolor(gold) mlcolor(orange_red) msize(large)
	|| scatter HR model if model==2 & clock==1, msymbol(O) mfcolor(gold) mlcolor(maroon) msize(large)
	|| scatter HR model if model==3 & clock==1, msymbol(O) mcolor(gold) mlcolor(purple) msize(large)
	|| scatter HR model if model==4 & clock==1, msymbol(O) mcolor(gold) mlcolor(navy) msize(large)
	|| scatter HR model if model==5 & clock==1, msymbol(O) mcolor(gold) mlcolor(dknavy) msize(large)

	xscale(range(.75 6.25)) 
	xlabel(0(1)5,valuelabels angle(50) labsize(medlarge))
	xtitle("")
	ytitle(Time-to-Event Effect-size (HR), margin(medium) size(medlarge))
	ylabel(.75(.25)2,labsize(medlarge) format(%9.2f))
	yline(1)
name(adjMortality, replace)
; #delimit cr
restore 
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/Figures/Mortality_adj.pdf" ,replace


//NAS Prevalence/ Incidence Figure 
preserve 
keep if Y==6 | Y==7
rename HR RR 
#delimit ; 
twoway rcap lb ub model if model==0 & clock==1, lcolor(orange_red) by(Y, cols(3) legend(off) note("")) scheme(plotplain)
	|| rcap lb ub model if model==2 & clock==1, lcolor(maroon)
	|| rcap lb ub model if model==3 & clock==1, lcolor(purple)
	|| rcap lb ub model if model==4 & clock==1, lcolor(navy)
	|| rcap lb ub model if model==5 & clock==1, lcolor(dknavy)
	|| scatter RR model if model==0 & clock==1, msymbol(O) mfcolor(gold) mlcolor(orange_red) msize(large)
	|| scatter RR model if model==2 & clock==1, msymbol(O) mfcolor(gold) mlcolor(maroon) msize(large)
	|| scatter RR model if model==3 & clock==1, msymbol(O) mcolor(gold) mlcolor(purple) msize(large)
	|| scatter RR model if model==4 & clock==1, msymbol(O) mcolor(gold) mlcolor(navy) msize(large)
	|| scatter RR model if model==5 & clock==1, msymbol(O) mcolor(gold) mlcolor(dknavy) msize(large)

	xscale(range(.75 6.25)) 
	xlabel(0(1)5,valuelabels angle(50) labsize(medlarge))
	xtitle("")
	ytitle(Effect-size (RR), margin(medium) size(medlarge))
	ylabel(.75(.25)1.5,labsize(medlarge) format(%9.2f))
	yline(1)
name(adjNASDisease, replace)
; #delimit cr
restore 
graph export "/Users/$db3275/OneDrive - cumc.columbia.edu/Projects/mPoA/DunedinPoAm45/eLife/R1/Figures/NASChronDx_adj.pdf" ,replace


