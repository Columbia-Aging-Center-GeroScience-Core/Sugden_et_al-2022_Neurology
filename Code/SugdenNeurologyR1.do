global db3275 "danielbelsky"

use "/Users/$db3275/OneDrive - cumc.columbia.edu/Framingham/Data/Backbone/FhamDNAmBackboneDWB.dta", clear 
destring dbgap_subject_id, replace force 
merge 1:1 dbgap_subject_id using "/Users/$db3275/OneDrive - cumc.columbia.edu/Framingham/Data/ConstructedVariables/FhamPheno_DunedinPACE.dta", nogen keepus(familyid max_parent_edu dmob rmob smoking_status smoking_quantity horvath hannum phenoage grimage dunedin_poam38 dunedin_poam45 ar_horvath zar_horvath ar_hannum zar_hannum ar_phenoage zar_phenoage ar_grimage zar_grimage zpoam zpace cd4t cd8t nk bcell mono gran cd8pcd28ncd45ran cd8naive cd4naive plasmablast smk cpd T_cvd T_chf T_chd T_afix T_dem dem T_stroke T_stroke_tia T_death death T_cvddeath T_chddeath B1 B2 dthrvwd idtype cvddeath chddeath datedth lastatt lastcon lastsoe chf cvd chd chddate chfdate cvddate dem_status dem_survdate ad_status vad_status stroke_tia stroketia_type sr_lac_abi stroketiadate stroke stroke_type strokedate)

//Export dataset to Gloria
save "/Users/danielbelsky/OneDrive - cumc.columbia.edu/Projects/GG/FraminghamSESMortality/FhamSurvivalVars/DatafromSugden22Neuroilogy&Belsky22eLife.dta", replace 

//N=2,479 w/ DNAm data
tab cohort DNA , miss 
//N=2,468 w/ Dementia Status
tab cohort dem , miss
//N=2,264 who are dementia free and have follow-up data (28) Dx'd w/ Dem prior to visit 8 DNAm collection 
count if dem!=. & date8<dem_survdate
//Confirm 28 Dem cases dx'd before DNAm collection
count if cohort == 1 & date8>dem_survdate
//1,227 Women and 1,037 men included in analysis
tab sex  if dem!=. & date8<dem_survdate
//2,230 White and 34 non-white 
tab race  if dem!=. & date8<dem_survdate
//Age range 40-92 at DNAm baseline M=66 (SD=9)
sum age8 if dem!=. & date8<dem_survdate
//N=151 individuals at risk Dx'd over follow-up
tab dem if dem!=. & date8<dem_survdate

//Of 2,264 at risk, 2,029 have education data; 50% college grads, 25% some college, 22% HS, 3% <HS
destring edyrs, replace force 
sum edyrs if dem!=. & date8<dem_survdate
tab educ if dem!=. & date8<dem_survdate

//Clock Var macro
global clocks "zpace zpoam zar_horvath zar_hannum zar_phenoage zar_grimage" 
//Cell control macro
global cells "cd4t cd8t nk bcell mono gran cd8pcd28ncd45ran cd8naive cd4naive plasmablast"



//Survival analysis -- survival + dementia-free survival 
foreach Y in death dem { 
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
}


