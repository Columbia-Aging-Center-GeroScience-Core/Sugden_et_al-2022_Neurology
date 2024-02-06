


capture drop pace3
egen pace3 = cut(zpace), group(3)
recode pace3 (0=1) (1=2) (2=3) 
	
stset dem_survdate, failure(dem=1) origin(date8) scale(365)

#delimit ;
sts graph, survival by(pace3) censored(single) censopts(lcolor(gs10%50) lwidth(thin))
risktable xlabel(0(3)15) 
risktable(, order(1 "Slow DunedinPACE" 2 "Average DunedinPACE" 3 "Fast DunedinPACE") failevents title(Number At Risk (Dementia), size(large)) size(medium))

legend(ring(0) pos(7) cols(1) symxsize(5) lab(1 "Slow") lab(2 "Average") lab(3 "Fast") size(medlarge)
	title(DunedinPACE, size(medlarge)) region(lcolor(white)) )
plot1opts(lwidth(thick) lcolor(blue))
plot2opts(lwidth(thick) lcolor(purple))
plot3opts(lwidth(thick) lcolor(red))

xtitle(Analysis Time (years), size(large))
ytitle(Survival, size(large))
ylabel(,angle(horiz) nogrid labsize(large) format(%9.2f))
xlabel(,labsize(large))

graphregion(color(white)) plotregion(color(white))
xsize(6) ysize(4)
name(fhamdemsurv,replace)
title("")
; #delimit cr
