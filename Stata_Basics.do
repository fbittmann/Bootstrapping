


*Felix Bittmann, 2019

clear all
version 15

sysuse nlsw88


********************************************************************************
*Compute descriptive statistics*
********************************************************************************

summarize wage, detail
return list

*Bootstrap the median*
bootstrap r(p50): summarize wage, detail
estat bootstrap

*More replications with BCa CIs*
bootstrap r(p50), reps(1000) bca nodots seed(123): summarize wage, detail
estat bootstrap, all

********************************************************************************
*Regression Models*
********************************************************************************

*Bootstrap Observations*
bootstrap, reps(500) nodots seed(123) saving(reg.dta, replace): regress wage c.ttl_exp i.south c.grade
estat bootstrap


*Bootstrap Residuals*
*https://www.statalist.org/forums/forum/general-stata-discussion/general/59938-trying-to-bootstrap-residuals

save "nlsw2.dta", replace
regress wage c.ttl_exp i.south c.grade
predict uhat, resid
keep uhat
save "residuals.dta", replace


cap program drop bootresiduals
program bootresiduals
version 15
use "residuals.dta", clear
bsample
merge using "nlsw2.dta"
regress wage c.ttl_exp i.south c.grade
predict xb
gen ystar = xb + uhat
regress ystar c.ttl_exp i.south c.grade
end

simulate _b, reps(500) dots(100) seed(123): bootresiduals
sum _b_ttl_exp _sim_3 _b_grade _b_cons
