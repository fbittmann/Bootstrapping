*Felix Bittmann, 2021
*Bootstrapping - An Integrated Approach with Python and Stata
*https://doi.org/10.1515/9783110693348
*This code is writte for Stata 16.1
*To run code step-by-step, highlight relevant parts and press Ctrl + D

********************************************************************************
*Example 1: Descriptive statistics
********************************************************************************
sysuse nlsw88, clear
summarize wage, detail
return list


bootstrap r(p50), seed(123) nowarn: summarize wage, detail
estat bootstrap


bootstrap r(p50), bca reps(500) seed(123): summarize wage, detail
estat bootstrap, all



********************************************************************************
*Example 2: Estimating models
********************************************************************************
regress wage ttl_exp south grade, vce(bootstrap, reps(500) nodots ///
	seed(123) saving(bdata.dta, replace))
estat bootstrap


use "bdata.dta", clear
summarize *



********************************************************************************
*Example 3: Writing own programs
********************************************************************************
capture program drop wage_ratio
program define wage_ratio, rclass
	syntax [if] [in]
	marksample touse
	summarize wage if collgrad == 0 & `touse' == 1
	local mean_nocoll = r(mean)
	summarize wage if collgrad == 1 & `touse' == 1
	local mean_coll = r(mean)
	return scalar ratio = `mean_nocoll' / `mean_coll'
end

sysuse nlsw88, clear
bootstrap r(ratio), reps(500) nodots bca seed(123) nowarn: wage_ratio
estat bootstrap, bca



********************************************************************************
*Example 4: Resampling residuals
********************************************************************************
capture program drop bs_residuals
program define bs_residuals
	syntax, command(str) ///
	depvar(varname) ///
	indepvar(varlist) ///
	reps(integer)

	tempfile originaldata residuals
	save `originaldata', replace
	`command' `depvar' `indepvar'
	predict r, residuals
	keep r
	save `residuals', replace

	local n_indepvar : list sizeof local(indepvar)
	local n_indepvar = `n_indepvar' + 1	//+1 for constant
	matrix outcomes = J(`reps', `n_indepvar', .)	//Create empty matrix
	forvalues rep = 1/`reps' {
		use `residuals', clear
		bsample
		merge using `originaldata'
		`command' `depvar' `indepvar'
		predict yhat, xb	//predicted values
		generate ystar = yhat + r
		`command' ystar `indepvar'
		foreach NUM of numlist 1/`n_indepvar' {
			matrix outcomes[`rep', `NUM'] = e(b)[1, `NUM']
		}
	}
	clear
	svmat outcomes
	summarize *
end


sysuse nlsw88, clear
set seed 123
bs_residuals, command(regress) depvar(wage) indepvar(ttl_exp south grade) reps(500)



********************************************************************************
*Example 5: Bagging
********************************************************************************
clear all
cap program drop bagg
program define bagg, rclass
	sysuse nlsw88, clear
	drop if missing(wage, union, collgrad, ttl_exp)
	bsample
	regress wage i.union i.collgrad c.ttl_exp
	margins, at(union = 1 collgrad = 1 ttl_exp = 10)
	return scalar prediction = r(table)[1, 1]
end

simulate p=r(prediction), seed(123) reps(999) dots(50): bagg
summarize p
centile p, centile(2.5 97.5)



********************************************************************************
*Example 6: Permutation tests
********************************************************************************
sysuse nlsw88, clear
drop if missing(union, hours)
tabstat hours, statistics(mean) by(union)

quiet ttest hours, by(union)
return list

permute union diff=(r(mu_1) - r(mu_2)), reps(9999) seed(123) dots(50) nowarn: ///
	ttest hours, by(union)
	
	
	
********************************************************************************
*Example 7: Longitudinal and cluster analyses
********************************************************************************
webuse pig, clear
graph bar (mean) weight, over(week) blabel(bar) scheme(plotplain)

xtset id week			//panelvar timevar
xtreg weight week, fe

quiet xtreg weight week, fe vce(bootstrap, reps(500) seed(123) nodots)
estat bootstrap, bc



********************************************************************************
*Example 8: User-written programs with longitudinal or clustered data
********************************************************************************
webuse fifeschool, clear
quiet mixed vrq || pid:
estat icc
return list

cap program drop booticc	//delete older versions
program define booticc, rclass
	syntax [if] [in], id(varname)
	marksample touse
	mixed vrq if `touse' == 1 || `id':
	estat icc
	return scalar icc = r(icc2)
end

booticc, id(pid)	//Test program

quiet bootstrap r(icc), cluster(pid) idcluster(new_id) seed(123): booticc, id(new_id)
estat bootstrap, bc



********************************************************************************
*Example 9: Storing and graphing results
********************************************************************************
webuse productivity, clear
tempname name
postfile `name' year theta lowerbound upperbound using "corrdata", replace

forvalues YEAR = 1970(2)1986 {
	quiet bootstrap r(rho), reps(999) seed(123): ///
	correlate public unemp if year == `YEAR'
	post `name' (`YEAR') (e(b)[1, 1]) (e(ci_bc)[1, 1]) (e(ci_bc)[2, 1])
}
postclose `name'

use "corrdata", clear
	twoway (scatter theta year) ///
	(rcap lowerbound upperbound year) ///
	, xlabel(1970(2)1986) scheme(plotplain) legend(off) ///
	ytitle("Correlation")
	
	
	
********************************************************************************
*Example 10: parallel
********************************************************************************
net install parallel, from(https://raw.github.com/gvegayon/parallel/master/) replace
mata mata mlib index

parallel initialize 4	//using 4 threads
sysuse nlsw88, clear
parallel bs, reps(500): logit union ttl_exp collgrad smsa
estat bootstrap, bc
