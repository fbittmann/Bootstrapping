

*Felix Bittmann, 2019

clear all
version 15

sysuse nlsw88

*First inspect where we find the stored results*
regress wage i.collgrad
margins collgrad
return list
matrix list r(table)


*Define Program*
capture program drop epercent
program define epercent, rclass
version 15
syntax [if]
marksample touse
quietly regress wage i.collgrad if `touse'
quietly margins collgrad
matrix output = r(table)
scalar wnocoll = output[1,1]
scalar wcoll = output[1,2]
return scalar result = (wnocoll / wcoll)
end

*Run bootstraps on program*
bootstrap r(result), reps(1000) dots(100) seed(123): epercent
estat bootstrap, all
