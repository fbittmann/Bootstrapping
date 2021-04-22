
#Felix Bittmann, 2019
#Requires Python 3.6 because of random.choices

import random
import math
import statistics as stats
import copy


def rational_approximation(t):
	#Source: https://www.johndcook.com/blog/python_phi_inverse/
    # Abramowitz and Stegun formula 26.2.23.
    # The absolute value of the error should be less than 4.5 e-4.
    c = [2.515517, 0.802853, 0.010328]
    d = [1.432788, 0.189269, 0.001308]
    numerator = (c[2]*t + c[1])*t + c[0]
    denominator = ((d[2]*t + d[1])*t + d[0])*t + 1.0
    return t - numerator / denominator


def normal_CDF_inverse(p):
	assert 0 < p < 1

	# See article above for explanation of this section.
	if p < 0.5:
		# F^-1(p) = - G^-1(p)
		return round(-rational_approximation( math.sqrt(-2.0*math.log(p)) ),4)
	else:
		# F^-1(p) = G^-1(1-p)
		return round(rational_approximation( math.sqrt(-2.0*math.log(1.0-p)) ),4)
        
        
def normal_CDF(x):
	"""Cumulative distribution function for the standard normal distribution"""
	return (1.0 + math.erf(x / math.sqrt(2.0))) / 2.0


def percentile(data, percent, key=lambda x:x):
	#Source: https://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/
	"""Computes the percentile of a list data"""
	assert 0.0000001 <= percent <= 1
	data = sorted(data)
	if not data:
		return None
	k = (len(data)-1) * percent
	f = math.floor(k)
	c = math.ceil(k)
	if f == c:
		return key(data[int(k)])
	d0 = key(data[int(f)]) * (c-k)
	d1 = key(data[int(c)]) * (k-f)
	return d0+d1
	

def jackknife(func, data):
	"""Calculates the jackknife result for a given list"""
	final = []
	for x in range(len(data)):
		inter = copy.deepcopy(data)
		inter.pop(x)
		res = func(inter)
		final.append(res)
	return final
	
	
def afunc(func, data):
	"""Calculates the acceleration coefficient for a given list"""
	jackknifes = jackknife(func, data)
	meanjack = stats.mean(jackknifes)
	#Enumerator
	ll = len(jackknifes)
	zsum = sum([(meanjack - jackknifes[i])**3 for i in range(ll)])
	#Denominator
	nsum = sum([(meanjack - jackknifes[i])**2 for i in range(ll)])
	return zsum / (6 * (nsum**1.5))
	
		
def allcis(data, func, reps1, reps2 = 0, alpha = 0.05, prec = 3):
	thetahat = func(data)
	thetastars = []
	tvalues = []
	laenge = len(data)	
	for x in range(reps1):
		sample = random.choices(data, k=laenge)
		thetastar = func(sample)
		thetastars.append(thetastar)
		if reps2 != 0:
			thetasinner = [func(random.choices(sample, k=laenge)) for i in range(reps2)]
			tvalues.append((stats.mean(thetasinner) - thetahat) / stats.stdev(thetasinner))
	seboot = stats.stdev(thetastars)
	meanboot = stats.mean(thetastars)
	
	###Normal Based Percentile###
	tcrit = abs(normal_CDF_inverse(1 - (alpha / 2)))
	cinormal = (round(thetahat - tcrit * seboot, prec), round(thetahat + tcrit * seboot, prec))
	
	###Percentile-t###
	cipercentile = (round(percentile(thetastars, alpha / 2), prec), round(percentile(thetastars, 1 - (alpha / 2)), prec))
	
	###BC###
	bias = meanboot - thetahat
	print("Bias / SE:", bias / seboot)		#Bias/SE should be smaller than 0.25
	x = 0
	for element in sorted(thetastars):
		if element <= thetahat:
			x += 1
		else:
			break
	pboot = x / len(thetastars)
	z0 = normal_CDF_inverse(pboot)
	yl = z0 + normal_CDF_inverse(alpha / 2)
	yu = z0 + normal_CDF_inverse(1 - (alpha / 2))
	Yl = normal_CDF(z0 + yl)
	Yu = normal_CDF(z0 + yu)
	cibc = (round(percentile(thetastars, Yl), prec), round(percentile(thetastars, Yu), prec))
	
	###BCa###
	acc = afunc(func, data)
	yl = (z0 + normal_CDF_inverse(alpha / 2) ) / (1 - acc * (z0 + normal_CDF_inverse(alpha / 2)) )
	yu = (z0 + normal_CDF_inverse(1 - (alpha / 2)) ) / (1 - acc * (z0 + normal_CDF_inverse(1 - (alpha / 2))) )
	Yl = normal_CDF(z0 + yl)
	Yu = normal_CDF(z0 + yu)
	cibca = (round(percentile(thetastars, Yl), prec), round(percentile(thetastars, Yu), prec))
	
	###Percentile-t###
	if reps2 != 0:
		tvalues = sorted(tvalues)
		perclow = percentile(tvalues, alpha / 2)
		perchigh = percentile(tvalues, 1 - (alpha / 2))
		cilow = round(thetahat - (perchigh * seboot), prec)
		cihigh = round(thetahat - (perclow * seboot), prec)
		cipercentilet = (cilow, cihigh)
		return (cinormal, cipercentile, cibc, cibca, cipercentilet)
	return (cinormal, cipercentile, cibc, cibca)



data2 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]


print(allcis(data2, stats.mean, 2500, reps2 = 200))


				
