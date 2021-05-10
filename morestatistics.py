#Felix Bittmann, 2021
#Bootstrapping - An Integrated Approach with Python and Stata
#https://doi.org/10.1515/9783110693348



import math
from statistics import mean

def rational_approximation(t):
	#Source: https://www.johndcook.com/blog/python_phi_inverse/
	# Abramowitz and Stegun formula 26.2.23.
	# The absolute value of the error should be less than 4.5 e-4.
	c = [2.515517, 0.802853, 0.010328]
	d = [1.432788, 0.189269, 0.001308]
	numerator = (c[2] * t + c[1]) * t + c[0]
	denominator = ((d[2] * t + d[1]) * t + d[0]) * t + 1
	return t - (numerator / denominator)


def inverse_normal_CDF(p):
	assert 0 < p < 1
	if p < 0.5:
		# F^-1(p) = - G^-1(p)
		return round(-rational_approximation(math.sqrt(-2 * math.log(p))), 4)
	else:
		# F^-1(p) = G^-1(1-p)
		return round(rational_approximation(math.sqrt(-2 * math.log(1 - p))), 4)


def normal_CDF(x):
	"""Cumulative distribution function for the standard normal distribution"""
	return (1 + math.erf(x / math.sqrt(2))) / 2


def jackknife(func, data):
	"""Calculates the jackknife coefficients for a given list"""
	output = []
	for element in data:
		copylist = data[:]		#Copy original list
		copylist.remove(element)
		coef = func(copylist)
		output.append(coef)
	return output


def acceleration_coefficient(func, data):
	"""Calculates the acceleration coefficient for a given list"""
	jackvalues = jackknife(func, data)
	mean_jackvalues = mean(jackvalues)
	nominator, denominator = 0, 0
	for element in jackvalues:
		nominator += (mean_jackvalues - element) ** 3
		denominator += (mean_jackvalues - element) ** 2
	return nominator / (6 * (denominator ** 1.5))


def percentile(data, percent, is_sorted=False):
	#Source: https://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/
	"""Computes the percentile of a given list"""
	assert 0 <= percent <= 100, "Percent must be in [0, 100]"
	if not is_sorted:
		data.sort()
	position = (len(data) - 1) * (percent / 100)
	if position.is_integer():
		return data[int(position)]
	else:
		lower = int(position)	#equivalent to math.floor(position)
		d0 = data[lower] * ((lower + 1) - position)
		d1 = data[lower + 1] * (position - lower)
		return d0 + d1
