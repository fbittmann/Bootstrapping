

#Felix Bittmann, 2019
#Requires Python 3.6


import random
import itertools
import statistics as stats
import time

"""To understand what itertools does in this code check out the documentation:
	https://docs.python.org/3.6/library/itertools.html
	"""

def exboot_naive(data, func, prec = 5):
	"""Naive algorithm, takes much longer"""
	if len(data) > 8:
		print("This will take a *very* long time!")
	thetastars = [func(element) for element in itertools.product(data, repeat = len(data))]
	return round(stats.stdev(thetastars), prec)
	

def exboot(data, func, prec = 5):
	"""Better algorithm that does not use every permutation, but selects out combinations
	and computes the number of permutations of combinations. Doing this and using a weighted
	algorithm for the standard deviation clearly improves speed"""
	length = len(data)
	if length > 9:
		print("This will take a *very* long time!")
	results = dict((element, len(set(itertools.permutations(element, length)))) for
		element in itertools.combinations_with_replacement(data, length))
	dict2 = {}
	for key, value in results.items():
		g = func(key)
		if g not in dict2:
			dict2[g] = value
		else:
			dict2[g] += value
	#Algorithm from Wikipedia: https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Frequency_weights
	wmean = weightmean(dict2)
	nniweights = sum([value for value in dict2.values()])
	numerator = sum([value * ((key - wmean)**2) for key, value in dict2.items()])
	denominator = nniweights - 1
	return round((numerator / denominator)**0.5, prec)


def weightmean(dictdata):
	"""Takes a Dictionary and calculates the weighted mean for the key - value pairs"""
	numerator = sum([key * value for key, value in dictdata.items()])
	denominator = sum([value for value in dictdata.values()])
	return numerator / denominator





testdata = [1,2,3,4,5,6,7]


print("Naive Algorithm")
start = time.time()
print(exboot_naive(testdata, stats.mean))
print(time.time() - start)
print("Better Algorithm")
start = time.time()
print(exboot(testdata, stats.mean))
print(time.time() - start)



"""Explanation: the naive algorithm takes the data and generates every possible
	permutation (n^n) and runs the function of interest for each. Then it sums 
	up the results. However, this is deadly slow and a lot of unnecessary work as
	many results are identical (for example, [1,2,3] and [3,2,1] will give the same
	result as both contain the same data, only in a different ordert, which
	usually does not matter).
	
	The better algorithm accounts for this and starts by generating only all
	possible combinations with replacement. For example, it will generate
	[1,2,3] and [1,1,2], but not [3,2,1] or [2,1,1] since these are only
	permutations of the already created ones. Then it is done, it uses each
	unique combination as a key in a dictionary and the value is the
	number of permutations for this combination. For example, for
	[1,1,1] there is only one permutation, but for [1,1,2] are 3 ones, so the
	value for this entry is 3. When done, it runs the function only for every
	key in the dictionary, which saves a massive amount of computations. After
	that, the weighted mean and weighted standard deviation is calculated using
	the algorithm from wikipedia.
	However, even with the better algorithm, the runtime explodes when more than 8
	data points are given, so in reality, exhaustive bootstrapping is
	rarely possible."""
