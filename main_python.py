

#Felix Bittmann, 2019
import random


def mean(data):
	"""Calculates the arithmetic mean of a given list"""
	return sum(data) / len(data)
	

def std(data):
	"""Calculates the standard deviation of a given list"""
	meandata = mean(data)
	sumdiff = 0
	for element in data:
		sumdiff = sumdiff + (element - meandata)**2
	return (sumdiff / (len(data) - 1))**0.5


def median(data):
	"""Returns the median of a given list"""
	data.sort()
	length = len(data)
	assert length >= 3, "Not enough items"
	if length % 2 == 1:		#odd number of items
		return data[length // 2]
	else:
		return (data[(length//2)-1] + data[(length//2)]) / 2
		
	
def bootsample(data):
	"""Creates one bootstrap sample from a given list"""
	bsample = []
	for x in range(0, len(data)):
		randomnum = random.randint(0, len(data) - 1)
		bsample.append(data[randomnum])
	return bsample
	
	
def bootse(data, func, reps):
	"""Calculates a bootstrapped standard error for a list and a given function"""
	assert reps > 10, "Not enough replications"
	thetastars = []
	for x in range(reps):
		sample = bootsample(data)
		thetastar = func(sample)
		thetastars.append(thetastar)
	return std(thetastars)
	

def normalci(data, func, reps, tcrit = 1.96, prec = 3):
	"""Calculates a normal-based confidence interval around a point estimate"""
	se = bootse(data, func, reps)
	point_estimate = func(data)
	ci_lower = round(point_estimate - tcrit * se, prec)
	ci_upper = round(point_estimate + tcrit * se, prec)
	return (ci_lower, ci_upper)
	
	
data1 = [1,2,3,4,5,6,7,8,9,10]

print(mean(data1))
print(std(data1))
print(bootse(data1, median, 5000))
print(normalci(data1, mean, 5000))
