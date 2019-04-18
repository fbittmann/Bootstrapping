
#Felix Bittmann, 2019
#Requires Python 3.6 or newer because of random.choices

import random
import statistics as stats

def pearsonsr(data, prec = 4):
	"""Unpacks a list of touples into two lists and calculates Pearson's R
	correlation coefficient"""
	a, b = [], []
	for element in data:
		a.append(element[0])
		b.append(element[1])
	assert len(a) == len(b)
	meana = stats.mean(a)
	meanb = stats.mean(b)
	numerator = sum([(a[x] - meana) * (b[x] - meanb) for x in range(len(a))])
	d1 = sum([(a[x] - meana)**2 for x in range(len(a))])
	d2 = sum([(b[x] - meanb)**2 for x in range(len(b))])
	return round(numerator / (d1 * d2)**0.5, prec)
	

def powerboot(func, data, reps, size, tcrit = 1.96, prec = 3):
	print("Size of each bootstrap resample:", size)
	theta = func(data)
	thetas = [func(random.choices(data, k = size)) for x in range(reps)]
	bootse = stats.stdev(thetas)
	ci_lower = round(theta - tcrit * bootse, prec)
	ci_upper = round(theta + tcrit * bootse, prec)
	return (ci_lower, ci_upper)


data1 = [(130, 24.626),(150, 22.869),(116, 26.38),(142, 23.216),(128, 21.526),
(120, 26.963),(126, 21.99),(134, 28.039),(138, 26.46),(122, 28.711),(175, 26.198),
(124, 25.336),(110, 22.512),(120, 18.673),(114, 28.403),(90, 23.363),(160, 26.228),
(118, 21.2373),(114, 26.009),(124, 30.247)]

print(pearsonsr(data1))
for x in (1,2,3,4,5,6,7):
	print(powerboot(pearsonsr, data1, 5000, 20 * x))
