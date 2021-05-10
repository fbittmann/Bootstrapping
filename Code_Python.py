#Felix Bittmann, 2021
#Bootstrapping - An Integrated Approach with Python and Stata
#https://doi.org/10.1515/9783110693348


#Example 1: Basic descriptive statistics in Python
print("Example 1")
from statistics import mean, stdev, median
data = [19, 29, 29, 30, 34, 36, 39, 47, 51, 52, 53, 60, 60, 64, 66, 68, 70]
print(mean(data))
print(stdev(data))
print(median(data))
print("\n" * 5)



#Example 2: Sampling with replacement
print("Example 2")
from random import choices
v = ["A", "B", "C", "D", "E", "F"]
print(choices(v, k=len(v)))
print("\n" * 5)



#Example 3: Basic bootstrap standard error
print("Example 3")
def bootstrap_se(data, func, reps):
	theta_hat_stars = []
	for i in range(reps):
		resample = choices(data, k=len(data))
		theta_hat_star = func(resample)
		theta_hat_stars.append(theta_hat_star)
	return stdev(theta_hat_stars)
print(bootstrap_se(data, mean, 50000))
print("\n" * 5)



#Example 4: Bootstrap normal CI
print("Example 4")
from random import choices
from statistics import mean, median, stdev
from morestatistics import *

def CI_normal(data, func, reps):
	theta_hat = func(data)
	theta_hat_stars = []
	for i in range(reps):
		resample = choices(data, k=len(data))
		theta_hat_star = func(resample)
		theta_hat_stars.append(theta_hat_star)
	SE_theta_hat = stdev(theta_hat_stars)
	CI_lower = theta_hat - 1.96 * SE_theta_hat
	CI_upper = theta_hat + 1.96 * SE_theta_hat
	return (CI_lower, CI_upper)
data = [19, 29, 29, 30, 34, 36, 39, 47, 51, 52, 53, 60, 60, 64, 66, 68, 70]
print(CI_normal(data, mean, 10000))
print("\n" * 5)



#Example 5: Bootstrap percentile CI
print("Example 5")
def CI_percentile(data, func, reps):
	theta_hat_stars = []
	for i in range(reps):
		resample = choices(data, k=len(data))
		theta_hat_star = func(resample)
		theta_hat_stars.append(theta_hat_star)
	CI_lower = percentile(theta_hat_stars, 2.5)
	CI_upper = percentile(theta_hat_stars, 97.5)
	return (CI_lower, CI_upper)
print(CI_percentile(data, mean, 10000))
print("\n" * 5)



#Example 6: Bootstrap BCa CI
print("Example 6")
def CI_BCa(data, func, reps, bca):
	theta_hat = func(data)
	theta_hat_stars = []
	for i in range(reps):
		resample = choices(data, k=len(data))
		theta_hat_star = func(resample)
		theta_hat_stars.append(theta_hat_star)
	n_smaller = sum([1 for theta_hat_star in theta_hat_stars if theta_hat_star <= theta_hat])
	share_smaller = n_smaller / len(theta_hat_stars)
	z0 = inverse_normal_CDF(share_smaller)	#0.975 --> 1.96
	if not bca:		#compute BC
		lower = normal_CDF(2 * z0 - 1.96)	#1.96 --> 0.975
		upper = normal_CDF(2 * z0 + 1.96)
	else:			#compute BCa
		a = acceleration_coefficient(func, data)
		lower = normal_CDF(z0 + ((z0 - 1.96) / (1 - a * (z0 - 1.96))))
		upper = normal_CDF(z0 + ((z0 + 1.96) / (1 - a * (z0 + 1.96))))
	CI_lower = percentile(theta_hat_stars, lower * 100)
	CI_upper = percentile(theta_hat_stars, upper * 100)
	return (CI_lower, CI_upper)
print(CI_BCa(data, mean, 10000, bca=False))
print(CI_BCa(data, mean, 10000, bca=True))
print("\n" * 5)



#Example 7: Bootstrap double CI
print("Example 7")
def CI_double(data, func, reps1, reps2):
	tvalues = []
	theta_hat_stars = []
	theta_hat = func(data)
	for r1 in range(reps1):
		bootsample = choices(data, k=len(data))
		theta_hat_star = func(bootsample)
		theta_hat_stars.append(theta_hat_star)
		innervalues = []
		for r2 in range(reps2):
			innersample = choices(bootsample, k=len(bootsample))
			innervalues.append(func(innersample))
		SE_theta_hat_star = stdev(innervalues)
		t = (theta_hat_star - theta_hat) / SE_theta_hat_star
		tvalues.append(t)
	SE_theta_hat = stdev(theta_hat_stars)
	lower = percentile(tvalues, 2.5)
	upper = percentile(tvalues, 97.5)
	CI_lower = theta_hat - SE_theta_hat * upper
	CI_upper = theta_hat - SE_theta_hat * lower
	return(CI_lower, CI_upper)
print(CI_double(data, mean, 10000, 100))
print("\n" * 5)



#Example 8: Permutation tests, unpaired
print("Example 8")
from random import choices, shuffle
from statistics import mean
from itertools import combinations, product

def permutationtest(func, data1, data2, reps, onesided):
	"""Implementing the unpaired permutation test"""
	assert len(data1) <= len(data2), "Enter the smaller group for data1"
	combined = data1 + data2
	tcrit = func(data1) - func(data2)
	more_extreme = 0

	if reps == 0:	#Exhaustive sampling
		ntotal = 0		#Counts the total number of combinations
		for permutation in combinations(combined, len(data1)):
			ntotal += 1
			copy = combined[:]		#Create a copy of combined
			for element in permutation:
				copy.remove(element)
			t = func(permutation) - func(copy)
			if ((onesided and tcrit <= 0 and t <= tcrit) or (onesided and tcrit > 0 and t >= tcrit) or (not onesided and abs(t) >= abs(tcrit))):
				more_extreme += 1
		return more_extreme / ntotal

	else:	#Random Sampling
		for i in range(reps):
			shuffle(combined)
			s1 = combined[:len(data1)]
			s2 = combined[len(data1):]
			t = func(s1) - func(s2)
			if ((onesided and tcrit <= 0 and t <= tcrit) or (onesided and tcrit > 0 and t >= tcrit) or (not onesided and abs(t) >= abs(tcrit))):
				more_extreme += 1
		return more_extreme / reps



#Example 9: Permutation tests, paried
print("Example 9")
def permutationtest_paired(func, data1, data2, reps, onesided):
	"""Implementing the paired permutation test"""
	assert len(data1) == len(data2), "Number of items must be equal in both lists"
	combined = data1 + data2
	diffdata = [e1 - e2 for e1, e2 in zip(data1, data2)]
	tcrit = func(diffdata)
	more_extreme = 0

	if reps == 0:						#Exhaustive sampling
		ntotal = 0		#Counts the total number of combinations
		for permutation in product([-1, 1], repeat=len(data1)):
			ntotal += 1
			t = func(value * sign for value, sign in zip(diffdata, permutation))
			if ((onesided and tcrit <= 0 and t <= tcrit) or	(onesided and tcrit > 0 and t >= tcrit) or (not onesided and abs(t) >= abs(tcrit))):
				more_extreme += 1
		return more_extreme / ntotal

	else:				#Random Sampling
		for i in range(reps):
			signs = choices([-1, 1], k=len(diffdata))
			t = func(value * sign for value, sign in zip(diffdata, signs))
			if ((onesided and tcrit <= 0 and t <= tcrit) or	(onesided and tcrit > 0 and t >= tcrit) or (not onesided and abs(t) >= abs(tcrit))):
				more_extreme += 1
		return more_extreme / reps
	
		
		
# Unpaired version, exhaustive, one sided
treatment = [94, 197, 16, 38, 99, 141, 23]
control = [52, 104, 146, 10, 51, 30, 40, 27, 46]
print(permutationtest(mean, treatment, control, reps=0, onesided=True))


# Paired version, 25,000 replications, two-sided
before= [6.37, 5.44, 5.58, 5.27, 5.11, 4.89, 4.70, 3.20]
after = [4.52, 5.69, 4.70, 3.81, 4.06, 3.22, 2.96, 3.5]
print(permutationtest_paired(mean, before, after, reps=25_000, onesided=False))
print("\n" * 5)



#Example 10: Multiprocessing in Python
print("Example 10")
import random
from statistics import stdev
from multiprocessing import Process, Queue

def thetas(data, func, reps, q):
	thetastars = [func(random.choices(data, k=len(data))) for x in range(reps)]
	q.put(thetastars)


def multiboot(data, func, reps, threads=4):
	"""Calculating Bootstrapped Standard Error using more than one thread"""
	q = Queue()
	threadlist = []
	for x in range(threads):
		threadlist.append(Process(target=thetas, args=(data, func, ((reps // threads) + 1), q)))
		threadlist[x].start()
	results = []
	for i in range(threads):
		results += q.get(True)
	seboot = stdev(results)
	for x in range(threads):
		threadlist[x].join()
	return round(seboot, 4)



#Example 11: Importing data from CSV
print("Example 11")
import csv
from statistics import mean

newdata = []
with open("testfile.csv", mode="r") as inputfile:
	reader = csv.reader(inputfile, delimiter=',')
	for row in reader:
		if row and row[0] == "ID":
			continue	#Skip first line with varnames
		elif row:
			individual = []
			for entry in row:
				individual.append(float(entry))
			newdata.append(individual)
print(newdata[:10])


#Extract wage information
wagedata = []
for entry in newdata:
	wagedata.append(entry[2])
print(len(wagedata))
print(wagedata[:10])
print(mean(wagedata))
