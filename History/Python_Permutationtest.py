
#Felix Bittmann, 2019
#Two-sided permutation test

import random
import statistics as stats
import itertools
import copy

def permutationtest(func, data1, data2, reps = 0, prec = 4):
	l1 = len(data1)
	l2 = len(data2)
	#If unequal number of elements per group, the first one should contain the lower number of items
	if l2 < l1:
		data1, data2 = data2, data1
		l1 = len(data1)
		l2 = len(data2)
	combined = data1 + data2
	empdiff = abs(func(data1) - func(data2))
		
	if reps == 0:						#All combinations
		if len(combined) > 16:
			print("This might take a while!")
		limit = 0
		length = 0		#Counts the total number of combinations
		for element in itertools.combinations(combined, l1):
			length += 1
			copylist = copy.deepcopy(combined)
			for subelement in element:
				copylist.remove(subelement)
			diff = abs(func(element) - func(copylist))
			if diff >= empdiff:
				limit += 1
		pvalue = limit / length
		return round(pvalue, prec)
		
	else:								#Random Sampling
		limit = 0
		for x in range(reps):
			random.shuffle(combined)
			s1 = combined[:l1]
			s2 = combined[l1:]
			diff = abs(func(s1) - func(s2))
			if diff >= empdiff:
				limit += 1
		pvalue = limit / reps
		return round(pvalue, prec)


t1 = [22,23,23,26,29]
t2 = [26,27,30,34,35]
print(permutationtest(stats.mean, t1, t2))
	

