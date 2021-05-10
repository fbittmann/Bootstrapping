
#Felix Bittmann, 2020
#Two-sided permutation test with multithreading

import math
import time
import random
import itertools
import statistics as stats
from multiprocessing import Process, Queue, Pool


def helper(res, queue):
	"""Helper Function for the multithreading processing"""
	
	reps = (res["reps"] // res["threads"]) + 1
	data = res["data1"] + res["data2"]
	alldiffs = []
	for i in range(reps):
		random.shuffle(data)
		s1 = data[:res["len1"]]
		s2 = data[res["len1"]:]
		diff = res["func"](s1) - res["func"](s2)
		alldiffs.append(diff)
	if res["empdiff"] < 0:
		overlimit = sum(1 for diff in alldiffs if diff <= res["empdiff"])
	else:
		overlimit = sum(1 for diff in alldiffs if diff >= res["empdiff"])
	results = (overlimit, reps)
	queue.put(results)
	
	
def helper_paired(res, queue):
	"""Helper Function for Multithreading for paired data"""
	
	reps = (res["reps"] // res["threads"]) + 1
	allthetas = []
	for i in range(reps):
		signs = random.choices([-1, 1], k=res["len1"])
		allthetas.append(
			res["func"](s * r for s, r in zip(signs, res["differences"]))
		)
	if res["empdiff"] > 0:
		more_extreme = sum(1 for theta in allthetas if theta >= res["empdiff"])
	else:
		more_extreme = sum(1 for theta in allthetas if theta <= res["empdiff"])
	results = (more_extreme, reps)
	queue.put(results)
	
	
def all_combos(res):
	"""This function generates all combinations for the exhaustive
	non-paired permutation"""
	
	alldata = res["data1"] + res["data2"]
	combos = itertools.combinations(alldata, res["len1"])
	for out1 in combos:
		out2 = alldata[:]	#Create copy of list with all items
		for element in out1:
			out2.remove(element)
		heap = (res, out1, out2)
		yield heap	#Packing into a single tuple since the consuming function is allowed only 1 argument
		
		
def more_extreme(heap):
	"""This function consumes the data produced by all_combos"""
	
	res, list1, list2 = heap
	diff = res["func"](list1) - res["func"](list2)
	if 0 < res["empdiff"] <= diff:
		return 1
	elif diff <= res["empdiff"] <= 0:
		return 1
	else:
		return 0
		
		
def all_combos_paired(res):
	"""This function generates all combinations for the exhaustive
	paired permutation"""
	
	signs = [-1, 1]
	for element in itertools.product(signs, repeat=res["len1"]):
		yield (res, element)
		
		
def more_extreme_paired(heap):
	"""This function consumes the data produced by all_combos_paired"""
	
	res, signs = heap
	newlist = [s * r for s, r in zip(res["differences"], signs)]
	theta = res["func"](newlist)
	if 0 < res["empdiff"] <= theta:
		return 1
	elif theta <= res["empdiff"] <= 0:
		return 1
	else:
		return 0
		
		
def display_results(res):
	"""Displays the final results of the computation"""
	
	display = {"func": "Function", "theta1": "Theta Group 1", "theta2": "Theta Group 2", "empdiff": "Group Difference",
		"len1": "N Group 1", "len2": "N Group 2", "reps": "Replications",
		"paired": "Paired Testing", "runtime": "Runtime", "threads": "Threads",
		"p_value": "P-Value",}
	p = res["prec"]
	for key, value in display.items():
		if isinstance(res[key], float):
			print(f"{value}: {res[key]:.{p}f}")
		else:
			print(f"{value}: {res[key]}")
			
			
def run_benchmark(res):
	"""Function to estimate the runtime of the desired computation"""
	
	testsize = 2000
	trash = []
	t_start = time.monotonic()
	if res["paired"]:
		if res["reps"] == 0:
			ntotal = 2 ** res["len1"]
		else:
			ntotal = res["reps"]		
		signs = [-1, 1]
		
		for i in range(testsize):
			allsigns = random.choices(signs, k=res["len1"])
			trash.append(res["func"](s * d for s, d in zip(allsigns, res["data1"])))
		trash2 = sum(1 for element in trash)
		t_end = time.monotonic()
		benchmark = t_end - t_start
		
	else:
		if res["reps"] == 0:
			ntotal = math.factorial(res["len1"] + res["len2"]) // ((math.factorial(res["len1"]) * math.factorial(res["len2"])))
		else:
			ntotal = res["reps"]
			
		data = res["data1"] + res["data2"]
		for i in range(testsize):
			random.shuffle(data)
			s1 = data[:res["len1"]]
			s2 = data[res["len1"]:]
			diff = res["func"](s1) - res["func"](s2)
			trash.append(diff)
		overlimit = sum(1 for x in trash if x < res["empdiff"])
		
		t_end = time.monotonic()
		benchmark = t_end - t_start
		
	estimated_runtime = ((ntotal / testsize) * benchmark) / res["threads"]
	estimated_runtime = estimated_runtime + (estimated_runtime * 0.15)		#Add buffer
	
	if 0 <= estimated_runtime < 5:
		print("Estimated runtime below 5 seconds")
	elif 5 <= estimated_runtime < 60:
		print(f"Estimated runtime about {estimated_runtime:.0f} seconds")
	elif 60 <= estimated_runtime < 3600:
		print(f"Estimated runtime about {estimated_runtime / 60:.1f} minutes")
	else:
		print(f"Estimated runtime about {estimated_runtime / 3600:.1f} hours")
		
		
def permutationtest(func, input1, input2, reps, prec=4, threads=4, paired=False, quiet=False, benchmark=False):
	"""Computes permutation test for paired and unpaired data
	
	Args:
		func (func): Function to compute the test with. Will be in most cases
			the arithmetic mean (statistics.mean)
		input1 (list): Containing data for group 1
		input2 (list): Containing data for group 2
		reps (int): Specifices how many random samples to take. When 0 is
			specified, an exhaustive test (testing *all* combinations) is run.
			Only feasible for small samples.
		prec (int): Specifies the display format of the results (default is 4)
		threads (int): Specifies the number of processes to create. The
			larger the faster the computation. Depends on your CPU (default is 4)
		paired (bool): Specifies whether to run a paired permutation test for
			dependent samples or a regular one. If a paired test is
			computed, the order of the elements in both lists is relevant
			(item 1 in input1 corresponds to item 1 in input2, etc...). If
			True, the number of elements in both given lists must be
			equal (default is False)
		quiet (bool): Specifies whether to display summary statistics at
			the end. If True, the results are only returned in a dict (default
			is False)
		benchmark (bool): Specifies whether to run a quick benchmark before
			the main compuation to estimate how long the computation will take.
			Gives a rough estimate of the expected runtime (default is False)
	
	Returns:
		res: a dict containing all relevant computed statistics and arguments
	"""
	
	assert threads > 0
	t_start = time.monotonic()
	res = locals()	#Collect all arguments in new dict
	#if unequal number of items, data1 should have fewer items
	alldata = sorted([input1, input2], key=len)
	res["data1"] = alldata[0]	#Shorter list
	res["data2"] = alldata[1]	#Longer list
	res["theta1"] = func(res["data1"])
	res["theta2"] = func(res["data2"])
	res["len1"] = len(res["data1"])
	res["len2"] = len(res["data2"])
	combined = res["data1"] + res["data2"]
	res["empdiff"] = res["theta1"] - res["theta2"]
	if benchmark:
		run_benchmark(res)
	
	if reps == 0 and res["len1"] + res["len2"] > 16 and not quiet:
		print("Warning, this operation (exhaustive) may take a *very* long time")
		print("Consider using random sampling instead")
	
	##### Not Paired - Exhaustive #####
	if not paired and reps == 0:
		with Pool(processes=res["threads"]) as pool:
			output = list(pool.imap_unordered(more_extreme, all_combos(res)))
			total = len(output)
			overlimit = sum(output)
			res["p_value"] = overlimit / total
			
			
	##### Not Paired - Random Sampling #####		
	if not paired and reps > 0:
		collections = Queue()
		threadlist = []
		for thread in range(res["threads"]):
			thread = Process(target=helper, args=(res, collections))
			thread.start()
			threadlist.append(thread)
		allres = []
		while len(allres) < res["threads"]:
			allres.append(collections.get(True))
			time.sleep(0.1)
		for thread in threadlist:
			thread.join()
			
		total_overlimit = sum(element[0] for element in allres)
		total_reps = sum(element[1] for element in allres)
		res["p_value"] = total_overlimit / total_reps
			
	
	if paired:
		if not res["len1"] == res["len2"]:
			raise AssertionError("Both groups must have the same number of elements in a paired test!")
			
		res["differences"] = [e1 - e2 for e1, e2 in zip(res["data1"], res["data2"])]
		res["empdiff"] = res["func"](res["differences"])
		
		##### Paired - Exhaustive #####
		if reps == 0:
			with Pool(processes=res["threads"]) as pool:
				output = list(pool.imap_unordered(more_extreme_paired, all_combos_paired(res)))
				total = len(output)
				overlimit = sum(output)
				res["p_value"] = overlimit / total
			
		##### Paired - Random Sampling #####
		if reps > 0:
			collections = Queue()
			threadlist = []
			for thread in range(res["threads"]):
				thread = Process(target=helper_paired, args=(res, collections))
				thread.start()
				threadlist.append(thread)
			
			allres = []
			while len(allres) < res["threads"]:
				allres.append(collections.get(True))
				time.sleep(0.1)
			for thread in threadlist:
				thread.join()
				
			assert len(allres) == res["threads"]
			total_more_extreme = sum(element[0] for element in allres)
			total_reps = sum(element[1] for element in allres)
			res["p_value"] = total_more_extreme / total_reps

	t_end = time.monotonic()
	res["runtime"] = round(t_end - t_start, 2)
	if not quiet:
		display_results(res)
	return res
		
		
		
		
		
		
		
		
		
		
if __name__ == '__main__':
	#Mouse Data (http://staff.ustc.edu.cn/~zwp/teach/nonpar/Permutation.pdf)
	#Result is about 0.14 (not significant, slide 13)
	treatment = [94, 197, 16, 38, 99, 141, 23]
	control = [52, 104, 146, 10, 51, 30, 40, 27, 46]


	#https://www.uvm.edu/~statdhtx/StatPages/ResamplingWithR/RandomMatchedSample/RandomMatchedSampleR.html
	#Result is about 0.016 (significant)
	before = [80.50,  84.90,  81.50,  82.60,  79.90,  88.70,  94.90,  76.30, 81.00,  80.50,
		85.00,  89.20,  81.30,  76.50,  70.00,  80.40,  83.30,  83.00,  87.70,  84.20,
		86.40,  76.50,  80.20,  87.80,  83.30,  79.70,  84.50,  80.80,  87.40]
	after = [82.20,  85.60,  81.40,  81.90,  76.40,  103.6,  98.40,  93.40,  73.40,  82.10,
		96.70,  95.30,  82.40,  72.50,  90.90,  71.30,  85.40,  81.60,  89.10,  83.90,
		82.70,  75.70,  82.60,  100.4,  85.20,  83.60,  84.60,  96.20,  86.70]


	#http://www.stat.uchicago.edu/~yibi/teaching/stat222/2017/Lectures/nonpar.pdf
	#p is 0.0156 (slide 50, one sided test)
	# ~ small1= (6.37, 5.44, 5.58, 5.27, 5.11, 4.89, 4.70, 3.20)
	# ~ small2 = (4.52, 5.69, 4.70, 3.81, 4.06, 3.22, 2.96, 3.53)

	
	# ~ print("Mouse data")
	# ~ res_mouse = permutationtest(stats.mean, treatment, control, reps=10**5, paired=False, benchmark=True)
	# ~ print("#" * 50)
	
	# ~ print("Mouse data (exhaustive)")
	# ~ res_mouse_ex = permutationtest(stats.mean, treatment, control, reps=0, paired=False, benchmark=True)
	# ~ print("#" * 50)
	
	# ~ print("Therapy data (paired)")
	# ~ res_therapy = permutationtest(stats.mean, before, after, reps=10**6, paired=True, benchmark=True)
	# ~ print("#" * 50)
	
	# ~ print("Coffee data (paired)")
	# ~ res_small = permutationtest(stats.mean, small1, small2, reps=10**5, paired=True, benchmark=True)
	# ~ print("#" * 50)
	
	# ~ print("Coffee data (paired, exhaustive)")
	# ~ res_small_ex = permutationtest(stats.mean, small1, small2, reps=0, paired=True, benchmark=True)
	# ~ print("#" * 50)
	
	
	
	
	#http://www.philender.com/courses/intro/notes/permute.html
	p1 = [55, 58, 60]
	p2 = [12, 22, 34]
	print(permutationtest(stats.mean, p1, p2, reps=0, paired=False, benchmark=False))


	

