
#Felix Bittmann, 2020
#Requires Python 3.6 because of random.choices

import time
import random
from statistics import mean, stdev
from multiprocessing import Queue, Process
from morestatistics import *


def run_benchmark(data, res):
	"""Computing a runtime estimate"""
	t_start = time.monotonic()
	ntotal = res["reps1"] * max(1, res["reps2"])
	testsize = 1000
	trash = [res["func"](random.choices(data, k=res["n"])) for i in range(testsize)]
	t_end = time.monotonic()
	time1 = t_end - t_start
	t_start = time.monotonic()
	trash = acceleration_coefficient(res["func"], data)
	t_end = time.monotonic()
	time2 = t_end - t_start
	estimated_runtime = ((ntotal / testsize) * time1) / res["threads"] + time2
	estimated_runtime = estimated_runtime + (estimated_runtime * 0.20)		#Add buffer
	
	if 0 <= estimated_runtime < 5:
		print("Estimated runtime below 5 seconds")
	elif 5 <= estimated_runtime < 60:
		print(f"Estimated runtime about {estimated_runtime:.0f} seconds")
	elif 60 <= estimated_runtime < 3600:
		print(f"Estimated runtime about {estimated_runtime / 60:.1f} minutes")
	else:
		print(f"Estimated runtime about {estimated_runtime / 3600:.1f} hours")
	

def multifunc(data, res, queue):
	"""Multihreading working function to generate bootstrap resamples"""
	theta_stars = []
	tvalues = []		#only needed for double bootstrap
	failed = 0			#bookkeeping
	failed_inner = 0	#bookkeeping
	reps1 = (res["reps1"] // res["threads"]) + 1

	for t1 in range(reps1):
		bootsample1 = random.choices(data, k=res["n"])
		try:
			theta_star = res["func"](bootsample1)
		except:
			failed += 1
		theta_stars.append(theta_star)
		
		if res["reps2"] > 0:
			innervalues = []
			for t2 in range(res["reps2"]):
				try:
					innervalues.append(res["func"](random.choices(bootsample1, k=len(bootsample1))))
				except:
					failed_inner += 1
			tvalues.append((theta_star - res["theta_hat"]) / stdev(innervalues))
	queue.put([theta_stars, tvalues, failed, failed_inner])
		
		
def bootstrap_ci(func, data, reps1, reps2=0, alpha=0.05, prec=3, threads=2, quiet=False, benchmark=False, seed=None):
	"""Computes Bootstrap Confidence Intervals for given data and function
	
	Args:
		func (func): Function to use to generate bootstrap results. This can be
			any Python function that works with the data given, for example the
			arithmetic mean, median, standard deviation, interquartile range, etc...
		data (list): Given data to compute confidence intervals for
		reps1 (int): Number of bootstrap replication samples. The larger the more
			precise the results. Must be greater than zero
		reps2 (int): If the double (iterated) bootstrap should be computed,
			a number must be given. Otherwise this type of interval will
			not be computed. Note that computation time rises in multiplicative
			fasion as the total number of samples to take is reps1 * reps2
			(default is 0)
		alpha (float): Nominal coverage of the CI is 1 - alpha. Commonly used
			are 0.05, 0.01 and 0.001. For a 95% CI, enter 0.05 (default is 0.05)
		prec (int): Number of decimal places to display in the results
			(default is 3)
		threads (int): Number of processes to run to speed up computation.
			Depends on your CPU (default is 2)
		quiet (bool): Specifies whether to display results in a nice fasion.
			If True, results are only returned as a dict (default is False)
		benchmark (book): Specifies whether to run a quick benchmark before
			the computation to give a coarse estimate of the expected runtime
			(default is False)
		seed (str): Input a seed for repeatable random draws (default is None)
		
	Returns:
		res: a dict containing all relevant computed statistics and arguments
	"""
		
	t_start = time.monotonic()	
	statistics = ("func","reps1", "reps2", "alpha", "prec", "theta_hat", "tvalues", "n", "se_boot", "mean_boot",
		"bias", "normal", "percentile", "bc", "bca", "double")
	res = {name: None for name in statistics}
	res["theta_hat"] = func(data)
	res["n"] = len(data)
	res["reps1"] = reps1
	res["reps2"] = reps2
	res["alpha"] = alpha
	res["prec"] = prec
	res["func"] = func
	res["failed"] = 0
	res["failed_inner"] = 0
	res["threads"] = threads

	
	if benchmark:
		run_benchmark(data, res)
	if seed:
		random.seed(seed)
	collection = Queue()	#collect all bootstrap results
	threadlist = []
	for i in range(threads):
		process = Process(target=multifunc, args=(data, res, collection))
		process.start()
		threadlist.append(process)
	
	tempdata = []
	while len(tempdata) < res["threads"]:
		tempdata.append(collection.get(True))
		time.sleep(0.5)
	
	theta_stars = []
	tvalues = []
	for element in tempdata:
		theta_stars += element[0]
		tvalues += element[1]
		res["failed"] += element[2]
		res["failed_inner"] += element[3]
		
	for thread in threadlist:				#Terminate all processes
		thread.join()
		
			
	res["mean_boot"] = mean(theta_stars)
	res["se_boot"] = stdev(theta_stars, res["mean_boot"])
	res["bias"] = res["mean_boot"] - res["theta_hat"]
	theta_stars.sort()		#Sort only once for all following computations

	
	### Normal Based ###
	tcrit = abs(inverse_normal_CDF(1 - (res["alpha"] / 2)))
	lower = res["theta_hat"] - tcrit * res["se_boot"]
	upper = res["theta_hat"] + tcrit * res["se_boot"]
	res["normal"] = (lower, upper)
	
	
	### Percentile ###
	lower = percentile(theta_stars, (res["alpha"] / 2) * 100, is_sorted=True)
	upper = percentile(theta_stars, (1 - (res["alpha"] / 2)) * 100, is_sorted=True)
	res["percentile"] = (lower, upper)
	
	
	### BC ###
	n_smaller = sum([1 for theta in theta_stars if theta < res["theta_hat"] ])
	share_smaller = n_smaller / len(theta_stars)
	z = inverse_normal_CDF(share_smaller)
	perclower = normal_CDF(2 * z - tcrit)
	percupper = normal_CDF(2 * z + tcrit)
	lower = percentile(theta_stars, perclower * 100, is_sorted=True)
	upper = percentile(theta_stars, percupper * 100, is_sorted=True)
	res["bc"] = (lower, upper)
	
	
	### BCa ###
	try:
		a = acceleration_coefficient(func, data)
		perclower = normal_CDF(z + ((z - tcrit) / (1 - a * (z - tcrit))))
		percupper = normal_CDF(z + ((z + tcrit) / (1 - a * (z + tcrit))))
		lower = percentile(theta_stars, perclower * 100, is_sorted=True)
		upper = percentile(theta_stars, percupper * 100, is_sorted=True)
		res["bca"] = (lower, upper)
	except:
		print("Computation of acceleration coefficient failed")
		res["bca"] = None
	
	
	### Double ###
	if reps2 > 0:
		tvalues.sort()
		perclower = percentile(tvalues, (res["alpha"] / 2) * 100, True)
		percupper = percentile(tvalues, (1 - (res["alpha"] / 2)) * 100, True)
		lower = res["theta_hat"] - res["se_boot"] * percupper
		upper = res["theta_hat"] - res["se_boot"] * perclower
		res["double"] = (lower, upper)
	
	res["runtime"] = time.monotonic() - t_start
	if not quiet:
		###Display results###
		for key, value in res.items():
			if isinstance(value, float):
				print(key, round(value, prec))
			elif isinstance(value, tuple):
				print(key, tuple((round(number, prec) for number in value)))
			else:
				print(key, value)
	return res
			



if __name__ == '__main__':
	# ~ testdata = [i + 1 for i in range(0, 200)]
	# ~ a = bootstrap_ci(mean, testdata, reps1=500, reps2=50, threads=4, benchmark=False)
	# ~ print(a)
	
	
	studentdata = [19, 29, 29, 30, 34, 36, 39, 47, 51, 52, 53, 60, 60, 64, 66, 68, 70]
	bootstrap_ci(kurtosis, studentdata, reps1=20_000, reps2=500, threads=5, benchmark=False)
	

	
