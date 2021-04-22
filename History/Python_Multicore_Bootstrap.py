

#Felix Bittmann, 2019
#Requires Python 3.6
#https://p16.praetorian.com/blog/multi-core-and-distributed-programming-in-python


import random
import statistics as stats
import multiprocessing as multi
import time



def thetas(data, func, reps, q):
	thetastars = [func(random.choices(data, k = len(data))) for x in range(reps)]
	q.put(thetastars)
	
	
def multiboot(data, func, reps, prec = 4, threads = 4):
	"""Calculating Bootstrapped Standard Error using more than one thread"""
	q = multi.Queue()
	threadlist = {}
	for x in range(threads):
		threadlist[x+1] = multi.Process(target=thetas, args=(data, func, (reps // threads + 1), q))
		threadlist[x+1].start()
	results = []
	for i in range(threads):
		results += q.get(True)
	seboot = stats.stdev(results)			#bottleneck, as this is not working in parallel
	for x in range(threads):
		threadlist[x+1].join()
	return round(seboot, prec)
	

def normalboot(data, func, reps, prec = 4):
	thetastars = [func(random.choices(data, k = len(data))) for x in range(reps)]
	return round(stats.stdev(thetastars), prec)
		


	
a = [x + 1 for x in range(2000)]

print("Non parallel SE calculation")
start = time.time()
print(normalboot(a, stats.mean, 30000))
print("Duration:", time.time() - start)

print("Parallel SE calculation")
start = time.time()
print(multiboot(a, stats.mean, 30000))
print("Duration:", time.time() - start)


"""Interestingly, the standard deviation calculation at the end in multiboot is the bottleneck.
	So, if you are have a large dataset and a comparatively low number of bootstrap replications,
	multiboot will be most efficient and use all four cores nicely. However, when the dataset is
	rather small and a high number of replications is desired, multiboot will be less efficient,
	however, still much faster than the single core implementation."""
