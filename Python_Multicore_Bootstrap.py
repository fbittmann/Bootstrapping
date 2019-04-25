

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
	
	
def multiboot(data, func, reps, prec = 4):
	"""Calculating Bootstrapped Standard Error using four cores"""
	q = multi.Queue()
	p1 = multi.Process(target=thetas, args=(data, func, (reps//4 + 1), q))
	p1.start()
	p2 = multi.Process(target=thetas, args=(data, func, (reps//4 + 1), q))
	p2.start()
	p3 = multi.Process(target=thetas, args=(data, func, (reps//4 + 1), q))
	p3.start()
	p4 = multi.Process(target=thetas, args=(data, func, (reps//4 + 1), q))
	p4.start()
	results = []
	for i in range(4):
		results += q.get(True)
	seboot = stats.stdev(results)			#bottleneck, as this is not working in parallel
	p1.join()
	p2.join()
	p3.join()
	p4.join()
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
