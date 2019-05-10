
#Felix Bittmann, 2019

import csv
import statistics as stats


newdata = []
with open("testfile.csv", mode="r") as inputfile:
	reader = csv.reader(inputfile, delimiter=',')
	for row in reader:
		if row[0] == "ID":
			continue
		else:
			individual = []
			for entry in row:
				"""Here we use a try / except construction to account for the fact that some
				data points might not be numeric. If this construct is not used and an error
				occurs, the script will stop"""
				try:
					individual.append(float(entry))
				except ValueError:
					print("Data point cannot be converted to int / float!")
					individual.append(entry)		#Add without converting
			newdata.append(individual)
print(newdata[:10])


"""Extracting information about wages"""
wagedata = []
for entry in newdata:
	wagedata.append(entry[2])
print(len(wagedata))
print(wagedata[:10])
print(stats.mean(wagedata))
