#!/usr/bin/env python

import numpy

#Question1
f = open("inflammation-01.csv", "r")
lines = f.readlines()
p5 = lines[4].split(",")

print(p5[0] + " " + p5[9] + " " + p5[-1])

Q2Avg = []

#Question 2 
for i in range(10):
	patient = lines[i].split(",")
	for x in range(40):
		patient[x] = int(patient[x])
	pAvg = numpy.average(patient)
	Q2Avg.append(pAvg)
	print(pAvg)

#Question 3
print(max(Q2Avg))
print(min(Q2Avg))

#Question 4
p1 = lines[0].split(",")

for i in range(40):
	p1[i] = int(p1[i]) 
	p5[i] = int(p5[i])
	pDiff = abs(p1[i] - p5[i])
	print(pDiff)

	

	

