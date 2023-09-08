#!/usr/bin/env python

import sys 

p1 = sys.argv[1]
p2 = sys.argv[2]
f = open(sys.argv[3])

def dInflammation(patient1, patient2, file):
	lines = file.readlines()
	patient1 = int(patient1)
	patient2 = int(patient2)
	patientData1 = lines[patient1].split(",")
	patientData2 = lines[patient2].split(",")

	for i in range(40):
		pDifference = abs(int(patientData2[i]) - int(patientData1[i]))
		print(pDifference)

dInflammation(p1, p2, f)



