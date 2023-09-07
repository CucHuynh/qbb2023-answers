#!/usr/bin/env python
import sys

fs = open(sys.argv[1])
lines = fs.readlines()

def findAverage(lines):
	nsum = 0
	for line in lines:
		num = int(line)
		#print(type(num))
		nsum += num

	return nsum/len(lines)

print(findAverage(lines))