#!/usr/bin/env python

nums = [3, 7, 9]

def findAverage(numset):
	nsum = 0
	for x in numset:
		nsum += x
	return nsum/len(numset)

print(findAverage(nums))