#!/usr/bin/env python

#get a starting frequency and a population size
#input parameters for function

#make a list to store allel frequencies 

#while our allele freq is between 0 and 1 (is not 0 or 1)
	#get the new allele frequency for next generation 
	#by drawing from the binomial distribution

	#convert number of sucesses into a frequency
	#store our allele frequency in the AF list

#return a list of allel frequency at each time point 
#number of generation to fixation is the length of the list

import numpy as np
import matplotlib.pyplot as plt

#Exercise 1 

def wrightfisher(allelefreq,pop):

	allelefreq_list = []

	while allelefreq > 0 and allelefreq < 1:
		success = np.random.binomial(2*pop, allelefreq)
		allelefreq = success/(2*pop)
		allelefreq_list.append(allelefreq)

	Generation = len(allelefreq_list)
	return [allelefreq_list, Generation]

results = wrightfisher(0.8, 6000)

fig, ax = plt.subplots()
'''
x_positions = range(0,results[1])
y_positions = results[0]
ax.plot(x_positions, y_positions)
ax.set_xlabel("Generation")
ax.set_ylabel("Allele Frequency")
ax.set_ylim(0,1)
ax.set_xlim(0)
ax.set_title("Wright Fisher Model")
fig.savefig( "WrightFisher1.png" )

plt.show()

#Exercise 2 
timefixation = []
for i in range(1000):
	results = wrightfisher(0.3, 500)
	timefixation.append(results[1])

ax.hist(timefixation)

ax.set_xlabel("Times to Fixation (generations)")
ax.set_ylabel("Number of Occurance")
ax.set_title("Histogram of Times to Fixation")
fig.savefig( "HistWrightFisher2.png" )
plt.show()
'''
#EXERCISE 2 RESUBMISSION 

timefixation = []
for i in range(1000):
	results = wrightfisher(0.3, 500)
	x_positions = range(0,results[1])
	y_positions = results[0]
	ax.plot(x_positions, y_positions)
	timefixation.append(results[1])

ax.hist(timefixation)

ax.set_xlabel("Times to Fixation (generations)")
ax.set_ylabel("Number of Occurance")
ax.set_title("Histogram of Times to Fixation")
fig.savefig( "HistWrightFisher2.png" )
plt.show()


#Exercise 3
'''
avgtimelist = []
populationsize = []
for x in range(5):
	randomnum = np.random.randint(50,200)
	populationsize.append(randomnum)
	timefixation = []
	for i in range(50):
		results = wrightfisher(0.5, randomnum)
		timefixation.append(results[1])

	avgtime = sum(timefixation)/len(timefixation)
	avgtimelist.append(avgtime)

ax.scatter(populationsize,avgtimelist)
ax.set_xlabel("Population Size")
ax.set_ylabel("Average Time to Fixation")
ax.set_title("The Effect of Population on Time to Fixation")
fig.savefig( "ChangingPopulation.png" )
'''
'''
avgtimelist = []
frequency = []
for i in range(1000):
	randomnum = np.random.uniform(0,1)
	frequency.append(randomnum)
	timefixation = []
	for i in range(50):
		results = wrightfisher(randomnum, 100)
		timefixation.append(results[1])

	avgtime = sum(timefixation)/len(timefixation)
	avgtimelist.append(avgtime)

ax.scatter(frequency,avgtimelist)
ax.set_xlabel("Allele Frequency")
ax.set_ylabel("Average Time to Fixation")
ax.set_title("The Effect of Allele Frequency on Time to Fixation")
fig.savefig( "ChangingFrequency.png" )

plt.show()
'''