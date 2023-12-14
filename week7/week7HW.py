#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np


#Part 3a & b

ONT = []
BI = []
for line in open("ONT.cpg.chr2.bedgraph"):
	col = line.rstrip().split()
	ONT.append(float(col[1]))

for line in open("bisulfite.cpg.chr2.bedgraph"):
	col = line.rstrip().split()
	BI.append(float(col[1]))


# Find reads that appear more than once in datasets
ONT_set = set()
for i in range(len(ONT)):
    if ONT[i] not in ONT_set:
        ONT_set.add(ONT[i])

BI_set = set()
for i in range(len(BI)):
    if BI[i] not in BI_set:
        BI_set.add(BI[i])

diff_ONT = ONT_set.difference(BI_set)
diff_BI = BI_set.difference(ONT_set)

BOTH = ONT_set.intersection(BI_set)


ONT_coverage = []
BI_coverage = []
for line in open("ONT.cpg.chr2.bedgraph"):
	col = line.rstrip().split()
	ONT_coverage.append(float(col[4]))

for line in open("bisulfite.cpg.chr2.bedgraph"):
	col = line.rstrip().split()
	BI_coverage.append(float(col[4]))

fig, ax1= plt.subplots()
ax1.hist(ONT_coverage, bins = 2000, color = 'grey', alpha = 0.5, label = "ONT coverage")
ax1.set_xlabel("Read Coverage")
ax1.set_ylabel("Frequency (site)")
ax1.hist(BI_coverage, bins = 2000, color = 'red', alpha = 0.5, label = "BI coverage")
ax1.set_title("ONT and Bisulfite Reads Coverage Distribution")
ax1.set_xlim(0,100)
ax1.legend()

fig.savefig("Read_coverage_for_BI_ONT.png")
fig.tight_layout()
#plt.show()

#Part 3c

ONT_methylation = []
BI_methylation = []
for line in open("ONT.cpg.chr2.bedgraph"):
	col = line.rstrip().split()
	ONT_methylation.append(float(col[3]))

for line in open("bisulfite.cpg.chr2.bedgraph"):
	col = line.rstrip().split()
	BI_methylation.append(float(col[3]))

ONT_score = []
BI_score = []

for i in range(len(ONT)):
	if ONT[i] in BOTH:
		ONT_score.append(ONT_methylation[i])

for j in range(len(BI)):
	if BI[j] in BOTH:
		BI_score.append(BI_methylation[j])

correlation = np.corrcoef(ONT_score, BI_score)

fig = plt.figure()
plt.imshow(np.log10(1+ np.histogram2d(ONT_score, BI_score, bins =100)[0]))
plt.xlabel("ONT Methylation Score")
plt.ylabel("BI Methylation Score")
plt.title(f'ONT and BI Methylation Scores (R= {correlation[0,1]:0.3f})')
plt.show()
fig.savefig("Correlation_Methylation_BI_ONT.png")


#Part 3d

ONT_normal = []
ONT_tumor = []
ONT_coverage = []
ONT_tumor_coverage = []
ONT_score = []
ONT_tumor_score = []

ontsite = []
ont1site = []
ontcoverage = []
ont1coverage = []
ontmeth = []
ont1meth = []

for line in open("normal.ONT.chr2.bedgraph"):
	col = line.rstrip().split()
	ONT_normal.append(float(col[1]))
	ONT_score.append(float(col[3]))
	ONT_coverage.append(int(col[4]))

for line in open("tumor.ONT.chr2.bedgraph"):
	col = line.rstrip().split()
	ONT_tumor.append(float(col[1]))
	ONT_tumor_score.append(float(col[3]))
	ONT_tumor_coverage.append(int(col[4]))

ONT_set = set()
for i in range(len(ONT_normal)):
    if ONT_normal[i] not in ONT_set:
        ONT_set.add(ONT_normal[i])


ONT_tumor_set = set()
for i in range(len(ONT_tumor)):
    if ONT_tumor[i] not in ONT_tumor_set:
        ONT_tumor_set.add(ONT_tumor[i])

bothsite = ONT_set.intersection(ONT_tumor_set)


ONT_methylation_normal = []
ONT_methylation_tumor = []
for i in range(len(ONT_normal)):
	if ONT_normal[i] in bothsite:
		ONT_methylation_normal.append(ONT_score[i])

for j in range(len(ONT_tumor)):
	if ONT_tumor[j] in bothsite:
		ONT_methylation_tumor.append(ONT_tumor_score[j])

change_ONT = []
for i in range(len(ONT_methylation_normal)):
	change = ONT_methylation_normal[i] - ONT_methylation_tumor[i]
	if change != 0:
		change_ONT.append(change)


BI_normal = []
BI_tumor = []
BI_coverage = []
BI_tumor_coverage = []
BI_score = []
BI_tumor_score = []


for line in open("normal.bisulfite.chr2.bedgraph"):
	col = line.rstrip().split()
	BI_normal.append(float(col[1]))
	BI_score.append(float(col[3]))
	BI_coverage.append(int(col[4]))

for line in open("tumor.bisulfite.chr2.bedgraph"):
	col = line.rstrip().split()
	BI_tumor.append(float(col[1]))
	BI_tumor_score.append(float(col[3]))
	BI_tumor_coverage.append(int(col[4]))

BI_set = set()
for i in range(len(BI_normal)):
    if BI_normal[i] not in BI_set:
        BI_set.add(BI_normal[i])


BI_tumor_set = set()
for i in range(len(BI_tumor)):
    if BI_tumor[i] not in BI_tumor_set:
        BI_tumor_set.add(BI_tumor[i])

bothsite = BI_set.intersection(BI_tumor_set)

BI_methylation_normal = []
BI_methylation_tumor = []
for i in range(len(BI_normal)):
	if BI_normal[i] in bothsite:
		BI_methylation_normal.append(BI_score[i])

for j in range(len(BI_tumor)):
	if BI_tumor[j] in bothsite:
		BI_methylation_tumor.append(BI_tumor_score[j])

change_BI = []
for i in range(len(BI_methylation_normal)):
	change_BI_stuff = BI_methylation_normal[i] - BI_methylation_tumor[i]
	if change_BI_stuff != 0:
		change_BI.append(change_BI_stuff)

fig, ax1= plt.subplots()
ax1.violinplot([change_ONT, change_BI])
ax1.set_xticks([1,2])
ax1.set_xticklabels(['ONT', 'Bisulfite'])
ax1.set_xlabel("Type of Sequencing")
ax1.set_ylabel("Methylation Change")
ax1.set_title('Distribution of Methylation Changes')
fig.savefig("violin_plot.png") 
fig.tight_layout()
plt.show()
plt.close