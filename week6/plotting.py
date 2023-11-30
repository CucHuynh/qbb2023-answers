#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 1.2
PC1 = []
PC2 = []

for line in open('gtPCA.eigenvec'):
	line = line.rstrip('\n').split()

	PC1.append(float(line[2]))
	PC2.append(float(line[3]))

fig, ax = plt.subplots()
ax.scatter(PC1, PC2)
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("PCs for genotypes File")
plt.show()
fig.savefig("gtPCs")

MAF = []

# 2.2 
f = open('plink.frq', 'r')
for line in f:
	line = line.split()
	if 'MAF' in line:
		continue

	data = line[4]
	MAF.append(data)

MAF = np.array(MAF).astype(float)

counts, bins = np.histogram(MAF, bins = 18)
MAFhistogram = plt.stairs(counts, bins, color = 'slategrey')
plt.title('Allele Frequencies')
plt.xlabel('Allele Frequency')
plt.savefig('gtMAF')


# Function for Step 3.2
def manhattanPlot(file):
	gwaRes = []
	info = []

	with open(file, 'r') as f:
		for line in f:
			line = line.split()
			if 'CHR' in line:
				info = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8],]

			if 'CHR' not in line:
				row = [int(line[0]), line[1], int(line[2]), line[3], line[4], line[5], float(line[6]), float(line[7]), float(line[8])]
				gwaRes.append(row)


	del gwaRes[0]
	df = pd.DataFrame(gwaRes, columns = info)
	del gwaRes

	# Create column of -log10 p values
	df['-log10P'] = -np.log10(df['P'])
	df.CHR = df.CHR.astype('category')
	df = df.sort_values('CHR', ascending = True)

	#
	df['idx'] = range(len(df))
	sorted_Gwas = df.groupby('CHR')

	filtDf = df[df['-log10P'] >= 4.0]
	filtDf_Grouped = filtDf.groupby('CHR')

	# Make Manhattan plot
	fig = plt.figure() # Set the figure size
	ax = fig.add_subplot(111)
	colors = ['slategrey', 'cornflowerblue', 'blue']
	x_labels = []
	x_labels_pos = []

	for n, (name, group) in enumerate(sorted_Gwas):
		group.plot(kind = 'scatter', x = 'idx', y = '-log10P', color = colors[n % len(colors)], ax = ax)
		x_labels.append(name)
		x_labels_pos.append((group['idx'].iloc[-1] - (group['idx'].iloc[-1] - group['idx'].iloc[0])/2))

	for n, (name, group) in enumerate(filtDf_Grouped):
		group.plot(kind = 'scatter', x = 'idx', y = '-log10P', color = 'crimson', ax = ax)

	ax.set_xticks(x_labels_pos)
	ax.set_xticklabels(x_labels)
	ax.set_xlim([0, len(df)])
	ax.set_ylim([0, 10])
	ax.set_xlabel('Chromosome')

	title = file.split('_GWAS.assoc.linear')[0] + ' Manhattan Plot'
	ax.set_title(title)
	fig.savefig(title)


# Function for Step 3.3
def plotSNPEffect(gwasAssoc, vcfF):
	gwaRes = []
	info = []

	with open(gwasAssoc, 'r') as file:
		for line in file:
			line = line.split()
			if 'CHR' in line:
				info = [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[8],]

			if 'CHR' not in line:
				row = [int(line[0]), line[1], int(line[2]), line[3], line[4], line[5], float(line[6]), float(line[7]), float(line[8])]
				gwaRes.append(row)

	df = pd.DataFrame(gwaRes, columns = info)
	df['-log10P'] = -np.log10(df['P'])
	df = df.sort_values('-log10P', ascending = True)
	maxPvalue = df.iloc[len(df) - 1]
	maxVar = maxPvalue.SNP
	maxVars = df[df['SNP'].str.contains(maxVar)]
	#print(maxVars)

	GTs = []
	SNPs = []

	for line in open(vcfF, 'r'):
		if line.startswith('#'):
			continue
		if maxVar in line:
			GTs = line.rstrip('\n').split('\t')
	#print(GTs)

	for i in range(9, len(GTs)):
		if GTs[i] == './.':
			SNPs.append('NA')
		if GTs[i] == '0/0':
			SNPs.append('AA')
		if GTs[i] == '0/1':
			SNPs.append('AT')
		if GTs[i] == '1/1':
			SNPs.append('TT')

	phenotype = pd.read_csv('GS451_IC50.txt', sep = '\t')
	phenotype['GT'] = SNPs
	#print(phenotype)

	AAdata = phenotype[phenotype['GT'] == 'AA']
	ATdata = phenotype[phenotype['GT'] == 'AT']
	#TTdata = phenotype[phenotype['GT'] == 'TT']
	# no TT SNPs in genotypes, left them out

	AA = AAdata['GS451_IC50'].values
	AT = ATdata['GS451_IC50'].values
	#TT = TTdata['GS451_IC50'].values

	AA = np.delete(AA, [16])

	labels = ['AA', 'AT']
	values = [AA, AT]
	print(values)
	plt.boxplot(values, vert = True, labels = labels)
	plt.title('Effect of SNP rs17113501 on Phenotype')
	plt.xlabel('Genotype')
	plt.ylabel('Phenotype')

	plt.savefig('snpEffectBoxplot')

# 3.2
manhattanPlot('GS451_IC50_GWAS.assoc.linear')
manhattanPlot('CB1908_IC50_GWAS.assoc.linear')

# 3.3
plotSNPEffect('GS451_IC50_GWAS.assoc.linear', 'genotypes.vcf')