#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt


#Initialize empty lists 

Read_depth_list = []
Genotype_quality_list = []
Allele_frequency_list = []
Effects = []

#Step 3.0
for line in open('annotated_Variation.vcf' , 'r'):
	if line.startswith('#'):
		continue
	fields = line.rstrip('\n').split('\t')
	header = fields[7].split(';')
	effectHeader = header[-1].split('|')
	Stats_Genotype = fields[14].split(':')
	#print(effectHeader)

	#3.1 
	Read_depth = header[7].split('DP=')[1]
	Read_depth_list.append(Read_depth)

	#3.2 
	Genotype_quality_list.append(Stats_Genotype[1])

	#3.3 
	Allele_frequency = header[3].split('AF=')[1]
	Allele_frequency_list.append(Allele_frequency) 

	#3.4
	effect = effectHeader[1]

	if '&' in effect:
		doubleEffect = effect.split('&')
		Effects.append(doubleEffect[0])
		Effects.append(doubleEffect[1])
	else:
		Effects.append(effect)

UGV = Effects.count('upstream_gene_variant')
DGV = Effects.count('downstream_gene_variant')
SV = Effects.count('synonymous_variant')
MV = Effects.count('missense_variant')
DII = Effects.count('disruptive_inframe_insertion')
SpRV = Effects.count('splice_region_variant')
StRV = Effects.count('stop_retained_variant')
SL = Effects.count('stop_lost')
SG = Effects.count('stop_gained')

EffectSummary = [UGV, DGV, SV, MV, DII, SpRV, StRV + SL + SG]
EffectLabels = ['Upstream\n Gene Variant', 'Downstream\n Gene Variant', 'Synonymous\n Variant', 'Missense\n Variant', 'Disruptive\n Inframe\n Insertion', 'Splice Region\n Variant', 'Stop\n Variants']

#3.4 Predicted effects
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
fig.suptitle('Exploratory Analysis of Variants')

ax1.hist(Read_depth_list, bins = 20)
ax1.set_title('Read Depth')

ax2.hist(Genotype_quality_list, bins = 10)
ax2.set_title('Genotype Quality')

ax3.hist(Allele_frequency_list, bins = 10)
ax3.set_title('Allele Frequency')

ax1.tick_params(axis = 'x', which = 'both', bottom = False, labelbottom = False)
ax2.tick_params(axis = 'x', which = 'both', bottom = False, labelbottom = False)
ax3.tick_params(axis = 'x', which = 'both', bottom = False, labelbottom = False)

ax4.bar(EffectLabels, EffectSummary)
effectGraph = ax4.bar(EffectLabels, EffectSummary)
ax4.bar_label(effectGraph, labels = EffectSummary)
ax4.tick_params(axis='x', labelsize = 7)
ax4.set_title('Predicted Variant Effects')

fig.savefig("Graph of the Variant")




