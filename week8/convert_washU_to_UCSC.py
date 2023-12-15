#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np

bait, washU, UCSC = sys.argv[1:]

interactions = []
baitchrom_coordinates = []
promo_names = []

with open(bait, 'r') as bf:
	for line in bf:
		data = line.split()
		promo_coord = data[1] + ':' + data[2]
		promo_name = data[4]
		baitchrom_coordinates.append(promo_coord)
		promo_names.append(promo_name)

promoter_dataframe = pd.DataFrame({'chrom_coordinates':baitchrom_coordinates, 'promoter':promo_names})

with open(washU, 'r') as wf:
	for line in wf:
		data = line.split()
		value = data[2]
		frag_1 = data[0].split(',')
		frag_2 = data[1].split(',')

		# bed fields
		interData = [frag_1[0], min(int(frag_1[1]), int(frag_2[1])), max(int(frag_1[2]), int(frag_2[2])), '.', 'score', value, '.', 0]
		chrom_coordinates = [frag_1[1], frag_1[2]]
		chrom_coordinates2 = [frag_2[1], frag_2[2]]

		coord_comb = frag_1[1] + ':' + frag_1[2]
		coord_comb2 = frag_2[1] + ':' + frag_2[2]

		data_from_bait = []

		if coord_comb in baitchrom_coordinates and coord_comb2 in baitchrom_coordinates:
			condition = promoter_dataframe['chrom_coordinates'] == coord_comb
			promo_name = promoter_dataframe.loc[condition]['promoter']
			promo_name = promo_name.to_string(header = False, index = False)

			condition2 = promoter_dataframe['chrom_coordinates'] == coord_comb2
			promo_name2 = promoter_dataframe.loc[condition2]['promoter']
			promo_name2 = promo_name2.to_string(header = False, index = False)

			data_from_bait = [interData[0], chrom_coordinates[0], chrom_coordinates[1], promo_name, '+', interData[0], chrom_coordinates2[0], chrom_coordinates2[1], promo_name2, '+']

		if coord_comb not in baitchrom_coordinates and coord_comb2 in baitchrom_coordinates:
			condition = promoter_dataframe['chrom_coordinates'] == coord_comb2
			promo_name = promoter_dataframe.loc[condition]['promoter']
			promo_name = promo_name.to_string(header = False, index = False)

			data_from_bait = [interData[0], chrom_coordinates[0], chrom_coordinates[1], promo_name, '+', interData[0], chrom_coordinates2[0], chrom_coordinates2[1], '.', '-']

		if coord_comb in baitchrom_coordinates and coord_comb2 not in baitchrom_coordinates:
			condition = promoter_dataframe['chrom_coordinates'] == coord_comb
			promo_name = promoter_dataframe.loc[condition]['promoter']
			promo_name = promo_name.to_string(header = False, index = False)

			data_from_bait = [interData[0], chrom_coordinates[0], chrom_coordinates[1], promo_name, '+', interData[0], chrom_coordinates2[0], chrom_coordinates2[1], '.', '-']

		fields = interData + data_from_bait
		
		interactions.append(fields)

headers = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 
		 'value', 'exp', 'color', 'sourceChrom', 'sourceStart', 
		 'sourceEnd', 'sourceName', 'sourceStrand', 'targetChrom', 
		 'targetStart', 'targetEnd', 'targetName', 'targetStrand']

dataframe = pd.DataFrame(interactions, columns = headers)

dataframe['value'] = dataframe['value'].astype(float)
max_value = max(dataframe['value'].value_counts())
dataframe['score'] = (dataframe['value'] / max_value) * 1000
dataframe['score'] = dataframe['score'].astype(int)
dataframe = dataframe.copy()

dataframe.to_csv(UCSC, sep = '\t', header = False, index = False)

header = 'track type=interact name="pCHIC" description="Chromatin interactions" useScore=on maxHeightPixels=200:100:50 visibility=full\n'

with open(UCSC, 'r') as f:
    ucsc_f = f.read()
    
header_ucsc = header + ucsc_f

with open(UCSC, 'w') as f:
    f.write(header_ucsc)

# 2.2
dataframe2 = pd.read_csv('chicago_output_UCSC.bed', names = headers, sep = '\t')

top_promoter_promoter = dataframe2[dataframe2['targetStrand'] == '+']
top_promoter_promoter = top_promoter_promoter.copy()
top_promoter_promoter.sort_values(by = 'score', ascending = False, inplace = True)

top_promoter_enhancer = dataframe2[dataframe2['targetStrand'] == '-']
top_promoter_enhancer = top_promoter_enhancer.copy()
top_promoter_enhancer.sort_values(by = 'score', ascending = False, inplace = True)

top_promoter_promoter.to_csv('max_promoter_promoter.txt', sep = ' ', header = True, index = False)
top_promoter_enhancer.to_csv('max_promoter_enhancer.txt', sep = ' ', header = True, index = False)





