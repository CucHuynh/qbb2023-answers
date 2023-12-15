#!/usr/bin/env python

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats import multitest
from pydeseq2 import preprocessing
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt
import csv

# Step 1.1 
counts_df = pd.read_csv("gtex_whole_blood_counts_formatted.txt", index_col = 0)

# read in metadata
metadata = pd.read_csv("gtex_metadata.txt", index_col = 0)

# Step 1.2 
counts_df_normed = preprocessing.deseq2_norm(counts_df)[0]
# log counts
counts_df_normed = np.log2(counts_df_normed + 1)

# Step 1.3 
full_design_df = pd.concat([counts_df_normed, metadata], axis = 1)

# Step 1.4 
model = smf.ols(formula = 'Q("DDX11L1") ~ SEX', data = full_design_df)
results = model.fit()

slope = results.params[1]
pval = results.pvalues[1]

# Step 1.5 
output_file = "./DE_Genes"

with open(output_file, mode = 'w', newline = '') as file:
	writer = csv.writer(file)

	writer.writerow(['Gene', 'Coefficient(Slope)', 'P-Value'])

	de_results = []
	for gene in counts_df_normed.columns:
		model = smf.ols(formula = f'Q("{gene}") ~ SEX', data = full_design_df)
		results = model.fit()

		slope = results.params[1]
		pval = results.pvalues[1]

		writer.writerow([gene, results.params['SEX'], results.pvalues['SEX']])

final_df = pd.read_csv("DE_Genes", index_col = 0)
pvals = final_df['P-Value'].to_numpy()

filtered_pvals = multitest.fdrcorrection(pvals, alpha = 0.1)
final_df['null_rejected'] = filtered_pvals[0]
FDR_df = final_df[final_df['null_rejected'] == True]
FDR_df.to_csv('DE_Genes_filtered', header = True)

# Exercise 2
# Load data into PyDESeq2 object
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors="SEX",
)

# Apply differential expression and extract results
dds.deseq2()
stat_res = DeseqStats(dds)
stat_res.summary()
results = stat_res.results_df
results.to_csv('DESeq2_results', header = True)

FDR_df = pd.read_csv("DE_Genes_filtered")
results = pd.read_csv("DESeq2_results")
results = results[~np.isnan(results['padj'])]

FDR_DESeq2 = results[results['padj'] < 0.1]
FDR_DESeq2.columns = ['Gene'] + list(results.columns[1:])
FDR_DESeq2.to_csv('DESeq2_filtered', header = True, index = True)

step2Genes = len(FDR_df)
DESeqGenes = len(FDR_DESeq2)

merge_df = pd.merge(FDR_df, FDR_DESeq2, how = 'inner', on = 'Gene')
mergeGenes = len(merge_df)

jaccard = ((mergeGenes)/(step2Genes + DESeqGenes)) * 100
print(jaccard)

# Exercise 3
colors = np.where(results['padj'] < 0.10, 'blue', 'gray')
plt.scatter(results['log2FoldChange'], -np.log10(results['padj']), color = colors, alpha = 0.5)
plt.axhline(-np.log10(0.10), color = 'black', linestyle = '--', label = 'Significance Threshold', alpha = 0.25)

plt.xlabel('Log2 Fold Change')
plt.ylabel('-log10(p-value)')
plt.title('FDR filtered (10%) DESeq2 Results')
plt.show()
plt.savefig('volcano_plot.png')



