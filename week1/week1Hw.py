#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as sps
import statsmodels.formula.api as smf
import statsmodels.api as sm


#step 1.1
dataset = pd.read_csv('./aau1043_dnm.csv')

#print(dataset)

#step 1.2
count = dataset.groupby(['Proband_id','Phase_combined']).size().reset_index(name = 'Count')
#print(count)
deNovoCount = {}
for i in range(0,len(count),2):
	deNovoCount[count.loc[i, "Proband_id"]] = [count.loc[i+1, "Count"],count.loc[i, "Count"]]
#print(deNovoCount)

#step 1.3
deNovoCountDF = pd.DataFrame.from_dict(deNovoCount, orient = 'index', columns = ['maternal_dnm', 'paternal_dnm'])
#print(deNovoCountDF)

#step 1.4
Age = pd.read_csv('./aau1043_parental_age.csv', index_col = 'Proband_id')
#print(Age)

#step 1.5
concat = pd.concat([deNovoCountDF,Age],axis = 1, join = 'outer')
#print(concat)


#step 2.1
fig, ax = plt.subplots()
ax.scatter(concat["Father_age"],concat["paternal_dnm"])
ax.set_xlabel("Father Age")
ax.set_ylabel("Paternal De Novo Mutations")
ax.set_title("Paternal DNMs versus Father Age")
fig.savefig("ex2_b.png")
plt.show()
fig, bx = plt.subplots()
bx.scatter(concat["Mother_age"],concat["maternal_dnm"], c = "red")
bx.set_xlabel("Mother Age")
bx.set_ylabel("Maternal De Novo Mutations")
bx.set_title("Maternal DNMs versus Mother Age")
fig.savefig("ex2_a.png")
plt.show()
