#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

#Exercise1

def simulate_coverage(coverage, genome_len, read_len, figname):

	coverage_arr = np.zeros(genome_len)
	
	num_reads = int(coverage * genome_len / read_len)

	low = 0
	high = genome_len - read_len 

	start_positions = np.random.randint(low=0, high=high + 1, size=num_reads) # high is exclusive

	for start in start_positions:
		coverage_arr[start: start+read_len] += 1

	x = np.arange(0, max(coverage_arr)+1)

	sim_0cov = genome_len - np.count_nonzero(coverage_arr)
	sim_0cov_pct = 100* sim_0cov / genome_len

	print(f'In the simulation, there are {sim_0cov} bases with 0 coverage')
	print(f'This is {sim_0cov_pct}% of the genome')

	# Get poisson distribution
	y_poisson = stats.poisson.pmf(x ,mu=coverage) * genome_len

	# Get normal distribution
	y_normal = stats.norm.pdf(x, loc=coverage, scale=np.sqrt(coverage)) * genome_len

	fig, ax = plt.subplots()
	ax.hist(coverage_arr, bins=x, align='left', label='Simulation')
	ax.plot(x, y_poisson, label='Poisson')
	ax.plot(x, y_normal, label='Normal')
	ax.set_xlabel('Coverage')
	ax.set_ylabel('Frequency (bp)')
	ax.legend()
	fig.tight_layout()
	fig.savefig(figname)

#step1.3
simulate_coverage(3, 1_000_000, 100, 'ex1_3x_cov.png')

#step1.4
simulate_coverage(10, 1_000_000, 100, 'ex1_10x_cov.png')


#step1.5
simulate_coverage(30, 1_000_000, 100, 'ex1_30x_cov.png')

#Exercise2

reads = ['ATTCA', 'ATTGA', 'CATTG', 'CTTAT', 'GATTG', 'TATTT', 'TCATT', 'TCTTA', 'TGATT', 'TTATT', 'TTCAT', 'TTCTT', 'TTGAT']

k = 3 
graph = set()

for read in reads:
	for i in range(len(read) - k - 1):
		kmer1 = read[i: i+k]
		kmer2 = read[i+1: i+1+k]
		edge= kmer1 + ' -> ' + kmer2
		graph.add(edge)

for edge in graph:
	print(edge)

out_file = 'edge.txt'2
out_fs = open(out_file, 'w')
out_fs.write('digraph {\n')
for edge in graph:
	out_fs.write(edge+ '\n')
out_fs.write('}\n')
out_fs.close()











