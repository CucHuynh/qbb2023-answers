#!/usr/bin/env python3

import sys

from livecoding import load_bedgraph, bin_array
import numpy
import scipy.stats
import matplotlib.pyplot as plt


def main():
    # Load file names and fragment width
    forward_fname, reverse_fname, out_fname = sys.argv[1:]

    # Define what genomic region we want to analyze
    target = "chr2R"
    chromstart = 10000000
    chromend =  12000000
    chromlen = chromend - chromstart

    # Load the sample bedgraph data, reusing the function we already wrote
    forwardsamp = load_bedgraph(forward_fname, target, chromstart, chromend)
    reversesamp = load_bedgraph(reverse_fname, target, chromstart, chromend)

    # Combine tag densities, shifting by our previously found fragment width
    frag_size = 198
    combinedsamp = numpy.zeros(chromlen, float)
    combinedsamp[frag_size//2:] += forwardsamp[:-frag_size//2]
    combinedsamp[:-frag_size//2] += reversesamp[frag_size//2:]

    # Load the control bedgraph data, reusing the function we already wrote
    forwardctrl = load_bedgraph('control.fwd.bg', target, chromstart, chromend)
    reversectrl = load_bedgraph('control.rev.bg', target, chromstart, chromend)

    # Combine tag densities
    combinedctrl = numpy.zeros(chromlen, float)
    combinedctrl[frag_size//2:] += forwardctrl[:-frag_size//2]
    combinedctrl[:-frag_size//2] += reversectrl[frag_size//2:]

    # Adjust the control to have the same coverage as our sample
    ctrlcoverage = numpy.sum(combinedctrl)
    sampcoverage = numpy.sum(combinedsamp)

    ctrladjusted = combinedctrl * (sampcoverage/ctrlcoverage)

    # Create a background mean using our previous binning function and a 1K window
    # Make sure to adjust to be the mean expected per base
    binsize = 1000
    backgroundscore = bin_array(ctrladjusted, binsize)
    globalscore = backgroundscore/1000

    # Find the mean tags/bp and make each background position the higher of the
    # the binned score and global background score
    meanscore = numpy.mean(ctrladjusted)
    higherscore = numpy.maximum(globalscore, meanscore)

    # Score the sample using a binsize that is twice our fragment size
    # We can reuse the binning function we already wrote
    samplescore = bin_array(combinedsamp,frag_size*2)

    # Find the p-value for each position (you can pass a whole array of values
    # and and array of means). Use scipy.stats.poisson for the distribution.
    # Remeber that we're looking for the probability of seeing a value this large
    # or larger
    # Also, don't forget that your background is per base, while your sample is
    # per 2 * width bases. You'll need to adjust your background
    pvaluematrix = 1-(scipy.stats.poisson.cdf(samplescore,(2*frag_size*higherscore)))

    # Transform the p-values into -log10
    # You will also need to set a minimum pvalue so you doen't get a divide by
    # zero error. I suggest using 1e-250
    logpvalue = -(numpy.log10(pvaluematrix + 1e-250))

    # Write p-values to a wiggle file
    # The file should start with the line
    # "fixedStep chrom=CHROM start=CHROMSTART step=1 span=1" where CHROM and
    # CHROMSTART are filled in from your target genomic region. Then you have
    # one value per line (in this case, representing a value for each basepair).
    # Note that wiggle files start coordinates at 1, not zero, so add 1 to your
    # chromstart. Also, the file should end in the suffix ".wig"
    write_wiggle(logpvalue,target,chromstart,f"{out_fname}.wig")
    
    # Write bed file with non-overlapping peaks defined by high-scoring regions 
    write_bed(logpvalue,target,chromstart,chromend,frag_size,f"{out_fname}.bed")

def write_wiggle(pvalues, chrom, chromstart, fname):
    output = open(fname, 'w')
    print(f"fixedStep chrom={chrom} start={chromstart + 1} step=1 span=1",
          file=output)
    for i in pvalues:
        print(i, file=output)
    output.close()

def write_bed(scores, chrom, chromstart, chromend, width, fname):
    chromlen = chromend - chromstart
    output = open(fname, 'w')
    while numpy.amax(scores) >= 10:
        pos = numpy.argmax(scores)
        start = pos
        while start > 0 and scores[start - 1] >= 10:
            start -= 1
        end = pos
        while end < chromlen - 1 and scores[end + 1] >= 10:
            end += 1
        end = min(chromlen, end + width - 1)
        print(f"{chrom}\t{start + chromstart}\t{end + chromstart}", file=output)
        scores[start:end] = 0
    output.close()


if __name__ == "__main__":
    main()




    