#!/bin/bash

#Exercise 1: Read alignment
#Step1.1
bwa index sacCer3.fa

#Step 1.2 and 1.3
for sample in *.fastq
do
	echo 'Aligning sample:' ${sample}
	bwa mem -t 4 -R "@RG\tID:${sample}\tSM:${sample}" sacCer3.fa ${sample} > ${sample}.sam 
		samtools sort ${sample}.sam > ${sample}.bam
		samtools index ${sample}.bam
done


#Exercise 2: Variant calling and annotation
#Step 2.1 
ls *.bam > list_of_bams
freebayes -f sacCer3.fa -p 4 -= -L list_of_bams > yeast_Variation.vcf

#Step 2.2
vcffilter yeast_Variation.vcf -f "QUAL > 20" > filtered_Variation.vcf

#Step 2.3
vcfallelicprimitives -k -g filtered_Variation.vcf > decomposed_Variation.vcf

#Step 2.4
snpEff download R64-1-1.105
snpEff ann R64-1-1.105 decomposed_Variation.vcf > annotated_Variation.vcf 
head -100 annotated_Variation.vcf > Small_Sample_Variation.vcf 