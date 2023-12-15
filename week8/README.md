
#STEP 1.1 
Rscript runChicago.R raw/PCHIC_Data/GM_rep1.chinput,raw/PCHIC_Data/GM_rep2.chinput,raw/PCHIC_Data/GM_rep3.chinput ./chicago_output -d ./raw/Design/ --en-feat-list ./raw/Features/featuresGM.txt -e washU_text

#STEP 1.2

These enrichments make sense because for CTCF, H3K4me1, H3K4me3, H3K27ac and H3K9me3 groups, the expected number (second colunm) is much lower than the results (first column) suggesting they're more associated with transcription activity. There is not significant difference in H3K27me3 group suggesting this group may not be interacting with transcriptional activity.


#STEP 2.2 

- Interactions between 2 promoters: 
chr20 44438565.0 44565593.0 . 1390.0 34.77 . 0.0 chr20 44562442.0 44565593.0 PCIF1 + chr20 44438565.0 44442365.0 UBE2C +
chr20 44438565.0 44607204.0 . 1371.0 34.29 . 0.0 chr20 44596299.0 44607204.0 FTLP1;ZNF335 + chr20 44438565.0 44442365.0 UBE2C +
chr21 26837918.0 26939577.0 . 1360.0 34.02 . 0.0 chr21 26837918.0 26842640.0 snoU13 + chr21 26926437.0 26939577.0 MIR155HG +
chr20 44452862.0 44565593.0 . 1355.0 33.89 . 0.0 chr20 44562442.0 44565593.0 PCIF1 + chr20 44452862.0 44471524.0 SNX21;TNNC2 +
chr20 17660712.0 17951709.0 . 1354.0 33.85 . 0.0 chr20 17946510.0 17951709.0 MGME1;SNX5 + chr20 17660712.0 17672229.0 RRBP1 +
chr20 24972345.0 25043735.0 . 1353.0 33.84 . 0.0 chr20 24972345.0 24985047.0 APMAP + chr20 25036380.0 25043735.0 ACSS1 +

- Interactions between promotor and enhancer: 
chr21 26797667.0 26939577.0 . 1325.0 33.13 . 0.0 chr21 26926437.0 26939577.0 MIR155HG + chr21 26797667.0 26799364.0 . -
chr20 55957140.0 56074932.0 . 1291.0 32.29 . 0.0 chr20 55957140.0 55973022.0 RBM38;RP4-800J21.3 + chr20 56067414.0 56074932.0 . -
chr21 26790966.0 26939577.0 . 1166.0 29.17 . 0.0 chr21 26926437.0 26939577.0 MIR155HG + chr21 26790966.0 26793953.0 . -
chr20 5585992.0 5628028.0 . 1155.0 28.88 . 0.0 chr20 5585992.0 5601172.0 GPCPD1 + chr20 5625693.0 5628028.0 . -
chr21 26793954.0 26939577.0 . 1049.0 26.23 . 0.0 chr21 26926437.0 26939577.0 MIR155HG + chr21 26793954.0 26795680.0 . -
chr20 5515866.0 5933156.0 . 1043.0 26.08 . 0.0 chr20 5929472.0 5933156.0 MCM8;TRMT6 + chr20 5515866.0 5523933.0 . -

#STEP 2.3
- chr20:5515866-5933156
The promoter and enhancer interaction in this case is not clear. For example, where the interaction displays, the enhancer in H3K4me1 group is not present but appear nearby which suggest that the way the enhancers bind may not detected. 

- chr21:26790966-2693957
The interaction between the promoter and enhancer makes sense because the enhancer region matches in each group and the transcription level is low while the transcription level at the promoter region is higher. 

 
