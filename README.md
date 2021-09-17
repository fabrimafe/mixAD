# Reference-free detection of genomic mixtures from low-coverage data
Used in Neandertal and Denisovan DNA from Pleistocene sediments, April 2017 Science 356(6338):eaam969
dx.doi.org/10.1126/science.aam9695

mix_noref implements a maximum likelihood estimate of the proportions of genomic mixtures, assessing wheter 1,2 or 3 or more genomes appear in the DNA samples and the phylogenetic relationships between the different genomic components (if present).
Note that the maximum likelihood formula implemented here makes the assumption of independent sites and neglects deamination. Hence we caution that regions dense in mutations, as well as mismappings, could lead to inflated likelihood estimates. For this reason in Slon et al.,2017 we examined separately transversion mutations (excluding deamination) and removed polymorphisms in close proximity. Such filters should be optimized depending on the specific dataset. In addition we suggest to consider as significant only mixtures strongly supported by low p-values/relative-likelihoods, e.g. <0.001.

In Slon et al.,2017, to distinguish endogenous mixtures from mixtures caused by the presence of modern contaminant DNA, we examined separately deaminated readsonly and all reads together. When one additional component was not detectable in the deaminated reads, but when all reads were analyzed, we considered the additional component as likely due to contamination.

The R code requires the installation of the packages optimx and tidyverse.

# Usage

A typical command line is:
./mix_noref.R -i inputfile -o outputfile
where inputfile is a 5 column tab separated file with position as first column and then the counts of A,G,C,T bases for each position, e.g.
1 0 0 10 12
2 0 20 0 0
....
16568 0 0 9 0
To extract such tables from a .bam file one can use samtools mpileup as in bam2tab.sh

The output file is a table structured with the models (1,2 or 3 genomes) represented as rows, from best to worst. The fields are:
-ngenomes: number of genomes in the model: 1,2 or 3/more
-p1,p2 and p3: the proportions of the genomes.
-perr: error rate. Note that error rates here are the probability that a base is mutated at random into one of the 4 bases.
-X1_2: for the model with 2 genomes, divergence between genome 1 and 2
-X1_1_2, X1_2_1, X1_2_2, X1_2_3: for the model with 3 genomes, the proportion of sites in which genomes 2 and 3 differ or are the same. For example, X1_1_2 indicates those for which the third genome differ from genome 1 and 2, X1_2_1 those for which the second differs from genome 1 and 3, and so on.
-value: log likelihood of the model
-xtime: time to compute the model
-degree.freedom
-AIC: Akaike information criterion for the model. Note that because of the composite-likelihood approximation these values should be considered meaningful only in relative terms. Used to sort models from best to worst and to compute relative likelihoods.
-lik_ratio_test_2g, lik_ratio_test_3g: p-value based on the likelihood ratio test comparing the model with 2genomes (or 3) versus the model in the given row
-relL_vsbest,relL_vs2ndbest: relative likelihood comparing the model in given row to the best and second best

For more information
./mix_noref.R --help






