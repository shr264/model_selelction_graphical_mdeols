# Inferring the Partial Correlation Structure of Allelic Effects and Incorporating it in Genome-wide Prediction

## Data description

There are two sets of simulated data, each one corresponds to scenarios 1 and 2 described in the paper and contain four replicates. Each replicate is a collection of three datasets containing: genotypic data (QTL genotypes), phenotypes along with additive genetic effects, residuals generation and other columns described below, and a file with individual QTL effects. Datasets simulated under scenario 1 contain 200 individuals, while those simulated under scenario 2 contain 380 individuals, the two scenarios considered 300 QTL, i.e., 300 regression variables. For further details on the differences between scenario 1 and 2 see the manuscript. On the other hand, the real dataset contains centered phenotypes corresponding to fertility, genotypes for 300 SNP that were found to be associated to fertility in Han and Peñagaricano (2016). This is a subset of the original one obtained after discarding half-sib families with less than 5 individuals. This editing process was necessary due to the fact that our frequentist methods require grouping individuals in families. 

## Data dictionary

Datasets Pop1_data_001,…, Pop1_data_004:
* Progeny: Individual ID.
* Sire: Individual’s sire ID.
* Dam: Individual’s dam ID.
* Sex: Individual’s gender, male (M) or female (F).
* G: Generation. 
* NMPrg: number of male progenies. 
* NFPrg: number of female progenies. 
* F: inbreeding coefficient
* Homo: percentage of homozygous loci.
* Phen: Phenotype, i.e., record.
* Res: Residual
* Polygene: It corresponds to a genetic effect different from the one corresponding to simulated QTL. In this simulation, this effect was null.
* QTL: Additive genetic effect coming from simulated QTL’s, in this simulation, the total additive genetic effect came from simulated QTL. It corresponds to the dot product between the vector containing coded genotypes (see manuscript for details) and the vector of individual QTL effects.

Datasets QTL_effects_001,…,QTL_effects_004: These files contain a single column corresponding to individual QTL effects. 
* V1: Individual QTL effects.

Datasets Pop1_qtl_001,…,Pop1_qtl_004: These files do not have column names. The firs column contains individual’s ID and it is followed by 600 numbers. Each consecutive pair corresponds to a QTL and codes the genotype; thus, 1 1 is a homozygous genotype for the first allele, 2 2 a homozygous genotype for allele 2, finally 1 2 and 2 1 indicate heterozygous genotypes. 

## Code description

The code used to perform the analyses described in the manuscript consists of seven R programs. There is a file to perform data analysis based on each one of our methods. Moreover, some files perform the required computations to obtain files that are used as inputs in other programs. Codes also include the computations used in model performance assessment, as described in the paper.
This set of programs permits to reproduce results shown in Table 1 and supplementary Table 1. Notice that these tables involve four replicates under two scenarios and that method Bayes DAG-Sel only works when there are more data than markers, that is, under scenario 2. Therefore, the corresponding code has to be run four times, while the remaining have to be executed eight times. For the sake of clarity, all codes are commented. 

* GLasso_EM: Implements the GLasso-EM method and computes some parameters measuring its performance. 
* CONCORD_EM: Implements the CONCORD-EM method and computes some parameters measuring its performance. 
* CSCS_EM: Implements the CSCS-EM method and computes some parameters measuring its performance.
* Model_Sel_Bayes_SSS: Implements the stochastic short-gun search algorithm to carry out Bayesian graphical model selection under a G-Wishart prior as described in the manuscript. 
* Model_Sel_Bayes_DAG_SSS: Implements the stochastic short-gun search algorithm to carry out Bayesian graphical model selection under a DAG-Wishart prior as described in the manuscript.
* Gibbs_Sampler_General_Graph: Performs Bayesian estimation of the precision matrix under a Gaussian concentration graph model once an undirected graph has been selected using Bayes G-Sel. 
* Gibbs_sampler_DAG_Wishart: Performs Bayesian estimation of the precision matrix under a Gaussian DAG model once a DAG has been selected using Bayes DAG-Sel. 




