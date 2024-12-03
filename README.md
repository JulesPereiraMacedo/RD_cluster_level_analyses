# Cluster-level analyses to estimate a risk difference in a cluster randomized trial with confounding individual-level covariates: a simulation study

This repository contains all the functions and programs to reconduct the simulation study of our paper:

Jules Antoine Pereira Macedo, Bruno Giraudeau, for the ESCIENT collaborators. Cluster-level analyses to estimate a risk difference in a cluster randomized trial with confounding individual-level covariates: a simulation study. (Journal name). (DATE). (DOI)

## Abstract
Cluster randomized trials (CRTs) may be analyzed using cluster-level analyses. For binary outcomes, proportions are estimated for each cluster, and a risk difference can be estimated. The confidence interval is estimated using a Student distribution. However, in doing so, individual-level characteristics are not adjusted for even though CRTs are known to be prone to recruitment/identification bias, possibly implying individual-level confounders. With a simulation study, we compared cluster-level analyses to estimate a risk difference for a two-arm parallel CRT with individual-level confounders and cluster-level covariates. We considered the unadjusted (UN) method, two two-stage procedure (TSP) methods considering a binomial or a Gaussian distribution, G-computation (GC), and targeted maximum likelihood estimation (TMLE) methods. As expected, the UN method was biased. TSP methods were also biased for scenarios with a treatment effect when the number of clusters per arm was small. GC and TMLE methods were unbiased. For these latter methods, adjustment on only individual-level covariates led to better performance measures (type I error rate, coverage rate and relative error of the standard error) than adjustment on both individual- and cluster-level covariates. In parallel, we observed convergence problems with GC methods. In the end, the TMLE method considering only individual-level covariates was the preferred approach to estimate a risk difference.

## Getting Started
### Dependencies
Data generation and analyses were run on a computing center of the Orléans-Tours university with R version 4.0, to parallelize faster but it can be run on a regular computer.
### Installing
To run all the files, start to run the package.R files to install all the package and to be able to lunch the simulation.
### Executing programs
This simulation study has 216 scenarios.

For each R files take care to change your file path where you want the simulation to be saved. It is named “Base_file =” in the R files (path in line 1675 in R file 1_Workspace_cluster_simu.R, for example).

We have 7 files, numeroted from 1 to 4. It needs to be run in the right order (1 → 4). And for the files 3a, 3b, 3c and 3d they can be run simultaneously

"1_Workspace_cluster_simu":
This R file create a R workspace. You will have to load it for each of the following R files. All the functions used are in this file, with comments, and also this file creates “list” variables "Scenarios_%n" from 1 to 216,  which correspond to the scenarios parameters for each scenarios.

"2_Data_cluster":
This R file creates and saves all the dataset (1000 dataset for each of the 216 scenarios).
The dataset of the scenarios were generated using the seed "10091958".
 
"3_a_Analyse_UW_TSP":
This R file does the analysis for the Unajusted method and the two stage procedure approach, with both Gaussian and binomial distribution. All methods are implemented in the Workspace file.

"3_b_Analyse_TMLE":
This R file does the analysis for the TMLE approach, adjusting on only individual covariates and on both level of covariates (individual and cluster levels). All methods are implemented in the Workspace file.

"3_c_Analyse_Gcomp_indv":
This R file does the analysis for the G-computation approach, adjusting on individual level covariates. This method is implemented in the Workspace file.

"3_d_Analyse_Gcomp_both":
This R file does the analysis for the G-computation approach, adjusting on both level of covariates (individual and cluster levels). This method is implemented in the Workspace file.

## Help
Users should contact pereiramacedo@univ-tours.fr (Pereira Macedo JA) if they have issues for running the programs.
## Authors
Jules Antoine Pereira Macedo: pereiramacedo@univ-tours.fr
