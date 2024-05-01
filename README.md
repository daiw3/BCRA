# BCRA

Ball Covariance Ranking and Aggregation (BCRA) hypothesis test

## Code

Required codes are located under ./**code/\***, including `Ball-1.3.12.tar.gz`, `0-tarv-transform-simu.R` and `BallDistanceVector.cpp`.

## Data

Data used for simulation demo are located under **./data/**\*. The genotype data used for simulations is too large to put here. It can be sent upon request

1)  *matrixB4.mat*: The \$\mathbf{B} as shown in \$Fig.2c
2)  *simu_data_demo*.RData: simulated data for single SNP-set simulation, generated from the code `0-single-snpset-generate-simu-data.R`
3)  *simu_data_multiset_demo*.RData: simulated data for single SNP-set simulation, generated from the code `0-single-snpset-generate-simu-data.R`
4)  *chr5_8_demo.RData*: make-up data for showing how to apply subsample-BCRA on UKB.

## Simulations 

### Step 1: Generate simulated Dataset

1.  **./demo/simulation/0-single-snpset-generate-simu-data.R** is the code to generate simulated dataset using real genotype dataset from UKB for *single SNP-set simulation*. It has multiple parameter inputs, including nonlinear settings (`non_linear_setting`), error correlation structures (`graph`), proportion of true signals (`p_snp`), proportion of subsets used in subsample-BCRA (`sub_p`), distance measurements (`dist_name`), noise level (`noise_level`)

    Given data privacy and capacity constraints, we did not put the original genotype data online. We generated a demo dataset `./data/simu_data_demo.RData` based on the this code for illustration purpose.

2.  **./demo/simulation/0-multi-snpset-generate-simu-data.R** is a demonstration code to generate synthetic data using real genotype dataset from UKB for *Multi SNP-set simulation*.

    Given data privacy and capacity constraints, we did not put the original genotype data online. We generated a demo dataset `./data/simu_data_multiset_demo.RData` based on the this code for illustration purpose.

### Step 2a: Run Single-SNP set simulation

1.  **./demo/simulation/1-simu_single_snpSet_subsample_BCRA.R** is a demonstration code to repeat the single SNP-set simulation for subsample-BCRA using the demo simulated dataset. It will generate result for 1 iteration. We put an example output under \`./results/single_SNPset/\`
2.  **./demo/simulation/1-simu_single_snpSet_BCRA_GWAS.R** is a demonstration code to repeat the single SNP-set simulation for BCRA and GWAS using the demo simulated dataset. It will generate result for 1 iteration. We put an example output under \`./results/single_SNPset/\`.

### Step 2b: Run Multi-SNP set simulation

1.  **./demo/simulation/1-simu_multi_snpSet_subsample_BCRA.R** is a demonstration code to repeat the single SNP-set simulation for both BCRA and subsample-BCRA using the demo simulated dataset. We put an example output under \`./results/multi_SNPset/\`.

### Step 3: Summarize simulation results and replicate Tables and Figures in the manuscript

1.  **./demo/simulation/2_reproduce_single_snpSet_simu.R** provides the code to calculate detection rate for single SNP-set simulation. We put one example output for the power simulation with $\pi$ = 0.01 (proportion of true signals) as in Fig.4 under the path `./results/Singleset_Figure4.RData`. This can reproduce the results of Fig.4 for $\pi$ = 0.01.

2.  **./demo/simulation/2_reproduce_multi_snpSet_simu.R** provides the code to calculate detection rate, SEN, SPE, PREC, NVP for Multi SNP-set simulation. We put one example output for the power simulation with the crossing case as in Table 2 under the path `./results/Multiset_Table2_3_S3_to_S8.RData`. This can reproduce the results of Table 2, 3, S3 to S8.

## Real Data Application (UKB)

1.  Calculate geodesic distance (other distance measures can be applied depending on your interest) of functional connectivity matrices across each pair of subjects.
2.  Partition Genotype into SNP-set: you can use gene/LD-block or physical locations (what we adopted) to do the partition.
3.  Genotype QC: we exclude: 1) subjects with more than 10% missing genotypes; 2) variants with missing genotype rate larger than 10%; and 3) variants that failed the Hardy-Weinberg test at $10^{-6}$ level.
4.  Run subsample-BCRA on a SNP-set using code located under `./demo/UKB/1-run-subsample-BCRA-SNPset.R`. Since we are unable to provide sensitive real data, we provided a demo dataset named as \*./data/chr5_8_demo.RData\* for illustration purpose to execuate the code. The output of the code will give selected SNPs (`selected_snps_int`) and permutation p-value for this SNP-set `pval_results_perm`.
=======
Ball Covariance Ranking and Aggregation (BCRA) hypothesis test
