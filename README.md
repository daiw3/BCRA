# BCRA

Ball Covariance Ranking and Aggregation (BCRA) hypothesis test

## Code

Required codes are located under ./**code/\***, including `Ball-1.3.12.tar.gz`, `0-tarv-transform-simu.R` and `BallDistanceVector.cpp`.

## Data

Data used for simulation demo are located under **./data/**\*.  

### Single SNP set simulation

1) *simu_data_demo*.RData: simulated data for single SNP-set simulation, generated from the code `./demo/0-single-snpset-generate-simu-data.R`. As the genotype data used for simulations are derived from UKB realdata, it has sensitive information to be shared. This dataset does not contain the full data we used for simulation but the code can be used to reproduce the full data we used for simulation once people get the access to UKB genotype data.

### Multi-SNP set simulation

1)  Multi-set \$\mathbf{B} as shown in\$Fig.2 are listed in files in the order: *matrixB3.mat*, *matrixB4.mat*, *emoji_smile_adverse.csv* and *emoji_make.xlsx*.
2)  *simu_data_multiset_demo*.RData: simulated data for multi SNP-set simulation, generated from the code `./demo/simulation/0-multi-snpset-generate-simu-data.R`. As the genotype data used for simulations are derived from UKB realdata, it has sensitive information to be shared. This dataset does not contain the full data we used for simulation but the code can be used to reproduce the full data we used for simulation once people get the access to UKB genotype data.

### UKB Real Data Demonstration
4)  *chr5_8_demo.RData*: a make-up data for showing how to apply subsample-BCRA on UKB as real genotype data contains sensitive information.

## Instructions to Perform Simulations 

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

2.  **./demo/simulation/2_reproduce_multi_snpSet_simu.R** provides the code to calculate detection rate, SEN, SPE, PREC, NVP for Multi SNP-set simulation. We put one example output for the power simulation with the crossing case as in Table 2 under the path `./results/Multiset_Table1_3_S3_to_S9.RData`. This can reproduce the results of Table 1 to 3, S3 to S9.

3. **./demo/reproduce_figure_tables/*.R** contain the code to replicate all figures and tables in the simulation studies, including Fig.3-4, Fig.S1-S3, Fig.S5, Table S1 to S2. 

## Replicate Tables and Figures for Simulations

1) Table 1 to 3, S3 to S9: **./demo/simulation/2_reproduce_multi_snpSet_simu.R**
2) Fig.3-4: **./demo/reproduce_figure_tables/Figure_3_4.R**
3) Fig.S1-S3, Fig.S5, Table S1 to S2: **./demo/reproduce_figure_tables/FigureS1_S5_Table_S1_S2.R**

## Real Data Application (UKB)
### Step 1: Calculate geodesic distance among functional connectivity matrices
1.  Calculate geodesic distance (other distance measures can be applied depending on your interest) of functional connectivity matrices across each pair of subjects.
### Step 2: Partition Genotype into SNP-set
2.  Partition Genotype into SNP-set: you can use gene/LD-block or physical locations (what we adopted) to do the partition. PLINK software can be used to fullfill the partition.
### Step 3: Genotype QC on each SNP-set
3.  Genotype QC: we exclude: 1) subjects with more than 10% missing genotypes; 2) variants with missing genotype rate larger than 10%; and 3) variants that failed the Hardy-Weinberg test at $10^{-6}$ level.
### Step 4: Run subsample-BCRA for each SNP-set
4.  Run subsample-BCRA on a SNP-set using code located under `./demo/UKB/1-run-subsample-BCRA-SNPset.R`. Since we are unable to provide sensitive real data, we provided a demo dataset named as \*./data/chr5_8_demo.RData\* for illustration purpose to execuate the code. The output of the code will give selected SNPs (`selected_snps_int`) and permutation p-value for this SNP-set `pval_results_perm`.
### Step 5: Summarize p-values associated with each SNP-set
5. We put a demonstration code `./demo/reproduce_figure_tables/FigureS1_S3_S5_Table_S1_S2.R` to plot p-values associated with each SNP-set (Figure S5) with the data file `./demo/reproduce_figure_tables/pval_UKB.RData` gives the SNPs selected with permutaed p-value for each SNP-set presented in the manuscript.
### Step 6: Visualize the results
6. The connectivity plots in the manuscript (Figure 6 and Figure S6-S1) can be generated using `chordDiagram` function in R.

## Running Time and Resource Requirements

### 1. Single SNP-set Simulation
- **Memory**: 32 GB per CPU
- **CPUs**: 4  
- **Time**: ~1 hour per iteration for one parameter set  
  - Example: Generating results for \(\pi = 0.01\) (proportion of true signals) as shown in Figure 4.

### 2. Multi SNP-set Simulation
- **Memory**: 32 GB per CPU
- **CPUs**: 4  
- **Time**: ~40 hours per iteration for one parameter set for BCRA and ~4 hours per iteration under the same resource conditions for subsample-BRCA.
  - Example: Generating results for \(\mathbf{B}\) in the "smile" case for BRCA.

### 3. UK Biobank Real Data Application
#### a. Distance Calculation
- **Memory**: At least 100 GB per CPU
- **CPUs**: 4  
- **Time**: ~1 week  
- **Note**: Additional memory is required to write and store distance matrices.

#### b. BCRA/Subsample Analysis
  - **Memory**: 32 GB per CPU 
  - **CPUs**: 4  
  - **Time**: BCRA: ~1 week per SNP-set with 5000 permutations. Subsample-BCRA: 2 days per SNP-set with 5000 permutations

=======
Ball Covariance Ranking and Aggregation (BCRA) hypothesis test
