## Mult-SNPset Simulation Data Generation"
## Wei Dai 2024/04/29

rm(list=ls())
library(snpStats)
library(Ball)
library(Matrix)
library(Rcpp)
library(huge)
library(R.matlab)

setwd('/home/wd278/BCRA') # change to the BCRA directory on your own PC
sourceCpp("./code/GeodesicDistance.cpp")

iter <- 1
error <- FALSE
num.snps.per.blk <- 500
num.blks <- 30
num_perm <- 1999
n <- 4000
noise_level <- 1

##load UKB Genotype data
# load("./data/multi_snpsets_geno.RData") # we're unable to provide this given privacy and capacity issues
load("~/project/brainFC_tarv/Simulation_RealData/UKB_chr16_30_blocks_500_SNPs_geno.RData")

geno <- do.call("cbind", geno_list)

## randomly select 4000 samples
geno <- geno[sample(1:nrow(geno), n), ]

## divided people into 2 sets
discovery.sample1 = sample(1:nrow(geno), nrow(geno)/2, replace = F)
discovery.sample2 = setdiff(1:nrow(geno), discovery.sample1)

snpsum <- as.data.frame(do.call("rbind", geno_snpsum))
##Transform geno data that 1 indicates rare variant
rm(list=setdiff(ls(), c("geno", "discovery.sample1", "discovery.sample2", "snpsum", 
                        "noise_level", "error", "num.blks", "num.snps.per.blk", "num_perm", "n")))

##Significant snps
sig_1=c(1,101,201,301,401,451)
sig_2=matrix(c(661,662,771,772,881,882),nrow=2)
sig_3=matrix(c(1011:1013,1221:1223, 1331:1333),nrow=3)
sigs=list(sig_1,sig_2,sig_3)

B_sig_1 = c(0.1, 0.3, 0.5, 0.5, 0.3, 0.1)
B_sig_2 = c(0.1, 0.3, 0.5)
B_sig_3 = c(0.1, 0.3, 0.5)
B_sigs=list(B_sig_1,B_sig_2,B_sig_3)

##Transform interactive snps
geno_sig_1= t(t(geno[,sig_1]) * B_sig_1)
geno_sig_2=NULL
for(i in 1:dim(sig_2)[2])
{
  geno_sig_2=cbind(geno_sig_2,as.integer(rowSums(geno[,sig_2[,i]])>0) * B_sig_2[i])
}
geno_sig_3=NULL
for(i in 1:dim(sig_3)[2])
{
  geno_sig_3=cbind(geno_sig_3,as.integer(rowSums(geno[,sig_3[,i]])>0) * B_sig_3[i])
}

transformed_geno=cbind(geno_sig_1,geno_sig_2,geno_sig_3)

rm(list=setdiff(ls(), c("geno","sigs", "B_sigs", 
                        "transformed_geno", "discovery.sample1", "discovery.sample2", 
                        "snpsum", "noise_level", "error", "num.blks", "num.snps.per.blk", "num_perm", "n")))

#### simulation true B #####
### case: smile
B_True = readMat("./data/matrixB4.mat")
B_True = B_True$matrixB4
B_True_vec <- as.vector(B_True)


rm(list=setdiff(ls(), c("geno","sigs","B_sigs", "transformed_geno","B_True","B_True_vec",
                        "discovery.sample1", "discovery.sample2", "snpsum",
                        "noise_level", "error", "num.blks", "num.snps.per.blk", "num_perm", "n")))

##Pheno
Y=matrix(NA,nrow=dim(transformed_geno)[1],ncol=length(B_True_vec))
for(i in 1:dim(transformed_geno)[1])
{
  Y[i,]=rowSums(B_True_vec %*% t(transformed_geno[i,]))
}
rm(B_True_vec)

annotation=as.character(c(1:dim(geno)[2]-1)%/%500) # 30 blocks

true_sig=unlist(sigs)

##Adding noise
ee <- matrix(rnorm(dim(Y)[1]*dim(Y)[2], mean=0, sd=noise_level),nrow=dim(Y)[1])
ee <- t(apply(ee, 1, function(ex){
  e_mat = matrix(ex, dim(B_True)[1], dim(B_True)[1])
  e_mat = as.matrix(forceSymmetric(e_mat))
  return(as.vector(e_mat))
}))
if(error){
  Y = matrix(0, nrow = dim(geno)[1], ncol = dim(B_True)[1]*dim(B_True)[1])
  Y_noise=Y+ee
  Y_noise[is.na(Y_noise)]=0
}else{
  Y_noise=Y+ee
  Y_noise[is.na(Y_noise)]=0
}

pheno_list = list()
for(i in 1:dim(Y_noise)[1]){
  pheno_list[[i]] = matrix(Y_noise[i, ],
                           nrow = sqrt(dim(Y_noise)[2]),
                           ncol = sqrt(dim(Y_noise)[2]))
}
Y_upper = sapply(pheno_list, function(x) x[upper.tri(x)])
Y_upper = as.data.frame(t(Y_upper))
rm(ee, Y)

dist_pheno_sample1 = 1 - cor(t(Y_upper[discovery.sample1,,drop = FALSE]),
                             method = "pearson")
dist_pheno_sample2 = 1 - cor(t(Y_upper[discovery.sample2,,drop = FALSE]),
                             method = "pearson")
cat("#######Data Simulation Finished. ########## \n")

save.image("./data/simu_data_multiset_demo.RData")

