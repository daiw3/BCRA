rm(list = ls())
setwd("/home/wd278/BCRA")
library(snpStats)
library(Ball)
library(Matrix)
library(Rcpp)
library(huge)

sub_p <- 0.1 ## control % of subset samples

n <- 2000 ## number of subjects 
num.snps.per.blk <- 500 ## number of SNPs per SNP-set
num.blks <- 1 ## number of SNP-sets

p_snp <- 0.01 ## % of SNPs with true signals in power simulation
error <- FALSE ## FALSE: power simulation; TRUE: Type-I Error

graph <- "identity" ## error structure of E
noise_level <- 1 ## \sigma for error
non_linear_setting <- "none" ## linear/non-linear between X and Y
dist_name <- "pearson" ## Euclidean/Pearson/Geodesic distance
weight <- 'constant' ## constant/probability/chi-squared for Bcov

iter <- 1 

#### Generate Simulated Data Set, save as simu_data_demo.RData ####

# ---- load part of chr19 data (0-3MB) of UKB ----
data.chr.blk <- read.plink("/gpfs/gibbs/project/zhang_heping/wd278/brainFC_tarv/Simulation_subsample/data/geno_simu")

## ---- extract geno type and conduct QC ----
geno <- data.chr.blk$genotypes
member <- data.chr.blk$fam$member

geno <- geno[sample(1:nrow(geno), n, replace = F),, drop = FALSE]
member <- as.numeric(rownames(geno))

# per-SNP quality control for call rate and HWE test 
snpsum <- col.summary(geno)
# remove SNPs with call rate less than 90% and MAF < 0.05
call.cut <- 0.9
minor <- 0.01
call.keep <- with(snpsum, ((!is.na(MAF)) & (MAF > minor)) & ((!is.na(Call.rate) & Call.rate>= call.cut)))
call.keep[is.na(call.keep)] <- FALSE
geno <- geno[, call.keep]
map.call.keep <- data.chr.blk$map[call.keep, ]
snpsum.call.keep <- snpsum[call.keep, ]

# remove SNPs with HWE test p-value less than 10-6
HWE.cut.off <- 10^-6
hwe.keep <- with(snpsum.call.keep, (!is.na(z.HWE) & (abs(z.HWE)<abs(qnorm(HWE.cut.off/2)))))
hwe.keep[is.na(hwe.keep)] <- FALSE
geno <- geno[, hwe.keep]
map <- map.call.keep[hwe.keep, ]
# snpsum <- snpsum.call.keep[hwe.keep, ]

rm(list = setdiff(ls(), c("non_linear_setting", "graph", "dist_name",
                          "geno", "map", "member", 
                          "n", "p_snp", "num.snps.per.blk", 
                          "noise_level", "num.blks", "error", "iter", "sub_p",
                          "weight")))

# Randomly sample 500 SNPs for a single block
snp.sample <- sample(ncol(geno), num.snps.per.blk * num.blks)
geno <- geno[, snp.sample]
map <- map[snp.sample, ]

# code checking
snpsum <- col.summary(geno)
geno <- t(as(geno, "numeric"))
coding.test <- snpsum[, 8] - snpsum[, 6]
geno[which(coding.test>0), ] <- 2 - geno[which(coding.test>0), ]
geno <- as.data.frame(t(geno))
rm(list = setdiff(ls(), c("non_linear_setting", "graph", "dist_name",
                          "geno", "map", "member", "snpsum", 
                          "num.snps.per.blk", "n", "p_snp", 
                          "noise_level", "num.blks", "error", "iter", "sub_p", "weight")))

# impute NAs with MAF
for(j in 1:ncol(geno)){
  geno_imp = sample(c(0,1,2),size = sum(is.na(geno[, j])), prob = unlist(c(snpsum[j,c(8,7,6)])), replace =T)
  geno[which(is.na(geno[,j])),j] = geno_imp
}


## divided people into 2 sets
discovery.sample1 = sample(1:nrow(geno), nrow(geno)/2, replace = F)
discovery.sample2 = setdiff(1:nrow(geno), discovery.sample1)

while(any(apply(geno[discovery.sample1,],2,function(x) length(unique(x))) == 1) || any(apply(geno[discovery.sample1,],2,function(x) length(unique(x))) == 1)){
  discovery.sample1 = sample(1:nrow(geno), nrow(geno)/2, replace = F)
  discovery.sample2 = setdiff(1:nrow(geno), discovery.sample1)
}


rm(list=setdiff(ls(), c("non_linear_setting", "graph", "dist_name",
                        "geno", "discovery.sample1", "discovery.sample2", 
                        "map", "member", "snpsum", "num.snps.per.blk", "num.blks",
                        "n", "p_snp", "noise_level", "error", "iter", "sub_p", "weight")))


# ---- B: four cases (0/1 matrix) ---- #####
B_True = as.matrix(read.csv("./data/emoji_smile_adverse.csv", header = T))
B_True[is.na(B_True)] = 0
dimnames(B_True) = NULL
B_True_vec <- as.vector(B_True)

rm(list=setdiff(ls(), c("non_linear_setting", "graph", "dist_name",
                        "geno", "map", "member","snpsum", "num.snps.per.blk", "num.blks", "error",
                        "B_True","B_True_vec", "sub_p",
                        "discovery.sample1", "discovery.sample2", "n", "p_snp", "noise_level", "iter", "weight")))



## Add beta onto B ~ U(0.4, 0.8) ##
sig = matrix(sample(ncol(geno), 
                    size = num.snps.per.blk * p_snp, 
                    replace = F),
             ncol = 5, byrow = T)

B_sig = runif(nrow(sig), 0.4, 0.8)
true_sig=unlist(c(sig))

## Convert X_ij to S_ig
transformed_geno=NULL
for(i in 1:dim(sig)[1])
{
  if(non_linear_setting == "none"){
    transformed_geno=cbind(transformed_geno,rowSums(geno[,sig[i,], drop = FALSE]) * B_sig[i])
  }
  
  if(non_linear_setting == "interaction_sum"){
    tmp <- cbind(geno[,sig[i,], drop = FALSE], geno[,sig[i,1], drop = FALSE] * geno[,sig[i,2], drop = FALSE] + geno[,sig[i,1], drop = FALSE] * geno[,sig[i,3], drop = FALSE] + geno[,sig[i,3], drop = FALSE] * geno[,sig[i,4], drop = FALSE] * geno[,sig[i,5], drop = FALSE])
    transformed_geno=cbind(transformed_geno,
                           B_sig[i] * rowSums(tmp))
  }
  
  if(non_linear_setting == "sin_cos"){
    tmp <- cbind(sin(geno[,sig[i,1:3], drop = FALSE]), 
                 cos(geno[,sig[i,4:5], drop = FALSE]),
                 sin(geno[,sig[i,1], drop = FALSE] * geno[,sig[i,2], drop = FALSE]) + 
                   cos(geno[,sig[i,1], drop = FALSE] * geno[,sig[i,3], drop = FALSE]) + 
                   sin(geno[,sig[i,3], drop = FALSE] * geno[,sig[i,4], drop = FALSE] * geno[,sig[i,5], drop = FALSE])
    )
    transformed_geno=cbind(transformed_geno,
                           B_sig[i] * rowSums(tmp))
  }
}

rm(list=setdiff(ls(), c("non_linear_setting",  "graph", "dist_name", "sub_p",
                        "geno","sig", "B_sig", "transformed_geno", "map", "member", "num.snps.per.blk", "num.blks", "error", "iter",
                        "discovery.sample1", "discovery.sample2", "snpsum", "n", 
                        "p_snp", "noise_level", "B_True","B_True_vec", "weight", "true_sig")))

# Y is generated based on Y = f(S) + E
Y=matrix(NA,nrow=dim(transformed_geno)[1],ncol=length(B_True_vec))
for(i in 1:dim(transformed_geno)[1])
{
  Y[i,]=colSums(t(transformed_geno[i,,drop = FALSE]) %*% t(as.matrix(B_True_vec)))
}
rm(B_True_vec)

## Simulate Omega for Error (E) (V, a symmetric positive definitive matrix in Wishart distribution)
if(graph == "identity"){
  omega = diag(1, dim(B_True)[2])
}
if(graph == "random"){
  set.seed(1195)
  omega.list <- huge.generator(n = 10, d = dim(B_True)[2], 
                               graph = "random", v = 1, u = 0.0001, prob = 0.1, vis = F)
  omega = omega.list$sigma
  theta = omega.list$theta
}
if(graph == "band"){
  set.seed(1195)
  omega.list <- huge.generator(n = 10, d = dim(B_True)[2], 
                               graph = "band", v = 1, u = 0.0001, g = 5, vis = F)
  omega = omega.list$sigma
  theta = omega.list$theta
}
if(graph == "hub"){
  set.seed(1195)
  omega.list <- huge.generator(n = 10, d = dim(B_True)[2], 
                               graph = "hub", v = 1, u = 0.0001, g = 5, vis = F)
  omega = omega.list$sigma
  theta = omega.list$theta
}

if(graph == "cluster"){
  set.seed(1195)
  omega.list <- huge.generator(n = 10, d = dim(B_True)[2], 
                               graph = "cluster", v = 1, u = 0.0001, g = 5, prob = 0.6, vis = F)
  omega = omega.list$sigma
  theta = omega.list$theta
}


## Adding E onto Y
ee <- rWishart(dim(geno)[1], dim(B_True)+10, omega)
ee <- aperm(ee, c(3,1,2))
normalized_ee <- t(apply(ee,1,function(x) scale(c(x), scale = T)))
if(error){
  Y = matrix(0, nrow = dim(geno)[1], ncol = dim(B_True)[1]*dim(B_True)[1])
}
Y_noise=Y+normalized_ee
Y_noise[is.na(Y_noise)]=0

pheno_list = list()
for(i in 1:dim(Y_noise)[1]){
  pheno_list[[i]] = matrix(Y_noise[i, ],
                           nrow = sqrt(dim(Y_noise)[2]),
                           ncol = sqrt(dim(Y_noise)[2]))
}
Y_upper = sapply(pheno_list, function(x) x[upper.tri(x)])
Y_upper = as.data.frame(t(Y_upper))
rm(ee, Y)

Y_eigen <- as.data.frame(t(sapply(pheno_list, function(x){
  tmp <- svd(x)
  tmp$d
})))

if(!dist_name %in% c("pearson", "geodesic")){
  dist_pheno_sample1 = dist(Y_upper[discovery.sample1,,drop = FALSE],
                            method = dist_name)
  dist_pheno_sample2 = dist(Y_upper[discovery.sample2,,drop = FALSE],
                            method = dist_name)
}

if(dist_name %in% c("pearson")){
  dist_pheno_sample1 = 1 - cor(t(Y_upper[discovery.sample1,,drop = FALSE]),
                               method = dist_name)
  dist_pheno_sample2 = 1 - cor(t(Y_upper[discovery.sample2,,drop = FALSE]),
                               method = dist_name)
}

if(dist_name %in% c("geodesic")){
  sourceCpp("./code/GeodesicDistance.cpp")
  start <- Sys.time()
  dist_pheno_sample1 = GeodesicDistance(pheno_list[discovery.sample1], 1e-3)
  dist_pheno_sample2 = GeodesicDistance(pheno_list[discovery.sample2], 1e-3)
  cat("Time to compute geodesic distace matrix: ", as.double(Sys.time() - start, "mins"), " mins.\n")
}

cat("#######Data Simulation Finished. ########## \n")

save.image("./data/simu_data_demo.RData")
#############################################################
################### End of Data Simulation ######################
#############################################################

