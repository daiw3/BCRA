rm(list = ls())
setwd("/home/wd278/BCRA")
source('./code/0-tarv-transform-simu.R')

library(snpStats)
library(Ball)
library(Matrix)
library(Rcpp)
library(huge)
load("./data/simu_data_demo.RData")



#############################################################
################### BCRA START ######################
#############################################################

cat("Starting BCRA ...\n")

#######################################
### Ranking
#######################################
cat("BCRA Ranking...\n")
start <- Sys.time()
tmp = try(bcorsis(x = geno[discovery.sample1,,drop = FALSE], 
                  y = dist_pheno_sample1,
                  distance = T, category = T))
end = Sys.time() -start
Z_BCRA = tmp$complete.info$statistic[,2]
cat("Done. BCRA Ranking. Time: ", as.double(end, "mins"), ".\n") # ~113.1078s for 15000 SNPs on 2000 samples


#######################################
### Transformation
#######################################
cat("BCRA Transformation...\n")
annotation=as.character(c(1:dim(geno)[2]-1)%/%num.snps.per.blk) # 1 block
geno.transed.BCRA <- TARV_transform(t(geno[discovery.sample1, ]), 
                                    annotation, 
                                    Z_BCRA, 
                                    direction="dominant") # transform into tarv required form
cat("Done. BCRA Transformation.\n")

#######################################
### Cut-off selection
#######################################
cat("BCRA Cut-off Selection...\n")

find.cutoff.BCRA.func <- function(var.idx){
  genotmp <- geno.transed.BCRA[, var.idx, drop = FALSE]
  cut.off.cand <- unique(genotmp[,1])
  if(length(cut.off.cand) == 1){
    next
  }
  cut.off.max <- max(cut.off.cand)
  cut.off.cand <- cut.off.cand[-which(cut.off.cand==cut.off.max)]
  
  geno.cut.df = NULL
  for(cut.off in cut.off.cand){
    geno.cut <- apply(genotmp, 2, function(x) as.numeric(x <= cut.off) )
    geno.cut.df = cbind(geno.cut.df, geno.cut)
  }
  tmp = try(bcorsis(x = geno.cut.df, 
                    y = dist_pheno_sample1,
                    category = T,
                    distance = T))
  bd_stat = cbind(cut.off.cand, tmp$complete.info$statistic)
  
  cutoff_df <- NULL
  for(weight in c("constant", "probability", "chisquare")){
    cut.best <- bd_stat[which.max(bd_stat[, paste0("bcor.", weight)]), 1]
    t.max <- bd_stat[which.max(bd_stat[,paste0("bcor.", weight)]),paste0("bcor.", weight)]
    cutoff_df = c(cutoff_df, cut.best, t.max)
  }
  return(c(var.list[var.idx], cutoff_df))
}

blk.size <- dim(geno)[2] / dim(geno.transed.BCRA)[2]
var.list <- colnames(geno.transed.BCRA)  
num.var <- dim(geno.transed.BCRA)[2]
start <- Sys.time()
cutoff.results.BCRA <- sapply(seq(num.var), find.cutoff.BCRA.func)
cutoff.results.BCRA <- as.data.frame(t(cutoff.results.BCRA))
colnames(cutoff.results.BCRA) <- c("super.var", 
                                   "cutoff_constant", "constant",
                                   "cutoff_probability", "probability",
                                   "cutoff_chisquare", "chisquare")
end <- Sys.time() - start
cat("Done. BCRA Cut-off Selection. Time: ", as.double(end, "mins"), ".\n")


#######################################
### Calculate final marginal p-values for BCRA
#######################################
blks = unique(as.numeric(annotation))
pvalues_calculation_BCRA <- function(weight){
  start <- Sys.time ()
  selected.snps = list()
  X.test  =  NULL
  
  for(blk in 1:length(blks)){
    blk_ind = (num.snps.per.blk * blks[blk]+1):(num.snps.per.blk*(blks[blk]+1))
    cut.best <- as.numeric(as.character(cutoff.results.BCRA[blk, paste0("cutoff_", weight)]))
    selected.snps[[blk]] = num.snps.per.blk*blks[blk] + order(Z_BCRA[blk_ind], decreasing = T)[1:cut.best]
    genotmp = as.numeric(apply(geno[discovery.sample2, selected.snps[[blk]], drop = FALSE] > 0, 1, any))
    X.test = cbind(X.test, genotmp)
  }
  colnames(X.test) = var.list
  
  pvalue_results = apply(X.test, 2, function(x){
    tmp = bcov.test(x = dist(x), 
                    y = dist_pheno_sample2,
                    num.permutations = 199,
                    distance = TRUE, 
                    weight = weight)
    return(c(tmp$p.value, as.double(Sys.time () - start, "mins")))
  })
  pvalue_results = t(pvalue_results)
  colnames(pvalue_results) = c("p", "time(mins)")
  
  return(list(X.test = X.test, 
              selected.snps = selected.snps, 
              pvalue_results = pvalue_results))
}
pvalue.results.BCRA = sapply(c("constant", "probability", "chisquare"), pvalues_calculation_BCRA)
cat("Final p-values for BCRA with weights (constant, probability, chisquare) are: ", sapply(pvalue.results.BCRA[c(3,6,9)], function(x) x[1]), ".\n")


#############################################################
################### mGWAS ######################
#############################################################
start.time <- Sys.time()
marginal.GWAS <- function(snp.id){
  m <- apply(Y_noise[discovery.sample1,  which(upper.tri(B_True)), drop = FALSE], 2, function(x){
    summary(lm(x ~ geno[discovery.sample1, snp.id]))$coefficients[2, ,drop = FALSE]
  })
  m <- t(m)
  m[, 4]
}
library(snowfall)
sfInit(parallel = T, cpus = 4)
sfExport("Y_noise", "B_True", "geno", "discovery.sample1")
start.time <- Sys.time()
gwas_pval <- sfSapply(1:ncol(geno), marginal.GWAS1)
# print(Sys.time() - start.time)
sfStop()
gwas.time <- as.double(Sys.time() - start.time, 'mins')

