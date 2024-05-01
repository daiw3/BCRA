rm(list = ls())
setwd("/home/wd278/BCRA")

library(snpStats)
library(Ball)
library(Matrix)
library(Rcpp)
library(huge)
source('./code/0-tarv-transform-simu.R')

load("./data/simu_data_multiset_demo.RData")

################################################
################ BRCA ################ 
################################################
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
annotation=as.character(c(1:dim(geno)[2]-1)%/%num.snps.per.blk) 
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
### Calculate permutation p-values for BCRA
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

#############################################
############### subsample-BCRA ##############
#############################################

cat("Starting subsample-BCRA ...\n")

dist_pheno_sample1 <- as.matrix(dist_pheno_sample1)
call_rate <- apply(geno[discovery.sample1,,drop = FALSE], 1, function(x) mean(x>0))
samp_ord <- order(call_rate, decreasing = T)
n_subset = (n/2) * 0.1
#######################################
### Ranking 
#######################################
cat("subsample-BCRA Ranking...\n")
start <- Sys.time()
tmp = try(bcorsis(x = geno[discovery.sample1[samp_ord[1:n_subset]],,drop = FALSE], 
                  y = dist_pheno_sample1[samp_ord[1:n_subset], samp_ord[1:n_subset]],
                  distance = T, category = T))
end = Sys.time() -start
Z_BCRA_sub = tmp$complete.info$statistic[,2] * (snpsum$MAF * (1-snpsum$MAF))
# plot(sort(Z_BCRA_sub, decreasing = T))

# unlist(sapply(true_sig, function(x) which(order(Z_BCRA_sub, decreasing = T) == x)))

#######################################
### Transformation
#######################################
cat("subsample-BCRA Transformation...\n")
annotation=as.character(c(1:dim(geno)[2]-1)%/%num.snps.per.blk) 
geno.transed.BCRA.sub <- TARV_transform(t(geno[discovery.sample1, ]), 
                                        annotation, 
                                        Z_BCRA_sub, 
                                        direction="dominant") # transform into tarv required form
cat("Done. subsample-BCRA Transformation.\n")

#######################################
### Cut-off selection
#######################################
cat("subsample-BCRA Cut-off Selection...\n")

find.cutoff.BCRA.sub.func <- function(var.idx){
  genotmp <- geno.transed.BCRA.sub[, var.idx, drop = FALSE]
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
  maf <- apply(geno.cut.df[samp_ord[1:n_subset],,drop = FALSE], 2, function(x){
    Ref_allele <- 2 * sum(x == 2, na.rm = T) + sum(x == 1, na.rm = T)
    alternate_allele <- 2 * sum(x == 0, na.rm = T) + sum(x == 1, na.rm = T)
    (2*Ref_allele)/(2*(Ref_allele + alternate_allele))
  })
  tmp = try(bcorsis(x = geno.cut.df[samp_ord[1:n_subset],,drop = FALSE], 
                    y = dist_pheno_sample1[samp_ord[1:n_subset],samp_ord[1:n_subset]],
                    category = T,
                    distance = T))
  bd_stat = cbind(cut.off.cand, tmp$complete.info$statistic * (maf * (1-maf)))
  
  cutoff_df <- NULL
  for(weight in c("constant", "probability", "chisquare")){
    cut.best <- bd_stat[which.max(bd_stat[, paste0("bcor.", weight)]), 1]
    t.max <- bd_stat[which.max(bd_stat[,paste0("bcor.", weight)]),paste0("bcor.", weight)]
    cutoff_df = c(cutoff_df, cut.best, t.max)
  }
  return(c(var.list[var.idx], cutoff_df))
}

blk.size <- dim(geno)[2] / dim(geno.transed.BCRA.sub)[2]
var.list <- colnames(geno.transed.BCRA.sub)  
num.var <- dim(geno.transed.BCRA.sub)[2]
start <- Sys.time()
cutoff.results.BCRA.sub <- sapply(seq(num.var), find.cutoff.BCRA.sub.func)
cutoff.results.BCRA.sub <- as.data.frame(t(cutoff.results.BCRA.sub))
colnames(cutoff.results.BCRA.sub) <- c("super.var", 
                                       "cutoff_constant", "constant",
                                       "cutoff_probability", "probability",
                                       "cutoff_chisquare", "chisquare")
end <- Sys.time() - start
cat("Done. subsample-BCRA Cut-off Selection. Time: ", as.double(end, "mins"), ".\n")


#######################################
### Calculate p-values for subsample-BCRA
#######################################
blks = unique(as.numeric(annotation))

dist_pheno_sample2 <- as.matrix(dist_pheno_sample2)
pvalues_calculation_subBCRA <- function(weight){
  # start <- Sys.time ()
  selected.snps = list()
  X.test = NULL
  
  for(blk in 1:length(blks)){
    blk_ind = (num.snps.per.blk * blks[blk]+1):(num.snps.per.blk*(blks[blk]+1))
    cut.best <- as.numeric(as.character(cutoff.results.BCRA.sub[blk, paste0("cutoff_", weight)]))
    selected.snps[[blk]] = num.snps.per.blk*blks[blk] + order(Z_BCRA_sub[blk_ind], decreasing = T)[1:cut.best]
    genotmp = as.numeric(apply(geno[discovery.sample2, selected.snps[[blk]], drop = FALSE] > 0, 1, any))
    X.test = cbind(X.test, genotmp)
  }
  colnames(X.test) = var.list
  
  maf <- apply(X.test, 2, function(x){
    Ref_allele <- 2 * sum(x == 2) + sum(x == 1)
    alternate_allele <- 2 * sum(x == 0) + sum(x == 1)
    (2*Ref_allele)/(2*(Ref_allele + alternate_allele))
  })
  
  suppressMessages({
    perms <- t(sapply(1:num_perm, function(ii){
      idx = sample(1: nrow(X.test), n_subset, replace = F)
      tmp = bcorsis(x = X.test[idx,,drop = FALSE], 
                    y = dist_pheno_sample2[idx,idx],
                    distance = TRUE, category = T)
      tmp$complete.info$statistic[,paste0("bcor.", weight)]
    }))
  })
  
  suppressMessages({obs = bcorsis(x = X.test, y = dist_pheno_sample2,
                                  distance = TRUE, category = T)})
  obs = obs$complete.info$statistic[,paste0("bcor.", weight)]
  obs = obs / (maf * (1-maf))
  # cat("Histogram of Permutations: \n")
  # hist(perms[,1])
  # abline(v = obs[1], col = 'red')
  p_subsample = sapply(1:num.blks, function(i){
    (1 + sum(perms[,i] >= obs[i])) / (1+nrow(perms))
  })
  return(list(X.test = X.test, 
              selected.snps = selected.snps, 
              pvalue_results = p_subsample))
}
pvalue.results.sub.BCRA = sapply(c("constant", "probability", "chisquare"), pvalues_calculation_subBCRA)
cat("Final p-values for subsample-BCRA with weights (constant, probability, chisquare) are: \n")
sapply(pvalue.results.sub.BCRA[c(3,6,9)], function(x) x)


