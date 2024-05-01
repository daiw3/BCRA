rm(list = ls())
setwd("/home/wd278/BCRA")
source('./code/0-tarv-transform-simu.R')

library(snpStats)
library(Ball)
library(Matrix)
library(Rcpp)
library(huge)
load("./data/simu_data_demo.RData")

############### subsample-BCRA #######
annotation=as.character(c(1:dim(geno)[2]-1)%/%num.snps.per.blk) # 1 block

cat("Starting subsample-BCRA ...\n")


dist_pheno_sample1 <- as.matrix(dist_pheno_sample1)
call_rate <- apply(geno[discovery.sample1,,drop = FALSE], 1, function(x) mean(x>0))
samp_ord <- order(call_rate, decreasing = T)
n_subset = (n/2) * sub_p
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
cat("TARV Transformation...\n")
geno.transed.BCRA.sub <- TARV_transform(t(geno[discovery.sample1, ]), 
                                        annotation, 
                                        Z_BCRA_sub, 
                                        direction="dominant") # transform into tarv required form
cat("Done. TARV Transformation.\n")

#######################################
### Cut-off selection
#######################################
cat("TARV Cut-off Selection...\n")

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
cat("Done. TARV Cut-off Selection. Time: ", as.double(end, "mins"), ".\n")


#######################################
### Calculate final marginal p-values for BCRA
#######################################
start.time <- Sys.time()

dist_pheno_sample2 <- as.matrix(dist_pheno_sample2)
start <- Sys.time ()
selected.snps = list()
X.test = NULL

blk = 1
blk_ind = 1:num.snps.per.blk
cut.best <- as.numeric(as.character(cutoff.results.BCRA.sub[blk, paste0("cutoff_", weight)]))
selected.snps[[blk]] = order(Z_BCRA_sub[blk_ind], decreasing = T)[1:cut.best]
genotmp = as.numeric(apply(geno[discovery.sample2, selected.snps[[blk]], drop = FALSE] > 0, 1, any))
X.test = cbind(X.test, genotmp)
colnames(X.test) = var.list

maf <- apply(X.test, 2, function(x){
  Ref_allele <- 2 * sum(x == 2) + sum(x == 1)
  alternate_allele <- 2 * sum(x == 0) + sum(x == 1)
  (2*Ref_allele)/(2*(Ref_allele + alternate_allele))
})

suppressMessages({
  perms <- t(sapply(1:199, function(ii){
    idx = sample(1: nrow(X.test), n_subset, replace = F)
    tmp = bcorsis(x = X.test[idx,,drop = FALSE], 
                  y = dist_pheno_sample2[idx,idx],
                  distance = TRUE, category = T)
    tmp$complete.info$statistic[1,paste0("bcor.", weight)]
  }))
})

suppressMessages({obs = bcorsis(x = X.test, y = dist_pheno_sample2,
                                distance = TRUE, category = T)})
obs = obs$complete.info$statistic[1,paste0("bcor.", weight)]
obs = obs / (maf * (1-maf))
cat("Observed Bcov: ", obs, "; p = ", (1 + sum(perms >= obs)) / (1+length(perms)), "\n")
cat("Histogram of Permutations: \n")
hist(perms)
abline(v = obs, col = 'red')

