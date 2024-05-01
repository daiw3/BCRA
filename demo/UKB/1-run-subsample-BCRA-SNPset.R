rm(list = ls())
library(snpStats) # snpStats for loading genetic data
library(Ball)
library(data.table)


iter = 1
chr = 5
blk = 8
# num_perm <- 999999
num_perm <- 99

setwd('/home/wd278/BCRA')
source("./code/0-BCRA-transformation-UKB.R") 

##################################################
#### BELOW CODE READ GENOTYPE, PHENOTYPE, IDs ####
### PLEASE UNCOOMENT THEM FOR YOUR OWN ANALYSIS ###

## load pheno data on discovery set
# load("./pheno.Data")
## load discovery set subject ID
# load("./UKB_discovery_id.RData")

# file_prefix <- paste0("chr", chr, "_", blk)

## Need the genotype arrary for chr5_8, which is after QC ###
# data.chr.blk <- read.plink(paste0("/vast/palmer/scratch/zhang_heping/wd278/brainFC_tarv_March/di-tarv-UKB_updated/UKB_updated-tarv-discover/divided/chr", chr, "/", file_prefix))

# geno.use <- data.chr.blk$genotypes
# map.finalmap.final <- data.chr.blk$map
# member <- data.chr.blk$fam$member
# 
# use.id = intersect(member, discovery_id)
# geno.use = geno.use[as.character(use.id), ]
# full_dist_mat_dis = full_dist_mat_dis[discovery_id %in% use.id, discovery_id %in% use.id]
# 
# snp.summary <- col.summary(geno.use)
# geno.use <- as(geno.use, "numeric")
# geno.use[,which(snp.summary[, 6] < snp.summary[, 8])] = 2 - geno.use[,which(snp.summary[, 6] < snp.summary[, 8])]
# 
# ## partition into 2 sets for future use
# set.seed(iter)
# id_set1 = sample(1:nrow(full_dist_mat_dis), round(nrow(full_dist_mat_dis)/2), replace = F)
# id_set2 = setdiff(1:nrow(full_dist_mat_dis), id_set1)

### Save Above Results into chr5_8_demo.RData
# save(geno.use, full_dist_mat_dis,
#      id_set1, id_set2, snp.summary, chr, blk, num_perm, iter,
#      file = '/home/wd278/BCRA/data/chr5_8_demo.RData')

rm(list = ls())
setwd('/home/wd278/BCRA')
source("./code/0-BCRA-transformation-UKB.R") 

load('./data/chr5_8_demo.RData')

cat("Starting subsample-BCRA Ranking...\n")
cat("There are ", ncol(geno.use), " SNPs waited for ranking. \n")

start.time = Sys.time()
call_rate <- apply(geno.use[id_set1,,drop = FALSE], 1, function(x) mean(x>0, na.rm = T))
samp_ord <- order(call_rate, decreasing = T)

n_subset = (nrow(geno.use)/2) * 0.1
if(n_subset < 2) n_subset = 2
tmp = suppressMessages({
  try(bcorsis(x = geno.use[id_set1[samp_ord[1:n_subset]], , drop = FALSE],
              y = full_dist_mat_dis[id_set1[samp_ord[1:n_subset]], id_set1[samp_ord[1:n_subset]], drop = FALSE],
              distance = T))
})
Z.local = tmp$complete.info$statistic[, 2] * (snp.summary$MAF * (1-snp.summary$MAF))
names(Z.local) <- rownames(snp.summary)
end = Sys.time() -start.time
cat("Ranking is Done. Time is ", as.double(end, "mins"), "\n")

# ----------  Transformation -----------------------------
cat("Starting subsample-BCRA Transformation ...\n")

# add annotation for transformation
file_prefix <- paste0("chr", chr, "_", blk)
annotation <- rep(file_prefix, dim(geno.use)[2])
geno.use1 <- geno.use[id_set1,, drop = FALSE]
geno.transed <- TARV_transform(t(geno.use1), 
                               annotation, 
                               Z.local, 
                               direction="dominant") # transform into tarv required form
cat("Transformation Done. \n")

# ------------Cut-off Selection ------------------------------------
cat("Starting Cut-off Selection ...\n")
find.cutoff.func <- function(var.idx){
  genotmp <- geno.transed
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
  tmp = suppressMessages({try(bcorsis(x = geno.cut.df[samp_ord[1:n_subset],,drop = FALSE], 
                                      y = full_dist_mat_dis[id_set1[samp_ord[1:n_subset]], id_set1[samp_ord[1:n_subset]]],
                                      distance = T))})
  maf <- apply(geno.cut.df[samp_ord[1:n_subset],,drop = FALSE], 2, function(x){
    Ref_allele <- 2 * sum(x == 2, na.rm = T) + sum(x == 1, na.rm = T)
    alternate_allele <- 2 * sum(x == 0, na.rm = T) + sum(x == 1, na.rm = T)
    (2*Ref_allele)/(2*(Ref_allele + alternate_allele))
  })
  bd_stat = cbind(cut.off.cand, tmp$complete.info$statistic * (maf * (1-maf)))
  
  cutoff_df <- NULL
  for(weight in c("constant", "probability", "chisquare")){
    cut.best <- bd_stat[which.max(bd_stat[, paste0("bcor.", weight)]), 1]
    t.max <- bd_stat[which.max(bd_stat[,paste0("bcor.", weight)]),paste0("bcor.", weight)]
    cutoff_df = c(cutoff_df, cut.best, t.max)
  }
  c(colnames(geno.transed), cutoff_df)
}
start <- Sys.time()
cutoff.results <- sapply(seq(ncol(geno.transed)), find.cutoff.func)
cutoff.results <- as.data.frame(t(cutoff.results))
colnames(cutoff.results) <- c("super.var", 
                              "cutoff_constant", "constant",
                              "cutoff_probability", "probability",
                              "cutoff_chisquare", "chisquare")
end <- Sys.time() - start
cat("Done. Cut-off Selection. Time: ", as.double(end, "mins"), ".\n")


# ------------extract information of contributing SNPs ------------------------------------
select_snps_func <- function(weight){
  z.ord <- order(Z.local, decreasing = T)
  # select the top SNPs according to the cut-off value
  snps.selected <- snp.summary[z.ord, ]
  snps.selected <- snps.selected[seq(as.numeric(as.character(cutoff.results[, paste0('cutoff_', weight)]))), ]
  snps.selected$cM <- cutoff.results$super.var
  colnames(snps.selected)[10] <- "super.var"
  snps.selected$snp.name <- rownames(snps.selected)
  row.names(snps.selected) <- NULL
  mag.effect.snps <- Z.local[snps.selected$snp.name]
  mag.effect.snps <- matrix(mag.effect.snps, ncol = 1)
  colnames(mag.effect.snps) <- c("sum_stat")
  snps.summary <- snp.summary[snps.selected$snp.name, ]
  snps.selected <- cbind(snps.selected, mag.effect.snps, snps.summary)
  row.names(snps.selected) = NULL
  return(snps.selected)
}
selected_snps_int = lapply(c("constant", "probability", "chisquare"), select_snps_func)


rm(list = setdiff(ls(),
                  c("id_set1", "id_set2", "n_subset",
                    "Z.local", # ranking results
                    "geno.transed", # transformed genotype for candidate cutoff
                    "cutoff.results", # cut-off results
                    "selected_snps_int", # selected SNPs for each super-variant
                    "snp.summary",
                    "iter", "chr", "blk", "num_perm",
                    "geno.use", "full_dist_mat_dis")))

# ----------  Marginal p-values ------------------------

cat("Start calculting marginal p-values...\n")
start.time <- Sys.time()
pvalues_calculation_BCRA_sub <- function(weight){
  start <- Sys.time ()
  
  if(weight == "constant"){
    i = 1
  }
  if(weight == "probability"){
    i = 2
  }
  if(weight == "chisquare"){
    i = 3
  }
  X.gene = as.numeric(apply(geno.use[id_set2, selected_snps_int[[i]]$snp.name, drop = FALSE] > 0, 1, any, na.rm = T))
  X.gene[is.na(X.gene)] <- 0
  
  X.gene <- data.frame(member=rownames(geno.use[id_set2,]), X.gene, stringsAsFactors = F)
  colnames(X.gene)[2] <- as.character(cutoff.results$super.var)
  
  suppressMessages({obs = bcorsis(x = X.gene[, 2, drop = FALSE],
                                  y = full_dist_mat_dis[id_set2, id_set2],
                                  distance = TRUE)})
  obs = obs$complete.info$statistic[,paste0("bcor.", weight)]
  maf <- apply(X.gene[, -1,drop = FALSE], 2, function(x){
    # mean(x > 0)
    Ref_allele <- 2 * sum(x == 2) + sum(x == 1)
    alternate_allele <- 2 * sum(x == 0) + sum(x == 1)
    (2*Ref_allele)/(2*(Ref_allele + alternate_allele))
  })
  obs = obs / (maf * (1-maf))
  
  suppressMessages({
    perms <- t(sapply(1:num_perm, function(ii){
      idx = sample(1: nrow(X.gene), log((nrow(geno.use)/2)), replace = F)
      tmp = bcorsis(x = X.gene[idx, -1, drop = FALSE],
                    y = full_dist_mat_dis[id_set2[idx], id_set2[idx]],
                    distance = TRUE)
      tmp$complete.info$statistic[,paste0("bcor.", weight)]
    }))
  })
  perms = c(perms)
  
  pvalue_results =(1 + sum(perms >= obs)) / (1+length(perms))
  
  
  pvalue_results = data.frame(super_variant = paste0("chr", chr, "_", blk),
                              p = pvalue_results,
                              time = as.double(Sys.time () - start, "mins"),
                              obs = obs,
                              mean_perms = mean(perms, na.rm = T),
                              sd_perms = sd(perms, na.rm = T)
  )
  
  colnames(pvalue_results) = c("super_variant", "p", "time(mins)", "obs", "mean_perm", "sd_perm")
  
  return(list(pvalue_results = pvalue_results, perms = perms, X.gene = X.gene))
}
pval_results_perm = lapply(c("constant", "probability", "chisquare"), pvalues_calculation_BCRA_sub)
cat("subsample-BCRA p-value Calculation Done. Time: ", as.double(Sys.time() - start.time, "mins"), ".\n")

cat("P-value for SNP-set: CHR-", chr, "-Set-", blk, ', is: \n')
sapply(pval_results_perm, function(x) x$pvalue_results)
