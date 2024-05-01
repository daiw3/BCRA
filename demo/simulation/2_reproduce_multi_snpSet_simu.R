rm(list = ls())
setwd('/home/wd278/BCRA')

library(Matrix)
error_B <- function(B, B_orig){
  
  perf <- rep(NA, times = 4)
  names(perf) <- c("SEN", "SPE", "PREC", "NPV")
  tp <- sum(B != 0 & B_orig != 0)
  tn <- sum(B == 0 & B_orig == 0)
  fp <- sum(B != 0 & B_orig == 0)
  fn <- sum(B == 0 & B_orig != 0)
  
  sen <- tp/(tp + fn)
  spe <- tn/(tn + fp)
  prec <- tp/(tp + fp)
  acc <- (tp + tn)/(tp + tn + fp + fn)
  fdr <- fp / (tp+fp)
  npv <- tn / (tn + fn)
  ppt <- tp/sum(B_orig != 0)
  
  perf["SEN"] <- sen
  perf["SPE"] <- spe
  perf["PREC"] <- prec
  perf["NPV"] <- npv
  # perf["PPT"] <- ppt
  return(perf)
}

error_B_1MB_window <- function(B, B_orig){
  
  perf <- rep(NA, times = 4)
  names(perf) <- c("SEN", "SPE", "PREC", "NPV")
  
  idx_window = sapply(genoBIM[which(B !=0), "position"], function(x){
    which((genoBIM$position >= x-50000) & (genoBIM$position <= x + 50000))
  })
  B[unique(unlist(idx_window))] = 1
  
  
  idx_window_orig = sapply(genoBIM[which(B_orig !=0), "position"], function(x){
    which((genoBIM$position >= x-50000) & (genoBIM$position <= x + 50000))
  })
  B_orig[unique(unlist(idx_window_orig))] = 1
  
  tp <- sum(B != 0 & B_orig != 0)
  tn <- sum(B == 0 & B_orig == 0)
  fp <- sum(B != 0 & B_orig == 0)
  fn <- sum(B == 0 & B_orig != 0)
  
  sen <- tp/(tp + fn)
  spe <- tn/(tn + fp)
  prec <- tp/(tp + fp)
  npv <- tn / (tn + fn)
  
  perf["SEN"] <- sen
  perf["SPE"] <- spe
  perf["PREC"] <- prec
  perf["NPV"] <- npv
  return(perf)
}

load('./data/simu_geno_ld.RData')

result.path <- './results/multi_SNPset/'
file.list <- list.files(result.path, full.names = F, pattern = "multi_snpSet*")


thr <- 0.05/30 ## significance level for BCRA/subsample-BCRA

detection_rate_iters <- NULL ## detection rate for each SNP-set (Table 2 & Table S3-S5)
SNP_four_metrics_iters <- NULL ## SEN, SPE, PREC, NVP at SNP-set level (Table 3 & Table S6-S8)
SNPset_four_metrics_iters <- NULL ## SEN, SPE, PREC, NVP at SNP level (Table 3 & Table S6-S8)
for(ifile in 1:length(file.list)){
  f = file.list[ifile]
  iter = ifile
  load(paste0(result.path, "/", f))
  
  ### SNPset-level detect rate ###
  detect_subsample <- as.data.frame(sapply(pvalue.results.sub.BCRA[c(3,6,9)], function(x) as.numeric(x < thr)))
  colnames(detect_subsample) <- c('constant', 'probability', 'chisquared')
  detect_subsample$SNP_set_index <- 1:num.blks
  detect_subsample$iteration <- iter
  detect_subsample$Method <- 'subsample-BCRA'
 
  detection_rate_iters <- rbind(detection_rate_iters, 
                                detect_subsample %>% filter(SNP_set_index %in% c(1,2,3)))
  
  ### SNPset-level sen + spe + prec + npv ###
  B_blk_true <- c(rep(1, 3), rep(0, 27))
  SNPset_four_metrics <- as.data.frame(t(sapply(1:3, function(i){
    error_B(detect_subsample[,i], B_blk_true)
  })))
  SNPset_four_metrics$weight <- c('constant', 'probability', 'chisquared')
  SNPset_four_metrics$iteration <- iter
  SNPset_four_metrics$Method <- 'subsample-BCRA'
  SNP_four_metrics_iters <- rbind(SNP_four_metrics_iters,
                                  SNPset_four_metrics)
  
  ### SNP-level sen + spe + prec + npv
  constant_idx = unlist(pvalue.results.sub.BCRA[[2]][which(detect_subsample$constant == 1)])
  probability_idx = unlist(pvalue.results.sub.BCRA[[5]][which(detect_subsample$probability == 1)])
  chisquared_idx = unlist(pvalue.results.sub.BCRA[[8]][which(detect_subsample$probability == 1)])
  
  B_snp_window <- matrix(0, nrow = 15000, ncol = 3)
  B_snp_window[constant_idx, 1] <- 1
  B_snp_window[probability_idx, 2] <- 1
  B_snp_window[chisquared_idx, 3] <- 1
  B_snp_window_true <- rep(0, 15000)
  B_snp_window_true[true_sig] = 1
  
  SNP_four_metrics <- as.data.frame(
    t(sapply(1:3, function(i) error_B_1MB_window(B_snp_window[, i], B_snp_window_true))
    )
  )
  SNP_four_metrics$weight <- c('constant', 'probability', 'chisquared')
  SNP_four_metrics$iteration <- iter
  SNP_four_metrics$Method <- 'subsample-BCRA'
  
  SNP_four_metrics_iters <- rbind(SNP_four_metrics_iters,
                                  SNP_four_metrics)
}

############################################################
## detection rate for each SNP-set (Table 2 & Table S3-S5) ##
detection_rate_table <- detection_rate_iters %>% group_by(Method, SNP_set_index) %>%
  summarise(detection_rate_constant = mean(constant), 
            detection_rate_probability = mean(probability), 
            detection_rate_chisquared = mean(chisquared))

######################################################################
## SNP-set and SNP-level SEN, SPE, PREC, NPV (Table 3 & Table S6-S8) ##

SNPset_table <- SNPset_four_metrics_iters %>% group_by(Method, weight) %>%
  summarise(mean(SEN), mean(SPE), mean(PREC), mean(NPV))

SNP_table <- SNP_four_metrics_iters %>% group_by(Method, weight) %>%
  summarise(mean(SEN), mean(SPE), mean(PREC), mean(NPV))
