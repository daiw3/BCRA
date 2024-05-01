rm(list = ls())
setwd('/home/wd278/BCRA')

library(dplyr)
library(ggplot2)
library(RColorBrewer)

error <- FALSE
p_snp <- 0.01 # proportion of true signals

result.path <- "./results/single_SNPset/"
files <- list.files(path = result.path, pattern = "*.RData")

snp_set = NULL
gwas_detection_rate = NULL

for(ifile in 1:length(files)){
  f = files[ifile]
  iter <- ifile
  
  load(paste0(result.path, f))
  
  
  d = data.frame(p = c(sapply(BCRA_results[[4]][c(3,6,9)], function(x) x[1]),
                       sapply(subsample_BCRA_results[[4]][c(3,6,9)], function(x) x[1])
  ),
  time = c(sapply(BCRA_results[[4]][c(3,6,9)], function(x) x[2]),
           sapply(subsample_BCRA_results[[4]][c(3,6,9)], function(x) x[2])
  )
  )
  d <- d %>% mutate(weight = rep(c("constant", "probability", "chisquared"), 2),
                    Method = rep(c("BCRA", "subsample-BCRA"), each = 3),
                    detect = as.numeric(p < 0.05),
                    iteration = iter)
  # gwas_detection_rate <- rbind(gwas_detection_rate,
  #                    c(mean(apply(gwas_pval[, true_sig], 2, min) < 0.05/(num.snps.per.blk*2016)), gwas.time)
  # )
  snp_set = rbind(snp_set, d)
  
}
rm(list = setdiff(ls(),
                  c("snp_set", "gwas_detection_rate",
                    "p_snp", "error_B", "error")))

### Average detect rate or error
load('./results/Singleset_Figure4.RData')
error <- FALSE
p_snp <- 0.01 # proportion of true signals

bind_rows(as.data.frame(snp_set %>% 
                          mutate(weight = factor(weight, levels = c("constant", "probability", "chisquared"))) %>%
                          group_by(Method, weight) %>%
                          summarise(average = mean(detect))),
          data.frame(Method  = "GWAS", weight = "GWAS", average  = mean(gwas_detection_rate[,1]))
) 


