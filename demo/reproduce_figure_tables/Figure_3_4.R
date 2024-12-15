rm(list = ls())
library(ggplot2)
library(dplyr)
setwd('/home/wd278/BCRA')
path_to_data <- "./demo/reproduce_figure_tables/data_to_replicate_figure_tables/"

## Figure 3

detect_rate_sparsity <- NULL
error_rate_sparsity <- NULL
for(p_ind in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5)){
  load(paste0(path_to_data, "single_SNPset_sparsity/single_SNPset_sparsity_", p_ind, "_error_FALSE.RData"))
  tmp <-  suppressMessages({bind_rows(as.data.frame(snp_set %>% 
                                                      mutate(weight = factor(weight, levels = c("constant", "probability", "chisquared"))) %>%
                                                      group_by(Method, weight) %>%
                                                      summarise(average = mean(detect))),
                                      data.frame(Method  = "GWAS", weight = "GWAS", average  = mean(gwas_rate[,1]))
  ) %>%
      mutate(p_snp = p_ind)
  })
  detect_rate_sparsity <- rbind(detect_rate_sparsity, tmp)
  
  load(paste0(path_to_data, "single_SNPset_sparsity/single_SNPset_sparsity_", p_ind, "_error_TRUE.RData"))
  tmp <- suppressMessages({
    bind_rows(as.data.frame(snp_set %>%
                              mutate(weight = factor(weight, levels = c("constant", "probability", "chisquared"))) %>%
                              group_by(Method, weight) %>%
                              summarise(average = mean(detect))),
              data.frame(Method  = "GWAS", weight = "GWAS", average  = mean(gwas_rate[,1]))
    ) %>%
      mutate(p_snp = p_ind)
  })
  error_rate_sparsity <- rbind(error_rate_sparsity, tmp)
}
rm(list = setdiff(ls(), c("detect_rate_sparsity", "error_rate_sparsity", "path_to_data")))

label0 = expand.grid(c("BCRA", "subsample-BCRA"), 
                     c("constant", "probability", "chisquared")) %>%
  arrange(Var1)
label0 = paste0(label0$Var1, ":", label0$Var2)

### Average detect rate  ###
detect_rate_sparsity %>% 
  mutate(label = paste0(Method, ":", weight), proportion = as.factor(p_snp)) %>%
  mutate(label = factor(label, levels = c(label0, "GWAS:GWAS"),
                        labels = c(label0, "GWAS"))) %>%
  ggplot(aes(x = proportion, y = average,
             color = label,
             group = label,
             shape = label)) +
  geom_jitter(size = 5, width = 0, height = 0) +
  geom_line(aes(linetype = label), linewidth = 1) +
  labs(x = "Proportion of True Signals",
       y = paste0("Average Detection Rate")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Method (Weight)" ,
                     values = c(brewer.pal(8, "Reds")[4:6],
                                brewer.pal(8, "Greens")[4:6],
                                brewer.pal(8, "Blues")[4]),
                     labels = c(label0, "GWAS")) + 
  scale_shape_manual(name = "Method (Weight)" ,
                     values = c(rep(c(19, 18, 17), 2), 15),
                     labels = c(label0, "GWAS")) +
  scale_linetype_manual(name = "Method (Weight)", 
                        values = c(rep(c(1, 2, 3), 2), 4),
                        labels = c(label0, "GWAS")) +
  scale_y_continuous(limits = c(0, 1)) 


### Type-I Error  ###
error_rate_sparsity %>% 
  mutate(label = paste0(Method, ":", weight), proportion = as.factor(p_snp)) %>%
  mutate(label = factor(label, levels = c(label0, "GWAS:GWAS"),
                        labels = c(label0, "GWAS"))) %>%
  ggplot(aes(x = proportion, y = average,
             color = label,
             group = label,
             shape = label)) +
  geom_jitter(size = 5, width = 0, height = 0) +
  geom_line(aes(linetype = label), linewidth = 1) +
  geom_hline(yintercept = 0.05, color = "grey") +
  labs(x = "Proportion of True Signals",
       y = paste0("Type-I Error")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Method (Weight)" ,
                     values = c(brewer.pal(8, "Reds")[4:6],
                                brewer.pal(8, "Greens")[4:6],
                                brewer.pal(8, "Blues")[4]),
                     labels = c(label0, "GWAS")) + 
  scale_shape_manual(name = "Method (Weight)" ,
                     values = c(rep(c(19, 18, 17), 2), 15),
                     labels = c(label0, "GWAS")) +
  scale_linetype_manual(name = "Method (Weight)", 
                        values = c(rep(c(1, 2, 3), 2), 4),
                        labels = c(label0, "GWAS")) +
  scale_y_continuous(limits = c(0, 0.1)) 

## Figure 4: P-value comparison
case <- 'rank2_crossing'
load(paste0(path_to_data, 'pvalue_summary_', case, '_additive_FALSE_error_FALSE.RData'))
colnames(subsample_BCRA_pvalues)[1:3] = c("subsample_BCRA_constant", "subsample_BCRA_prob", "subsample_BCRA_chisquared")
colnames(BCRA_pvalues)[1:3] = c("BCRA_constant", "BCRA_prob", "BCRA_chisquared")

par(mfrow=c(2,2), mar = c(1,1,1,1), pty="s")
plot(sort(BCRA_pvalues[,1]), sort(subsample_BCRA_pvalues[,1]), 
     col = "grey", 
     xlab= "BCRA: p-value", ylab = "subsample-BCRA: p-value", 
     main = "weight: constant",
     xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, col = "red", lty=2, lwd = 2)

plot(sort(BCRA_pvalues[,2]), sort(subsample_BCRA_pvalues[,2]), col = "grey", 
     xlab= "BCRA: p-value", ylab = "subsample-BCRA: p-value", 
     main = "weight: probability",
     xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, col = "red", lty=2, lwd = 2)

plot(sort(BCRA_pvalues[,3]), sort(subsample_BCRA_pvalues[,3]), col = "grey", 
     xlab= "BCRA: p-value", ylab = "subsample-BCRA: p-value", 
     xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, col = "red", lty=2, lwd = 2)
dev.off()


