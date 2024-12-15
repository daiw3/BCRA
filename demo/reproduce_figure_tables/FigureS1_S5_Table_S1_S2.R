rm(list = ls())
library(ggplot2)
library(dplyr)
library(RColorBrewer)

setwd('/home/wd278/BCRA')
path_to_data <- "./demo/reproduce_figure_tables/data_to_replicate_figure_tables/"

## Table S1 - nonlinear
p_snp = 0.05
error = FALSE
detect_rate <- data.frame(Method = c(rep(c("BCRA", "subsample-BCRA"), each = 3), "GWAS"),
                          weight = c(rep(c("constant", "probability", "chisquared"), 2), "GWAS"))

for(nonlinear in c("none", "interaction_sum", "sin_cos", "exponential", "BD", "PS")){
  if(nonlinear == "PS"){
    load(paste0(path_to_data, "single_two_nonlinear/two_SNPset_prop_", p_snp, "_nonlinear_", nonlinear, "_error_", error, ".RData"))
  }else{
    load(paste0(path_to_data, "single_two_nonlinear/single_SNPset_prop_", p_snp, "_nonlinear_", nonlinear, "_error_", error, ".RData"))
  }
  tmp <-suppressMessages({bind_rows(as.data.frame(snp_set %>% 
                                                    mutate(weight = factor(weight, levels = c("constant", "probability", "chisquared"))) %>%
                                                    group_by(Method, weight) %>%
                                                    summarise(average = mean(detect))),
                                    data.frame(Method  = "GWAS", weight = "GWAS", average  = mean(gwas_rate[,1]))
  ) 
  })
  detect_rate <- cbind(detect_rate, tmp %>% select(average)) 
}
colnames(detect_rate)[3:8] <- c("linear", "interaction_sum", "sin_cos", "exponential", "BD", "PS")

## Table S2 - correlated error

detect_rate <- data.frame(Method = c(rep(c("BCRA", "subsample-BCRA"), each = 3), "GWAS"),
                          weight = c(rep(c("constant", "probability", "chisquared"), 2), "GWAS"))
p_snp = 0.05
for(graph in c("identity", "random", "band", "hub", "cluster")){
  load(paste0(path_to_data, "single_SNPset_wishart_prop/single_SNPset_wishart_prop_", p_snp, "_graph_",graph, ".RData"))
  tmp <-suppressMessages({bind_rows(as.data.frame(snp_set %>% 
                                                    mutate(weight = factor(weight, levels = c("constant", "probability", "chisquared"))) %>%
                                                    group_by(Method, weight) %>%
                                                    summarise(average = mean(detect))),
                                    data.frame(Method  = "GWAS", weight = "GWAS", average  = mean(gwas_rate[,1]))
  ) 
  })
  detect_rate <- cbind(detect_rate, tmp %>% select(average)) 
}
colnames(detect_rate)[3:7] <- c("identity", "random", "band", "hub", "cluster")

## Figure S1 - distance
p_snp <- 0.05
d <- read.csv(paste0(path_to_data, "single_SNPset_prop_", p_snp, "_distance_detect_rate.csv"), header = T)
### Average detect rate  ###

d %>% mutate(label0 = paste0(Method, ":", weight)) %>% 
  reshape(varying = c("euclidean", "pearson", "geodesic"),
          v.names = "Average_detect_rate",
          idvar = "label0", direction = "long", 
          timevar = "distance_measure", 
          times = c("euclidean", "pearson", "geodesic")) %>%
  mutate(label = factor(label0, levels = c(label0, "GWAS:GWAS"),
                        labels = c(label0, "GWAS"))) %>%
  ggplot(aes(x = distance_measure, y = Average_detect_rate,
             color = label,
             fill = label,
             group = label,
             shape = label)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x = "3 Distance Measurements",
       y = paste0("Average Detection Rate")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name = "Method (Weight)" ,
                     values = c(brewer.pal(8, "Reds")[4:6],
                                brewer.pal(8, "Greens")[4:6],
                                brewer.pal(8, "Blues")[4])) + 
  scale_fill_manual(name = "Method (Weight)" ,
                    values = c(brewer.pal(8, "Reds")[4:6],
                               brewer.pal(8, "Greens")[4:6],
                               brewer.pal(8, "Blues")[4])) + 
  scale_shape_manual(name = "Method (Weight)" ,
                     values = c(rep(c(19, 18, 17), 2), 15)) +
  scale_linetype_manual(name = "Method (Weight)", 
                        values = c(rep(c(1, 2, 3), 2), 4)) +
  scale_y_continuous(limits = c(0, 1)) 

## Figure S3

snp_set_all <- NULL

for(n in c(500, 600, 700, 800, 900, 1000)){
  for(nSNPs in c(50, 300, 500, 650, 800, 1000)){
    for(sub_p in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)){
      load(paste0(path_to_data, "single_SNPset_nsub/single_SNPset_sample_", n, "_nSNPs_", nSNPs, "_psub_", sub_p, "_results.RData"))
      snp_set_all <- suppressMessages({rbind(snp_set_all, snp_set %>% 
                                               mutate(weight = factor(weight, levels = c("constant", "probability", "chisquared"))) %>%
                                               filter(Method == "subsample-BCRA") %>%
                                               group_by(n, nSNPs, proportion, Method, weight) %>%
                                               summarise(average = mean(detect))
      )
      })
    }
  }
}

label = expand.grid(c("subsample-BCRA"), c("consant", "probability", "chisquared")) %>%
  arrange(Var1)
label = paste0(label$Var1, ":", label$Var2)

### Average detect rate  ###
snp_set_all %>% 
  mutate(p_n_ratio = nSNPs / n) %>%
  group_by(proportion, p_n_ratio, weight) %>%
  summarise(avg = mean(average)) %>%
  arrange(desc(p_n_ratio)) %>%
  ggplot(aes(x = p_n_ratio, y= avg, 
             color = weight)) +
  geom_point(position =  position_dodge(width = 0.1), size = 3, aes(shape = weight), alpha = 0.8) + 
  scale_color_manual(name = "Weight" ,
                     values = c(brewer.pal(3, "Set1")),
                     label = label) +
  scale_shape_manual(name = "Weight" ,
                     values = c(15,16,17),
                     label = label) +
  facet_wrap(~proportion, nrow = 4) +
  theme_classic() +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0, 2), 
                     breaks = c(0.1, 0.2, 0.3, 0.5, 1, 1.5, 2)) +
  labs(x = "#SNPs / # subjects ",
       y = paste0("Detect Rate")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45))

## Figure S5: UKB P-values
load(path_to_data, 'pval_UBK.RData')
all_pvals_sp %>% group_by(chr, blk, iter, super_var) %>%
  summarise(p = min(min(constant, probability), chisquared)) %>%
  mutate(p = ifelse(p < 1e-6, 1e-6, p)) %>%
  ggplot(aes(x = iter, y = -log10(p), color = super_var)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05/2723), linetype = 'dashed', color = 'grey40') +
  facet_wrap(~super_var, nrow = 5) +
  theme_bw()


