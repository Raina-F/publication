
##### MR-pqtl-GAL1-interval #####
library(MRPRESSO) 
setwd("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/241106/GAL1-interval")
exposure <- read.csv("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/exposure_prot-a-1726.csv")

list_outcomes <- c("ebi-a-GCST006908")
for (i in 1:length(list_outcomes)){
  outcome_dat <- extract_outcome_data(snps = exposure$SNP, outcomes = (list_outcomes[i]), proxies = FALSE)  # Extract outcome data
  write.csv(outcome_dat, paste0("outcome_", list_outcomes[i], ".csv"))
  dat <- harmonise_data(exposure, outcome_dat, action = 3)   # Align alleles, handle palindromic sequences
  write.csv(dat, paste0("dat_", list_outcomes[i], ".csv"))
  
  dat$EAF2 <- (1 - dat$eaf.exposure)
  dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)
  dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.outcome)
  dat$FSTAT <- ((dat$samplesize.outcome - 1 - 1) / 1) * (dat$PVE / (1 - dat$PVE))  # F-statistic
  write.csv(dat, paste0("datF_", list_outcomes[i], ".csv"))
  results <- mr(dat)  # MR computation
  write.csv(results, paste0("Results_", list_outcomes[i], ".csv"))
  OR <- generate_odds_ratios(results)   # Odds Ratios  
  write.csv(OR, paste0("OR_", list_outcomes[i], ".csv"))
  
  # 1. Horizontal pleiotropy  2. Identify outlier SNPs  3. Compare results with and without outlier SNPs
  # mr_presso(BetaOutcome = 'beta.outcome',
  #           BetaExposure = 'beta.exposure', 
  #           SdOutcome = 'se.outcome', 
  #           SdExposure = 'se.exposure', 
  #           data = dat, OUTLIERtest = TRUE, 
  #           DISTORTIONtest = TRUE, SignifThreshold = 0.05, NbDistribution = 5000, seed = NULL)
  
  # Visualization
  library(TwoSampleMR)
  single <- mr_leaveoneout(dat)
  mr_leaveoneout_plot(single)
  ggsave(paste0("leaveoneout_", list_outcomes[i], ".pdf"), width = 6, height = 6)
  mr_scatter_plot(results, dat)
  ggsave(paste0("mr_scatter_", list_outcomes[i], ".pdf"), width = 6, height = 6)
  results_single <- mr_singlesnp(dat)
  mr_forest_plot(results_single)
  ggsave(paste0("forestPlot_", list_outcomes[i], ".pdf"), width = 6, height = 6)
  mr_funnel_plot(results_single)
  ggsave(paste0("funnelPlot_", list_outcomes[i], ".pdf"), width = 6, height = 6)
  het <- mr_heterogeneity(dat)
  write.csv(het, paste0("het_", list_outcomes[i], ".csv"))
  pleio <- mr_pleiotropy_test(dat)
  write.csv(pleio, paste0("pleio_", list_outcomes[i], ".csv"))
}

#####MR-pqtl-CD44-interval#####
library(MRPRESSO) 
setwd("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/241106/CD44-interval")
exposure <- read.csv("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/exposure_prot-a-449.csv")
list_outcomes <- c("ebi-a-GCST006908")
for (i in 1:length(list_outcomes)){
  outcome_dat <- extract_outcome_data(snps = exposure$SNP, outcomes = (list_outcomes[i]), proxies = TRUE)  # Extract outcome data
  write.csv(outcome_dat, paste0("outcome_", list_outcomes[i], ".csv"))
  dat <- harmonise_data(exposure, outcome_dat, action = 3)   # Align alleles, palindromic sequences
  write.csv(dat, paste0("dat_", list_outcomes[i], ".csv"))
  
  dat$EAF2 <- (1 - dat$eaf.exposure)
  dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)
  dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.outcome)
  dat$FSTAT <- ((dat$samplesize.outcome - 1 - 1)/1) * (dat$PVE/(1 - dat$PVE))  # F-statistic
  write.csv(dat, paste0("datF_", list_outcomes[i], ".csv"))
  results <- mr(dat)  # MR calculation
  write.csv(results, paste0("Results_", list_outcomes[i], ".csv"))
  OR <- generate_odds_ratios(results)   # OR (Odds Ratio)
  write.csv(OR, paste0("OR_", list_outcomes[i], ".csv"))
  
  # 1. Horizontal pleiotropy  2. Detect outlier SNPs  3. Check whether the results differ with or without outliers
  # mr_presso(BetaOutcome = 'beta.outcome',
  #           BetaExposure = 'beta.exposure', 
  #           SdOutcome = 'se.outcome', 
  #           SdExposure = 'se.exposure', 
  #           data = dat, OUTLIERtest = TRUE, 
  #           DISTORTIONtest = TRUE, SignifThreshold = 0.05, NbDistribution = 5000, seed = NULL)
  
  # Visualization
  library(TwoSampleMR)
  single <- mr_leaveoneout(dat)
  mr_leaveoneout_plot(single)
  ggsave(paste0("leaveoneout_", list_outcomes[i], ".pdf"), width = 5, height = 5)
  mr_scatter_plot(results, dat)
  ggsave(paste0("mr_scatter_", list_outcomes[i], ".pdf"), width = 4, height = 4)
  results_single <- mr_singlesnp(dat)
  mr_forest_plot(results_single)
  ggsave(paste0("forestPlot_", list_outcomes[i], ".pdf"), width = 5, height = 5)
  mr_funnel_plot(results_single)
  ggsave(paste0("funnelPlot_", list_outcomes[i], ".pdf"), width = 4, height = 4)
  het <- mr_heterogeneity(dat)
  write.csv(het, paste0("het_", list_outcomes[i], ".csv"))
  pleio <- mr_pleiotropy_test(dat)
  write.csv(pleio, paste0("pleio_", list_outcomes[i], ".csv"))
}




#####MR-pqtl-IL6R-interval#####
library(MRPRESSO) 
setwd("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/241106/IL6R-interval")
exposure <- read.csv("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/exposure_prot-a-1540.csv")
list_outcomes <- c("ebi-a-GCST006908")
for (i in 1:length(list_outcomes)){
  outcome_dat <- extract_outcome_data(snps = exposure$SNP, outcomes = (list_outcomes[i]), proxies = TRUE)  # Extract outcome data
  write.csv(outcome_dat, paste0("outcome_", list_outcomes[i], ".csv"))
  dat <- harmonise_data(exposure, outcome_dat, action = 3)   # Align alleles, palindromic sequences
  write.csv(dat, paste0("dat_", list_outcomes[i], ".csv"))
  
  dat$EAF2 <- (1 - dat$eaf.exposure)
  dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)
  dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.outcome)
  dat$FSTAT <- ((dat$samplesize.outcome - 1 - 1)/1) * (dat$PVE/(1 - dat$PVE))  # F-statistic
  write.csv(dat, paste0("datF_", list_outcomes[i], ".csv"))
  results <- mr(dat)  # MR calculation
  write.csv(results, paste0("Results_", list_outcomes[i], ".csv"))
  OR <- generate_odds_ratios(results)   # OR (Odds Ratio)
  write.csv(OR, paste0("OR_", list_outcomes[i], ".csv"))
  
  # 1. Horizontal pleiotropy  2. Detect outlier SNPs  3. Check whether the results differ with or without outliers
  # mr_presso(BetaOutcome = 'beta.outcome',
  #           BetaExposure = 'beta.exposure', 
  #           SdOutcome = 'se.outcome', 
  #           SdExposure = 'se.exposure', 
  #           data = dat, OUTLIERtest = TRUE, 
  #           DISTORTIONtest = TRUE, SignifThreshold = 0.05, NbDistribution = 5000, seed = NULL)
  
  # Visualization
  library(TwoSampleMR)
  single <- mr_leaveoneout(dat)
  mr_leaveoneout_plot(single)
  ggsave(paste0("leaveoneout_", list_outcomes[i], ".pdf"), width = 5, height = 5)
  mr_scatter_plot(results, dat)
  ggsave(paste0("mr_scatter_", list_outcomes[i], ".pdf"), width = 4, height = 4)
  results_single <- mr_singlesnp(dat)
  mr_forest_plot(results_single)
  ggsave(paste0("forestPlot_", list_outcomes[i], ".pdf"), width = 5, height = 5)
  mr_funnel_plot(results_single)
  ggsave(paste0("funnelPlot_", list_outcomes[i], ".pdf"), width = 4, height = 4)
  het <- mr_heterogeneity(dat)
  write.csv(het, paste0("het_", list_outcomes[i], ".csv"))
  pleio <- mr_pleiotropy_test(dat)
  write.csv(pleio, paste0("pleio_", list_outcomes[i], ".csv"))
}


######UKB-PPP Exposure Processing.TREM2######
# GENPOS is the Hg38 version, and the numbers in ID are from the Hg19 version
library(dplyr)
# TREM2
trem2_raw <- read.table("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/TREM2_Q9NZC2_OID20731_v1_Inflammation/discovery_chr6_TREM2:Q9NZC2:OID20731:v1:Inflammation.txt",
                        sep = " ",
                        header = T)
genpos_lower <- 41158506 - 500000
genpos_upper <- 41163186 + 500000
# Filter data
trem2_cis_pqtl <- trem2_raw %>%
  filter(GENPOS >= genpos_lower & GENPOS <= genpos_upper & LOG10P > 5)
chrom6map <- read_tsv("/Users/yuyao/Desktop/Bioinfo-Big-Files/olink_ukb/olink_rsid_map/olink_rsid_map_mac5_info03_b0_7_chr6_patched_v2.tsv.gz")
trem2_cis_pqtl <- trem2_cis_pqtl %>%
  left_join(chrom6map %>% select(ID, rsid), by = "ID") %>%
  rename(SNP = rsid)  # Rename the extracted rsid as SNP

head(trem2_cis_pqtl)
write.csv(trem2_cis_pqtl,"/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/trem2_cis_pqtl_chorom6_500kb_logp5.csv")

exposure <- TwoSampleMR::format_data(trem2_cis_pqtl,type='exposure',snp_col = "SNP",
                                     beta_col = "BETA",
                                     se_col = "SE",
                                     eaf_col = "A1FREQ",
                                     effect_allele_col = "ALLELE1",
                                     other_allele_col = "ALLELE0",
                                     pval_col = "LOG10P",
                                     samplesize_col = "N",
                                     log_pval = T)
write.csv(exposure,"/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/exposure_trem2_cis_pqtl_chorom6_500kb_logp5.csv")
exposure_clump <- clump_data(exposure,clump_kb = 10000, clump_r2 = 0.1)
write.csv(exposure_clump,file = "/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/exposure_trem2_cis_pqtl_chorom6_500kb_logp5_clump_0.1_10000.csv")
if (T) {
  setwd("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/smad3_stroke_screen")
  # exposure <- read_csv(list[i])
  Traitinfo <-extract_outcome_data(snps = exposure_clump$SNP,outcomes = list_stroke_GWASID,proxies=T)# Extract outcome data
  mr_result <- getMR_Atlas_modify(exposure,Traitinfo)# 5min
  write.csv(mr_result,paste0("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/","TREM2_AS_fast.csv"))
}

######UKB-PPP Exposure Processing.SMAD3######
# GENPOS is the Hg38 version, and the numbers in ID are from the Hg19 version
library(dplyr)
# SMAD3
smad3_raw <- read.table("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/SMAD3_P84022_OID31418_v1_Oncology_II/discovery_chr15_SMAD3:P84022:OID31418:v1:Oncology_II.txt",
                        sep = " ",
                        header = T)
genpos_lower <- 67063763 - 500000
genpos_upper <- 67195173 + 500000
# Filter data
smad3_cis_pqtl <- smad3_raw %>%
  filter(GENPOS >= genpos_lower & GENPOS <= genpos_upper & LOG10P > 5)
chrom15map <- read_tsv("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/olink_rsid_map_mac5_info03_b0_7_chr15_patched_v2.tsv")
smad3_cis_pqtl <- smad3_cis_pqtl %>%
  left_join(chrom15map %>% select(ID, rsid), by = "ID") %>%
  rename(SNP = rsid)  # Rename the extracted rsid as SNP

head(smad3_cis_pqtl)
write.csv(smad3_cis_pqtl,"/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/smad3_cis_pqtl_chorom15_500kb_logp5.csv")

exposure <- TwoSampleMR::format_data(smad3_cis_pqtl,type='exposure',snp_col = "SNP",
                                     beta_col = "BETA",
                                     se_col = "SE",
                                     eaf_col = "A1FREQ",
                                     effect_allele_col = "ALLELE1",
                                     other_allele_col = "ALLELE0",
                                     pval_col = "LOG10P",
                                     samplesize_col = "N",
                                     log_pval = T)
write.csv(exposure,"/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/exposure_smad3_cis_pqtl_chorom15_500kb_logp5.csv")
exposure_clump <- clump_data(exposure,clump_kb = 10000, clump_r2 = 0.1)
write.csv(exposure_clump,file = "/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/exposure_smad3_cis_pqtl_chorom15_500kb_logp5_clump_0.1_10000.csv")

if (T) {
  setwd("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/smad3_stroke_screen")
  # exposure <- read_csv(list[i])
  Traitinfo <-extract_outcome_data(snps = exposure_clump$SNP,outcomes = list_stroke_GWASID,proxies=T)# Extract outcome data
  mr_result <- getMR_Atlas_modify(exposure,Traitinfo)# 5min
  write.csv(mr_result,paste0("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/smad3_stroke_screen/","SMAD3_AS_fast.csv"))
}

#####MR-pqtl-SMAD3-ukb-ppp#####
library(MRPRESSO) 
setwd("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/241106/SMAD3-ukb-ppp")
exposure <- read.csv("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/exposure_smad3_cis_pqtl_chorom15_500kb_logp5_clump_0.1_10000.csv")
list_outcomes <- c("ebi-a-GCST006908")
for (i in 1:length(list_outcomes)){
  outcome_dat <- extract_outcome_data(snps = exposure$SNP, outcomes = (list_outcomes[i]), proxies = T)  # Extract outcome data
  write.csv(outcome_dat, paste0("outcome_", list_outcomes[i], ".csv"))
  dat <- harmonise_data(exposure, outcome_dat, action = 2)   # Align alleles, handle palindromic sequences
  write.csv(dat, paste0("dat_", list_outcomes[i], ".csv"))
  
  dat$EAF2 <- (1 - dat$eaf.exposure)
  dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)
  dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.outcome)
  dat$FSTAT <- ((dat$samplesize.outcome - 1 - 1) / 1) * (dat$PVE / (1 - dat$PVE))  # F-statistic
  write.csv(dat, paste0("datF_", list_outcomes[i], ".csv"))
  results <- mr(dat)  # MR calculation
  write.csv(results, paste0("Results_", list_outcomes[i], ".csv"))
  OR <- generate_odds_ratios(results)   # Odds Ratios  
  write.csv(OR, paste0("OR_", list_outcomes[i], ".csv"))
  
  # 1. Horizontal pleiotropy 
  # 2. Identify outlier SNPs 
  # 3. Compare results with and without outlier SNPs
  # mr_presso(BetaOutcome = 'beta.outcome',
  #           BetaExposure = 'beta.exposure', 
  #           SdOutcome = 'se.outcome', 
  #           SdExposure = 'se.exposure', 
  #           data = dat, OUTLIERtest = TRUE, 
  #           DISTORTIONtest = TRUE, SignifThreshold = 0.05, NbDistribution = 5000, seed = NULL)
  
  # Visualization
  library(TwoSampleMR)
  single <- mr_leaveoneout(dat)
  mr_leaveoneout_plot(single)
  ggsave(paste0("leaveoneout_", list_outcomes[i], ".pdf"), width = 5, height = 5)
  mr_scatter_plot(results, dat)
  ggsave(paste0("mr_scatter_", list_outcomes[i], ".pdf"), width = 4, height = 4)
  results_single <- mr_singlesnp(dat)
  mr_forest_plot(results_single)
  ggsave(paste0("forestPlot_", list_outcomes[i], ".pdf"), width = 5, height = 5)
  mr_funnel_plot(results_single)
  ggsave(paste0("funnelPlot_", list_outcomes[i], ".pdf"), width = 4, height = 4)
  het <- mr_heterogeneity(dat)
  write.csv(het, paste0("het_", list_outcomes[i], ".csv"))
  pleio <- mr_pleiotropy_test(dat)
  write.csv(pleio, paste0("pleio_", list_outcomes[i], ".csv"))
}

#####MR-pqtl-TREM2-ukb-ppp#####
library(MRPRESSO) 
setwd("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/241106/TREM2-ukb-ppp")
exposure <- read.csv("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/UKB-PPP-pqtl/ukb-ppp-pqtl-exposure/exposure_trem2_cis_pqtl_chorom6_500kb_logp5_clump_0.1_10000.csv")
list_outcomes <- c("ebi-a-GCST006908")
for (i in 1:length(list_outcomes)){
  outcome_dat <- extract_outcome_data(snps = exposure$SNP, outcomes = (list_outcomes[i]), proxies = T)  # Extract outcome data
  write.csv(outcome_dat, paste0("outcome_", list_outcomes[i], ".csv"))
  dat <- harmonise_data(exposure, outcome_dat, action = 3)   # Align alleles, handle palindromic sequences
  write.csv(dat, paste0("dat_", list_outcomes[i], ".csv"))
  
  dat$EAF2 <- (1 - dat$eaf.exposure)
  dat$MAF <- pmin(dat$eaf.exposure, dat$EAF2)
  dat$PVE <- mapply(PVEfx, dat$beta.exposure, dat$MAF, dat$se.exposure, N = dat$samplesize.outcome)
  dat$FSTAT <- ((dat$samplesize.outcome - 1 - 1) / 1) * (dat$PVE / (1 - dat$PVE))  # F-statistic
  write.csv(dat, paste0("datF_", list_outcomes[i], ".csv"))
  results <- mr(dat)  # MR calculation
  write.csv(results, paste0("Results_", list_outcomes[i], ".csv"))
  OR <- generate_odds_ratios(results)   # Odds Ratios  
  write.csv(OR, paste0("OR_", list_outcomes[i], ".csv"))
  
  # 1. Horizontal pleiotropy 
  # 2. Identify outlier SNPs 
  # 3. Compare results with and without outlier SNPs
  # mr_presso(BetaOutcome = 'beta.outcome',
  #           BetaExposure = 'beta.exposure', 
  #           SdOutcome = 'se.outcome', 
  #           SdExposure = 'se.exposure', 
  #           data = dat, OUTLIERtest = TRUE, 
  #           DISTORTIONtest = TRUE, SignifThreshold = 0.05, NbDistribution = 5000, seed = NULL)
  
  # Visualization
  library(TwoSampleMR)
  single <- mr_leaveoneout(dat)
  mr_leaveoneout_plot(single)
  ggsave(paste0("leaveoneout_", list_outcomes[i], ".pdf"), width = 5, height = 5)
  mr_scatter_plot(results, dat)
  ggsave(paste0("mr_scatter_", list_outcomes[i], ".pdf"), width = 4, height = 4)
  results_single <- mr_singlesnp(dat)
  mr_forest_plot(results_single)
  ggsave(paste0("forestPlot_", list_outcomes[i], ".pdf"), width = 5, height = 5)
  mr_funnel_plot(results_single)
  ggsave(paste0("funnelPlot_", list_outcomes[i], ".pdf"), width = 4, height = 4)
  het <- mr_heterogeneity(dat)
  write.csv(het, paste0("het_", list_outcomes[i], ".csv"))
  pleio <- mr_pleiotropy_test(dat)
  write.csv(pleio, paste0("pleio_", list_outcomes[i], ".csv"))
}

#####MR Forest Plot - V2 - 241107######
library(TwoSampleMR)
library(MRPRESSO)
library(tidyverse)
library(data.table)
library(ggplot2)

data <- read_csv("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/exposure/forest_data_stroke_241107.csv")
# Define color for each protein using hex codes
Protein_colors <- c('SMAD3' =   '#FA8072', 
                    'IL6R' = '#FFA500',  
                    'CD44' = '#008080', 
                    'TREM2' = '#91D1C2FF', 
                    'GAL-1' = '#000080')

# Create forest plot
ggplot_object <- ggplot(data, aes(y = Protein, x = Mean, xmin = Lower, xmax = Upper, color = Protein)) +
  geom_point() +
  geom_errorbarh(aes(height = 0.1)) +  # Draw error bars
  scale_color_manual(values = Protein_colors) +  # Assign colors
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_blank(),  # Remove background
    axis.line = element_line(colour = "black"),  # Add axis lines
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Add border
  ) +
  labs(x = "OR", y = "Proteins") +
  facet_wrap(~Method, scales = "free_x", ncol = 2) +  # 分面显示
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  # 在x=0处添加竖直虚线
  theme(strip.background = element_blank())  # 移除面标签的背景


print(ggplot_object)
ggsave("/Users/yuyao/Desktop/RESEARCH/Plaque proteomics_AHA_WCN/AHA-MR/fig5b_forest_MR_stroke_241107.pdf",
       height = 2.5, width = 7)