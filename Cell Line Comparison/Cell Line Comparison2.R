setwd("E:/Simon/HSP90DIA_may2024/Results/Cell Line Comparison")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(BSDA)
# 
# #NDNB 1 comparison
# #import processed data for each cell line
# cell_10A <- read.table("../10A_may20/WTvsNDNBProcessed2.tsv", sep = "\t", header = T)
# cell_231 <- read.table("../231_may23/WTvsNDNBProcessed2.tsv", sep = "\t", header = T)
# cell_468 <- read.table("../468_may23/WTvsNDNBProcessed2.tsv", sep = "\t", header = T)
# 
# 
# #set filtering value
# combined_cutoff_pVal <- 0.1
# #join data frames on proteins
# cell_10Avs231 <- full_join(cell_10A, cell_231, by = c("Accession", "proteinName", "gene"),
#                        suffix = c(".10A", ".231"))
# #filter by cutoff
# cell_10Avs231 <- cell_10Avs231[cell_10Avs231$adj.P.Val.10A < combined_cutoff_pVal |
#                                  cell_10Avs231$adj.P.Val.231 < combined_cutoff_pVal,]
# 
# #same as above for different comparisons
# cell_10Avs468 <- full_join(cell_10A, cell_468, by = c("Accession", "proteinName", "gene"),
#                            suffix = c(".10A", ".468"))
# 
# cell_10Avs468 <- cell_10Avs468[cell_10Avs468$adj.P.Val.10A < combined_cutoff_pVal |
#                                  cell_10Avs468$adj.P.Val.468 < combined_cutoff_pVal,]
# 
# 
# 
# cell_231vs468 <- full_join(cell_231, cell_468, by = c("Accession", "proteinName", "gene"),
#                            suffix = c(".231", ".468"))
# 
# cell_231vs468 <- cell_231vs468[cell_231vs468$adj.P.Val.231 < combined_cutoff_pVal |
#                                  cell_231vs468$adj.P.Val.468 < combined_cutoff_pVal,]
# 
# #cutoffs for coloring and labelling

# #plot
# ggplot(cell_10Avs231) +
#   geom_point(aes(x = logFC.10A, y = logFC.231,
#                  color = abs(logFC.10A - logFC.231) > colorCutoff),
#              size = 1.2) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   geom_abline(slope = 1) +
#   geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
#   geom_label_repel(data = cell_10Avs231[abs(cell_10Avs231$logFC.10A - cell_10Avs231$logFC.231) > labelCutoff,],
#                    aes(x = logFC.10A, y = logFC.231, label = gene)) +
#   labs(x = "Log Fold Change MCF-10A",
#        y = "Log Fold Change MDA-MB-231",
#        title = "NDNB1")
# #save
# ggsave("plots2/NDNB1-10A vs 231.png", width = 10, height = 10)
# 
# 
# 
# 
# 
# colorCutoff <- 2
# labelCutoff <- 2
# ggplot(cell_10Avs468) +
#   geom_point(aes(x = logFC.10A, y = logFC.468,
#                  color = abs(logFC.10A - logFC.468) > colorCutoff),
#              size = 1.2) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   geom_abline(slope = 1) +
#   geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
#   geom_label_repel(data = cell_10Avs468[abs(cell_10Avs468$logFC.10A - cell_10Avs468$logFC.468) > labelCutoff,],
#                    aes(x = logFC.10A, y = logFC.468, label = gene)) +
#   labs(x = "Log Fold Change MCF-10A",
#        y = "Log Fold Change MDA-MB-486",
#        title = "NDNB1")
# 
# ggsave("plots2/NDNB1-10A vs 486.png", width = 10, height = 10)
# 
# 
# 
# 
# colorCutoff <- 2
# labelCutoff <- 2.5
# ggplot(cell_231vs468) +
#   geom_point(aes(x = logFC.231, y = logFC.468,
#                  color = abs(logFC.231 - logFC.468) > colorCutoff),
#              size = 1.2) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   geom_abline(slope = 1) +
#   geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
#   geom_label_repel(data = cell_231vs468[abs(cell_231vs468$logFC.231 - cell_231vs468$logFC.468) > labelCutoff,],
#                    aes(x = logFC.231, y = logFC.468, label = gene)) +
#   labs(x = "Log Fold Change MDA-MB-231",
#        y = "Log Fold Change MDA-MB-486",
#        title = "NDNB1")
# 
# ggsave("plots2/NDNB1-231 vs 486.png", width = 10, height = 10)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #NDNB1182 comparison
# cell_10A <- read.table("../10A_may20/WTvs1182Processed2.tsv", sep = "\t", header = T)
# cell_231 <- read.table("../231_may23/WTvs1182Processed2.tsv", sep = "\t", header = T)
# cell_468 <- read.table("../468_may23/WTvs1182Processed2.tsv", sep = "\t", header = T)
# 
# 
# 
# combined_cutoff_pVal <- 0.1
# cell_10Avs231 <- full_join(cell_10A, cell_231, by = c("Accession", "proteinName", "gene"),
#                            suffix = c(".10A", ".231"))
# 
# cell_10Avs231 <- cell_10Avs231[cell_10Avs231$adj.P.Val.10A < combined_cutoff_pVal |
#                                  cell_10Avs231$adj.P.Val.231 < combined_cutoff_pVal,]
# 
# 
# cell_10Avs468 <- full_join(cell_10A, cell_468, by = c("Accession", "proteinName", "gene"),
#                            suffix = c(".10A", ".468"))
# 
# cell_10Avs468 <- cell_10Avs468[cell_10Avs468$adj.P.Val.10A < combined_cutoff_pVal |
#                                  cell_10Avs468$adj.P.Val.468 < combined_cutoff_pVal,]
# 
# 
# 
# cell_231vs468 <- full_join(cell_231, cell_468, by = c("Accession", "proteinName", "gene"),
#                            suffix = c(".231", ".468"))
# 
# cell_231vs468 <- cell_231vs468[cell_231vs468$adj.P.Val.231 < combined_cutoff_pVal |
#                                  cell_231vs468$adj.P.Val.468 < combined_cutoff_pVal,]
# 
# 
# colorCutoff <- 2
# labelCutoff <- 3
# ggplot(cell_10Avs231) +
#   geom_point(aes(x = logFC.10A, y = logFC.231,
#                  color = abs(logFC.10A - logFC.231) > colorCutoff),
#              size = 1.2) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   geom_abline(slope = 1) +
#   geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
#   geom_label_repel(data = cell_10Avs231[abs(cell_10Avs231$logFC.10A - cell_10Avs231$logFC.231) > labelCutoff,],
#                    aes(x = logFC.10A, y = logFC.231, label = gene)) +
#   labs(x = "Log Fold Change MCF-10A",
#        y = "Log Fold Change MDA-MB-231",
#        title = "NDNB1182")
# 
# ggsave("plots2/NDNB1182-10A vs 231.png", width = 10, height = 10)
# 
# 
# 
# 
# 
# colorCutoff <- 2
# labelCutoff <- 2
# ggplot(cell_10Avs468) +
#   geom_point(aes(x = logFC.10A, y = logFC.468,
#                  color = abs(logFC.10A - logFC.468) > colorCutoff),
#              size = 1.2) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   geom_abline(slope = 1) +
#   geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
#   geom_label_repel(data = cell_10Avs468[abs(cell_10Avs468$logFC.10A - cell_10Avs468$logFC.468) > labelCutoff,],
#                    aes(x = logFC.10A, y = logFC.468, label = gene)) +
#   labs(x = "Log Fold Change MCF-10A",
#        y = "Log Fold Change MDA-MB-486",
#        title = "NDNB1182")
# 
# ggsave("plots2/NDNB1182-10A vs 486.png", width = 10, height = 10)
# 
# 
# 
# 
# colorCutoff <- 2
# labelCutoff <- 2.5
# ggplot(cell_231vs468) +
#   geom_point(aes(x = logFC.231, y = logFC.468,
#                  color = abs(logFC.231 - logFC.468) > colorCutoff),
#              size = 1.2) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   geom_abline(slope = 1) +
#   geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
#   geom_label_repel(data = cell_231vs468[abs(cell_231vs468$logFC.231 - cell_231vs468$logFC.468) > labelCutoff,],
#                    aes(x = logFC.231, y = logFC.468, label = gene)) +
#   labs(x = "Log Fold Change MDA-MB-231",
#        y = "Log Fold Change MDA-MB-486",
#        title = "NDNB1182")
# 
# ggsave("plots2/NDNB1182-231 vs 486.png", width = 10, height = 10)
# 
# 

colorCutoff <- 2
labelCutoff <- 3

#statistical comparison
pvalueCutoff <- 0.05
FCdiffCutoffForTests <- 1

cell_10A <- read.table("../10A_may20/WTvsNDNBProcessed2.tsv", sep = "\t", header = T)
cell_231 <- read.table("../231_may23/WTvsNDNBProcessed2.tsv", sep = "\t", header = T)
cell_468 <- read.table("../468_may23/WTvsNDNBProcessed2.tsv", sep = "\t", header = T)


sampleSize <- 5
cell_10A$SD <- ((cell_10A$CI.R - cell_10A$logFC) * sqrt(sampleSize)) / 1.96

cell_231$SD <- ((cell_231$CI.R - cell_231$logFC) * sqrt(sampleSize)) / 1.96

cell_468$SD <- ((cell_468$CI.R - cell_468$logFC) * sqrt(sampleSize)) / 1.96

cell_10Avs231 <- full_join(cell_10A, cell_231, by = c("Accession", "proteinName", "gene"),
                           suffix = c(".10A", ".231"))

cell_10Avs468 <- full_join(cell_10A, cell_468, by = c("Accession", "proteinName", "gene"),
                           suffix = c(".10A", ".468"))


cell_231vs468 <- full_join(cell_231, cell_468, by = c("Accession", "proteinName", "gene"),
                           suffix = c(".231", ".468"))


cell_10Avs231$pvalofFCdiff <- NA
for (i in 1:nrow(cell_10Avs231)) {
  if (!is.na(cell_10Avs231$logFC.10A[i]) &
      !is.na(cell_10Avs231$logFC.231[i]) &
      !is.na(cell_10Avs231$SD.10A[i]) &
      !is.na(cell_10Avs231$SD.231[i])) {
    cell_10Avs231$pvalofFCdiff[i] <- tsum.test(mean.x=cell_10Avs231$logFC.10A[i], s.x=cell_10Avs231$SD.10A[i], n.x=sampleSize,
                                               mean.y=cell_10Avs231$logFC.231[i], s.y=cell_10Avs231$SD.231[i], n.y=sampleSize)[[3]]
  }
}

cell_10Avs231$FCdiff <- cell_10Avs231$logFC.231 - cell_10Avs231$logFC.10A



cell_10Avs231$pvalofFCdiffNEG10LOG <- -log(cell_10Avs231$pvalofFCdiff, base = 10)




#BHcorrection
pvalueCutoff <- 0.05
cell_10Avs231 <- cell_10Avs231[order(cell_10Avs231$pvalofFCdiff),]
validTests <- cell_10Avs231[!is.na(cell_10Avs231$pvalofFCdiff) &
                              abs(cell_10Avs231$FCdiff) > FCdiffCutoffForTests,]
validTestsN <- nrow(validTests)


cell_10Avs231$rank <- 1:nrow(cell_10Avs231)
cell_10Avs231$critValue <- NA
cell_10Avs231$BHSIG <- F
for(i in 1:validTestsN) {
  cell_10Avs231$critValue[i] <- (cell_10Avs231$rank[i] / validTestsN) * pvalueCutoff
  cell_10Avs231$BHSIG[i] <- cell_10Avs231$pvalofFCdiff[i] <= cell_10Avs231$critValue[i]
}


sigvals <- length(cell_10Avs231$pvalofFCdiff[cell_10Avs231$BHSIG])
BHcutoffPvalue <- max(cell_10Avs231$pvalofFCdiff[cell_10Avs231$BHSIG])


blackline <- -log(pvalueCutoff, base = 10)
redline <- -log(BHcutoffPvalue, base = 10)

ggplot(cell_10Avs231) +
  geom_point(aes(x = FCdiff, y = pvalofFCdiffNEG10LOG), color = "grey40") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Log Fold Change Difference (MDA-MB-231 / MCF-10A)",
       y = "Significance(-log10 p-value)",
       title = "NDNB1")+
  geom_point(data = cell_10Avs231[cell_10Avs231$pvalofFCdiff <= pvalueCutoff &
                                    abs(cell_10Avs231$FCdiff) > FCdiffCutoffForTests,],
             aes(x = FCdiff, y = pvalofFCdiffNEG10LOG, fill = FCdiff < 0),
             shape = 21, size = 3, stroke = 2) +
  geom_label_repel(data = cell_10Avs231[cell_10Avs231$pvalofFCdiff <= pvalueCutoff &
                                         abs(cell_10Avs231$FCdiff) > FCdiffCutoffForTests,],
                  aes(x = FCdiff, y = pvalofFCdiffNEG10LOG,
                      label = gene), size = 3, alpha = 0.7) +
  geom_vline(xintercept = c(FCdiffCutoffForTests, -FCdiffCutoffForTests), linetype = 2) +
  geom_hline(yintercept = blackline, linetype = 2) +
  geom_hline(yintercept = redline, linetype = 2, color = "red")


# +
#   geom_label(x = 0, y = blackline - 0.05, label = "p-value: 0.05", size = 3) +
#   geom_label(x = 0, y = redline + 0.05, label = paste0("B-H 0.05 (p-value: ", round(BHcutoffPvalue, 4), ")"), size = 3,
#              color = "red")
ggsave("plots2/NDNB1-10A vs 231 VOLCANO.png", width = 10, height = 10)


ggplot(cell_10Avs231) +
  geom_point(aes(x = logFC.10A, y = logFC.231,
                 color = abs(logFC.10A - logFC.231) > colorCutoff),
             size = 1.2) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
  geom_label_repel(data = cell_10Avs231[abs(cell_10Avs231$logFC.10A - cell_10Avs231$logFC.231) > labelCutoff,],
                   aes(x = logFC.10A, y = logFC.231, label = gene)) +
  labs(x = "Log Fold Change MCF-10A",
       y = "Log Fold Change MDA-MB-231",
       title = "NDNB1")
#save
ggsave("plots2/NDNB1-10A vs 231.png", width = 10, height = 10)










cell_10Avs468$pvalofFCdiff <- NA
for (i in 1:nrow(cell_10Avs468)) {
  if (!is.na(cell_10Avs468$logFC.10A[i]) &
      !is.na(cell_10Avs468$logFC.468[i]) &
      !is.na(cell_10Avs468$SD.10A[i]) &
      !is.na(cell_10Avs468$SD.468[i])) {
    cell_10Avs468$pvalofFCdiff[i] <- tsum.test(mean.x=cell_10Avs468$logFC.10A[i], s.x=cell_10Avs468$SD.10A[i], n.x=sampleSize,
                                               mean.y=cell_10Avs468$logFC.468[i], s.y=cell_10Avs468$SD.468[i], n.y=sampleSize)[[3]]
  }
}

cell_10Avs468$FCdiff <- cell_10Avs468$logFC.468 - cell_10Avs468$logFC.10A



cell_10Avs468$pvalofFCdiffNEG10LOG <- -log(cell_10Avs468$pvalofFCdiff, base = 10)

#BHcorrection
pvalueCutoff <- 0.05
cell_10Avs468 <- cell_10Avs468[order(cell_10Avs468$pvalofFCdiff),]
validTests <- cell_10Avs468[!is.na(cell_10Avs468$pvalofFCdiff) &
                              abs(cell_10Avs468$FCdiff) > FCdiffCutoffForTests,]
validTestsN <- nrow(validTests)


cell_10Avs468$rank <- 1:nrow(cell_10Avs468)
cell_10Avs468$critValue <- NA
cell_10Avs468$BHSIG <- F
for(i in 1:validTestsN) {
  cell_10Avs468$critValue[i] <- (cell_10Avs468$rank[i] / validTestsN) * pvalueCutoff
  cell_10Avs468$BHSIG[i] <- cell_10Avs468$pvalofFCdiff[i] <= cell_10Avs468$critValue[i]
}


sigvals <- length(cell_10Avs468$pvalofFCdiff[cell_10Avs468$BHSIG])
BHcutoffPvalue <- max(cell_10Avs468$pvalofFCdiff[cell_10Avs468$BHSIG])


ggplot(cell_10Avs468) +
  geom_point(aes(x = FCdiff, y = pvalofFCdiffNEG10LOG), color = "grey40") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Log Fold Change Difference (MDA-MB-468 / MCF-10A)",
       y = "Significance(-log10 p-value)",
       title = "NDNB1") +
  geom_point(data = cell_10Avs468[cell_10Avs468$pvalofFCdiff <= pvalueCutoff &
                                    abs(cell_10Avs468$FCdiff) > FCdiffCutoffForTests,],
             aes(x = FCdiff, y = pvalofFCdiffNEG10LOG, fill = FCdiff < 0),
             shape = 21, size = 3, stroke = 2) +
  geom_label_repel(data = cell_10Avs468[cell_10Avs468$pvalofFCdiff <= pvalueCutoff &
                                          abs(cell_10Avs468$FCdiff) > FCdiffCutoffForTests,],
                   aes(x = FCdiff, y = pvalofFCdiffNEG10LOG,
                       label = gene), size = 3, alpha = 0.7) +
  geom_vline(xintercept = c(FCdiffCutoffForTests, -FCdiffCutoffForTests), linetype = 2) +
  geom_hline(yintercept = -log(pvalueCutoff, base = 10), linetype = 2) +
  geom_hline(yintercept = -log(BHcutoffPvalue, base = 10), linetype = 2, color = "red")

# +
#   geom_label(x = 0, y = -log(pvalueCutoff, base = 10) - 0.05, label = "p-value: 0.05", size = 3) +
#   geom_label(x = 0, y = -log(BHcutoffPvalue, base = 10) + 0.05, label = paste0("B-H 0.05 (p-value: ", round(BHcutoffPvalue, 4), ")"), size = 3,
#              color = "red")
ggsave("plots2/NDNB1-10A vs 468 VOLCANO.png", width = 10, height = 10)


ggplot(cell_10Avs468) +
  geom_point(aes(x = logFC.10A, y = logFC.468,
                 color = abs(logFC.10A - logFC.468) > colorCutoff),
             size = 1.2) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
  geom_label_repel(data = cell_10Avs468[abs(cell_10Avs468$logFC.10A - cell_10Avs468$logFC.468) > labelCutoff,],
                   aes(x = logFC.10A, y = logFC.468, label = gene)) +
  labs(x = "Log Fold Change MCF-10A",
       y = "Log Fold Change MDA-MB-468",
       title = "NDNB1")
#save
ggsave("plots2/NDNB1-10A vs 468.png", width = 10, height = 10)
  






cell_231vs468$pvalofFCdiff <- NA
for (i in 1:nrow(cell_231vs468)) {
  if (!is.na(cell_231vs468$logFC.231[i]) &
      !is.na(cell_231vs468$logFC.468[i]) &
      !is.na(cell_231vs468$SD.231[i]) &
      !is.na(cell_231vs468$SD.468[i])) {
    cell_231vs468$pvalofFCdiff[i] <- tsum.test(mean.x=cell_231vs468$logFC.231[i], s.x=cell_231vs468$SD.231[i], n.x=sampleSize,
                                               mean.y=cell_231vs468$logFC.468[i], s.y=cell_231vs468$SD.468[i], n.y=sampleSize)[[3]]
  }
}

cell_231vs468$FCdiff <- cell_231vs468$logFC.468 - cell_231vs468$logFC.231



cell_231vs468$pvalofFCdiffNEG10LOG <- -log(cell_231vs468$pvalofFCdiff, base = 10)



#BHcorrection
pvalueCutoff <- 0.05
cell_231vs468 <- cell_231vs468[order(cell_231vs468$pvalofFCdiff),]
validTests <- cell_231vs468[!is.na(cell_231vs468$pvalofFCdiff) &
                              abs(cell_231vs468$FCdiff) > FCdiffCutoffForTests,]
validTestsN <- nrow(validTests)


cell_231vs468$rank <- 1:nrow(cell_231vs468)
cell_231vs468$critValue <- NA
cell_231vs468$BHSIG <- F
for(i in 1:validTestsN) {
  cell_231vs468$critValue[i] <- (cell_231vs468$rank[i] / validTestsN) * pvalueCutoff
  cell_231vs468$BHSIG[i] <- cell_231vs468$pvalofFCdiff[i] <= cell_231vs468$critValue[i]
}

sigvals <- length(cell_231vs468$pvalofFCdiff[cell_231vs468$BHSIG])
BHcutoffPvalue <- max(cell_231vs468$pvalofFCdiff[cell_231vs468$BHSIG])


ggplot(cell_231vs468) +
  geom_point(aes(x = FCdiff, y = pvalofFCdiffNEG10LOG), color = "grey40") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Log Fold Change Difference (MDA-MB-468 / MDA-MB-231)",
       y = "Significance(-log10 p-value)",
       title = "NDNB1") +
  geom_point(data = cell_231vs468[cell_231vs468$pvalofFCdiff <= pvalueCutoff &
                                    abs(cell_231vs468$FCdiff) > FCdiffCutoffForTests,],
             aes(x = FCdiff, y = pvalofFCdiffNEG10LOG, fill = FCdiff < 0),
             shape = 21, size = 3, stroke = 2) +
  geom_label_repel(data = cell_231vs468[cell_231vs468$pvalofFCdiff <= pvalueCutoff &
                                          abs(cell_231vs468$FCdiff) > FCdiffCutoffForTests,],
                   aes(x = FCdiff, y = pvalofFCdiffNEG10LOG,
                       label = gene), size = 3, alpha = 0.7) +
  geom_vline(xintercept = c(FCdiffCutoffForTests, -FCdiffCutoffForTests), linetype = 2) +
  geom_hline(yintercept = -log(pvalueCutoff, base = 10), linetype = 2) +
  geom_hline(yintercept = -log(BHcutoffPvalue, base = 10), linetype = 2, color = "red")

# +
#   geom_label(x = 0, y = -log(pvalueCutoff, base = 10) - 0.05, label = "p-value: 0.05", size = 3) +
#   geom_label(x = 0, y = -log(BHcutoffPvalue, base = 10) + 0.05, label = paste0("B-H 0.05 (p-value: ", round(BHcutoffPvalue, 4), ")"), size = 3,
#              color = "red")
ggsave("plots2/NDNB1-231 vs 468 VOLCANO.png", width = 10, height = 10)


ggplot(cell_231vs468) +
  geom_point(aes(x = logFC.231, y = logFC.468,
                 color = abs(logFC.231 - logFC.468) > colorCutoff),
             size = 1.2) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
  geom_label_repel(data = cell_231vs468[abs(cell_231vs468$logFC.231 - cell_231vs468$logFC.468) > labelCutoff,],
                   aes(x = logFC.231, y = logFC.468, label = gene)) +
  labs(x = "Log Fold Change MDA-MB-231",
       y = "Log Fold Change MDA-MB-468",
       title = "NDNB1")
#save
ggsave("plots2/NDNB1-231 vs 468.png", width = 10, height = 10)



#NDNB1182 comparison

pvalueCutoff <- 0.05
FCdiffCutoffForTests <- 1
cell_10A <- read.table("../10A_may20/WTvs1182Processed2.tsv", sep = "\t", header = T)
cell_231 <- read.table("../231_may23/WTvs1182Processed2.tsv", sep = "\t", header = T)
cell_468 <- read.table("../468_may23/WTvs1182Processed2.tsv", sep = "\t", header = T)



sampleSize <- 5
cell_10A$SD <- ((cell_10A$CI.R - cell_10A$logFC) * sqrt(sampleSize)) / 1.96

cell_231$SD <- ((cell_231$CI.R - cell_231$logFC) * sqrt(sampleSize)) / 1.96

cell_468$SD <- ((cell_468$CI.R - cell_468$logFC) * sqrt(sampleSize)) / 1.96



cell_10Avs231 <- full_join(cell_10A, cell_231, by = c("Accession", "proteinName", "gene"),
                           suffix = c(".10A", ".231"))

cell_10Avs468 <- full_join(cell_10A, cell_468, by = c("Accession", "proteinName", "gene"),
                           suffix = c(".10A", ".468"))


cell_231vs468 <- full_join(cell_231, cell_468, by = c("Accession", "proteinName", "gene"),
                           suffix = c(".231", ".468"))


cell_10Avs231$pvalofFCdiff <- NA
for (i in 1:nrow(cell_10Avs231)) {
  if (!is.na(cell_10Avs231$logFC.10A[i]) &
      !is.na(cell_10Avs231$logFC.231[i]) &
      !is.na(cell_10Avs231$SD.10A[i]) &
      !is.na(cell_10Avs231$SD.231[i])) {
    cell_10Avs231$pvalofFCdiff[i] <- tsum.test(mean.x=cell_10Avs231$logFC.10A[i], s.x=cell_10Avs231$SD.10A[i], n.x=sampleSize,
                                               mean.y=cell_10Avs231$logFC.231[i], s.y=cell_10Avs231$SD.231[i], n.y=sampleSize)[[3]]
  }
}

cell_10Avs231$FCdiff <- cell_10Avs231$logFC.231 - cell_10Avs231$logFC.10A



cell_10Avs231$pvalofFCdiffNEG10LOG <- -log(cell_10Avs231$pvalofFCdiff, base = 10)




#BHcorrection
pvalueCutoff <- 0.05
cell_10Avs231 <- cell_10Avs231[order(cell_10Avs231$pvalofFCdiff),]
validTests <- cell_10Avs231[!is.na(cell_10Avs231$pvalofFCdiff) &
                              abs(cell_10Avs231$FCdiff) > FCdiffCutoffForTests,]
validTestsN <- nrow(validTests)


cell_10Avs231$rank <- 1:nrow(cell_10Avs231)
cell_10Avs231$critValue <- NA
cell_10Avs231$BHSIG <- F
for(i in 1:validTestsN) {
  cell_10Avs231$critValue[i] <- (cell_10Avs231$rank[i] / validTestsN) * pvalueCutoff
  cell_10Avs231$BHSIG[i] <- cell_10Avs231$pvalofFCdiff[i] <= cell_10Avs231$critValue[i]
}


sigvals <- length(cell_10Avs231$pvalofFCdiff[cell_10Avs231$BHSIG])
BHcutoffPvalue <- max(cell_10Avs231$pvalofFCdiff[cell_10Avs231$BHSIG])



ggplot(cell_10Avs231) +
  geom_point(aes(x = FCdiff, y = pvalofFCdiffNEG10LOG), color = "grey40") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Log Fold Change Difference (MDA-MB-231 / MCF-10A)",
       y = "Significance(-log10 p-value)",
       title = "NDNB1182") +
  geom_point(data = cell_10Avs231[cell_10Avs231$pvalofFCdiff <= pvalueCutoff &
                                    abs(cell_10Avs231$FCdiff) > FCdiffCutoffForTests,],
             aes(x = FCdiff, y = pvalofFCdiffNEG10LOG, fill = FCdiff < 0),
             shape = 21, size = 3, stroke = 2) +
  geom_label_repel(data = cell_10Avs231[cell_10Avs231$pvalofFCdiff <= pvalueCutoff &
                                          abs(cell_10Avs231$FCdiff) > FCdiffCutoffForTests,],
                   aes(x = FCdiff, y = pvalofFCdiffNEG10LOG,
                       label = gene), size = 3, alpha = 0.7) +
  geom_vline(xintercept = c(FCdiffCutoffForTests, -FCdiffCutoffForTests), linetype = 2) +
  geom_hline(yintercept = -log(pvalueCutoff, base = 10), linetype = 2) +
  geom_hline(yintercept = -log(BHcutoffPvalue, base = 10), linetype = 2, color = "red")


# +
#   geom_label(x = 0, y = -log(pvalueCutoff, base = 10) - 0.05, label = "p-value: 0.05", size = 3) +
#   geom_label(x = 0, y = -log(BHcutoffPvalue, base = 10) + 0.05, label = paste0("B-H 0.05 (p-value: ", round(BHcutoffPvalue, 4), ")"), size = 3,
#              color = "red")
ggsave("plots2/NDNB1182-10A vs 231 VOLCANO.png", width = 10, height = 10)


ggplot(cell_10Avs231) +
  geom_point(aes(x = logFC.10A, y = logFC.231,
                 color = abs(logFC.10A - logFC.231) > colorCutoff),
             size = 1.2) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
  geom_label_repel(data = cell_10Avs231[abs(cell_10Avs231$logFC.10A - cell_10Avs231$logFC.231) > labelCutoff,],
                   aes(x = logFC.10A, y = logFC.231, label = gene)) +
  labs(x = "Log Fold Change MCF-10A",
       y = "Log Fold Change MDA-MB-231",
       title = "NDNB1182")
#save
ggsave("plots2/NDNB1182-10A vs 231.png", width = 10, height = 10)










cell_10Avs468$pvalofFCdiff <- NA
for (i in 1:nrow(cell_10Avs468)) {
  if (!is.na(cell_10Avs468$logFC.10A[i]) &
      !is.na(cell_10Avs468$logFC.468[i]) &
      !is.na(cell_10Avs468$SD.10A[i]) &
      !is.na(cell_10Avs468$SD.468[i])) {
    cell_10Avs468$pvalofFCdiff[i] <- tsum.test(mean.x=cell_10Avs468$logFC.10A[i], s.x=cell_10Avs468$SD.10A[i], n.x=sampleSize,
                                               mean.y=cell_10Avs468$logFC.468[i], s.y=cell_10Avs468$SD.468[i], n.y=sampleSize)[[3]]
  }
}

cell_10Avs468$FCdiff <- cell_10Avs468$logFC.468 - cell_10Avs468$logFC.10A



cell_10Avs468$pvalofFCdiffNEG10LOG <- -log(cell_10Avs468$pvalofFCdiff, base = 10)

#BHcorrection
pvalueCutoff <- 0.05
cell_10Avs468 <- cell_10Avs468[order(cell_10Avs468$pvalofFCdiff),]
validTests <- cell_10Avs468[!is.na(cell_10Avs468$pvalofFCdiff) &
                              abs(cell_10Avs468$FCdiff) > FCdiffCutoffForTests,]
validTestsN <- nrow(validTests)


cell_10Avs468$rank <- 1:nrow(cell_10Avs468)
cell_10Avs468$critValue <- NA
cell_10Avs468$BHSIG <- F
for(i in 1:validTestsN) {
  cell_10Avs468$critValue[i] <- (cell_10Avs468$rank[i] / validTestsN) * pvalueCutoff
  cell_10Avs468$BHSIG[i] <- cell_10Avs468$pvalofFCdiff[i] <= cell_10Avs468$critValue[i]
}


sigvals <- length(cell_10Avs468$pvalofFCdiff[cell_10Avs468$BHSIG])
BHcutoffPvalue <- max(cell_10Avs468$pvalofFCdiff[cell_10Avs468$BHSIG])


ggplot(cell_10Avs468) +
  geom_point(aes(x = FCdiff, y = pvalofFCdiffNEG10LOG), color = "grey40") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Log Fold Change Difference (MDA-MB-468 / MCF-10A)",
       y = "Significance(-log10 p-value)",
       title = "NDNB1182") +
  geom_point(data = cell_10Avs468[cell_10Avs468$pvalofFCdiff <= pvalueCutoff &
                                    abs(cell_10Avs468$FCdiff) > FCdiffCutoffForTests,],
             aes(x = FCdiff, y = pvalofFCdiffNEG10LOG, fill = FCdiff < 0),
             shape = 21, size = 3, stroke = 2) +
  geom_label_repel(data = cell_10Avs468[cell_10Avs468$pvalofFCdiff <= pvalueCutoff &
                                          abs(cell_10Avs468$FCdiff) > FCdiffCutoffForTests,],
                   aes(x = FCdiff, y = pvalofFCdiffNEG10LOG,
                       label = gene), size = 3, alpha = 0.7) +
  geom_vline(xintercept = c(FCdiffCutoffForTests, -FCdiffCutoffForTests), linetype = 2) +
  geom_hline(yintercept = -log(pvalueCutoff, base = 10), linetype = 2) +
  geom_hline(yintercept = -log(BHcutoffPvalue, base = 10), linetype = 2, color = "red")

# 
# +
#   geom_label(x = 0, y = -log(pvalueCutoff, base = 10) - 0.05, label = "p-value: 0.05", size = 3) +
#   geom_label(x = 0, y = -log(BHcutoffPvalue, base = 10) + 0.05, label = paste0("B-H 0.05 (p-value: ", round(BHcutoffPvalue, 4), ")"), size = 3,
#              color = "red")
ggsave("plots2/NDNB1182-10A vs 468 VOLCANO.png", width = 10, height = 10)


ggplot(cell_10Avs468) +
  geom_point(aes(x = logFC.10A, y = logFC.468,
                 color = abs(logFC.10A - logFC.468) > colorCutoff),
             size = 1.2) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
  geom_label_repel(data = cell_10Avs468[abs(cell_10Avs468$logFC.10A - cell_10Avs468$logFC.468) > labelCutoff,],
                   aes(x = logFC.10A, y = logFC.468, label = gene)) +
  labs(x = "Log Fold Change MCF-10A",
       y = "Log Fold Change MDA-MB-468",
       title = "NDNB1182")
#save
ggsave("plots2/NDNB1182-10A vs 468.png", width = 10, height = 10)







cell_231vs468$pvalofFCdiff <- NA
for (i in 1:nrow(cell_231vs468)) {
  if (!is.na(cell_231vs468$logFC.231[i]) &
      !is.na(cell_231vs468$logFC.468[i]) &
      !is.na(cell_231vs468$SD.231[i]) &
      !is.na(cell_231vs468$SD.468[i])) {
    cell_231vs468$pvalofFCdiff[i] <- tsum.test(mean.x=cell_231vs468$logFC.231[i], s.x=cell_231vs468$SD.231[i], n.x=sampleSize,
                                               mean.y=cell_231vs468$logFC.468[i], s.y=cell_231vs468$SD.468[i], n.y=sampleSize)[[3]]
  }
}

cell_231vs468$FCdiff <- cell_231vs468$logFC.468 - cell_231vs468$logFC.231



cell_231vs468$pvalofFCdiffNEG10LOG <- -log(cell_231vs468$pvalofFCdiff, base = 10)



#BHcorrection
pvalueCutoff <- 0.05
cell_231vs468 <- cell_231vs468[order(cell_231vs468$pvalofFCdiff),]
validTests <- cell_231vs468[!is.na(cell_231vs468$pvalofFCdiff) &
                              abs(cell_231vs468$FCdiff) > FCdiffCutoffForTests,]
validTestsN <- nrow(validTests)


cell_231vs468$rank <- 1:nrow(cell_231vs468)
cell_231vs468$critValue <- NA
cell_231vs468$BHSIG <- F
for(i in 1:validTestsN) {
  cell_231vs468$critValue[i] <- (cell_231vs468$rank[i] / validTestsN) * pvalueCutoff
  cell_231vs468$BHSIG[i] <- cell_231vs468$pvalofFCdiff[i] <= cell_231vs468$critValue[i]
}


sigvals <- length(cell_231vs468$pvalofFCdiff[cell_231vs468$BHSIG])
BHcutoffPvalue <- max(cell_231vs468$pvalofFCdiff[cell_231vs468$BHSIG])


ggplot(cell_231vs468) +
  geom_point(aes(x = FCdiff, y = pvalofFCdiffNEG10LOG), color = "grey40") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Log Fold Change Difference (MDA-MB-468 / MDA-MB-231)",
       y = "Significance(-log10 p-value)",
       title = "NDNB1182") +
  geom_point(data = cell_231vs468[cell_231vs468$pvalofFCdiff <= pvalueCutoff &
                                    abs(cell_231vs468$FCdiff) > FCdiffCutoffForTests,],
             aes(x = FCdiff, y = pvalofFCdiffNEG10LOG, fill = FCdiff < 0),
             shape = 21, size = 3, stroke = 2) +
  geom_label_repel(data = cell_231vs468[cell_231vs468$pvalofFCdiff <= pvalueCutoff &
                                          abs(cell_231vs468$FCdiff) > FCdiffCutoffForTests,],
                   aes(x = FCdiff, y = pvalofFCdiffNEG10LOG,
                       label = gene), size = 3, alpha = 0.7) +
  geom_vline(xintercept = c(FCdiffCutoffForTests, -FCdiffCutoffForTests), linetype = 2) +
  geom_hline(yintercept = -log(pvalueCutoff, base = 10), linetype = 2) +
  geom_hline(yintercept = -log(BHcutoffPvalue, base = 10), linetype = 2, color = "red")
# 
# +
#   geom_label(x = 0, y = -log(pvalueCutoff, base = 10) - 0.05, label = "p-value: 0.05", size = 3) +
#   geom_label(x = 0, y = -log(BHcutoffPvalue, base = 10) + 0.05, label = paste0("B-H 0.05 (p-value: ", round(BHcutoffPvalue, 4), ")"), size = 3,
#              color = "red")
ggsave("plots2/NDNB1182-231 vs 468 VOLCANO.png", width = 10, height = 10)


ggplot(cell_231vs468) +
  geom_point(aes(x = logFC.231, y = logFC.468,
                 color = abs(logFC.231 - logFC.468) > colorCutoff),
             size = 1.2) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
  geom_label_repel(data = cell_231vs468[abs(cell_231vs468$logFC.231 - cell_231vs468$logFC.468) > labelCutoff,],
                   aes(x = logFC.231, y = logFC.468, label = gene)) +
  labs(x = "Log Fold Change MDA-MB-231",
       y = "Log Fold Change MDA-MB-468",
       title = "NDNB1182")
#save
ggsave("plots2/NDNB1182-231 vs 468.png", width = 10, height = 10)


