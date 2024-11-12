#set working directory
setwd("E:/Simon/HSP90DIA_may2024/Results/10A_may20/")

#import libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggrepel)
library(QFeatures)
library(corrplot)
library(limma)
library(tibble)

#import protein results
proteins <- read.table("20240523_141927_10A_may20_Report-proteins-pivot-ND-withprecursorIDcount.tsv",
                       sep = "\t",
                       header = T)

proteinIDcolN <- 3 #first three columns contain protein ID info instead of quant info


#simplify names of injections
names(proteins)[(proteinIDcolN + 1):ncol(proteins)] <- str_remove(str_remove((str_extract(names(proteins)[(proteinIDcolN + 1):ncol(proteins)],
                                                                                          "[^_]*_BR[0-9]+.*[0-9]{4}.*")), "_Slot[0-9]+.[0-9]+"), "[0-9]+.d.PG.")
names(proteins)
precursorsMeasuredDF <- proteins[,which(grepl("NrOfPrecursors", names(proteins)) | grepl("^PG.", names(proteins)))]
proteins <- proteins[,which(grepl("Quantity", names(proteins)) | grepl("^PG.", names(proteins)))]
proteins[precursorsMeasuredDF <= 1] <- NaN


names(proteins) <- str_remove_all(names(proteins), "_Quantity")



#dictionaries for keeping track of linked genes, protein names, and accessions
geneDict <- as.list(proteins$PG.Genes)
names(geneDict) <- proteins$PG.ProteinGroups
protNameDict <- as.list(proteins$PG.ProteinNames)
names(protNameDict) <- proteins$PG.ProteinGroups

#add proteins as Qprot object
Qprot <- readQFeatures(table = proteins,
                       ecol = which(grepl("_BR", names(proteins))),
                       fnames = "PG.ProteinGroups",
                       name = "raw_proteins")

#print to check and add metadata to columns
colnames(Qprot[[1]])
Qprot$bioRep <- str_extract(colnames(Qprot[[1]]), "BR[0-9]+")
Qprot$techRep <- str_extract(colnames(Qprot[[1]]), "[0-9]+$")
Qprot$condition <- str_extract(colnames(Qprot[[1]]), "^[^_]*")
Qprot$Inj <- colnames(Qprot[[1]])
#print to check
colData(Qprot)

#log transform, base 2 so that later FC is in base 2
Qprot <- addAssay(Qprot,
                  logTransform((Qprot[["raw_proteins"]]),
                               base = 2),
                  name = "log_proteins")

#normalize
Qprot <- addAssay(Qprot,
                  normalize(Qprot[["log_proteins"]],
                            method = "center.median"),
                  name = "norm_proteins")

#save figure of log and normaliztions
png("plots2/LogNormalization.png", units = "px", width = 1500, height = 1000)
par(mfrow = c(1,3))
limma::plotDensities(assay(Qprot[[1]]), legend = F)
limma::plotDensities(assay(Qprot[[2]]), legend = F)
limma::plotDensities(assay(Qprot[[3]]), legend = F)
dev.off()


#add in old column data to normalized data
colData(Qprot[["norm_proteins"]]) <- colData(Qprot)

#create a PCA plot
prot_pca <- Qprot[["norm_proteins"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp(scale = TRUE, center = TRUE)
pca_df <- as.data.frame(prot_pca$x)
pca_df$condition <- Qprot$condition
pca_df$br <- Qprot$bioRep

pca_df %>%
  ggplot(aes(x = PC1, y = PC2, colour = condition, shape = br)) +
  geom_point(size = 4) +
  labs(colour = element_blank(),shape = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "PC1", y = "PC2") +
  ggtitle("Protein-level PCA plot")+
  theme_bw()

ggsave("plots2/proteinPCA.png", width = 8, height = 5)



#pull out normalized proteins as a dataframe to plot correlation
prot_df <- Qprot[["norm_proteins"]] %>% assay() %>% as.data.frame()
#create correlation matrix
corr_matrixP <- cor(prot_df,
                    method = "pearson",
                    use = "pairwise.complete.obs")
layout(matrix(c(1)), 1)
png("plots2/corrplot.png", width = 1000, height = 1000)
corrplot(corr_matrixP,
         method = "color",
         tl.col = "black",
         number.font = 1,
         tl.cex = 0.65,
         type = "upper")
dev.off()





#differential analysis

#pull out normalized values into object. not necessary, but make syntax easier
brSummary <- Qprot[["norm_proteins"]] %>% assay() %>% data.frame() %>% rownames_to_column("Accession")
#pivot longer, one obs per row
brSummaryL <- brSummary %>% pivot_longer(cols = which(!grepl("Accession", names(brSummary))))
brSummaryL$name <- str_replace_all(brSummaryL$name, "X1182", "1182") #bug fix for column starting with a number
#meta data
brSummaryL$bioRep <- str_extract(brSummaryL$name, "BR[0-9]+")
brSummaryL$techRep <- str_extract(brSummaryL$name, "[0-9]+$")
brSummaryL$condition <- str_extract(brSummaryL$name, "^[^_]*")

#group bioreps together and average
brSummaryComb <- brSummaryL %>% group_by(bioRep, condition, Accession) %>%
  summarise(Area = mean(value, na.rm = T))


quad_vectorized <- function(y_vec, a, b, c1) {
  answers <- sapply(y_vec, function(y) {
    c <- c1 - y
    answer <- (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
  })
  return(answers)
}

NDNB1norm <- brSummaryComb[brSummaryComb$condition == "NDNB" |
                             brSummaryComb$condition == "DMSO",] %>% group_by(condition, Accession) %>%
  summarise(AreaMean = mean(Area, na.rm = T))

NDNB1normWide <- NDNB1norm %>% pivot_wider(names_from = condition, values_from = AreaMean)

NDNB1normWide <- NDNB1normWide[!is.na(NDNB1normWide$DMSO) & !is.na(NDNB1normWide$NDNB),]

fitQuad <- lm(NDNB1normWide$NDNB ~ poly(NDNB1normWide$DMSO, 2, raw = T))
coQuat <- coef(fitQuad)
a <- coQuat[[3]]
b <- coQuat[[2]]
c <- coQuat[[1]]

quadDF <- data.frame(x = seq(-10,10, 0.1))
quadDF$y <- ((quadDF$x)^2 * a) + (quadDF$x) * b + c

ggplot(NDNB1normWide) +
  geom_point(aes(x = DMSO, y = NDNB)) +
  geom_abline(slope = 1, color = "red", linetype = 2, size = 1) +
  geom_line(data = quadDF, aes(x,y), size = 1 ,color = "green2", linetype = 2) +
  theme_bw()

ggsave("plots2/preQuadNormNDNB1.png", width = 8, height = 8)

ggplot(NDNB1normWide) +
  geom_point(aes(x = DMSO, y = quad_vectorized(NDNB, a, b, c))) +
  geom_abline(slope = 1, color = "red", linetype = 2, size = 1) +
  theme_bw()


ggsave("plots2/postQuadNormNDNB1.png", width = 8, height = 8)

brSummaryComb$Area[brSummaryComb$condition == 'NDNB'] <- quad_vectorized(brSummaryComb$Area[brSummaryComb$condition == 'NDNB'],
                                                                         a,b,c)

NDNB1182norm <- brSummaryComb[brSummaryComb$condition == "1182" |
                                brSummaryComb$condition == "DMSO",] %>% group_by(condition, Accession) %>%
  summarise(AreaMean = mean(Area, na.rm = T))

NDNB1182normWide <- NDNB1182norm %>% pivot_wider(names_from = condition, values_from = AreaMean)

NDNB1182normWide <- NDNB1182normWide[!is.na(NDNB1182normWide$DMSO) & !is.na(NDNB1182normWide[["1182"]]),]

fitQuad <- lm(NDNB1182normWide[["1182"]] ~ poly(NDNB1182normWide$DMSO, 2, raw = T))
coQuat <- coef(fitQuad)
a <- coQuat[[3]]
b <- coQuat[[2]]
c <- coQuat[[1]]

quadDF <- data.frame(x = seq(-10,10, 0.1))
quadDF$y <- ((quadDF$x)^2 * a) + (quadDF$x) * b + c

ggplot(NDNB1182normWide) +
  geom_point(aes(x = DMSO, y = `1182`)) +
  geom_abline(slope = 1, color = "red", linetype = 2, size = 1) +
  geom_line(data = quadDF, aes(x,y), size = 1 ,color = "green2", linetype = 2) +
  theme_bw()

ggsave("plots2/preQuadNormNDNB1182.png", width = 8, height = 8)

ggplot(NDNB1182normWide) +
  geom_point(aes(x = DMSO, y = `1182` %>% quad_vectorized(a, b, c))) +
  geom_abline(slope = 1, color = "red", linetype = 2, size = 1) +
  theme_bw()

ggsave("plots2/postQuadNormNDNB1182.png", width = 8, height = 8)


brSummaryComb$Area[brSummaryComb$condition == '1182'] <- quad_vectorized(brSummaryComb$Area[brSummaryComb$condition == '1182'],
                                                                         a,b,c)


#pivot back to wide format for qfeat
brSummaryWide <- brSummaryComb %>% pivot_wider(names_from = c(condition, bioRep),
                                               values_from = Area,
                                               names_sep = "_",
                                               names_sort = T)

#reimport to qfeatures
QprotBR <- readQFeatures(table = brSummaryWide,
                         ecol = 2:ncol(brSummaryWide),
                         fnames = "Accession",
                         name = "norm_proteins_comb")
colnames(QprotBR[[1]])

#metadata
QprotBR$bioRep <- str_extract(colnames(QprotBR[[1]]), "BR[0-9]+")
QprotBR$condition <- str_extract(colnames(QprotBR[[1]]), "^[^_]*")
QprotBR$ID <- colnames(QprotBR[[1]])

colData(QprotBR)


#pull out as object
proteinsBR <- QprotBR[["norm_proteins_comb"]]
#set variables as factor, and clarify the reference as DMSO
proteinsBR$condition <- factor(QprotBR$condition)
proteinsBR$condition <- relevel(proteinsBR$condition, ref = "DMSO")

#create a linear matrix model using model.matrix to assign which columns go with which coefficients
model_design <- model.matrix(~proteinsBR$condition)

#fit proteins to linear model with lmFit and the matrix we just created
fitted_lm <- proteinsBR %>%
  assay() %>%
  lmFit(design = model_design) %>%
  eBayes() 

#makes the result easy to read in a table, sorted by most significant genes
#select each coefficient to see it related to DMSO
#adjust by benjamini hochburg for adj.P.values for multiple hypothesis testing
DE_1182 <- topTable(fit = fitted_lm,
                    adjust.method = "BH",
                    number = Inf,
                    coef  = "proteinsBR$condition1182",
                    sort.by = "p", confint = T)
DE_1182$Accession <- rownames(DE_1182)

#calculate corrected z score based on B-H adujsted p-value. Altman, https://doi.org/10.1136/bmj.d2090
DE_1182$corrected_z <- sqrt(0.743 - (2.404*log(DE_1182$adj.P.Val))) - 0.862
#calculate the one directional 95% confidence interval
DE_1182$error <-  abs(DE_1182$logFC /DE_1182$corrected_z) * 1.96

#arbitrary plotting cutoffs
adjSigCutoff <- 0.05
LFCcutoff <- 1
labelLFCcutoff <- 1

#pull out interesting proteins to label
proteins_sig <- na.omit(unique(DE_1182$Accession[DE_1182$adj.P.Val < adjSigCutoff]))
proteins_sig2 <- na.omit(unique(DE_1182$Accession[DE_1182$adj.P.Val < adjSigCutoff &
                                                    abs(DE_1182$logFC) > LFCcutoff]))
proteins_label_pos <- na.omit(unique(DE_1182$Accession[DE_1182$adj.P.Val < adjSigCutoff &
                                                         DE_1182$logFC > labelLFCcutoff]))
proteins_label_neg <- na.omit(unique(DE_1182$Accession[DE_1182$adj.P.Val < adjSigCutoff &
                                                         DE_1182$logFC < -labelLFCcutoff]))

#add in gene and protein name info from saved dictionaries
DE_1182$gene <- sapply(DE_1182$Accession, function(Accession) geneDict[[Accession]])
DE_1182$proteinName <- sapply(DE_1182$Accession, function(Accession) protNameDict[[Accession]])

#plot
ggplot() +
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(y = "Log2 Fold Change", x = "Normalized Average Expression (Log)",
       title = "MCF-10A | NDNB1182 vs DMSO") +
  geom_point(data = DE_1182[!DE_1182$Accession %in% proteins_sig,],
             aes(x = AveExpr, y = logFC), color = "grey50", alpha = 0.8, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 1) + #add line at LFC = 0
  geom_hline(yintercept = c(-1,1), linetype = 2) + #add lines at LFC = -1, 1
  geom_point(data = DE_1182[DE_1182$Accession %in% proteins_sig,],
             aes(x = AveExpr, y = logFC, color = logFC > 0), size = 0.8) +
  geom_segment(data = DE_1182[DE_1182$Accession %in% proteins_sig2,],
               aes(x = AveExpr, y = logFC - error, yend = logFC + error), lwd = 1) + #add 95% CI lines
  geom_segment(data = DE_1182[DE_1182$Accession %in% proteins_sig2,],
               aes(x = AveExpr, y = logFC - error, yend = logFC + error,
                   color = logFC > 0), lwd = 0.5) +
  geom_point(data = DE_1182[DE_1182$Accession %in% proteins_sig2,],
             aes(x = AveExpr, y = logFC, fill = logFC > 0), size = 3,
             shape = 21, stroke = 1.5) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_label_repel(data = DE_1182[DE_1182$Accession %in% proteins_label_pos,],
                   aes(x = AveExpr, y = logFC, label = gene), nudge_y = 1, nudge_x = -1.2) +
  geom_label_repel(data = DE_1182[DE_1182$Accession %in% proteins_label_neg,],
                   aes(x = AveExpr, y = logFC, label = gene), nudge_y = -1, nudge_x = -1.2)


ggsave("plots2/DEplot_1182vsWT.png", width = 12, height = 8)





ggplot() +
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(y = "B-H adj p-value (-log)", x = "Log2 Fold Change",
       title = "MCF-10A | NDNB1182 vs DMSO") +
  geom_point(data = DE_1182[!DE_1182$Accession %in% proteins_sig,],
             aes(x = logFC, y = -log(adj.P.Val, base = 10)), color = "grey50", alpha = 0.8, size = 0.8) +
  geom_hline(yintercept = -log(adjSigCutoff, base = 10), linetype = 2) + #add line at LFC = 0
  geom_vline(xintercept = c(-1,1), linetype = 2) + #add lines at LFC = -1, 1
  geom_point(data = DE_1182[DE_1182$Accession %in% proteins_sig,],
             aes(x = logFC, y = -log(adj.P.Val, base = 10), color = logFC > 0), size = 0.8) +
  geom_segment(data = DE_1182[DE_1182$Accession %in% proteins_sig2,],
               aes(x = logFC - error, xend = logFC + error, y = -log(adj.P.Val, base = 10)), lwd = 1) + #add 95% CI lines
  geom_segment(data = DE_1182[DE_1182$Accession %in% proteins_sig2,],
               aes(x = logFC - error, xend = logFC + error, -log(adj.P.Val, base = 10),
                   color = logFC > 0), lwd = 0.5) +
  geom_point(data = DE_1182[DE_1182$Accession %in% proteins_sig2,],
             aes(x = logFC, y = -log(adj.P.Val, base = 10), fill = logFC > 0), size = 3,
             shape = 21, stroke = 1.5) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_label_repel(data = DE_1182[DE_1182$Accession %in% proteins_label_pos,],
                   aes(x = logFC, y = -log(adj.P.Val, base = 10), label = gene), nudge_y = 1, nudge_x = 0.5) +
  geom_label_repel(data = DE_1182[DE_1182$Accession %in% proteins_label_neg,],
                   aes(x = logFC, y = -log(adj.P.Val, base = 10), label = gene), nudge_y = 1, nudge_x = -0.5)

ggsave("plots2/WT vs 1182 volcano.png",width = 12, height = 12)

#save data
write.table(DE_1182, "WTvs1182Processed2.tsv", sep = "\t",
            row.names = F, quote = F)


#SAME AS ABOVE BUT FOR OTHER TREATMENT
DE_NDNB <- topTable(fit = fitted_lm,
                    adjust.method = "BH",
                    number = Inf,
                    coef  = "proteinsBR$conditionNDNB",
                    sort.by = "p", confint = T)
DE_NDNB$Accession <- rownames(DE_NDNB)


DE_NDNB$corrected_z <- sqrt(0.743 - (2.404*log(DE_NDNB$adj.P.Val))) - 0.862
DE_NDNB$error <-  abs(DE_NDNB$logFC /DE_NDNB$corrected_z) * 1.96

adjSigCutoff <- 0.05
LFCcutoff <- 1
labelLFCcutoff <- 1
proteins_sig <- na.omit(unique(DE_NDNB$Accession[DE_NDNB$adj.P.Val < adjSigCutoff]))
proteins_sig2 <- na.omit(unique(DE_NDNB$Accession[DE_NDNB$adj.P.Val < adjSigCutoff &
                                                    abs(DE_NDNB$logFC) > LFCcutoff]))

proteins_label_pos <- na.omit(unique(DE_NDNB$Accession[DE_NDNB$adj.P.Val < adjSigCutoff &
                                                         DE_NDNB$logFC > labelLFCcutoff]))
proteins_label_neg <- na.omit(unique(DE_NDNB$Accession[DE_NDNB$adj.P.Val < adjSigCutoff &
                                                         DE_NDNB$logFC < -labelLFCcutoff]))

DE_NDNB$gene <- sapply(DE_NDNB$Accession, function(Accession) geneDict[[Accession]])
DE_NDNB$proteinName <- sapply(DE_NDNB$Accession, function(Accession) protNameDict[[Accession]])


ggplot() +
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(y = "Log2 Fold Change", x = "Normalized Average Expression (Log)",
       title = "MCF-10A | NDNB1 vs DMSO") +
  geom_point(data = DE_NDNB[!DE_NDNB$Accession %in% proteins_sig,],
             aes(x = AveExpr, y = logFC), color = "grey50", alpha = 0.8, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 1) + #add line at LFC = 0
  geom_hline(yintercept = c(-1,1), linetype = 2) + #add lines at LFC = -1, 1
  geom_point(data = DE_NDNB[DE_NDNB$Accession %in% proteins_sig,],
             aes(x = AveExpr, y = logFC, color = logFC > 0), size = 0.8) +
  geom_segment(data = DE_NDNB[DE_NDNB$Accession %in% proteins_sig2,],
               aes(x = AveExpr, y = logFC - error, yend = logFC + error), lwd = 1) + #add 95% CI lines
  geom_segment(data = DE_NDNB[DE_NDNB$Accession %in% proteins_sig2,],
               aes(x = AveExpr, y = logFC - error, yend = logFC + error,
                   color = logFC > 0), lwd = 0.5) +
  geom_point(data = DE_NDNB[DE_NDNB$Accession %in% proteins_sig2,],
             aes(x = AveExpr, y = logFC, fill = logFC > 0), size = 3,
             shape = 21, stroke = 1.5) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_label_repel(data = DE_NDNB[DE_NDNB$Accession %in% proteins_label_pos,],
                   aes(x = AveExpr, y = logFC, label = gene), nudge_y = 1, nudge_x = -1.2) +
  geom_label_repel(data = DE_NDNB[DE_NDNB$Accession %in% proteins_label_neg,],
                   aes(x = AveExpr, y = logFC, label = gene), nudge_y = -1, nudge_x = -1.2)


ggsave("plots2/DEplot_NDNBvsWT.png", width = 12, height = 8)




ggplot() +
  theme_bw(base_size = 25) +
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  labs(y = "B-H adj p-value (-log)", x = "Log2 Fold Change",
       title = "MCF-10A | NDNB1 vs DMSO") +
  geom_point(data = DE_NDNB[!DE_NDNB$Accession %in% proteins_sig,],
             aes(x = logFC, y = -log(adj.P.Val, base = 10)), color = "grey50", alpha = 0.8, size = 0.8) +
  geom_hline(yintercept = -log(adjSigCutoff, base = 10), linetype = 2) + #add line at LFC = 0
  geom_vline(xintercept = c(-1,1), linetype = 2) + #add lines at LFC = -1, 1
  geom_point(data = DE_NDNB[DE_NDNB$Accession %in% proteins_sig,],
             aes(x = logFC, y = -log(adj.P.Val, base = 10), color = logFC > 0), size = 0.8) +
  geom_segment(data = DE_NDNB[DE_NDNB$Accession %in% proteins_sig2,],
               aes(x = logFC - error, xend = logFC + error, y = -log(adj.P.Val, base = 10)), lwd = 1) + #add 95% CI lines
  geom_segment(data = DE_NDNB[DE_NDNB$Accession %in% proteins_sig2,],
               aes(x = logFC - error, xend = logFC + error, -log(adj.P.Val, base = 10),
                   color = logFC > 0), lwd = 0.5) +
  geom_point(data = DE_NDNB[DE_NDNB$Accession %in% proteins_sig2,],
             aes(x = logFC, y = -log(adj.P.Val, base = 10), fill = logFC > 0), size = 3,
             shape = 21, stroke = 1.5) +
  scale_color_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_label_repel(data = DE_NDNB[DE_NDNB$Accession %in% proteins_label_pos,],
                   aes(x = logFC, y = -log(adj.P.Val, base = 10), label = gene), nudge_y = 1, nudge_x = 0.5) +
  geom_label_repel(data = DE_NDNB[DE_NDNB$Accession %in% proteins_label_neg,],
                   aes(x = logFC, y = -log(adj.P.Val, base = 10), label = gene), nudge_y = 1, nudge_x = -0.5)

ggsave("plots2/WT vs NDNB volcano.png",width = 12, height = 12)



write.table(DE_NDNB, "WTvsNDNBProcessed2.tsv", sep = "\t",
            row.names = F, quote = F)



combined_cutoff_pVal <- 0.1

DE_joined <- full_join(DE_NDNB, DE_1182, by = c("Accession", "proteinName", "gene"),
                       suffix = c(".NDNB", ".1182"))

DE_joined <- DE_joined[DE_joined$adj.P.Val.1182 < combined_cutoff_pVal |
                         DE_joined$adj.P.Val.NDNB < combined_cutoff_pVal,]


colorCutoff <- 1
labelCutoff <- 1
ggplot(DE_joined) +
  geom_point(aes(x = logFC.NDNB, y = logFC.1182,
                 color = abs(logFC.NDNB - logFC.1182) > colorCutoff),
             size = 1.2) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_abline(slope = 1) +
  geom_abline(slope = 1, intercept = c(colorCutoff, -colorCutoff), linetype = 2) +
  geom_label_repel(data = DE_joined[abs(DE_joined$logFC.NDNB - DE_joined$logFC.1182) > labelCutoff,],
                   aes(x = logFC.NDNB, y = logFC.1182, label = gene)) +
  labs(x = "Log Fold Change NDNB1",
       y = "Log Fold Change NDNB1182",
       title = "MCF-10A")

ggsave("plots2/NDNB vs 1182 scatter.png", width = 8, height = 8)

