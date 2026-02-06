
#Fig 4B

setwd("C:/Users/zhang/Desktop/R/DESeq2+VIP/")
# limma 
library(tidyverse)
library(limma)

expr <- read.csv("../a.OPLS-DA/1.metabolites_exp_soil-C4.csv", header=T, row.names = 1) 
group <- read.csv("group.csv",row.names = 1)
group <- group[colnames(expr), , drop=F]
expr_log <- log2(expr)

# Build a design matrix
group$group <- relevel(factor(group$group), ref = "S_CK") 

design <- model.matrix(~ group, data = group)
design

fit <- lmFit(expr_log, design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "groupS_CMA", number = Inf, adjust.method = "fdr")#logFC——代表基因表达的对数2倍变化（log2 fold change）

sig.meta <- subset(res, adj.P.Val < 0.05 & abs(logFC) >= 2.5) 
diff_gene_C4 <- sig.meta
write.csv(sig.meta,file = "1.limma-diff_soil-C4.csv")


#limma+VIP
VIP.c4 <- read.csv("VIP_General List_root.txt", sep = "\t")
VIP.C4.p <- read.csv("Descriptive Statistics - M2.txt", sep = "\t")
VIP.C4.p <- VIP.C4.p[,1:2]
colnames(VIP.c4)[1] <- "meta"
colnames(VIP.c4)[2] <- "vip"
colnames(VIP.C4.p)[1] <- "meta"
VIP.c4 <- merge(VIP.c4, VIP.C4.p, by.x = "meta", by.y = "meta", all.x = T)

VIP.c4 <- VIP.c4[VIP.c4$vip > 1.8 & VIP.c4$Probability < 0.05,] 
rownames(VIP.c4) <- VIP.c4$meta

commonrows <- intersect(rownames(diff_gene_C4), rownames(VIP.c4))
diff_gene_C4 <- diff_gene_C4[commonrows, ]
VIP.c4 <- VIP.c4[commonrows, ]
DV.c4 <- cbind(diff_gene_C4, VIP.c4)

#relative.abundance
re.ab <- read.csv("metabolites.csv", header=T, row.names = 1)
library(dplyr)
re.ab <- re.ab / colSums(re.ab)
re.ab <- re.ab[commonrows,]
re.ab <- re.ab %>% mutate(mean4 = rowMeans(across(1:6)), meanC = rowMeans(across(7:12)),
                          sum = rowSums(across(1:12)))
DV.c4 <- cbind(DV.c4, re.ab)

metaname <- read.csv("metabolites_full_table.csv")
metaname <- metaname[,2:3]
rownames(metaname) <- metaname$meta
metaname <- metaname[commonrows,]
DV.c4 <- cbind(DV.c4, metaname)
DV.c4 <- DV.c4[order(-DV.c4$sum),]

DV.c4[which(DV.c4$log2FoldChange > 0), "group"] <- "#6888F5"
DV.c4[which(DV.c4$log2FoldChange < 0), "group"] <- "#D77071"
write.csv(DV.c4, "soil.csv")


#plot
DV.c4.plot <- data.frame(group = c(rep("CK", nrow(DV.c4)), rep("CMA", nrow(DV.c4))),
                         meta = factor(c(DV.c4$meta, DV.c4$meta), levels = unique(DV.c4$meta)),
                         reab = c(DV.c4$meanC, DV.c4$mean4),
                         order = c(1:nrow(DV.c4), 1:nrow(DV.c4)))
DV.c4.label <- DV.c4.plot %>% group_by(group) %>%
  arrange(group, desc(order)) %>%
  mutate(y.pos = cumsum(reab) - 0.5 * reab) %>% ungroup()

library(ggplot2)
library(ggalluvial)
library(see)
library(ggrepel)
pdf("soil.pdf", width = 8, height = 10)
ggplot(DV.c4.label, aes(x = group, y = reab, fill = meta,
                        stratum = meta, alluvium = meta, label = meta)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_y_continuous(expand = c(0, 0), trans = "sqrt") + 
  scale_x_discrete(limits = c("CMA", "CK")) +
  geom_alluvium(width = 0.5, alpha = 0.4, colour = "black", curve_type = "arctan") +
  geom_stratum(width = 0.5, size = 0.4) + 
  scale_fill_manual(values = DV.c4$group) + 
  theme_classic() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
dev.off()
