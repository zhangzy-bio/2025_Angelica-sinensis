
###==========================differential ASVs analysis========================###

setwd("C:/Users/zhang/Desktop/R/DESeq2/")

library(tidyverse)
library(DESeq2)

ASV_table <- read.csv("asv.csv", row.names = 1) 

#group information
conditionC1 <- factor(c(rep("PC",6),rep("CK",6)),levels = c("PC","CK"))
conditionC2 <- factor(c(rep("BM",6),rep("CK",6)),levels = c("BM","CK"))
conditionC3 <- factor(c(rep("BA",6),rep("CK",6)),levels = c("BA","CK"))
conditionC4 <- factor(c(rep("CMA",6),rep("CK",6)),levels = c("CMA","CK"))
condition41 <- factor(c(rep("PC",6),rep("CMA",6)),levels = c("PC","CMA"))
condition42 <- factor(c(rep("BM",6),rep("CMA",6)),levels = c("BM","CMA"))
condition43 <- factor(c(rep("BA",6),rep("CMA",6)),levels = c("BA","CMA"))
colDataC1 <- data.frame(row.names = colnames(ASV_table[,c(1:6,25:30)]),conditionC1)
colDataC2 <- data.frame(row.names = colnames(ASV_table[,c(7:12,25:30)]),conditionC2)
colDataC3 <- data.frame(row.names = colnames(ASV_table[,c(13:18,25:30)]),conditionC3)
colDataC4 <- data.frame(row.names = colnames(ASV_table[,c(19:24,25:30)]),conditionC4)
colData41 <- data.frame(row.names = colnames(ASV_table[,c(1:6,19:24)]),condition41)
colData42 <- data.frame(row.names = colnames(ASV_table[,c(7:12,19:24)]),condition42)
colData43 <- data.frame(row.names = colnames(ASV_table[,c(13:18,19:24)]),condition43)

#dds matrix
ddsC1 <- DESeqDataSetFromMatrix(ASV_table[,c(1:6,25:30)],colDataC1,design = ~conditionC1)
ddsC2 <- DESeqDataSetFromMatrix(ASV_table[,c(7:12,25:30)],colDataC2,design = ~conditionC2)
ddsC3 <- DESeqDataSetFromMatrix(ASV_table[,c(13:18,25:30)],colDataC3,design = ~conditionC3)
ddsC4 <- DESeqDataSetFromMatrix(ASV_table[,c(19:24,25:30)],colDataC4,design = ~conditionC4)
dds41 <- DESeqDataSetFromMatrix(ASV_table[,c(1:6,19:24)],colData41,design = ~condition41)
dds42 <- DESeqDataSetFromMatrix(ASV_table[,c(7:12,19:24)],colData42,design = ~condition42)
dds43 <- DESeqDataSetFromMatrix(ASV_table[,c(13:18,19:24)],colData43,design = ~condition43)

#difference analysis
ddsC1 <- DESeq(ddsC1)
ddsC2 <- DESeq(ddsC2)
ddsC3 <- DESeq(ddsC3)
ddsC4 <- DESeq(ddsC4)
dds41 <- DESeq(dds41)
dds42 <- DESeq(dds42)
dds43 <- DESeq(dds43)

#identification of differentially abundant ASVs
resC1 <- results(ddsC1,contrast = c("conditionC1","PC","CK"))
resC1 <- resC1[order(resC1$pvalue),]
resC1data <- merge(as.data.frame(resC1), as.data.frame(counts(ddsC1, normalized = T)), by = "row.names", sort = F)
diff_ASV_C1 <- subset(resC1data, padj<0.05 & abs(log2FoldChange)>2)

resC2 <- results(ddsC2,contrast = c("conditionC2","BM","CK"))
resC2 <- resC2[order(resC2$pvalue),]
resC2data <- merge(as.data.frame(resC2), as.data.frame(counts(ddsC2, normalized = T)), by = "row.names", sort = F)
diff_ASV_C2 <- subset(resC2data, padj<0.05 & abs(log2FoldChange)>2)

resC3 <- results(ddsC3,contrast = c("conditionC3","BA","CK"))
resC3 <- resC3[order(resC3$pvalue),]
resC3data <- merge(as.data.frame(resC3), as.data.frame(counts(ddsC3, normalized = T)), by = "row.names", sort = F)
diff_ASV_C3 <- subset(resC3data, padj<0.05 & abs(log2FoldChange)>2)

resC4 <- results(ddsC4,contrast = c("conditionC4","CMA","CK"))
resC4 <- resC4[order(resC4$pvalue),]
resC4data <- merge(as.data.frame(resC4), as.data.frame(counts(ddsC4, normalized = T)), by = "row.names", sort = F)
diff_ASV_C4 <- subset(resC4data, padj<0.05 & abs(log2FoldChange)>2)

res41 <- results(dds41,contrast = c("condition41","PC","CMA"))
res41 <- res41[order(res41$pvalue),]
res41data <- merge(as.data.frame(res41), as.data.frame(counts(dds41, normalized = T)), by = "row.names", sort = F)
diff_ASV_41 <- subset(res41data, padj<0.05 & abs(log2FoldChange)>2)

res42 <- results(dds42,contrast = c("condition42","BM","CMA"))
res42 <- res42[order(res42$pvalue),]
res42data <- merge(as.data.frame(res42), as.data.frame(counts(dds42, normalized = T)), by = "row.names", sort = F)
diff_ASV_42 <- subset(res42data, padj<0.05 & abs(log2FoldChange)>2)

res43 <- results(dds43,contrast = c("condition43","BA","CMA"))
res43 <- res43[order(res43$pvalue),]
res43data <- merge(as.data.frame(res43), as.data.frame(counts(dds43, normalized = T)), by = "row.names", sort = F)
diff_ASV_43 <- subset(res43data, padj<0.05 & abs(log2FoldChange)>2)


#tax information
tax <- read.csv("tax.txt", row.names = 1, sep='\t')

taxC1 <- subset(tax, rownames(tax) %in% rownames(diff_ASV_C1))
taxC1 <- taxC1[match(rownames(diff_ASV_C1), rownames(taxC1)),]
diff_ASV_C1.bubble <- cbind(diff_ASV_C1, taxC1)
diff_ASV_C1.bubble <- diff_ASV_C1.bubble[,c("Phylum","log2FoldChange")]#Phylum,Order,Genus
diff_ASV_C1.bubble.c <- diff_ASV_C1.bubble %>% 
  group_by(Phylum) %>%
  summarise(positive = sum(log2FoldChange > 0),
            negative = sum(log2FoldChange < 0))

taxC2 <- subset(tax, rownames(tax) %in% rownames(diff_ASV_C2))
taxC2 <- taxC2[match(rownames(diff_ASV_C2), rownames(taxC2)),]
diff_ASV_C2.bubble <- cbind(diff_ASV_C2, taxC2)
diff_ASV_C2.bubble <- diff_ASV_C2.bubble[,c("Phylum","log2FoldChange")]
diff_ASV_C2.bubble.c <- diff_ASV_C2.bubble %>% 
  group_by(Phylum) %>%
  summarise(positive = sum(log2FoldChange > 0),
            negative = sum(log2FoldChange < 0))

taxC3 <- subset(tax, rownames(tax) %in% rownames(diff_ASV_C3))
taxC3 <- taxC3[match(rownames(diff_ASV_C3), rownames(taxC3)),]
diff_ASV_C3.bubble <- cbind(diff_ASV_C3, taxC3)
diff_ASV_C3.bubble <- diff_ASV_C3.bubble[,c("Phylum","log2FoldChange")]
diff_ASV_C3.bubble.c <- diff_ASV_C3.bubble %>% 
  group_by(Phylum) %>%
  summarise(positive = sum(log2FoldChange > 0),
            negative = sum(log2FoldChange < 0))

taxC4 <- subset(tax, rownames(tax) %in% rownames(diff_ASV_C4))
taxC4 <- taxC4[match(rownames(diff_ASV_C4), rownames(taxC4)),]
diff_ASV_C4.bubble <- cbind(diff_ASV_C4, taxC4)
diff_ASV_C4.bubble <- diff_ASV_C4.bubble[,c("Phylum","log2FoldChange")]
diff_ASV_C4.bubble.c <- diff_ASV_C4.bubble %>% 
  group_by(Phylum) %>%
  summarise(positive = sum(log2FoldChange > 0),
            negative = sum(log2FoldChange < 0))

tax41 <- subset(tax, rownames(tax) %in% rownames(diff_ASV_41))
tax41 <- tax41[match(rownames(diff_ASV_41), rownames(tax41)),]
diff_ASV_41.bubble <- cbind(diff_ASV_41, tax41)
diff_ASV_41.bubble <- diff_ASV_41.bubble[,c("Phylum","log2FoldChange")]
diff_ASV_41.bubble.c <- diff_ASV_41.bubble %>% 
  group_by(Phylum) %>%
  summarise(positive = sum(log2FoldChange > 0),
            negative = sum(log2FoldChange < 0))

tax42 <- subset(tax, rownames(tax) %in% rownames(diff_ASV_42))
tax42 <- tax42[match(rownames(diff_ASV_42), rownames(tax42)),]
diff_ASV_42.bubble <- cbind(diff_ASV_42, tax42)
diff_ASV_42.bubble <- diff_ASV_42.bubble[,c("Phylum","log2FoldChange")]
diff_ASV_42.bubble.c <- diff_ASV_42.bubble %>% 
  group_by(Phylum) %>%
  summarise(positive = sum(log2FoldChange > 0),
            negative = sum(log2FoldChange < 0))

tax43 <- subset(tax, rownames(tax) %in% rownames(diff_ASV_43))
tax43 <- tax43[match(rownames(diff_ASV_43), rownames(tax43)),]
diff_ASV_43.bubble <- cbind(diff_ASV_43, tax43)
diff_ASV_43.bubble <- diff_ASV_43.bubble[,c("Phylum","log2FoldChange")]
diff_ASV_43.bubble.c <- diff_ASV_43.bubble %>% 
  group_by(Phylum) %>%
  summarise(positive = sum(log2FoldChange > 0),
            negative = sum(log2FoldChange < 0))


diff_ASV_C1.bubble.c <- as.data.frame(diff_ASV_C1.bubble.c)
row.names(diff_ASV_C1.bubble.c) <- diff_ASV_C1.bubble.c$Phylum
colnames(diff_ASV_C1.bubble.c) <- c("Phy","1PC_pos","1PC_neg")

diff_ASV_C2.bubble.c <- as.data.frame(diff_ASV_C2.bubble.c)
row.names(diff_ASV_C2.bubble.c) <- diff_ASV_C2.bubble.c$Phylum
colnames(diff_ASV_C2.bubble.c) <- c("Phy","2BM_pos","2BM_neg")

diff_ASV_C3.bubble.c <- as.data.frame(diff_ASV_C3.bubble.c)
row.names(diff_ASV_C3.bubble.c) <- diff_ASV_C3.bubble.c$Phylum
colnames(diff_ASV_C3.bubble.c) <- c("Phy","3BA_pos","3BA_neg")

diff_ASV_C4.bubble.c <- as.data.frame(diff_ASV_C4.bubble.c)
row.names(diff_ASV_C4.bubble.c) <- diff_ASV_C4.bubble.c$Phylum
colnames(diff_ASV_C4.bubble.c) <- c("Phy","4CMA_pos","4CMA_neg")

diff_ASV_41.bubble.c <- as.data.frame(diff_ASV_41.bubble.c)
row.names(diff_ASV_41.bubble.c) <- diff_ASV_41.bubble.c$Phylum
colnames(diff_ASV_41.bubble.c) <- c("Phy","5PC_pos","5PC_neg")

diff_ASV_42.bubble.c <- as.data.frame(diff_ASV_42.bubble.c)
row.names(diff_ASV_42.bubble.c) <- diff_ASV_42.bubble.c$Phylum
colnames(diff_ASV_42.bubble.c) <- c("Phy","6BM_pos","6BM_neg")

diff_ASV_43.bubble.c <- as.data.frame(diff_ASV_43.bubble.c)
row.names(diff_ASV_43.bubble.c) <- diff_ASV_43.bubble.c$Phylum
colnames(diff_ASV_43.bubble.c) <- c("Phy","7BA_pos","7BA_neg")


dfs <- list(diff_ASV_C1.bubble.c, diff_ASV_C2.bubble.c, diff_ASV_C3.bubble.c, diff_ASV_C4.bubble.c,
            diff_ASV_41.bubble.c, diff_ASV_42.bubble.c, diff_ASV_43.bubble.c)
dfs.merge <- dfs[[1]]
for (i in 2:length(dfs)) {
  dfs.merge <- merge(dfs.merge, dfs[[i]], by="Phy", all=TRUE)
}
dfs.merge[is.na(dfs.merge)] <- 0
dfs.merge <- dfs.merge[rowSums(dfs.merge[,2:ncol(dfs.merge)]) > 3,]
rownames(dfs.merge) <- dfs.merge$Phy


#bubble plot
library(tidyr)
dfs.merge.long <- pivot_longer(dfs.merge, cols = -"Phy", 
                               names_to = "group", values_to = "counts") 
dfs.merge.long <- dfs.merge.long %>%
  mutate(groupp = case_when(
    group == "1PC_pos" ~ "PC+", group == "1PC_neg" ~ "PC-", 
    group == "2BM_pos" ~ "BM+", group == "2BM_neg" ~ "BM-", 
    group == "3BA_pos" ~ "BA+", group == "3BA_neg" ~ "BA-", 
    group == "4CMA_pos" ~ "CMA+", group == "4CMA_neg" ~ "CMA-", 
    group == "5PC_pos" ~ "PC+", group == "5PC_neg" ~ "PC-", 
    group == "6BM_pos" ~ "BM+", group == "6BM_neg" ~ "BM-", 
    group == "7BA_pos" ~ "BA+", group == "7BA_neg" ~ "BA-"
  ))
dfs.merge.long$group <- factor(dfs.merge.long$group,levels = c("1PC_pos","1PC_neg","2BM_pos","2BM_neg",
                                                               "3BA_pos","3BA_neg","4CMA_pos","4CMA_neg",
                                                               "5PC_pos","5PC_neg","6BM_pos","6BM_neg",
                                                               "7BA_pos","7BA_neg"))
pdf(file = "bubble.pdf", width = 7, height = 4)
ggplot(dfs.merge.long, aes(x = group, y = Phy)) +
  geom_point(aes(size = counts, color = groupp), alpha = 1) +
  scale_color_manual(values = c("#FF7F00","#CD8162","#1E90FF","#4682B4",
                                "#8470FF","#483D8B","#4F4F4F","#9C9C9C"),
                     limits = c("PC+","PC-","BM+","BM-","BA+","BA-","CMA+","CMA-")) +
  scale_x_discrete(labels = c("1PC_pos"="PC+","1PC_neg"="PC-","2BM_pos"="BM+","2BM_neg"="BM-",
                              "3BA_pos"="BA+","3BA_neg"="BA-","4CMA_pos"="CMA+","4CMA_neg"="CMA-",
                              "5PC_pos"="PC+","5PC_neg"="PC-","6BM_pos"="BM+","6BM_neg"="BM-",
                              "7BA_pos"="BA+","7BA_neg"="BA-")) +
  scale_size_continuous(range = c(-1, 4)) +
  geom_text(aes(label = counts), size = 3) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#EAEAEA"),
        panel.border = element_rect(colour = "black",fill = NA),
        axis.text.x = element_text(angle = 45, vjust = 0.7)) +
  geom_vline(xintercept = 8.5)
dev.off()






#functional prediction analysis

setwd("C:/Users/zhang/Desktop/R/predict-function/")
library(tidyverse)
library(DESeq2)
#KEGG pathway
KEGG.raw <- read.csv("KO1-4.txt", sep = "\t")
KEGG.raw <- KEGG.raw[,c(1,4)]

KEGG <- read.csv("pred_metagenome_unstrat_16s-root.tsv", sep = "\t")
KEGG.m <- merge(KEGG.raw, KEGG, by = 1, all.x = T)
KEGG.m <- unique(KEGG.m)
KEGG.m <- na.omit(KEGG.m)
KEGG.m <- KEGG.m[,-1]
KEGG.m <- aggregate(KEGG.m[,-1], by=list(KEGG.m$Pathway),sum)
rownames(KEGG.m) <- KEGG.m$Group.1
KEGG.m <- KEGG.m[,-1]
KEGG.m <- round(KEGG.m,0)

#differential pathway analysis, DESeq2
conditionC1.k <- factor(c(rep("PC",6),rep("CK",6)), levels = c("PC","CK"))
colDataC1.k <- data.frame(row.names = colnames(KEGG.m)[c(1:6,25:30)], conditionC1.k)#1:6,7:12,13:18,19:24
# colDataC1.k <- data.frame(row.names = colnames(KEGG.m)[c(7:12,25:30)], conditionC1.k)
# colDataC1.k <- data.frame(row.names = colnames(KEGG.m)[c(13:18,25:30)], conditionC1.k)
# colDataC1.k <- data.frame(row.names = colnames(KEGG.m)[c(19:24,25:30)], conditionC1.k)
ddsC1.k <- DESeqDataSetFromMatrix(KEGG.m[,c(1:6,25:30)],colDataC1.k,design = ~ conditionC1.k)
# ddsC1.k <- DESeqDataSetFromMatrix(KEGG.m[,c(7:12,25:30)],colDataC1.k,design = ~ conditionC1.k)
# ddsC1.k <- DESeqDataSetFromMatrix(KEGG.m[,c(13:18,25:30)],colDataC1.k,design = ~ conditionC1.k)
# ddsC1.k <- DESeqDataSetFromMatrix(KEGG.m[,c(19:24,25:30)],colDataC1.k,design = ~ conditionC1.k)
ddsC1.k <- DESeq(ddsC1.k)
resC1.k <- results(ddsC1.k,contrast = c("conditionC1.k","PC","CK"))
KEGG.m$log2fc <- resC1.k$log2FoldChange
KEGG.m$padj <- resC1.k$padj
resC1.k <- resC1.k[order(resC1.k$pvalue),]
resC1data.k <- merge(as.data.frame(resC1.k), as.data.frame(counts(ddsC1.k, normalized = T)), by = "row.names", sort = F)
diff_ASV_C1.k <- subset(resC1data.k, padj<0.05 & abs(log2FoldChange)>1)
rownames(diff_ASV_C1.k) <- diff_ASV_C1.k$Row.names

row.c1 <- rownames(diff_ASV_C1.k)
row.c2 <- rownames(diff_ASV_C1.k)
row.c3 <- rownames(diff_ASV_C1.k)
row.c4 <- rownames(diff_ASV_C1.k)

ko.row.all <- unique(c(as.character(row.c1),as.character(row.c2),
                       as.character(row.c3),as.character(row.c4)))

KEGG.m <- KEGG.m[rownames(KEGG.m) %in% ko.row.all, ]
KEGG.m <- KEGG.m[,-c(ncol(KEGG.m)-1,ncol(KEGG.m))] 

#heatmap
KEGG.m$PC <- rowMeans(KEGG.m[,1:6])
KEGG.m$BM <- rowMeans(KEGG.m[,7:12])
KEGG.m$BA <- rowMeans(KEGG.m[,13:18])
KEGG.m$CMA <- rowMeans(KEGG.m[,19:24])
KEGG.m$CK <- rowMeans(KEGG.m[,25:30])
KEGG.m <- KEGG.m[,c((ncol(KEGG.m)-4):ncol(KEGG.m))]
max(KEGG.m)
min(KEGG.m)

#data standardization
KEGG.m.long <- gather(KEGG.m, key = "Variable", value = "Value")
ggplot(KEGG.m.long, aes(x = Value)) +
  geom_histogram(binwidth = 100, fill = "orange", color = "black") +
  theme_minimal()

KEGG.m.log <- log10(KEGG.m)
KEGG.m.long <- gather(KEGG.m.log, key = "Variable", value = "Value")
ggplot(KEGG.m.long, aes(x = Value)) +
  geom_histogram(binwidth = 0.1, fill = "orange", color = "black") +
  theme_minimal()

max(KEGG.m.log)
min(KEGG.m.log)

#plot
KEGG.m.log$C1 <- KEGG.m.log[,1] - KEGG.m.log[,5]
KEGG.m.log$C2 <- KEGG.m.log[,2] - KEGG.m.log[,5]
KEGG.m.log$C3 <- KEGG.m.log[,3] - KEGG.m.log[,5]
KEGG.m.log$C4 <- KEGG.m.log[,4] - KEGG.m.log[,5]

KEGG.m.log <- KEGG.m.log[order(KEGG.m.log[[ncol(KEGG.m.log)]], decreasing = T),]

KEGG.m.log$signC1 <- NA_character_
matched.rowC1 <- rownames(KEGG.m.log) %in% row.c1
KEGG.m.log[matched.rowC1 & KEGG.m.log[,6] > 0, "signC1"] <- "+"
KEGG.m.log[matched.rowC1 & KEGG.m.log[,6] < 0, "signC1"] <- "-"

KEGG.m.log$signC2 <- NA_character_
matched.rowC2 <- rownames(KEGG.m.log) %in% row.c2
KEGG.m.log[matched.rowC2 & KEGG.m.log[,7] > 0, "signC2"] <- "+"
KEGG.m.log[matched.rowC2 & KEGG.m.log[,7] < 0, "signC2"] <- "-"

KEGG.m.log$signC3 <- NA_character_
matched.rowC3 <- rownames(KEGG.m.log) %in% row.c3
KEGG.m.log[matched.rowC3 & KEGG.m.log[,8] > 0, "signC3"] <- "+"
KEGG.m.log[matched.rowC3 & KEGG.m.log[,8] < 0, "signC3"] <- "-"

KEGG.m.log$signC4 <- NA_character_
matched.rowC4 <- rownames(KEGG.m.log) %in% row.c4
KEGG.m.log[matched.rowC4 & KEGG.m.log[,9] > 0, "signC4"] <- "+"
KEGG.m.log[matched.rowC4 & KEGG.m.log[,9] < 0, "signC4"] <- "-"

KEGG.m.log.sig <- KEGG.m.log[,-c(1:9)]
KEGG.m.log.sig$signC5 <- NA_character_
colnames(KEGG.m.log.sig) <- colnames(KEGG.m.log)[1:5]
KEGG.m.log.sig[is.na(KEGG.m.log.sig)] <- ""

library(pheatmap)
ph <- pheatmap(KEGG.m.log[,1:5], 
               cluster_cols = F, cluster_rows = F,
               color = colorRampPalette(colors = c("blue","#9090FF","white","#FFAEAE","red"))(100),
               angle_col = 45, cellwidth = 16, cellheight = 13,
               display_numbers = KEGG.m.log.sig)

save_pheatmap_pdf <- function(x, filename, width=15, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(ph, "heatmap.pdf")













###==========================differential metabolite analysis========================###


setwd("C:/Users/zhang/Desktop/R/DESeq2+VIP/")
#DESeq2 
meta <- read.csv("metabolites.csv", header=T, row.names = 1) 
meta <- round(meta)

conditionC4 <- factor(c(rep("CK",6),rep("CMA",6)),
                      levels = c("CK","CMA"))
colDataC4 <- data.frame(row.names = colnames(meta[,c(7:12,1:6)]),conditionC4)
colDataC4

library(DESeq2)
ddsC4 <- DESeqDataSetFromMatrix(meta[,c(7:12,1:6)],colDataC4,design = ~conditionC4)
ddsC4 <- DESeq(ddsC4)

resC4 <- results(ddsC4,contrast = c("conditionC4","CMA","CK"))
resC4 <- resC4[order(resC4$pvalue),]
head(resC4)
summary(resC4)

diff_gene_C4 <- subset(resC4, padj<0.001 & abs(log2FoldChange)>2.5)
dim(diff_gene_C4)
count(diff_gene_C4$log2FoldChange>0)
count(diff_gene_C4$log2FoldChange<0)
diff_gene_C4 <- as.data.frame(diff_gene_C4)
write.csv(diff_gene_C4,file = "DESeq2-diff.csv")

#DESeq2+VIP
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










###==========================co-network analysis========================###

setwd("C:/Users/zhang/Desktop/R/Network/")

library(microeco)
library(meconetcomp)
library(magrittr)
library(igraph) 
library(ggplot2)
library(ape)
library(dplyr)
library(WGCNA)
library(rgexf)

ASV_table <- read.csv("asv.csv", row.names = 1) 
ASV_table <- ASV_table[,1:24]
ASV_table <- t(ASV_table)

#ASV data filtering
prev_prop_df <- apply(X = ASV_table, MARGIN = 2, 
                      FUN = function(x){sum(x > 0)/dim(ASV_table)[1]})
is.above.prev <- prev_prop_df >= 0.2
table(is.above.prev)

total_relative_abundance <- colSums(ASV_table)/sum(colSums(ASV_table)) 
is.above.abundance <- total_relative_abundance > 0.0002
table(is.above.abundance)

table(is.above.prev & is.above.abundance)
target_ASV <- colnames(ASV_table)[is.above.prev & is.above.abundance]

#relative abundance
ASV_table.prop <- as.data.frame(t(apply(ASV_table, 1, function(ASV) ASV/sum(ASV))))
target_ASV_table.prop = ASV_table.prop[, target_ASV]
print(rowSums(target_ASV_table.prop))
target_ASV_table.prop <- as.data.frame(t(target_ASV_table.prop))

#group and taxonomy
group <- read.csv("group.txt",row.names = 1,sep = "\t")
group$type2 <- c(rep(c("GBI"),24),rep(c("CK"),6))
group <- group[1:24,]
phylo <- read.tree("tree.nwk", tree.names = NULL)
taxonomy <- read.csv("taxonomy.txt", row.names = 1, sep='\t') 
taxonomy <- subset(taxonomy, row.names(taxonomy) %in% row.names(target_ASV_table.prop))

#creat a microtable object
soil_amp <- microtable$new(otu_table = target_ASV_table.prop, sample_table = group, phylo_tree = phylo, tax_table = taxonomy)
soil_amp

soil_amp$tax_table %<>% subset(Kingdom %in% c("d__Bacteria","k__Fungi")) 
soil_amp$filter_pollution(taxa = c("mitochondria", "chloroplast"))
soil_amp$tidy_dataset()

#correlation
t1 <- trans_network$new(dataset = soil_amp, 
                        cal_cor = "WGCNA", 
                        taxa_level = "OTU",
                        cor_method = "spearman") 

#build modular network
t1$cal_network(p_thres = 0.05, 
               COR_cut = 0.6,
               add_taxa_name = c("Kingdom","Phylum"),
               delete_unlinked_nodes = TRUE,
               COR_optimization = FALSE)

t1$cal_module(
  method = "cluster_fast_greedy",
  module_name_prefix = "M")

module_colors <- c("M1" = "#291173",
                   "M2" = "#356529",
                   "M3" = "#611A44",
                   "M4" = "#629586",
                   "M5" = "#A5C8DD",
                   "M6" = "#EC6B2D",
                   "M7" = "#6D8EF7",
                   "M8" = "#7935F6")
other_color <- "grey"

V(t1$res_network)$color <- ifelse(V(t1$res_network)$module %in% names(module_colors), 
                                  module_colors[V(t1$res_network)$module], 
                                  other_color)
table(V(t1$res_network)$module)

t1$save_network("n_module_16s+its-soil1234.gexf")







#network stability analysis

library(microeco) 
library(meconetcomp) 
library(magrittr)
library(igraph)
library(ggplot2)
library(ape)
library(dplyr)

soil_amp_network <- list()

tmp <- clone(soil_amp) 
tmp$sample_table %<>% subset(type2 == "GBI") 
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0002) 
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.6) 
soil_amp_network$GBI <- tmp

tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(type == "PC")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0002)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
soil_amp_network$PC <- tmp

tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(type == "BM")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0002)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
soil_amp_network$BM <- tmp 

tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(type == "BA")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0002)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
soil_amp_network$BA <- tmp

tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(type == "CMA")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0002)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
soil_amp_network$CMA <- tmp 

tmp <- clone(soil_amp)
tmp$sample_table %<>% subset(type == "CK")
tmp$tidy_dataset()
tmp <- trans_network$new(dataset = tmp, cor_method = "spearman", filter_thres = 0.0002)
tmp$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)
soil_amp_network$CK <- tmp 


#robustness
tmp <- robustness$new(soil_amp_network, remove_strategy = c("edge_rand", "edge_strong", "node_rand", "node_degree_high"), 
                      remove_ratio = seq(0, 0.99, 0.1), measure = c("Eff", "Eigen", "Pcr"), run = 100)

write.csv(tmp$res_table,"robustness_detail.csv")
write.csv(tmp$res_summary,"robustness_summary.csv")

tmp$plot(linewidth = 1)
ggsave("r_pos_robustness.pdf", tmp$plot(linewidth = 1), width = 10, height = 7)

#cohesion
t1 <- cohesionclass$new(soil_amp_network)
write.csv(t1$res_list$sample,"Cohesion.csv")
t1$cal_diff(method = "anova")
g1<-t1$plot(measure = "r_pos")
g1
ggsave("r_pos_cohesion.pdf", g1, width = 7, height = 7)







#keystone species 

library(igraph)
library(ggplot2)
library(dplyr)
igraph <- read_graph("network.graphml",format = "graphml")
wtc <- cluster_louvain(igraph, NA)
modularity(wtc)
V(igraph)$module <- membership(wtc)

set.seed(123)
#z-score
within_module_deg_z_score <- function(g, A=NULL, weighted=F) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse = F, names = T, attr = "weight")
    } else {
      A <- as_adj(g, sparse = F, names = T)
    }
  }
  memb <- vertex_attr(g, "module")
  N <- max(memb)
  nS <- tabulate(memb)
  z <- Ki <- rep.int(0, dim(A)[1L])
  Ksi <- sigKsi <- rep.int(0, N)
  names(z) <- names(Ki) <- rownames(A)
  for (S in seq_len(N)) {
    x <- rowSums(A[memb == S, memb == S])
    Ki[memb == S] <- x
    Ksi[S] <- sum(x) / nS[S]
    sigKsi[S] <- sqrt(sum((x - Ksi[S])^2) / (nS[S]-1))
  }
  z <- (Ki - Ksi[memb]) / sigKsi[memb]
  z[is.infinite(z)] <- 0
  df <- data.frame(Ki,z,row.names = names(Ki))
  return(df)
}

part_coeff <- function(g, A=NULL, weighted=F) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse = F, attr = "weight")
    } else {
      A <- as_adj(g, sparse = F)
    }
  }
  memb <- vertex_attr(g, "module")
  Ki <- colSums(A)
  N <- max(memb)
  Kis <- t(rowsum(A, memb))
  pi <- 1 - ((1 / Ki^2) * rowSums(Kis^2))
  names(pi) <- rownames(A)
  return(pi)
}

zi <- within_module_deg_z_score(igraph)
pi <- part_coeff(igraph)
zi_pi <- data.frame(zi, pi)
zi_pi <- na.omit(zi_pi)

zi_pi[which(zi_pi$z < 2.5 & zi_pi$pi < 0.62), "type"] <- "Peripherals"
zi_pi[which(zi_pi$z < 2.5 & zi_pi$pi > 0.62), "type"] <- "Connectors"
zi_pi[which(zi_pi$z > 2.5 & zi_pi$pi < 0.62), "type"] <- "Module hubs"
zi_pi[which(zi_pi$z > 2.5 & zi_pi$pi > 0.62), "type"] <- "Network hubs"

tax <- read.csv("taxonomy.txt", row.names = 1, sep='\t')
zi_pi.plot <- merge(zi_pi, tax, by = "row.names", all.x = T)
zi_pi.plot$group <- c(rep("Soil",nrow(zi_pi.plot)))
write.csv(zi_pi.plot, "zi_pi.plot.csv")

ggplot(zi_pi.plot, aes(pi, z, shape = Kingdom)) +
  geom_point(aes(color = type), alpha = 0.5, size = 4) +
  scale_color_manual(values = c("gray","red","blue","purple"),
                     limits = c("Peripherals","Connectors","Module hubs","Network hubs")) +
  scale_shape_manual(values = c("d__Bacteria" = 16, "k__Fungi" = 17)) +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.key = element_blank(),
        panel.border = element_rect(fill=NA,color="black",linetype="solid")) +
  labs(x = "pi", y = "zi", color = "") +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)

zi_pi.hub <- subset(zi_pi, type %in% c("Connectors", "Module hubs", "Network hubs"))
tax.hub <- tax[rownames(zi_pi.hub),]
zi_pi.hub <- merge(zi_pi.hub, tax.hub, by = "row.names", all.x = T)
rownames(zi_pi.hub) <- zi_pi.hub$Row.names

write.csv(zi_pi.hub, "zi_pi.hub.csv")
write.csv(zi_pi, "zi_pi.csv")

#plot
plot.root <- read.csv("zi_pi.plot_root.csv", header = T, row.names = 1)
plot.soil <- read.csv("zi_pi.plot_soil.csv", header = T, row.names = 1)
plot.all <- rbind(plot.root, plot.soil)
plot.all[which(plot.all$type == "Peripherals"), "color"] <- "gray"
plot.all[which(plot.all$type == "Connectors" & plot.all$group == "Root"), "color"] <- "red"
plot.all[which(plot.all$type == "Module hubs" & plot.all$group == "Root"), "color"] <- "red"
plot.all[which(plot.all$type == "Network hubs" & plot.all$group == "Root"), "color"] <- "red"
plot.all[which(plot.all$type == "Connectors" & plot.all$group == "Soil"), "color"] <- "blue"
plot.all[which(plot.all$type == "Module hubs" & plot.all$group == "Soil"), "color"] <- "blue"
plot.all[which(plot.all$type == "Network hubs" & plot.all$group == "Soil"), "color"] <- "blue"

pdf(file = "zi-pi_plot_total.pdf", width = 6, height = 6)
ggplot(plot.all, aes(pi, z, shape = Kingdom, group = group)) +
  geom_point(aes(color = color), size = 5) +
  scale_color_manual(values = c("gray","#3CB371","#EEB422"),
                     limits = c("gray","red","blue")) +
  scale_shape_manual(values = c("d__Bacteria" = 16, "k__Fungi" = 17)) +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.key = element_blank(),
        panel.border = element_rect(fill=NA,color="black",linetype="solid",size = 1),
        axis.text = element_text(size = 12)) +
  labs(x = "pi", y = "zi", color = "", size = 30) +
  geom_vline(xintercept = 0.62, size = 1) +
  geom_hline(yintercept = 2.5, size = 1)
dev.off()




