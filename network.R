
#Fig3A

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







#Fig3B and 3C

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
tmp$tidy_dataset() # 清理对象中的所有文件
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







#Fig3D

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

