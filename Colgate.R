setwd("~/OneDrive - University of California, San Diego Health/Projects/Colgate_Wender") #set you own working directory

library(tidyverse)
library(mixOmics)
library(ggpubr)
library(CoDaSeq)
library(vegan)
library(pheatmap)

### Read data ###

data <- read_csv("teeth_quant_table.csv") %>% arrange(Sample)
metadata <- read_csv("Metadata.csv") %>% arrange(Sample)

all(data$Sample == metadata$Sample)

### Generate PCA ###

# PCA raw 
PCA_whole <- mixOmics::pca(column_to_rownames(data, var = "Sample"), ncomp = 2, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X, metadata)
PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "Attribute_type", alpha = 0.6, 
            title = paste("PCA RAW -", "Attribute_type", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores %>% group_by(Attribute_type) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = Attribute_type), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PERMANOVA raw
dist_metabolites <- vegdist(column_to_rownames(data, var = "Sample"), method = "euclidean")
permanova <- adonis2(dist_metabolites ~ PCA_whole_scores$Attribute_type, PCA_whole_scores, na.action = na.omit)

# PCA CLR
data_clr <- codaSeq.clr(column_to_rownames(data, var = "Sample") + 1) %>% as.data.frame()
PCA_whole_clr <- mixOmics::pca(data_clr, ncomp = 2, scale = FALSE)
PCA_whole_clr_scores <- data.frame(PCA_whole_clr$variates$X, metadata)
PCA_clr_plot <- PCA_whole_clr_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = "Attribute_type", alpha = 0.6, 
            title = paste("PCA CLR -", "Attribute_type", sep = " "),
            xlab = paste("PC1 (", round(PCA_whole_clr$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole_clr$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_clr_scores %>% group_by(Attribute_type) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = Attribute_type), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PERMANOVA CLR
dist_metabolites_clr <- vegdist(data_clr, method = "euclidean")
permanova_clr <- adonis2(dist_metabolites_clr ~ PCA_whole_clr_scores$Attribute_type, PCA_whole_clr_scores, na.action = na.omit)

#ggsave(plot = PCA_clr_plot, filename = "PCA_CLR.svg", device = "svg", dpi = "retina")

### Generate pairwise PLS-DA models ###

# PLS_DA - Control vs MPS
data_clr_CT_MPS <- data_clr %>% dplyr::filter(metadata$Attribute_type != "H2O2")
metadata_clr_CT_MPS <- metadata %>% dplyr::filter(metadata$Attribute_type != "H2O2")

PLSDA_CT_MPS <- mixOmics::plsda(data_clr_CT_MPS, metadata_clr_CT_MPS$Attribute_type, ncomp = 3, scale = TRUE)
PLSDA_CT_MPS_scores <- data.frame(PLSDA_CT_MPS$variates$X, metadata_clr_CT_MPS)

PLSDA_CT_MPS_plot <- PLSDA_CT_MPS_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Attribute_type", alpha = 0.6, title = "PLS-DA CLR - CT v MPS", 
            xlab = paste("Component 1 (", round(PLSDA_CT_MPS$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_CT_MPS$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = c("#F8766D", "#619CFF")) +
  geom_point(data = PLSDA_CT_MPS_scores %>% group_by(Attribute_type) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Attribute_type), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_CT_MPS <- plotLoadings(PLSDA_CT_MPS, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib, importance)

perf_plsda_CT_MPS <- perf(PLSDA_CT_MPS, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_plsda_CT_MPS, legend = FALSE)

VIPs_CT_MPS <- as.data.frame(mixOmics::vip(PLSDA_CT_MPS))
VIPs_CT_MPS_filter <- dplyr::filter(VIPs_CT_MPS, VIPs_CT_MPS$comp1 > 1)
VIPs_CT_MPS_filter$ID <- rownames(VIPs_CT_MPS_filter)
VIPs_CT_MPS_select <- VIPs_CT_MPS_filter %>% dplyr::select(ID, comp1)
VIPs_CT_MPS_Load <- VIPs_CT_MPS_select %>% left_join(Loadings_CT_MPS, by = c("ID" = "rowname")) %>% arrange(desc(comp1))
VIPs_CT_MPS_Load_Top25 <- VIPs_CT_MPS_Load %>% head(., 25)

Loadings_CT_MPS_plot <- plotLoadings(PLSDA_CT_MPS, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_CT_MPS_Load_Top25$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "ascending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLS-DA - CT v MPS") + geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

#write_csv(VIPs_CT_MPS_Load, "PLSDA_VIP_CT_MPS_Full.csv")
#ggsave(plot = PLSDA_CT_MPS_plot, filename = "PLSDA_CT_MPS_plot.svg", device = "svg", dpi = "retina")


# PLS_DA - Control vs H2O2
data_clr_CT_H2O2 <- data_clr %>% dplyr::filter(metadata$Attribute_type != "MPS")
metadata_clr_CT_H2O2 <- metadata %>% dplyr::filter(metadata$Attribute_type != "MPS")

PLSDA_CT_H2O2 <- mixOmics::plsda(data_clr_CT_H2O2, metadata_clr_CT_H2O2$Attribute_type, ncomp = 3, scale = TRUE)
PLSDA_CT_H2O2_scores <- data.frame(PLSDA_CT_H2O2$variates$X, metadata_clr_CT_H2O2)

PLSDA_CT_H2O2_plot <- PLSDA_CT_H2O2_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Attribute_type", alpha = 0.6, title = "PLS-DA CLR - CT v H2O2",
            xlab = paste("Component 1 (", round(PLSDA_CT_H2O2$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_CT_H2O2$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = c("#F8766D", "#00BA38")) +
  geom_point(data = PLSDA_CT_H2O2_scores %>% group_by(Attribute_type) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Attribute_type), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_CT_H2O2 <- plotLoadings(PLSDA_CT_H2O2, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib, importance)

perf_plsda_CT_H2O2 <- perf(PLSDA_CT_H2O2, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_plsda_CT_H2O2, legend = TRUE)

VIPs_CT_H2O2 <- as.data.frame(mixOmics::vip(PLSDA_CT_H2O2))
VIPs_CT_H2O2_filter <- dplyr::filter(VIPs_CT_H2O2, VIPs_CT_H2O2$comp1 > 1)
VIPs_CT_H2O2_filter$ID <- rownames(VIPs_CT_H2O2_filter)
VIPs_CT_H2O2_select <- VIPs_CT_H2O2_filter %>% dplyr::select(ID, comp1)
VIPs_CT_H2O2_Load <- VIPs_CT_H2O2_select %>% left_join(Loadings_CT_H2O2, by = c("ID" = "rowname")) %>% arrange(desc(comp1))
VIPs_CT_H2O2_Load_Top25 <- VIPs_CT_H2O2_Load %>% head(., 25)

Loadings_CT_H2O2_plot <- plotLoadings(PLSDA_CT_H2O2, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_CT_H2O2_Load_Top25$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "ascending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLS-DA - CT v H2O2") + geom_hline(yintercept = 0, linetype = 2, color = "lightgray")


#write_csv(VIPs_CT_H2O2_Load, "PLSDA_VIP_CT_H2O2_Full.csv")
#ggsave(plot = PLSDA_CT_H2O2_plot, filename = "PLSDA_CT_H2O2_plot.svg", device = "svg", dpi = "retina")


# PLS_DA - MPS vs H2O2
data_clr_MPS_H2O2 <- data_clr %>% dplyr::filter(metadata$Attribute_type != "Control")
metadata_clr_MPS_H2O2 <- metadata %>% dplyr::filter(metadata$Attribute_type != "Control")

PLSDA_MPS_H2O2 <- mixOmics::plsda(data_clr_MPS_H2O2, metadata_clr_MPS_H2O2$Attribute_type, ncomp = 3, scale = TRUE)
PLSDA_MPS_H2O2_scores <- data.frame(PLSDA_MPS_H2O2$variates$X, metadata_clr_MPS_H2O2)

PLSDA_MPS_H2O2_plot <- PLSDA_MPS_H2O2_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Attribute_type", alpha = 0.6, title = "PLS-DA CLR - MPS v H2O2",
            xlab = paste("Component 1 (", round(PLSDA_MPS_H2O2$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_MPS_H2O2$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = c("#00BA38", "#619CFF")) +
  geom_point(data = PLSDA_MPS_H2O2_scores %>% group_by(Attribute_type) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Attribute_type), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_MPS_H2O2 <- plotLoadings(PLSDA_MPS_H2O2, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib, importance)

perf_plsda_MPS_H2O2 <- perf(PLSDA_MPS_H2O2, validation = "Mfold", progressBar = TRUE, auc = TRUE, folds = 4, nrepeat = 100) 
plot(perf_plsda_MPS_H2O2, legend = TRUE)

VIPs_MPS_H2O2 <- as.data.frame(mixOmics::vip(PLSDA_MPS_H2O2))
VIPs_MPS_H2O2_filter <- dplyr::filter(VIPs_MPS_H2O2, VIPs_MPS_H2O2$comp1 > 1)
VIPs_MPS_H2O2_filter$ID <- rownames(VIPs_MPS_H2O2_filter)
VIPs_MPS_H2O2_select <- VIPs_MPS_H2O2_filter %>% dplyr::select(ID, comp1)
VIPs_MPS_H2O2_Load <- VIPs_MPS_H2O2_select %>% left_join(Loadings_MPS_H2O2, by = c("ID" = "rowname")) %>% arrange(desc(comp1))
VIPs_MPS_H2O2_Load_Top25 <- VIPs_MPS_H2O2_Load %>% head(., 25)

Loadings_MPS_H2O2_plot <- plotLoadings(PLSDA_MPS_H2O2, plot = FALSE, contrib = "max") %>%
  arrange(GroupContrib) %>%
  rownames_to_column(var = "Feature") %>%
  dplyr::filter(Feature %in% VIPs_MPS_H2O2_Load_Top25$ID) %>%
  ggdotchart(x = "Feature", y = "importance",
             color = "GroupContrib", sorting = "ascending",  
             dot.size = 2, add = "segments",                  
             add.params = list(color = "lightgray", size = 0.1),
             title = "Loadings PLS-DA - MPS v H2O2") + geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

#write_csv(VIPs_MPS_H2O2_Load, "PLSDA_VIP_MPS_H2O2_Full.csv")
#ggsave(plot = PLSDA_MPS_H2O2_plot, filename = "PLSDA_MPS_H2O2_plot.svg", device = "svg", dpi = "retina", width = 4, height = 4)


### Check shared top features between Control vs H2O2 and Control vs MPS

VIP_shared <- VIPs_CT_H2O2_Load_Top25 %>% inner_join(VIPs_CT_MPS_Load_Top25, by = "ID")


### Generate heatmap on the top 25 discriminating metabolites obtained from PLS-DA models ###

data_filter <- data_clr %>% dplyr::select(VIPs_CT_H2O2_Load_Top25$ID | VIPs_CT_MPS_Load_Top25$ID)

#list_wender <- colnames(data_filter) %>% as.data.frame()
#write_csv(list_wender, "list_wender.csv")

# Import CANOPUS output obtained from list of feature of interest
canopus <- read_csv("canopus_compound_summary.csv")
colnames(data_filter) <- canopus$`NPC#class`

hmap <- pheatmap(data_filter, labels_row = metadata$Attribute_type, fontsize_row = 6, fontsize_col = 6, kmeans_k = 3)

#ggsave(filename = "hetmap_class.svg", plot = hmap, device = "svg", dpi = "retina")


# Check total overlapping significant features obtained from PLS-DA models
full_vips <- rbind(VIPs_MPS_H2O2_Load, VIPs_CT_H2O2_Load, VIPs_CT_MPS_Load)
full_vips_unique <- full_vips %>% dplyr::distinct(ID, .keep_all = TRUE)
