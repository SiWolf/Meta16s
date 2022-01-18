# --------------------------------------------------------------------------------------------------------
# Title: Meta16s.R
# Author: Silver A. Wolf
# Last Modified: Tue, 18.01.2022
# Version: 0.2.4
# --------------------------------------------------------------------------------------------------------

# Libraries

library("circlize")
library("ComplexHeatmap")
library("ggpubr")
library("metagMisc")
library("microbiome")
library("microbiomeExplorer")
library("randomcoloR")
library("openxlsx")
library("phyloseq")
library("stringr")

# --------------------------------------------------------------------------------------------------------

# [1] Import and preprocess analysis data

# BIOM File
data.biom <- import_biom("metadata/final_otu_table_clean.biom", parseFunction = parse_taxonomy_default)

# Metadata
meta <- read.csv("metadata/16s_Horses_Overview_Reordered.csv", sep = "\t")
meta.filtered <- meta[meta$Microbiome == "Gut" & meta$Type == "16s" & meta$Comment != "EXCLUDED",]

# Group order
groups.order <- c("SSG", "5DG", "CONTROL")

# --------------------------------------------------------------------------------------------------------

# [2] Diversity Estimations

# Sample Depth
sample_depth <- sort(sample_sums(data.biom))

# Alpha Diversity (Raw)
data.alpha <- microbiome::alpha(data.biom)

# Alpha diversity (Rarefy)
data.rarefy <- rarefy_even_depth(data.biom, rngseed = 1, sample.size = 10000, replace = FALSE, trimOTUs = FALSE)
data.alpha.rarefy <- microbiome::alpha(data.rarefy)
#data.rarefy <- aggregate_top_taxa(data.rarefy, 22, "Rank2")

# Beta Diversity (Bray-Curtis distance)
braycurtis <- phyloseq::distance(data.rarefy, method = "bray")
data.bray <- as.matrix(braycurtis)

# Beta Diversity (PCoA)
braycurtis.pcoa <- ordinate(physeq = data.rarefy, method = "PCoA", distance = "bray")
data.pcoa <- as.data.frame(braycurtis.pcoa$vectors, row.names = NULL, optional = FALSE, cut.names = FALSE, col.names = names(braycurtis.pcoa$vectors), fix.empty.names = TRUE, stringsAsFactors = default.stringsAsFactors())

# Add Metadata
meta.sorted = meta.filtered[match(rownames(data.alpha), meta.filtered$SampleID),]
data.alpha$HORSE = meta.sorted$HorseID
data.alpha$AB_GROUP = factor(meta.sorted$AB_Group, levels = groups.order)
data.alpha$TIMEPOINT = meta.sorted$Timepoint
data.alpha.rarefy$HORSE = meta.sorted[meta.sorted$SampleID %in% rownames(data.alpha.rarefy),]$HorseID
data.alpha.rarefy$AB_GROUP = factor(meta.sorted[meta.sorted$SampleID %in% rownames(data.alpha.rarefy),]$AB_Group, levels = groups.order)
data.alpha.rarefy$TIMEPOINT = meta.sorted[meta.sorted$SampleID %in% rownames(data.alpha.rarefy),]$Timepoint
data.pcoa$HORSE = meta.sorted[meta.sorted$SampleID %in% rownames(data.alpha.rarefy),]$HorseID
data.pcoa$AB_GROUP = factor(meta.sorted[meta.sorted$SampleID %in% rownames(data.alpha.rarefy),]$AB_Group, levels = groups.order)
data.pcoa$TIMEPOINT = meta.sorted[meta.sorted$SampleID %in% rownames(data.alpha.rarefy),]$Timepoint

# Export Diversities
write.csv(data.alpha, file = "results/tab_div_alpha_raw.csv", quote = FALSE)
write.csv(data.alpha.rarefy, file = "results/tab_div_alpha_rarefy.csv", quote = FALSE)
write.csv(data.bray, file = "results/tab_div_beta_distance.csv", quote = FALSE)
write.csv(data.pcoa, file = "results/tab_div_beta_pcoa.csv", quote = FALSE)

# Export human-readable OTU table
data.otu <- phyloseq_to_df(data.rarefy)
write.csv(data.otu, file = "results/tab_otu.csv", row.names = FALSE, quote = FALSE)

# PCA
colours.days = c("t0" = "#00BA38",
                 "t1" = "#F8766D",
                 "t2" = "#619CFF"
                 )

colours.groups = c("SSG" = "#00ff7f",
                   "5DG" = "#ffa500",
                   "CONTROL" = "#00bfff"
                   )

eigenvalue_pc1 = round(braycurtis.pcoa$values$Relative_eig[1]*100, 1)
eigenvalue_pc2 = round(braycurtis.pcoa$values$Relative_eig[2]*100, 1)

# Timepoints - All Samples
png("results/div_pca_time_all.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa,
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.55, 0.25),
          ylim = c(-0.4, 0.5),
          color = "TIMEPOINT",
          shape = "TIMEPOINT",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.days,
          title = "All Samples - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Timepoints - SSG
png("results/div_pca_time_ssg.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$AB_GROUP == "SSG", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.55, 0.25),
          ylim = c(-0.4, 0.5),
          color = "TIMEPOINT",
          shape = "TIMEPOINT",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.days,
          title = "SSG - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Timepoints - 5DG
png("results/div_pca_time_5dg.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$AB_GROUP == "5DG", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.55, 0.25),
          ylim = c(-0.4, 0.5),
          color = "TIMEPOINT",
          shape = "TIMEPOINT",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.days,
          title = "5DG - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Timepoints - Control
png("results/div_pca_time_control.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$AB_GROUP == "CONTROL", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.55, 0.25),
          ylim = c(-0.4, 0.5),
          color = "TIMEPOINT",
          shape = "TIMEPOINT",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.days,
          title = "CONTROL - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Groups - All Samples
png("results/div_pca_group_all.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa,
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.55, 0.25),
          ylim = c(-0.4, 0.5),
          color = "AB_GROUP",
          shape = "AB_GROUP",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.groups,
          title = "All Samples - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Groups - t0
png("results/div_pca_group_t0.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$TIMEPOINT == "t0", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.55, 0.25),
          ylim = c(-0.4, 0.5),
          color = "AB_GROUP",
          shape = "AB_GROUP",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.groups,
          title = "t0 - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Groups - t1
png("results/div_pca_group_t1.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$TIMEPOINT == "t1", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.55, 0.25),
          ylim = c(-0.4, 0.5),
          color = "AB_GROUP",
          shape = "AB_GROUP",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.groups,
          title = "t1 - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Groups - t2
png("results/div_pca_group_t2.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$TIMEPOINT == "t2", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.55, 0.25),
          ylim = c(-0.4, 0.5),
          color = "AB_GROUP",
          shape = "AB_GROUP",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.groups,
          title = "t2 - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Boxplots
boxplot.comparisions <- list(c("t0", "t1"), c("t1", "t2"), c("t0","t2"))

png("results/div_box.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(data.alpha.rarefy, aes(x = TIMEPOINT, y = diversity_shannon, fill = TIMEPOINT)) +
  geom_boxplot() +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(4.5, 7.4)) +
  scale_y_continuous(breaks = c(5, 6, 7)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.comparisions,
                     method = "wilcox.test",
                     label.y = c(6.9, 7.1, 7.3),
                     size = 3,
                     paired = TRUE)
dev.off()

png("results/div_box_even.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(data.alpha.rarefy, aes(x = TIMEPOINT, y = evenness_simpson, fill = TIMEPOINT)) +
  geom_boxplot() +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(0.02, 0.27)) +
  scale_y_continuous(breaks = c(0.05, 0.15, 0.25)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.comparisions,
                     method = "wilcox.test",
                     label.y = c(0.22, 0.24, 0.26),
                     size = 3,
                     paired = TRUE)
dev.off()

# Additional tests

# Significant Differences between SSG and Control (t0) -> yes
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t0" & data.alpha.rarefy$AB_GROUP != "5DG",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# Significant Differences between SSG and Control (t1) -> yes
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t1" & data.alpha.rarefy$AB_GROUP != "5DG",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# Significant Differences between 5DG and Control (t0) -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t0" & data.alpha.rarefy$AB_GROUP != "SSG",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# Significant Differences between 5DG and Control (t1) -> yes
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t1" & data.alpha.rarefy$AB_GROUP != "SSG",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, p.adjust.method = "BH")
stat.res

# Barplots
bar_data_aggregated <- aggregate_top_taxa(data.rarefy, 9, "Rank2")
bar_data_melted <- psmelt(bar_data_aggregated)
num_taxa <- length(unique(bar_data_melted$OTU))
palette <- distinctColorPalette(num_taxa)

colours.phyla = c("Actinobacteria" = "#e2f6f3",
                  "Bacteroidetes" = "#414770",
                  "Candidatus Saccharibacteria" = "#b4d6b4",
                  "Firmicutes" = "#7d97ca",
                  "Other" = "#95f8d6",
                  "Proteobacteria" = "#a7ccf1",
                  "Spirochaetes" = "#f7edc0",
                  "SR1" = "#48a0a1",
                  "Unclassified" = "#d4d7db",
                  "Verrucomicrobia" = "#e4a080"
                  )

bar_ext_horse = c()
bar_ext_group = c()
bar_ext_time = c()
i = 1

for(e in bar_data_melted$Sample){
  current_line = data.alpha.rarefy[rownames(data.alpha.rarefy) == e,]
  bar_ext_horse[i] <- current_line$HORSE
  bar_ext_group[i] <- as.character(current_line$AB_GROUP)
  bar_ext_time[i] <- current_line$TIMEPOINT
  i = i + 1
}

bar_data_melted$HORSE <- bar_ext_horse
bar_data_melted$AB_GROUP <- factor(bar_ext_group, levels = groups.order)
bar_data_melted$TIMEPOINT <- bar_ext_time

# Rename OTUs accordingly
bar_data_melted$OTU[bar_data_melted$OTU == "p__"] <- "Unclassified"
bar_data_melted$OTU <- gsub("p__", "", bar_data_melted$OTU)
#bar_data_melted$OTU <- gsub(" 1", "", bar_data_melted$OTU)
bar_data_melted$OTU <- bar_data_melted$OTU

# Individual Horses - 5DG
bar_data_5dg <- bar_data_melted[bar_data_melted$AB_GROUP == "5DG",]

png("results/tax_bar_horses_5dg.png", width = 30, height = 15, units = "cm", res = 500)
ggplot(bar_data_5dg, aes(fill = OTU, y = Abundance, x = TIMEPOINT)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colours.phyla) +
  labs(title = "Relative Abundance (5DG)",
       x = "Horses",
       y = "Relative Abundance (%)",
       fill = "Top 10 Phyla"
       ) +
  facet_grid(~ HORSE) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        )
dev.off()

# Individual Horses - SSG
bar_data_ssg <- bar_data_melted[bar_data_melted$AB_GROUP == "SSG",]

png("results/tax_bar_horses_ssg.png", width = 30, height = 15, units = "cm", res = 500)
ggplot(bar_data_ssg, aes(fill = OTU, y = Abundance, x = TIMEPOINT)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colours.phyla) +
  labs(title = "Relative Abundance (SSG)",
       x = "Horses",
       y = "Relative Abundance (%)",
       fill = "Top 10 Phyla"
       ) +
  facet_grid(~ HORSE) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        )
dev.off()

# Individual Horses - Control
bar_data_control <- bar_data_melted[bar_data_melted$AB_GROUP == "CONTROL",]

png("results/tax_bar_horses_control.png", width = 30, height = 15, units = "cm", res = 500)
ggplot(bar_data_control, aes(fill = OTU, y = Abundance, x = TIMEPOINT)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colours.phyla) +
  labs(title = "Relative Abundance (Control)",
       x = "Horses",
       y = "Relative Abundance (%)",
       fill = "Top 10 Phyla"
       ) +
  facet_grid(~ HORSE) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        )
dev.off()

# Summarized groups

png("results/tax_bar_sum.png", width = 15, height = 15, units = "cm", res = 500)
ggplot(bar_data_melted, aes(fill = OTU, y = Abundance, x = TIMEPOINT)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colours.phyla) +
  labs(title = "Mean Abundance (Groups)",
       x = "AB Groups",
       y = "Relative Abundance (%)",
       fill = "Top 10 Phyla"
       ) +
  facet_grid(~ AB_GROUP) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        )
dev.off()

# Percentages of individual taxa
per_taxa_total = sum(bar_data_melted$Abundance)/108

per_bact_total = sum(bar_data_melted[bar_data_melted$OTU == "Bacteroidetes", ]$Abundance)/108
per_bact_norm = (per_bact_total/per_taxa_total) * 100

per_firm_total = sum(bar_data_melted[bar_data_melted$OTU == "Firmicutes", ]$Abundance)/108
per_firm_norm = (per_firm_total/per_taxa_total) * 100

per_prot_total = sum(bar_data_melted[bar_data_melted$OTU == "Proteobacteria", ]$Abundance)/108
per_prot_norm = (per_prot_total/per_taxa_total) * 100

per_spir_total = sum(bar_data_melted[bar_data_melted$OTU == "Spirochaetes", ]$Abundance)/108
per_spir_norm = (per_spir_total/per_taxa_total) * 100

per_verr_total = sum(bar_data_melted[bar_data_melted$OTU == "Verrucomicrobia", ]$Abundance)/108
per_verr_norm = (per_verr_total/per_taxa_total) * 100

per_bact_norm
per_firm_norm
per_prot_norm
per_spir_norm
per_verr_norm

# Abundance heatmap

abundance_matrix <- matrix(0, 9, length(unique(bar_data_melted$OTU)))
abundance_columns <- unique(bar_data_melted$OTU)
abundance_rows <- c("SSG_t0", "5DG_t0", "CONTROL_t0", "SSG_t1", "5DG_t1", "CONTROL_t1", "SSG_t2", "5DG_t2", "CONTROL_t2")

j = 1

for (c in abundance_columns){
  i = 1
  for (r in abundance_rows){
    group <- strsplit(r, "_")[[1]][1]
    timepoint <- strsplit(r, "_")[[1]][2]
    abundance_matrix[i, j] <- mean(bar_data_melted[bar_data_melted$TIMEPOINT == timepoint & bar_data_melted$AB_GROUP == group & bar_data_melted$OTU == c,]$Abundance)
    i = i + 1
  }
  j = j + 1
}

abundance_matrix <- t(abundance_matrix)

rownames(abundance_matrix) <- abundance_columns
colnames(abundance_matrix) <- c("SSG", "5DG", "CONTROL", "SSG", "5DG", "CONTROL", "SSG", "5DG", "CONTROL")

png("results/tax_heatmap.png", width = 30, height = 20, units = "cm", res = 500)
Heatmap(log2(abundance_matrix + 1),
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        #col = c("black", "darkred", "red", "orange", "yellow"),
        col = c("black", "red"),
        column_split = c(rep("t0", 3), rep("t1", 3), rep("t2", 3)),
        row_title = "Top 10 Phyla",
        name = "log2(Abundance)",
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10)
        )
dev.off()

# Boxplots for abundance of specific families

boxplot.families.biom <- aggregate_top_taxa(data.rarefy, 22, "Rank5")
boxplot.families.melted <- psmelt(boxplot.families.biom)

box.ext.horse = c()
box.ext.group = c()
box.ext.time = c()
j = 1

for(f in boxplot.families.melted$Sample){
  current.line = data.alpha.rarefy[rownames(data.alpha.rarefy) == f,]
  box.ext.horse[j] <- current.line$HORSE
  box.ext.group[j] <- as.character(current.line$AB_GROUP)
  box.ext.time[j] <- current.line$TIMEPOINT
  j = j + 1
}

boxplot.families.melted$HORSE <- box.ext.horse
boxplot.families.melted$AB_GROUP <- factor(box.ext.group, levels = groups.order)
boxplot.families.melted$TIMEPOINT <- box.ext.time

G1 <- boxplot.families.melted[boxplot.families.melted$Rank5 == "f__Ruminococcaceae", ]
png("results/tax_box_f_ruminococcaceae.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(G1, aes(x = TIMEPOINT, y = log2(Abundance + 1), fill = TIMEPOINT)) +
  geom_boxplot() +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(0, 15.5)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.comparisions,
                     method = "wilcox.test",
                     label.y = c(12.5, 13.5, 14.5),
                     size = 3,
                     paired = TRUE,
                     method.args = list(exact = FALSE)) +
  ggtitle("Abundance - Ruminococcaceae (Family)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

G2 <- boxplot.families.melted[boxplot.families.melted$Rank5 == "f__Enterobacteriaceae", ]
png("results/tax_box_f_enterobacteriaceae.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(G2, aes(x = TIMEPOINT, y = log2(Abundance + 1), fill = TIMEPOINT)) +
  geom_boxplot() +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(0, 15.5)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.comparisions,
                     method = "wilcox.test",
                     label.y = c(12.5, 13.5, 14.5),
                     size = 3,
                     paired = TRUE,
                     method.args = list(exact = FALSE)) +
  ggtitle("Abundance - Enterobacteriaceae (Family)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Calculate differences between timepoints

# Alpha diversity
alpha.5dg.t0 = median(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t0", ]$diversity_shannon)
alpha.5dg.t1 = median(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t1", ]$diversity_shannon)
alpha.5dg.t2 = median(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t2", ]$diversity_shannon)

alpha.5dg.diff.t0.t1 = round(alpha.5dg.t1 - alpha.5dg.t0, 2)
alpha.5dg.diff.t1.t2 = round(alpha.5dg.t2 - alpha.5dg.t1, 2)

alpha.ssg.t0 = median(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t0", ]$diversity_shannon)
alpha.ssg.t1 = median(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t1", ]$diversity_shannon)
alpha.ssg.t2 = median(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t2", ]$diversity_shannon)

alpha.ssg.diff.t0.t1 = round(alpha.ssg.t1 - alpha.ssg.t0, 2)
alpha.ssg.diff.t1.t2 = round(alpha.ssg.t2 - alpha.ssg.t1, 2)

# Abundance of Enterobacteriaceae
entero.5dg.t0 = log2(median(G2[G2$AB_GROUP == "5DG" & G2$TIMEPOINT == "t0", ]$Abundance))
entero.5dg.t1 = log2(median(G2[G2$AB_GROUP == "5DG" & G2$TIMEPOINT == "t1", ]$Abundance))
entero.5dg.t2 = log2(median(G2[G2$AB_GROUP == "5DG" & G2$TIMEPOINT == "t2", ]$Abundance))

entero.5dg.diff.t0.t1 = round(entero.5dg.t1 - entero.5dg.t0, 2)
entero.5dg.diff.t1.t2 = round(entero.5dg.t2 - entero.5dg.t1, 2)

entero.ssg.t0 = log2(median(G2[G2$AB_GROUP == "SSG" & G2$TIMEPOINT == "t0", ]$Abundance))
entero.ssg.t1 = log2(median(G2[G2$AB_GROUP == "SSG" & G2$TIMEPOINT == "t1", ]$Abundance))
entero.ssg.t2 = log2(median(G2[G2$AB_GROUP == "SSG" & G2$TIMEPOINT == "t2", ]$Abundance))

entero.ssg.diff.t0.t1 = round(entero.ssg.t1 - entero.ssg.t0, 2)
entero.ssg.diff.t1.t2 = round(entero.ssg.t2 - entero.ssg.t1, 2)

# --------------------------------------------------------------------------------------------------------

# [03] Run MicrobiomeExplorer (16s/WGS)

converted_biom <- readData(filepath = "metadata/final_otu_table_clean.biom", type = "BIOM")
saveRDS(converted_biom, "results/kraken2.rds")

#runMicrobiomeExplorer()