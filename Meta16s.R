# --------------------------------------------------------------------------------------------------------
# Title: Meta16s.R
# Author: Silver A. Wolf
# Last Modified: Fri, 10.06.2022
# Version: 0.3.5
# --------------------------------------------------------------------------------------------------------

# Libraries

library("circlize")
library("ComplexHeatmap")
library("dplyr")
library("ggplot2")
library("ggpubr")
library("kableExtra")
library("knitr")
library("metagMisc")
library("microbiome")
library("microbiomeExplorer")
library("micropower")
library("randomcoloR")
library("openxlsx")
library("parallel")
library("phyloseq")
library("sp")
library("stringr")
library("vegan")

# --------------------------------------------------------------------------------------------------------

# [1] Import and preprocess analysis data

# BIOM File
data.biom <- import_biom("output/final_otu_table_clean.biom", parseFunction = parse_taxonomy_default)

# Metadata
meta <- read.xlsx("metadata/22_02_Horses_Overview.xlsx", sheet = 1)
meta.filtered <- meta[meta$AB_Group != "SWITCHED" & !(meta$AB_Group == "CONTROL" & meta$Timepoint == "t2"), ]

# Group order
groups.order <- c("SSG", "5DG", "CONTROL")

# --------------------------------------------------------------------------------------------------------

# [2] Diversity Estimations

# Sample Depth
sample_depth <- sort(sample_sums(data.biom))
min(sample_depth)
max(sample_depth)
median(sample_depth)
sum(sample_depth)

# Alpha Diversity (Raw)
data.alpha <- microbiome::alpha(data.biom)

# Plot Rarefaction Curve
png("results/div_rarefaction_curve.png", width = 16, height = 16, units = "cm", res = 500)
rarecurve(t(otu_table(data.biom)), step = 50, cex = 0.5)
dev.off()

# Alpha diversity (Rarefy)
data.rarefy <- rarefy_even_depth(data.biom, rngseed = 1, sample.size = min(sample_depth), replace = FALSE, trimOTUs = FALSE)
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
data.otu <- phyloseq_to_df(data.biom)
write.csv(data.otu, file = "results/tab_otu.csv", row.names = FALSE, quote = FALSE)
data.otu.rarefy <- phyloseq_to_df(data.rarefy)
write.csv(data.otu.rarefy, file = "results/tab_otu_rarefy.csv", row.names = FALSE, quote = FALSE)

# Sort OTU table by metadata
otu.tmp.1 <- data.otu.rarefy[,1:7]
otu.tmp.2 <- data.otu.rarefy[,9:ncol(data.otu.rarefy)]
otu.tmp.3 <- meta.sorted[with(meta.sorted, order(AB_Group, HorseID, Day)), ]
otu.tmp.4 <- otu.tmp.2[,otu.tmp.3$SampleID]
otu.tmp.5 <- data.frame(otu.tmp.1, otu.tmp.4)

otu.row.1 <- c(rep("", 7), otu.tmp.3$HorseID)
otu.row.2 <- c(rep("", 7), otu.tmp.3$Name)
otu.row.3 <- c(rep("", 7), otu.tmp.3$AB_Group)
otu.row.4 <- c(rep("", 7), otu.tmp.3$Timepoint)

otu.tmp.6 <- rbind(colnames(otu.tmp.5), otu.tmp.5)
otu.tmp.6 <- rbind(otu.row.4, otu.tmp.6)
otu.tmp.6 <- rbind(otu.row.3, otu.tmp.6)
otu.tmp.6 <- rbind(otu.row.2, otu.tmp.6)
otu.tmp.6 <- rbind(otu.row.1, otu.tmp.6)

write.table(otu.tmp.6, file = "results/tab_otu_rarefy_reordered.csv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

# Total OTUs
nrow(data.otu)

# Unique Genus Level Taxa
length(unique(data.otu$Rank6))

# Diversity values per group
max(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)
min(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)
median(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)
mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)

mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t1",]$diversity_shannon)
mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "SSG" & data.alpha.rarefy$TIMEPOINT == "t2",]$diversity_shannon)

max(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)
min(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)
median(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)
mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)

mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t1",]$diversity_shannon)
mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "5DG" & data.alpha.rarefy$TIMEPOINT == "t2",]$diversity_shannon)

max(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "CONTROL" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)
min(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "CONTROL" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)
median(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "CONTROL" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)
mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "CONTROL" & data.alpha.rarefy$TIMEPOINT == "t0",]$diversity_shannon)

mean(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP == "CONTROL" & data.alpha.rarefy$TIMEPOINT == "t1",]$diversity_shannon)

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
          xlim = c(-0.5, 0.25),
          ylim = c(-0.25, 0.4),
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

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t0", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t1", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t2", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

# Timepoints - SSG
png("results/div_pca_time_ssg.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$AB_GROUP == "SSG", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.5, 0.25),
          ylim = c(-0.25, 0.4),
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

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t0" & data.pcoa$AB_GROUP == "SSG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t1" & data.pcoa$AB_GROUP == "SSG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t2" & data.pcoa$AB_GROUP == "SSG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

# Timepoints - 5DG
png("results/div_pca_time_5dg.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$AB_GROUP == "5DG", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.5, 0.25),
          ylim = c(-0.25, 0.4),
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

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t0" & data.pcoa$AB_GROUP == "5DG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t1" & data.pcoa$AB_GROUP == "5DG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t2" & data.pcoa$AB_GROUP == "5DG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

# Timepoints - Control
png("results/div_pca_time_control.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$AB_GROUP == "CONTROL", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.5, 0.25),
          ylim = c(-0.25, 0.4),
          color = "TIMEPOINT",
          shape = "TIMEPOINT",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.days[1:2],
          title = "CONTROL - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
dev.off()

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t0" & data.pcoa$AB_GROUP == "CONTROL", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t1" & data.pcoa$AB_GROUP == "CONTROL", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

# Groups - All Samples
png("results/div_pca_group_all.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa,
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.5, 0.25),
          ylim = c(-0.25, 0.4),
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

chull.raw <- data.pcoa[data.pcoa$AB_GROUP == "SSG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$AB_GROUP == "5DG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$AB_GROUP == "CONTROL", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

# Groups - t0
png("results/div_pca_group_t0.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$TIMEPOINT == "t0", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.5, 0.25),
          ylim = c(-0.25, 0.4),
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

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t0" & data.pcoa$AB_GROUP == "SSG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t0" & data.pcoa$AB_GROUP == "5DG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t0" & data.pcoa$AB_GROUP == "CONTROL", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

# Groups - t1
png("results/div_pca_group_t1.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$TIMEPOINT == "t1", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.5, 0.25),
          ylim = c(-0.25, 0.4),
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

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t1" & data.pcoa$AB_GROUP == "SSG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t1" & data.pcoa$AB_GROUP == "5DG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t1" & data.pcoa$AB_GROUP == "CONTROL", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

# Groups - t2
png("results/div_pca_group_t2.png", width = 16, height = 16, units = "cm", res = 500)
ggscatter(data.pcoa[data.pcoa$TIMEPOINT == "t2", ],
          x = "Axis.1",
          y = "Axis.2",
          xlab = paste("PC1 (", eigenvalue_pc1, "%)", sep = ""),
          ylab = paste("PC2 (", eigenvalue_pc2, "%)", sep = ""),
          xlim = c(-0.5, 0.25),
          ylim = c(-0.25, 0.4),
          color = "AB_GROUP",
          shape = "AB_GROUP",
          #star.plot = TRUE,
          #mean.point = TRUE,
          ellipse = TRUE,
          ellipse.alpha = 0.3,
          ellipse.border.remove = TRUE,
          ellipse.type = "convex",
          palette = colours.groups[1:2],
          title = "t2 - Beta Diversity PCA (Bray–Curtis Dissimilarity)"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
dev.off()

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t2" & data.pcoa$AB_GROUP == "SSG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t2" & data.pcoa$AB_GROUP == "5DG", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

chull.raw <- data.pcoa[data.pcoa$TIMEPOINT == "t2" & data.pcoa$AB_GROUP == "CONTROL", 1:2]
chull.hpts <- chull(x = chull.raw$Axis.1, y = chull.raw$Axis.2)
chull.hpts <- c(chull.hpts, chull.hpts[1])
chull.coords <- chull.raw[chull.hpts,]
chull.poly <- Polygon(chull.coords, hole = F)
chull.area <- chull.poly@area
chull.area

# Boxplots (Groupwise)
boxplot.groups <- list(c("SSG", "5DG"))

png("results/div_box_group_shan.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP != "CONTROL",], aes(x = AB_GROUP, y = diversity_shannon, fill = AB_GROUP)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~TIMEPOINT, scale = "free") +
  coord_cartesian(ylim = c(4.6, 7.4)) +
  scale_y_continuous(breaks = c(5, 6, 7)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.groups) +
  stat_compare_means(comparisons = boxplot.groups,
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(6.9, 7.1, 7.3),
                     size = 3,
                     paired = FALSE)
dev.off()

png("results/div_box_group_even.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP != "CONTROL",], aes(x = AB_GROUP, y = evenness_simpson, fill = AB_GROUP)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~TIMEPOINT, scale = "free") +
  coord_cartesian(ylim = c(0.02, 0.27)) +
  scale_y_continuous(breaks = c(0.05, 0.15, 0.25)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.groups) +
  stat_compare_means(comparisons = boxplot.groups,
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(0.22, 0.24, 0.26),
                     size = 3,
                     paired = FALSE)
dev.off()

# Boxplots (Timewise)
boxplot.timepoints <- list(c("t0", "t1"), c("t1", "t2"), c("t0","t2"))

png("results/div_box_time_shan_group.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP != "CONTROL",], aes(x = TIMEPOINT, y = diversity_shannon, fill = TIMEPOINT)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(4.6, 7.4)) +
  scale_y_continuous(breaks = c(5, 6, 7)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.timepoints,
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(6.9, 7.1, 7.3),
                     size = 3,
                     paired = TRUE)
dev.off()

png("results/div_box_time_shan_control.png", width = 10, height = 10, units = "cm", res = 500)
ggplot(data.alpha.rarefy[(data.alpha.rarefy$AB_GROUP == "CONTROL") & (data.alpha.rarefy$TIMEPOINT != "t2"),], aes(x = TIMEPOINT, y = diversity_shannon, fill = TIMEPOINT)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(4.6, 7.4)) +
  scale_y_continuous(breaks = c(5, 6, 7)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.timepoints[1],
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(6.9, 7.1, 7.3),
                     size = 3,
                     paired = TRUE)
dev.off()

png("results/div_box_time_even_group.png", width = 20, height = 10, units = "cm", res = 500)
ggplot(data.alpha.rarefy[data.alpha.rarefy$AB_GROUP != "CONTROL",], aes(x = TIMEPOINT, y = evenness_simpson, fill = TIMEPOINT)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(0.02, 0.27)) +
  scale_y_continuous(breaks = c(0.05, 0.15, 0.25)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.timepoints,
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(0.22, 0.24, 0.26),
                     size = 3,
                     paired = TRUE)
dev.off()

png("results/div_box_time_even_control.png", width = 10, height = 10, units = "cm", res = 500)
ggplot(data.alpha.rarefy[(data.alpha.rarefy$AB_GROUP == "CONTROL") & (data.alpha.rarefy$TIMEPOINT != "t2"),], aes(x = TIMEPOINT, y = evenness_simpson, fill = TIMEPOINT)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(0.02, 0.27)) +
  scale_y_continuous(breaks = c(0.05, 0.15, 0.25)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.timepoints[1],
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(0.22, 0.24, 0.26),
                     size = 3,
                     paired = TRUE)
dev.off()

# Additional tests

# Significant Differences between SSG and Control (t0) -> yes
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t0" & data.alpha.rarefy$AB_GROUP != "5DG",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# Significant Differences between SSG and Control (t1) -> yes
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t1" & data.alpha.rarefy$AB_GROUP != "5DG",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# Significant Differences between 5DG and Control (t0) -> no
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t0" & data.alpha.rarefy$AB_GROUP != "SSG",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# Significant Differences between 5DG and Control (t1) -> yes
stat.data <- data.alpha.rarefy[data.alpha.rarefy$TIMEPOINT == "t1" & data.alpha.rarefy$AB_GROUP != "SSG",]
stat.df <- data.frame(GROUP = stat.data$AB_GROUP, AMR = stat.data$diversity_shannon)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$AMR), stat.df$GROUP, alternative = "two.sided", p.adjust.method = "BH")
stat.res

# Barplots
bar_data_aggregated <- aggregate_top_taxa(data.rarefy, 9, "Rank2")
bar_data_melted <- psmelt(bar_data_aggregated)
num_taxa <- length(unique(bar_data_melted$OTU))
palette <- distinctColorPalette(num_taxa)

colours.phyla = c("Actinomycetota" = "#e2f6f3",
                  "Armatimonadota" = "#48a0a1",
                  "Bacillota" = "#7d97ca",
                  "Bacteroidota" = "#414770",
                  "Candidatus Saccharibacteria" = "#b4d6b4",
                  "Other" = "#95f8d6",
                  "Pseudomonadota" = "#a7ccf1",
                  "Spirochaetota" = "#f7edc0",
                  "Unclassified" = "#d4d7db",
                  "Verrucomicrobiota" = "#e4a080"
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

# Summarized groups (5DG/SSG)

png("results/tax_bar_sum_5dg_ssg.png", width = 15, height = 15, units = "cm", res = 500)
ggplot(bar_data_melted[bar_data_melted$AB_GROUP != "CONTROL", ], aes(fill = OTU, y = Abundance, x = TIMEPOINT)) + 
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

# Summarized groups (CONTROL)

png("results/tax_bar_sum_control.png", width = 10, height = 15, units = "cm", res = 500)
ggplot(bar_data_melted[bar_data_melted$AB_GROUP == "CONTROL", ], aes(fill = OTU, y = Abundance, x = TIMEPOINT)) + 
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
per_taxa_total = sum(bar_data_melted$Abundance)/length(sample_depth)

per_bact_total = sum(bar_data_melted[bar_data_melted$OTU == "Bacteroidota", ]$Abundance)/length(sample_depth)
per_bact_norm = (per_bact_total/per_taxa_total) * 100

per_firm_total = sum(bar_data_melted[bar_data_melted$OTU == "Bacillota", ]$Abundance)/length(sample_depth)
per_firm_norm = (per_firm_total/per_taxa_total) * 100

per_prot_total = sum(bar_data_melted[bar_data_melted$OTU == "Pseudomonadota", ]$Abundance)/length(sample_depth)
per_prot_norm = (per_prot_total/per_taxa_total) * 100

per_spir_total = sum(bar_data_melted[bar_data_melted$OTU == "Spirochaetota", ]$Abundance)/length(sample_depth)
per_spir_norm = (per_spir_total/per_taxa_total) * 100

per_verr_total = sum(bar_data_melted[bar_data_melted$OTU == "Verrucomicrobiota", ]$Abundance)/length(sample_depth)
per_verr_norm = (per_verr_total/per_taxa_total) * 100

per_bact_norm
per_firm_norm
per_prot_norm
per_spir_norm
per_verr_norm

# Abundance heatmap

abundance_matrix <- matrix(0, 8, length(unique(bar_data_melted$OTU)))
abundance_columns <- unique(bar_data_melted$OTU)
abundance_rows <- c("SSG_t0", "5DG_t0", "SSG_t1", "5DG_t1", "SSG_t2", "5DG_t2", "CONTROL_t0", "CONTROL_t1")

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
colnames(abundance_matrix) <- c("SSG", "5DG", "SSG", "5DG", "SSG", "5DG", "t0", "t1")

png("results/tax_heatmap.png", width = 30, height = 20, units = "cm", res = 500)
Heatmap(log2(abundance_matrix + 1),
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        #col = c("black", "darkred", "red", "orange", "yellow"),
        #col = c("black", "red"),
        col = c("grey", "orange", "red", "darkred"),
        column_split = c(rep("t0", 2), rep("t1", 2), rep("t2", 2), rep("CONTROL", 2)),
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

png("results/tax_box_f_ruminococcaceae_time_group.png", width = 18, height = 10, units = "cm", res = 500)
ggplot(G1[G1$AB_GROUP != "CONTROL",], aes(x = TIMEPOINT, y = log2(Abundance + 1), fill = TIMEPOINT)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(0, 15.5)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.timepoints,
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(12.5, 13.5, 14.5),
                     size = 3,
                     paired = TRUE,
                     method.args = list(exact = FALSE)) +
  ggtitle("Abundance - Ruminococcaceae (Family)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("results/tax_box_f_ruminococcaceae_time_control.png", width = 11, height = 10, units = "cm", res = 500)
ggplot(G1[G1$AB_GROUP == "CONTROL",], aes(x = TIMEPOINT, y = log2(Abundance + 1), fill = TIMEPOINT)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(0, 15.5)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.timepoints[1],
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(12.5, 13.5, 14.5),
                     size = 3,
                     paired = TRUE,
                     method.args = list(exact = FALSE)) +
  ggtitle("Abundance - Ruminococcaceae (Family)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("results/tax_box_f_ruminococcaceae_groups.png", width = 18, height = 10, units = "cm", res = 500)
ggplot(G1[G1$AB_GROUP != "CONTROL",], aes(x = AB_GROUP, y = log2(Abundance + 1), fill = AB_GROUP)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~TIMEPOINT, scale = "free") +
  coord_cartesian(ylim = c(0, 15.5)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.groups) +
  stat_compare_means(comparisons = boxplot.groups,
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(12.5, 13.5, 14.5),
                     size = 3,
                     paired = FALSE,
                     method.args = list(exact = FALSE)) +
  ggtitle("Abundance - Ruminococcaceae (Family)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

G2 <- boxplot.families.melted[boxplot.families.melted$Rank5 == "f__Enterobacteriaceae", ]

png("results/tax_box_f_enterobacteriaceae_time_group.png", width = 18, height = 10, units = "cm", res = 500)
ggplot(G2[G2$AB_GROUP != "CONTROL",], aes(x = TIMEPOINT, y = log2(Abundance + 1), fill = TIMEPOINT)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(0, 15.5)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.timepoints,
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(12.5, 13.5, 14.5),
                     size = 3,
                     paired = TRUE,
                     method.args = list(exact = FALSE)) +
  ggtitle("Abundance - Enterobacteriaceae (Family)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("results/tax_box_f_enterobacteriaceae_time_control.png", width = 11, height = 10, units = "cm", res = 500)
ggplot(G2[G2$AB_GROUP == "CONTROL",], aes(x = TIMEPOINT, y = log2(Abundance + 1), fill = TIMEPOINT)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~AB_GROUP, scale = "free") +
  coord_cartesian(ylim = c(0, 15.5)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.days) +
  stat_compare_means(comparisons = boxplot.timepoints[1],
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(12.5, 13.5, 14.5),
                     size = 3,
                     paired = TRUE,
                     method.args = list(exact = FALSE)) +
  ggtitle("Abundance - Enterobacteriaceae (Family)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("results/tax_box_f_enterobacteriaceae_groups.png", width = 18, height = 10, units = "cm", res = 500)
ggplot(G2[G2$AB_GROUP != "CONTROL",], aes(x = AB_GROUP, y = log2(Abundance + 1), fill = AB_GROUP)) +
  geom_boxplot(alpha = 0.9) +
  facet_wrap(~TIMEPOINT, scale = "free") +
  coord_cartesian(ylim = c(0, 15.5)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15)) +
  stat_boxplot(geom = "errorbar", width = 0.5) +
  scale_fill_manual(values = colours.groups) +
  stat_compare_means(comparisons = boxplot.groups,
                     alternative = "two.sided",
                     method = "wilcox.test",
                     label.y = c(12.5, 13.5, 14.5),
                     size = 3,
                     paired = FALSE,
                     method.args = list(exact = FALSE)) +
  ggtitle("Abundance - Enterobacteriaceae (Family)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Analyses for abundance of genus

boxplot.genus.biom <- aggregate_top_taxa(data.rarefy, 22, "Rank6")
boxplot.genus.melted <- psmelt(boxplot.genus.biom)

box.ext.horse = c()
box.ext.group = c()
box.ext.time = c()
j = 1

for(f in boxplot.genus.melted$Sample){
  current.line = data.alpha.rarefy[rownames(data.alpha.rarefy) == f,]
  box.ext.horse[j] <- current.line$HORSE
  box.ext.group[j] <- as.character(current.line$AB_GROUP)
  box.ext.time[j] <- current.line$TIMEPOINT
  j = j + 1
}

boxplot.genus.melted$HORSE <- box.ext.horse
boxplot.genus.melted$AB_GROUP <- factor(box.ext.group, levels = groups.order)
boxplot.genus.melted$TIMEPOINT <- box.ext.time

G3 <- boxplot.genus.melted[boxplot.genus.melted$Rank6 == "g__Escherichia/Shigella", ]

stat.data <- G3[G3$TIMEPOINT != "t2" & G3$AB_GROUP == "5DG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, ABU = stat.data$Abundance)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$ABU), stat.df$GROUP, alternative = "greater", p.adjust.method = "BH", paired = TRUE, exact = FALSE)
stat.res

stat.data <- G3[G3$TIMEPOINT != "t0" & G3$AB_GROUP == "5DG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, ABU = stat.data$Abundance)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$ABU), stat.df$GROUP, alternative = "less", p.adjust.method = "BH", paired = TRUE, exact = FALSE)
stat.res

stat.data <- G3[G3$TIMEPOINT != "t2" & G3$AB_GROUP == "SSG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, ABU = stat.data$Abundance)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$ABU), stat.df$GROUP, alternative = "greater", p.adjust.method = "BH", paired = TRUE, exact = FALSE)
stat.res

stat.data <- G3[G3$TIMEPOINT != "t0" & G3$AB_GROUP == "SSG",]
stat.df <- data.frame(GROUP = stat.data$TIMEPOINT, ABU = stat.data$Abundance)
stat.res <- pairwise.wilcox.test(as.numeric(stat.df$ABU), stat.df$GROUP, alternative = "less", p.adjust.method = "BH", paired = TRUE, exact = FALSE)
stat.res

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

converted_biom <- readData(filepath = "output/final_otu_table_clean.biom", type = "BIOM")
saveRDS(converted_biom, "results/kraken2.rds")

#runMicrobiomeExplorer()

# --------------------------------------------------------------------------------------------------------

# [04] Power Calculation (16s)

#power.cores = detectCores()
#power.depth = 1000
#power.matrix <- as(otu_table(data.biom), "matrix")

#jaccard <- calcWJstudy(power.matrix)
#m <- mean(jaccard)
#s <- sd(jaccard)

#weight_hm <- hashMean(rare_levels = runif(100, 0, 1), rep_per_level = 10, otu_number = power.depth, sequence_depth = power.depth)
#means <- mclapply(weight_hm, function(x) (mean(lowerTriDM(calcWJstudy(x)))), mc.cores = power.cores)
#names(means) <- sapply(strsplit(names(means), "_"), FUN = function(x) {x[[1]]})
#mean_df <- data.frame(subsampling = as.numeric(names(means)), means = as.numeric(means), target = abs(as.numeric(means) - m))
#subsample <- mean_df[which(mean_df$target == min(mean_df$target)), ]$subsampling
#png("results/power_mean_subsample.png", width = 16, height = 16, units = "cm", res = 500)
#qplot(data = mean_df, x = subsampling, y = means)
#dev.off()

#weight_hsd <- hashSD(rare_depth = subsample, otu_number_range = 10^runif(n = 100, min = 1, max = 5), sim_number = 30, sequence_depth = power.depth)
#sds <- mclapply(weight_hsd, function(x) (sd(lowerTriDM(calcWJstudy(x)))), mc.cores = power.cores)
#names(sds) <- sapply(names(sds), function(x) substring(x, 4))
#sds_df <- data.frame(otunum = as.numeric(names(sds)), sd = as.numeric(sds), target = abs(as.numeric(sds) - s))
#otunum <- sds_df[which(sds_df$target == min(sds_df$target)), ]$otunum
#png("results/power_sd_subsample.png", width = 16, height = 16, units = "cm", res = 500)
#qplot(data = sds_df, x = otunum, y = sd) + geom_smooth(method = "lm") + scale_x_log10() + scale_y_log10() + xlab("log10(otunum)") + ylab("log10(sd)")
#dev.off()

#sp <- simPower(group_size_vector = c(10, 10), otu_number = otunum, rare_depth = subsample, sequence_depth = power.depth, effect_range = seq(0, 0.3, length.out = 100))
#wj <- mclapply(sp, function(x) calcWJstudy(x), mc.cores = power.cores)
#mean(as.numeric(mclapply(wj, mean, mc.cores = power.cores)))
#mean(as.numeric(mclapply(wj, sd, mc.cores = power.cores)))

#bp30 <- bootPower(wj, boot_number = 10, subject_group_vector = c(10, 10), alpha = 0.05)
#bp30_model <- subset(bp30, power < 0.95 & power > 0.2)
#bp30_model <- data.frame(log_omega2 = log10(bp30_model$simulated_omega2), log_power = log10(bp30_model$power))
#bp30_model <- subset(bp30_model, log_omega2 > -Inf)
#bp30_lm <- lm(log_omega2 ~ log_power, data = bp30_model)
#power80 <- 10^predict(bp30_lm, newdata = data.frame(log_power = log10(0.8)))
#power90 <- 10^predict(bp30_lm, newdata = data.frame(log_power = log10(0.9)))
#png("results/power_linear_model.png", width = 16, height = 16, units = "cm", res = 500)
#ggplot2::qplot(data = bp30_model, x = log_power, y = log_omega2) +
#  geom_smooth(method = "lm") +
#  ggtitle("log10(Omega2) to log10(Power) with 10 Subjects/Group") +
#  xlab("log10(Power)") +
#  ylab("log10(Omega2)")
#dev.off()

# --------------------------------------------------------------------------------------------------------

# [05] Calculation of specific taxa abundances

samples.ssg.t0 <- meta.sorted[meta.sorted$AB_Group == "SSG" & meta.sorted$Timepoint == "t0", ]$SampleID
samples.ssg.t1 <- meta.sorted[meta.sorted$AB_Group == "SSG" & meta.sorted$Timepoint == "t1", ]$SampleID
samples.ssg.t2 <- meta.sorted[meta.sorted$AB_Group == "SSG" & meta.sorted$Timepoint == "t2", ]$SampleID

samples.5dg.t0 <- meta.sorted[meta.sorted$AB_Group == "5DG" & meta.sorted$Timepoint == "t0", ]$SampleID
samples.5dg.t1 <- meta.sorted[meta.sorted$AB_Group == "5DG" & meta.sorted$Timepoint == "t1", ]$SampleID
samples.5dg.t2 <- meta.sorted[meta.sorted$AB_Group == "5DG" & meta.sorted$Timepoint == "t2", ]$SampleID

samples.control.t0 <- meta.sorted[meta.sorted$AB_Group == "CONTROL" & meta.sorted$Timepoint == "t0", ]$SampleID
samples.control.t1 <- meta.sorted[meta.sorted$AB_Group == "CONTROL" & meta.sorted$Timepoint == "t1", ]$SampleID

samples.taxa <- c("p__Actinomycetota",
                  "p__Armatimonadota",
                  "p__Bacillota",
                  "p__Bacteroidota",
                  "p__Candidatus Saccharibacteria",
                  "p__Pseudomonadota",
                  "p__Spirochaetota",
                  "p__SR1",
                  "p__Verrucomicrobiota",
                  "c__Subdivision5",
                  "g__Allisonella",
                  "g__Escherichia/Shigella",
                  "g__Fusobacterium",
                  "g__Lactobacillus",
                  "g__Phascolarctobacterium",
                  "g__Rhodococcus",
                  "g__Roseburia",
                  "g__Ruminococcus",
                  "g__Salmonella",
                  "g__Staphylococcus",
                  "g__Streptococcus",
                  "f__Acidaminococcaceae",
                  "f__Bacteroidaceae",
                  "f__Enterobacteriaceae",
                  "f__Lachnospiraceae",
                  "f__Moraxellaceae",
                  "f__Planococcaceae",
                  "f__Ruminococcaceae",
                  "f__Veillonellaceae"
                  )

c1 = c()
c2 = c()
c3 = c()
c4 = c()
c5 = c()
c6 = c()
c7 = c()
c8 = c()
x = 1

for (t in samples.taxa){
  t1 <- data.otu.rarefy[data.otu.rarefy$Rank1 == t |
                        data.otu.rarefy$Rank2 == t |
                        data.otu.rarefy$Rank3 == t |
                        data.otu.rarefy$Rank4 == t |
                        data.otu.rarefy$Rank5 == t |
                        data.otu.rarefy$Rank6 == t |
                        data.otu.rarefy$Rank7 == t, ]
  
  t2 <- t1[samples.ssg.t0]
  t3 <- colSums(t2)
  t4 <- sum(t3)/length(samples.ssg.t0)
  t5 <- (t4/min(sample_depth))*100
  t6 <- round(t5, 2)
  c1[x] <- t6
  
  t2 <- t1[samples.ssg.t1]
  t3 <- colSums(t2)
  t4 <- sum(t3)/length(samples.ssg.t1)
  t5 <- (t4/min(sample_depth))*100
  t6 <- round(t5, 2)
  c2[x] <- t6
  
  t2 <- t1[samples.ssg.t2]
  t3 <- colSums(t2)
  t4 <- sum(t3)/length(samples.ssg.t2)
  t5 <- (t4/min(sample_depth))*100
  t6 <- round(t5, 2)
  c3[x] <- t6
  
  t2 <- t1[samples.5dg.t0]
  t3 <- colSums(t2)
  t4 <- sum(t3)/length(samples.5dg.t0)
  t5 <- (t4/min(sample_depth))*100
  t6 <- round(t5, 2)
  c4[x] <- t6
  
  t2 <- t1[samples.5dg.t1]
  t3 <- colSums(t2)
  t4 <- sum(t3)/length(samples.5dg.t1)
  t5 <- (t4/min(sample_depth))*100
  t6 <- round(t5, 2)
  c5[x] <- t6
  
  t2 <- t1[samples.5dg.t2]
  t3 <- colSums(t2)
  t4 <- sum(t3)/length(samples.5dg.t2)
  t5 <- (t4/min(sample_depth))*100
  t6 <- round(t5, 2)
  c6[x] <- t6
  
  t2 <- t1[samples.control.t0]
  t3 <- colSums(t2)
  t4 <- sum(t3)/length(samples.control.t0)
  t5 <- (t4/min(sample_depth))*100
  t6 <- round(t5, 2)
  c7[x] <- t6
  
  t2 <- t1[samples.control.t1]
  t3 <- colSums(t2)
  t4 <- sum(t3)/length(samples.control.t1)
  t5 <- (t4/min(sample_depth))*100
  t6 <- round(t5, 2)
  c8[x] <- t6

  x = x + 1
}

taxa.abundancies <- data.frame(TAXA = samples.taxa,
                               SSG_t0 = c1, SSG_t1 = c2, SSG_t2 = c3,
                               FDG_t0 = c4, FDG_t1 = c5, FDG_t2 = c6,
                               CON_t0 = c7, CON_t1 = c8
                               )

write.csv(taxa.abundancies, file = "results/tab_otu_species_mean.csv", quote = FALSE, row.names = FALSE)

# --------------------------------------------------------------------------------------------------------

# [06] Assessing differences in beta diversities

beta.horses <- unique(meta.sorted$HorseID)
beta.intra <- c()
beta.inter <- c()

for (h in beta.horses){
  t <- meta.sorted[meta.sorted$HorseID == h, ]$SampleID
  r <- data.bray[t[1], ]
  # We select t0 as a starting point and query its ID to receive t1 and t2 distances
  beta.intra <- c(beta.intra, r[(names(r) %in% t)])
  # Afterwards, we assess the distances of this t0 sample to all other samples
  beta.inter <- c(beta.inter, r[!(names(r) %in% t)])
}

beta.intra <- beta.intra[!(beta.intra %in% c(0))]
beta.df <- data.frame(GROUP = c(rep("Intra", length(beta.intra)), rep("Inter", length(beta.inter))), VAL = c(beta.intra, beta.inter))
beta.res <- pairwise.wilcox.test(beta.df$VAL, beta.df$GROUP, paired = FALSE, alternative = "less", p.adjust.method = "BH")

# --------------------------------------------------------------------------------------------------------

# [07] Calculating log2FC on family level

df_ssg_t0 <- NULL
df_ssg_t1 <- NULL
df_ssg_t2 <- NULL

df_5dg_t0 <- NULL
df_5dg_t1 <- NULL
df_5dg_t2 <- NULL

taxa_list <- unique(data.otu.rarefy$Rank5)

for (t in taxa_list){
  t1 <- data.otu.rarefy[data.otu.rarefy$Rank1 == t |
                        data.otu.rarefy$Rank2 == t |
                        data.otu.rarefy$Rank3 == t |
                        data.otu.rarefy$Rank4 == t |
                        data.otu.rarefy$Rank5 == t |
                        data.otu.rarefy$Rank6 == t |
                        data.otu.rarefy$Rank7 == t, ]
  
  # Control
  c1 <- t1[samples.control.t0]
  c2 <- colSums(c1)
  c3 <- median(c2)
  
  if (c3 > 0){
    # SSG
    s1 <- t1[samples.ssg.t0]
    s2 <- colSums(s1)
    s3 <- (s2 - c3)/c3
    s4 <- sign(s3)*round(log(abs(s3), 2), 2)
    s5 <- c(t, s4)
    rbind(df_ssg_t0, s5) -> df_ssg_t0
    
    s1 <- t1[samples.ssg.t1]
    s2 <- colSums(s1)
    s3 <- (s2 - c3)/c3
    s4 <- sign(s3)*round(log(abs(s3), 2), 2)
    s5 <- c(t, s4)
    rbind(df_ssg_t1, s5) -> df_ssg_t1
    
    s1 <- t1[samples.ssg.t2]
    s2 <- colSums(s1)
    s3 <- (s2 - c3)/c3
    s4 <- sign(s3)*round(log(abs(s3), 2), 2)
    s5 <- c(t, s4)
    rbind(df_ssg_t2, s5) -> df_ssg_t2
    
    # 5DG
    s1 <- t1[samples.5dg.t0]
    s2 <- colSums(s1)
    s3 <- (s2 - c3)/c3
    s4 <- sign(s3)*round(log(abs(s3), 2), 2)
    s5 <- c(t, s4)
    rbind(df_5dg_t0, s5) -> df_5dg_t0
    
    s1 <- t1[samples.5dg.t1]
    s2 <- colSums(s1)
    s3 <- (s2 - c3)/c3
    s4 <- sign(s3)*round(log(abs(s3), 2), 2)
    s5 <- c(t, s4)
    rbind(df_5dg_t1, s5) -> df_5dg_t1
    
    s1 <- t1[samples.5dg.t2]
    s2 <- colSums(s1)
    s3 <- (s2 - c3)/c3
    s4 <- sign(s3)*round(log(abs(s3), 2), 2)
    s5 <- c(t, s4)
    rbind(df_5dg_t2, s5) -> df_5dg_t2
  }
}

df_ssg_t0[df_ssg_t0 == "NaN"] <- "0"
df_ssg_t1[df_ssg_t1 == "NaN"] <- "0"
df_ssg_t2[df_ssg_t2 == "NaN"] <- "0"

df_5dg_t0[df_5dg_t0 == "NaN"] <- "0"
df_5dg_t1[df_5dg_t1 == "NaN"] <- "0"
df_5dg_t2[df_5dg_t2 == "NaN"] <- "0"

df_ssg_t0 <- as.data.frame(df_ssg_t0)
df_ssg_t1 <- as.data.frame(df_ssg_t1)
df_ssg_t2 <- as.data.frame(df_ssg_t2)

df_5dg_t0 <- as.data.frame(df_5dg_t0)
df_5dg_t1 <- as.data.frame(df_5dg_t1)
df_5dg_t2 <- as.data.frame(df_5dg_t2)

new.names.old = meta.sorted[match(colnames(df_ssg_t0)[c(-1)], meta.sorted$SampleID),]
colnames(df_ssg_t0) <- c("Family", new.names.old$HorseID)
df_ssg_t0 <- df_ssg_t0[,order(colnames(df_ssg_t0))]
df_ssg_t0 <- df_ssg_t0[order(df_ssg_t0$Family),]

new.names.old = meta.sorted[match(colnames(df_ssg_t1)[c(-1)], meta.sorted$SampleID),]
colnames(df_ssg_t1) <- c("Family", new.names.old$HorseID)
df_ssg_t1 <- df_ssg_t1[,order(colnames(df_ssg_t1))]
df_ssg_t1 <- df_ssg_t1[order(df_ssg_t1$Family),]

new.names.old = meta.sorted[match(colnames(df_ssg_t2)[c(-1)], meta.sorted$SampleID),]
colnames(df_ssg_t2) <- c("Family", new.names.old$HorseID)
df_ssg_t2 <- df_ssg_t2[,order(colnames(df_ssg_t2))]
df_ssg_t2 <- df_ssg_t2[order(df_ssg_t2$Family),]

new.names.old = meta.sorted[match(colnames(df_5dg_t0)[c(-1)], meta.sorted$SampleID),]
colnames(df_5dg_t0) <- c("Family", new.names.old$HorseID)
df_5dg_t0 <- df_5dg_t0[,order(colnames(df_5dg_t0))]
df_5dg_t0 <- df_5dg_t0[order(df_5dg_t0$Family),]

new.names.old = meta.sorted[match(colnames(df_5dg_t1)[c(-1)], meta.sorted$SampleID),]
colnames(df_5dg_t1) <- c("Family", new.names.old$HorseID)
df_5dg_t1 <- df_5dg_t1[,order(colnames(df_5dg_t1))]
df_5dg_t1 <- df_5dg_t1[order(df_5dg_t1$Family),]

new.names.old = meta.sorted[match(colnames(df_5dg_t2)[c(-1)], meta.sorted$SampleID),]
colnames(df_5dg_t2) <- c("Family", new.names.old$HorseID)
df_5dg_t2 <- df_5dg_t2[,order(colnames(df_5dg_t2))]
df_5dg_t2 <- df_5dg_t2[order(df_5dg_t2$Family),]

list_of_datasets <- list("SSG_t0" = df_ssg_t0,
                         "SSG_t1" = df_ssg_t1,
                         "SSG_t2" = df_ssg_t2,
                         "5DG_t0" = df_5dg_t0,
                         "5DG_t1" = df_5dg_t1,
                         "5DG_t2" = df_5dg_t2
                         )

write.xlsx(list_of_datasets, file = "results/tab_otu_family_fc.xlsx")