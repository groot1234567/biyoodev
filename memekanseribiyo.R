
# ============================================
# Meme Kanseri Gen D0fadesi Analizi (GSE70947)
# Bioconductor TabanlD1 R Scripti
# ============================================

# Gerekli paketleri yükle (sadece ilk sefer için gerekliyse)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "Biobase", "pheatmap", "EnhancedVolcano"))
install.packages("ggplot2")

# Kütüphaneleri çağır
library(GEOquery)
library(Biobase)
library(limma)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)

# GEO veri setini indir
gse <- getGEO("GSE70947", GSEMatrix = TRUE)
eset <- gse[[1]]

# Log2 
exprs_data <- exprs(eset)
if (max(exprs_data, na.rm = TRUE) > 100) {
  exprs(eset) <- log2(exprs_data + 1)
}

# Quantile normalizasyon
exprs(eset) <- normalizeBetweenArrays(exprs(eset), method = "quantile")

# Grup isimlerini düzelt (geçerli R ismi formatD1 haline getir)
group_raw <- pData(eset)$`source_name_ch1`
group <- factor(make.names(group_raw))  

# Design matrisi oluştur
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Contrast matrisi tanımla (grup isimleri düzeltildiği için artık)
contrast.matrix <- makeContrasts(
  TumorVsNormal = breast.adenocarcinoma - breast.tissue.adjacent.to.tumor,
  levels = design
)

# Modeli kur ve genleri bul
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# SonuçlarD1 al
results <- topTable(fit2, adjust = "fdr", number = Inf)
significant_genes <- subset(results, adj.P.Val < 0.05 & abs(logFC) > 1) 

# PCA için NA, NaN, Inf içeren satırlar?? (genleri) kaldır
exprs_clean <- exprs(eset)
exprs_clean <- exprs_clean[complete.cases(exprs_clean), ] 
exprs_clean <- exprs_clean[apply(exprs_clean, 1, function(x) all(is.finite(x))), ]  

# PCA analizi (sorunsuz)
pca <- prcomp(t(exprs_clean), scale. = TRUE)

# Görselleötirme için veri 
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Group = group)

# PCA plot çiz
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")


# PCA Plot
pca <- prcomp(t(exprs(eset)), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Group = group)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot", x = "PC1", y = "PC2")

# Volcano Plot
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano Plot',
                subtitle = "Tumor vs Normal",
                legendPosition = 'right') 

# Anlamlı genlerin ekspresyon matrisini al
heat_data <- exprs(eset)[rownames(significant_genes), ]

# NA / Inf içeren satırları kaldır
heat_data_clean <- heat_data[complete.cases(heat_data), ]
heat_data_clean <- heat_data_clean[apply(heat_data_clean, 1, function(x) all(is.finite(x))), ]

# Z-score standardizasyonu
heat_data_scaled <- t(scale(t(heat_data_clean)))

# Grupları renklemek için anotasyon
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(heat_data_scaled)

# Heatmap çizimi
pheatmap(heat_data_scaled,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         main = "Heatmap of Significant Genes")


# Heatmap (anlamlD1 genlerle)
heat_data <- exprs(eset)[rownames(significant_genes), ]
heat_data_scaled <- t(scale(t(heat_data)))  # z-score dC6nC<EC<mC<
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(heat_data_scaled)

pheatmap(heat_data_scaled,
         show_rownames = FALSE,
         annotation_col = annotation_col,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         main = "Heatmap of Significant Genes")
