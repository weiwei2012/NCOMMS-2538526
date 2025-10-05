library(Seurat)
library(Matrix)
library(SeuratObject)
library(harmony)
library(devtools)
library(clustree)
library(cowplot)
scRNAlist1 <- readRDS("~/8w_harmony.rds")
scRNAlist2 <- readRDS("~/10w_harmony.rds")
scRNAlist3 <- readRDS("~/12w_harmony.rds")
scRNAlist4 <- readRDS("~/15w_harmony.rds")
#点击安装embryo.RData
scRNAlist5 <- embryo
scRNAlist6 <- readRDS("~/Pulp_harmony.rds")
scRNA_harmony <- merge(scRNAlist1, y=c(scRNAlist2, scRNAlist3, scRNAlist4, scRNAlist5, scRNAlist6))
scRNA_harmony[["percent.mt"]] <- PercentageFeatureSet(scRNA_harmony, pattern = "^MT-") #鼠源的需要换成mt
scRNA_harmony <- subset(scRNA_harmony, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)   
scRNA_harmony <- NormalizeData(scRNA_harmony, normalization.method = "LogNormalize", scale.factors=10000)
scRNA_harmony <- FindVariableFeatures(scRNA_harmony, selection.method = "vst", nfeatures = 2000)
scRNA_harmony <- ScaleData(scRNA_harmony, vars.to.regress = "percent.mt")
scRNA_harmony <- RunPCA(scRNA_harmony, verbose=FALSE)
scRNA_harmony@meta.data[["orig.ident"]]
unique_ident <- unique(scRNA_harmony$orig.ident)
scRNA_harmony$sample <- sub("-.*", "", scRNA_harmony@meta.data[["orig.ident"]])
scRNA_harmony@meta.data[["sample"]]
unique_samples <- unique(scRNA_harmony$sample)

system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:10)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.4)
table(scRNA_harmony$seurat_clusters)
saveRDS(scRNA_harmony, '24new_human_all.rds')

all.markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, group.by = "seurat_clusters")
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

named_clusters <- c("Dental Follicle", "Dental Papilla", "Dental Follicle","Dental Epithelium Cell", "Dental Papilla", 
                    "Endothelial Cell", "Proliferating Cell", "Monocyte", "Lymphocyte", "Erythroblast", 
                    "Fibroblast", "Vascular Smooth Muscle Cell", "Neuronal Cell")
scRNA_harmony$CellType <- factor(scRNA_harmony$seurat_clusters, labels = named_clusters)

#大亚群UMAP图
custom_colors <- c("#A589B8", "#92C1FF", "#D16969", 
                   "#F0B2B2", "#DCE9FF", "#C8C8C8", 
                   "#888888", "#F5D3D3", "#EBD3FF", 
                   "#C171A1", "#EBCF98")
plot1 = DimPlot(scRNA_harmony, 
                reduction = "umap", 
                group.by='CellType',
                cols = custom_colors,) 
ggsave("大亚群umap.png", plot = plot1, width = 8, height = 6)

#各时间点样本umap图
custom_colors <- c("#E6E6D8", "#D9D9C7", "#CDCDB6", "#C0C0A5", "#B3B394", "#A6A683", "#999972")

plot2 = DimPlot(scRNA_harmony, 
                    split.by = "sample", 
                    group.by='sample',
                    reduction = "umap", 
                    cols = custom_colors)
ggsave("各时间点umap.png", plot = plot2, width = 20, height = 5)

#大亚群高表达基因特征展示
pdf(file = "各大亚群Top3基因热图.pdf")
GroupHeatmap(
  srt = scRNA_harmony,
  features = c(
               "KRT19", "KRT5", "KRT14",
               "DLX6-AS1", "FGF3", "RSPO4",
               "COL12A1", "POSTN", "OGN",
               "TOP2A", "NUSAP1", "UBE2C",
               "ESAM", "CLDN5","VWF",
               "PDGFRB", "RGS5", "BGN",
               "TYROBP", "AIF1", "CD14",
               "TRBC2", "CD3E", "CD48",
               "ALAS2", "HBG1", "HBA1",
               "ACTC1", "SGCA", "TNNT1",
               "S100B", "MPZ", "PLP1"),
  group.by = c("CellType"),
  heatmap_palette = "YlOrRd",
  cell_annotation_palette = c("Dark2", "Paired", "Paired"),
  show_row_names = TRUE, row_names_side = "left",
  add_dot = T, dot_size = unit(6, "mm"),
  add_reticle = F
)
dev.off()

# 对基因做umap映射图
pdf(file = "各大亚群Top基因umap1.pdf")
FeatureDimPlot(
  srt = scRNA_harmony,
  features = c("COL12A1", "POSTN", "OGN",
               "DLX6-AS1", "FGF3", "RSPO4",
               "DSPP", "DMP1"),
  reduction = "umap", 
  theme_use = "theme_blank"
)
dev.off()

pdf(file = "各大亚群Top基因umap2.pdf")
FeatureDimPlot(
  srt = scRNA_harmony,
  features = c(
    "KRT19", "KRT5", "KRT14",
    "ESAM", "CLDN5","WNT6",
    "WNT7B", "WNT10B", "JAG2"),
  reduction = "umap", 
  theme_use = "theme_blank"
)
dev.off()

pdf(file = "各大亚群Top基因umap3.pdf")
FeatureDimPlot(
  srt = scRNA_harmony,
  features = c(
    "TOP2A", "NUSAP1", "UBE2C",
    "TYROBP", "AIF1", "CD14",
    "TRBC2", "CD3E", "CD48"
  ),
  reduction = "umap", 
  theme_use = "theme_blank"
)
dev.off()

pdf(file = "各大亚群Top基因umap4.pdf")
FeatureDimPlot(
  srt = scRNA_harmony,
  features = c(
    "ALAS2", "HBG1", "HBA1",
    "PDGFRB", "RGS5", "BGN",
    "ACTC1", "SGCA", "TNNT1",
    "S100B", "MPZ", "PLP1"),
  reduction = "umap", 
  theme_use = "theme_blank"
)
dev.off()

#上皮牙乳头映射
library(dplyr)
scRNA_harmony <- SCTransform(scRNA_harmony, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
# 加载单张ST数据
tooth_molar_15w <- readRDS("~/tooth_molar_15w.rds")
tooth_molar_15w <- RunPCA(tooth_molar_15w, assay = "SCT", verbose = FALSE)
tooth_molar_15w <- FindNeighbors(tooth_molar_15w, reduction = "pca", dims = 1:30)
tooth_molar_15w <- FindClusters(tooth_molar_15w, resolution = 1.0, verbose = FALSE)
tooth_molar_15w <- RunUMAP(tooth_molar_15w, reduction = "pca", dims = 1:30)
tooth_molar_15w <- SCTransform(tooth_molar_15w, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
anchors <- FindTransferAnchors(reference = scRNA_harmony, query = tooth_molar_15w, normalization.method = "SCT",
                               k.anchor = 10, # 可以调整k.anchor为更小的值
                               k.filter = NA)
predictions.assay <- TransferData(anchorset = anchors, refdata = scRNA_harmony$CellType, prediction.assay = TRUE, 
                                  weight.reduction = tooth_molar_15w[["pca"]], dims = 1:30)
tooth_molar_15w[["predictions"]] <- predictions.assay
DefaultAssay(tooth_molar_15w) <- "predictions"

pdf("大群映射Celltype_15.pdf")
SpatialFeaturePlot(tooth_molar_15w, features = c("Dental Papilla", "Epithelium Cell", "Dental Follicle"),
                   pt.size.factor = 1, ncol = 2, crop = F, alpha = c(0.1, 1)
)
dev.off()


tooth_molar_10w <- readRDS("~/tooth_molar_10w.rds")
tooth_molar_10w <- RunPCA(tooth_molar_10w, assay = "SCT", verbose = FALSE)
tooth_molar_10w <- FindNeighbors(tooth_molar_10w, reduction = "pca", dims = 1:30)
tooth_molar_10w <- FindClusters(tooth_molar_10w, resolution = 2.0, verbose = FALSE)
tooth_molar_10w <- RunUMAP(tooth_molar_10w, reduction = "pca", dims = 1:30)
tooth_molar_10w <- SCTransform(tooth_molar_10w, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
anchors <- FindTransferAnchors(reference = scRNA_harmony, 
                               query = tooth_molar_10w, 
                               normalization.method = "SCT", 
                               k.anchor = 5, # 可以调整k.anchor为更小的值
                               k.filter = NA) # 可以设置k.filter为NA，或者为一个较小的值
predictions.assay <- TransferData(anchorset = anchors, refdata = scRNA_harmony$CellType, prediction.assay = TRUE, 
                                  weight.reduction = tooth_molar_10w[["pca"]], dims = 1:30)
tooth_molar_10w[["predictions"]] <- predictions.assay
DefaultAssay(tooth_molar_10w) <- "predictions"
pdf("三大群映射Celltype_10w.pdf")
SpatialFeaturePlot(tooth_molar_10w, features = c("Dental Papilla", "Epithelium Cell", "Dental Follicle"),
                   pt.size.factor = 1, ncol = 2, crop = F, alpha = c(0.1, 1)
)
dev.off()

tooth_molar_12w <- readRDS("~/tooth_molar_12w.rds")
tooth_molar_12w <- RunPCA(tooth_molar_12w, assay = "SCT", verbose = FALSE)
tooth_molar_12w <- FindNeighbors(tooth_molar_12w, reduction = "pca", dims = 1:30)
tooth_molar_12w <- FindClusters(tooth_molar_12w, resolution = 3.0, verbose = FALSE)
tooth_molar_12w <- RunUMAP(tooth_molar_12w, reduction = "pca", dims = 1:30)
tooth_molar_12w <- SCTransform(tooth_molar_12w, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
# FindAnchor,得到了每个班每个点的预测分数
anchors <- FindTransferAnchors(reference = scRNA_harmony, 
                               query = tooth_molar_12w, 
                               normalization.method = "SCT", 
                               k.anchor = 5, # 可以调整k.anchor为更小的值
                               k.filter = NA) # 可以设置k.filter为NA，或者为一个较小的值
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = scRNA_harmony$CellType, 
                                  prediction.assay = TRUE, 
                                  weight.reduction = tooth_molar_12w[["pca"]], 
                                  dims = 1:30,
                                  k.weight = 5) # 确保 k.weight 小于锚点数

tooth_molar_12w[["predictions"]] <- predictions.assay
DefaultAssay(tooth_molar_12w) <- "predictions"
pdf("三大群映射Celltype_12w.pdf")
SpatialFeaturePlot(tooth_molar_12w, features = c("Dental Papilla", "Epithelium Cell", "Dental Follicle"),
                   pt.size.factor = 1, ncol = 2, crop = F, alpha = c(0.1, 1)
)
dev.off()

# 点击加载spatial.RDara
tooth_molar_17w <- spatial
tooth_molar_17w <- RunPCA(tooth_molar_17w, assay = "SCT", verbose = FALSE)
tooth_molar_17w <- FindNeighbors(tooth_molar_17w, reduction = "pca", dims = 1:30)
tooth_molar_17w <- FindClusters(tooth_molar_17w, resolution = 1.0, verbose = FALSE)
tooth_molar_17w <- RunUMAP(tooth_molar_17w, reduction = "pca", dims = 1:30)
tooth_molar_17w <- SCTransform(tooth_molar_17w, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
# FindAnchor,得到了每个班每个点的预测分数
anchors <- FindTransferAnchors(reference = scRNA_harmony, 
                               query = tooth_molar_17w, 
                               normalization.method = "SCT", 
                               k.anchor = 5, # 可以调整k.anchor为更小的值
                               k.filter = NA) # 可以设置k.filter为NA，或者为一个较小的值
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = scRNA_harmony$CellType, 
                                  prediction.assay = TRUE, 
                                  weight.reduction = tooth_molar_17w[["pca"]], 
                                  dims = 1:30,
                                  k.weight = 5) # 确保 k.weight 小于锚点数
tooth_molar_17w[["predictions"]] <- predictions.assay
DefaultAssay(tooth_molar_17w) <- "predictions"
pdf("三大群映射Celltype_17w.pdf")
SpatialFeaturePlot(tooth_molar_17w, features = c("Dental Papilla", "Epithelium Cell", "Dental Follicle"),
                   pt.size.factor = 1, ncol = 2, crop = F, alpha = c(0.1, 1)
)
dev.off()

#各大发育通路在上皮牙乳头的表达情况
Idents(scRNA_harmony) <- "CellType"
selected_DEDP <- c("Epithelium Cell", "Dental Papilla")
sub_DEDP <- subset(scRNA_harmony, idents = selected_DEDP)
sub_DEDP <- NormalizeData(sub_DEDP, normalization.method = "LogNormalize", scale.factors=10000)
sub_DEDP@meta.data[["sample"]]
unique_sample<- unique(sub_DEDP$sample)
print(unique_sample)
sub_DEDP$sample_cluster <- paste(sub_DEDP@meta.data[["sample"]], sub_DEDP@meta.data[["CellType"]], sep = "-")
sub_DEDP@meta.data[["sample_cluster"]]
unique_sample_cluster<- unique(sub_DEDP$sample_cluster)
genes_of_interest <- c(wnt_genes, bmp_genes, fgf_genes, shh_genes, eda_genes, notch_genes, tgf_genes)
expr_data <- FetchData(sub_DEDP_subset, vars = genes_of_interest)
pathway_expr <- data.frame(
  Pathway = c(rep("WNT", length(wnt_genes)), 
              rep("BMP", length(bmp_genes)), 
              rep("FGF", length(fgf_genes)), 
              rep("SHH", length(shh_genes)), 
              rep("EDA", length(eda_genes)), 
              rep("NOTCH", length(notch_genes)),
              rep("TGF", length(tgf_genes))),
  Expression = colSums(expr_data)
)
pathway_summary <- aggregate(Expression ~ Pathway, data = pathway_expr, sum)

data_long <- tidyr::pivot_longer(data, cols = -Pathway, names_to = "Week", values_to = "Value")
data_long <- data_long %>%
  group_by(Week) %>%
  mutate(Percentage = Value / sum(Value) * 100)
data_long$Pathway <- factor(data_long$Pathway, levels = c("WNT", "NOTCH", "SHH", "EDA", "BMP", "TGF", "FGF"))
ggplot(data_long, aes(x = Week, y = Percentage, fill = Pathway, alluvium = Pathway, stratum = Pathway)) +
  geom_col(position = 'stack', width = 0.7)+
  geom_alluvium(aes(stratum = Pathway), width = 0.7, alpha = 0.2,color="black",linewidth = 0.5,curve_type = "sigmoid") +
  geom_stratum(width = 0.7, alpha = 0.5,color="black", linewidth = 0.5) +
  scale_fill_manual(values = c("#FFE082", "#FF8F00","#B39DDB", "#C8C8C8", "#80CBC4", "#006400", "#D4A890")) +
  scale_y_continuous(expand = c(0,0)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  labs(title = "Pathway Percentage Over Time_Legend", y = "Percentage", x = "Week")

#pathway 上皮间质表达占比饼图
# 提取两个细胞群
cluster1_cells <- WhichCells(sub_DEDP, ident = "Epithelium Cell")
cluster2_cells <- WhichCells(sub_DEDP, ident = "Dental Papilla")
genes_of_interest <- c(wnt_genes, bmp_genes, fgf_genes, shh_genes, eda_genes, notch_genes, tgf_genes)
cluster1_expr <- FetchData(sub_DEDP, vars = genes_of_interest, cells = cluster1_cells)
cluster1_expr <- data.frame(
  Pathway = c(rep("WNT", length(wnt_genes)), 
              rep("BMP", length(bmp_genes)), 
              rep("FGF", length(fgf_genes)), 
              rep("SHH", length(shh_genes)), 
              rep("EDA", length(eda_genes)), 
              rep("NOTCH", length(notch_genes)),
              rep("TGF", length(tgf_genes))),
  Expression = colSums(cluster1_expr)
)
cluster2_expr <- FetchData(sub_DEDP, vars = genes_of_interest, cells = cluster2_cells)
cluster2_expr <- data.frame(
  Pathway = c(rep("WNT", length(wnt_genes)), 
              rep("BMP", length(bmp_genes)), 
              rep("FGF", length(fgf_genes)), 
              rep("SHH", length(shh_genes)), 
              rep("EDA", length(eda_genes)), 
              rep("NOTCH", length(notch_genes)),
              rep("TGF", length(tgf_genes))),
  Expression = colSums(cluster2_expr)
)

# 合并数据
cluster1_expr$Cluster <- "Cluster1"
cluster2_expr$Cluster <- "Cluster2"
expr_data <- rbind(cluster1_expr, cluster2_expr)
gene_cluster_summary_wide$total <- gene_cluster_summary_wide$Cluster1 + gene_cluster_summary_wide$Cluster2
gene_cluster_summary_wide$Cluster1_percent <- gene_cluster_summary_wide$Cluster1 / gene_cluster_summary_wide$total
gene_cluster_summary_wide$Cluster2_percent <- gene_cluster_summary_wide$Cluster2 / gene_cluster_summary_wide$total

# 绘制单个通路的饼图函数
draw_pie <- function(gene_name) {
  data <- data.frame(cluster = c("Cluster1", "Cluster2"),
                     percent = c(gene_cluster_summary_wide[gene_cluster_summary_wide$Pathway  == gene_name, ]$Cluster1_percent,
                                 gene_cluster_summary_wide[gene_cluster_summary_wide$Pathway  == gene_name, ]$Cluster2_percent))
  ggplot(data, aes(x = "", y = percent, fill = cluster)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    labs(title = paste("Pathway:", gene_name)) +
    scale_fill_manual(values = c("#CE5E5E", "#91BAFF")) +
    theme_void()
}

# 获取所有基因名称
genes <- unique(gene_cluster_summary_wide$Pathway)

# 创建一个空的列表来存储每个基因的饼图
pie_plots <- list()

# 为每个基因绘制饼图并添加到列表中
for (gene in genes) {
  pie_plots[[gene]] <- draw_pie(gene)
}

tryCatch({
  # 将所有基因的饼图排列在一个图形中
  combined_plot <- grid.arrange(grobs = pie_plots, ncol = ceiling(length(genes)/5))
  # 保存图形为 PDF 文件
  ggsave("Fig1E_pathway_receptors_pie_charts.pdf", combined_plot, width = 10, height = 10)
}, error = function(e) {
  cat("Error occurred while arranging plots:", conditionMessage(e), "\n")
})
