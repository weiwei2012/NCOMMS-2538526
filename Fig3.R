library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)

# 基因映射
tooth_molar_15w <- readRDS("~/tooth_molar_15w.rds")
tooth_molar_15w <- RunPCA(tooth_molar_15w, assay = "SCT", verbose = FALSE)
tooth_molar_15w <- FindNeighbors(tooth_molar_15w, reduction = "pca", dims = 1:30)
tooth_molar_15w <- FindClusters(tooth_molar_15w, resolution = 1.0, verbose = FALSE)
tooth_molar_15w <- RunUMAP(tooth_molar_15w, reduction = "pca", dims = 1:30)
DP <- subset(tooth_molar_15w, idents = c(4, 12))

tooth_molar_17w <- spatial
tooth_molar_17w <- SCTransform(tooth_molar_17w, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
tooth_molar_17w <- RunPCA(tooth_molar_17w, assay = "SCT", verbose = FALSE)
tooth_molar_17w <- FindNeighbors(tooth_molar_17w, reduction = "pca", dims = 1:30)
tooth_molar_17w <- FindClusters(tooth_molar_17w, resolution = 1.0, verbose = FALSE)
tooth_molar_17w <- RunUMAP(tooth_molar_17w, reduction = "pca", dims = 1:30)
DP <- subset(tooth_molar_17w, idents = c(10, 3))

SpatialFeaturePlot(DP, features = c("DLX6-AS1","LEF1","FZD2", "HEY1", "NOTCH2"),
                   pt.size.factor = 1, min.cutoff = 0, ncol = 3, 
                   crop = F, alpha = c(0,1)) + scale_alpha()


# DP分群umap及小提琴图
custom_colors <- c("#0066CC", "#99CCFF", "#33CC99")
plot2 = DimPlot(sub_DP, 
                reduction = "umap", 
                group.by='Cell_Type',
                cols = custom_colors,) 

pdf(file = "Fig2C_DP三亚群基因小提琴图.pdf")
FeatureStatPlot(
  srt = sub_DP, 
  group.by = "Cell_Type4",
  stat.by = c('DLX6-AS1',
              "FGF3",
              "KIF5C",
              "PTCH1",
              "S100A6",
              
              "NR2F1",
              "TWIST1",
              "SOX11",
              "CCND2",
              "CCN1",
              
              "MT-RNR2",
              "DSTN",
              "GPX3", 
              "TF",
              "CYP1B1"
              
  ), 
  legend.position = "top", legend.direction = "horizontal",
  stack = T,
  plot_type = "violin",
  palcolor = c("#0066CC", "#99CCFF", "#33CC99")
)
dev.off()


#MA
# Convert to factor
sub_DP$sample <- as.factor(sub_DP$sample)
Idents(sub_DP) <- "sample"
# Get differential expression between clusters
de_result <- FindMarkers(sub_DP, only.pos = FALSE, ident.1 = "17w", ident.2 = "15w", logfc.threshold = 0.25, test.use = "wilcox")

genes_for_de <- rownames(de_result)
# 在计算平均 log2 表达量时，只考虑这些基因
expr_matrix <- as.matrix(GetAssayData(sub_DP))
sample_groups <- unique(sub_DP@meta.data$sample)
avg_log2_expr <- numeric(length(genes_for_de))
for (i in seq_along(genes_for_de)) {
  gene_idx <- which(rownames(expr_matrix) == genes_for_de[i])
  gene_expr <- expr_matrix[gene_idx, ]
  avg_expr_per_group <- sapply(sample_groups, function(group) {
    mean(gene_expr[sub_DP@meta.data$sample == group])
  })
  avg_log2_expr[i] <- log2(mean(avg_expr_per_group))
}

de_result$avg_log2_expr <- avg_log2_expr
ggplot(de_result, aes(x = avg_log2_expr, y = avg_log2FC)) +
  geom_point(aes(color = p_val_adj < 0.05), size = 2, alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "#8B1A1A")) +
  labs(x = "Average Log2 Expression", y = "Log2 Fold Change", color = "Significant") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank()) 

# GO
de_genes <- rownames(de_result)[de_result$p_val_adj < 0.05 & abs(de_result$avg_log2FC) > 1]
ego_mf <- enrichGO(gene = de_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego_mf_df <- as.data.frame(ego_mf)
ego_mf_df$significant <- ifelse(ego_mf_df$p.adjust < 0.05, "Significant", "Not Significant")
num_rows <- nrow(ego_mf_df)
gray_palette <- colorRampPalette(c("white", "gray90", "gray80", "gray70", "gray60", "gray50", "gray40", "gray30", "gray20", "gray10"))(num_rows + 10)
ego_mf_df$color <- ifelse(1:nrow(ego_mf_df) <= length(gray_palette), 
                          gray_palette[as.integer((ego_mf_df$p.adjust / 0.05) * length(gray_palette))], 
                          gray_palette[length(gray_palette)])

p2 <- ggplot(ego_mf_df, aes(x = reorder(Description, -p.adjust), y = Count, fill = color)) +
  geom_bar(stat = "identity") +
  labs(x = "GO Term", y = "Count", title = "GO Enrichment (MF)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_identity(name = "p.adjust")
ego_bp <- enrichGO(gene = de_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego_bp_df <- as.data.frame(ego_bp)
ego_bp_df$significant <- ifelse(ego_bp_df$p.adjust < 0.05, "Significant", "Not Significant")
ego_bp_df_top20 <- ego_bp_df[1:30, ]
ego_bp_df_top20$color <- gray_palette[pmax(1, pmin(as.integer((ego_bp_df_top20$p.adjust / 0.05) * 100), length(gray_palette)))]
p3 <- ggplot(ego_bp_df_top20, aes(x = reorder(Description, -p.adjust), y = Count, fill = color)) +
  geom_bar(stat = "identity") +
  labs(x = "GO Term", y = "Count", title = "GO Enrichment (BP)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_identity(name = "p.adjust")

ego_cc <- enrichGO(gene = de_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ego_cc_df <- as.data.frame(ego_cc)
ego_cc_df$significant <- ifelse(ego_cc_df$p.adjust < 0.05, "Significant", "Not Significant")
p4 <- ggplot(ego_cc_df, aes(x = reorder(Description, -p.adjust), y = Count, fill = color)) +
  geom_bar(stat = "identity") +
  labs(x = "GO Term", y = "Count", title = "GO Enrichment (CC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_identity(name = "p.adjust")


##dlx6as1与runx2 bmp等基因的相关性图##############
data <- sub_DP@assays$RNA@data
# 提取相关基因的表达值
dlx6_as1 <- data["DLX6-AS1",]
dmp1 <- data["DMP1",]
dspp <- data["DSPP",]
runx2 <- data["RUNX2",]
bmp4 <- data["BMP4",]
dlx1 <- data["DLX1",]
lef1 <- data["LEF1",]
notch2 <- data["NOTCH2",]
hes1 <- data["HES1",]
lrp5 <- data["LRP5",]

cor_data <- data.frame(DLX6_AS1 = dlx6_as1, DMP1 = dmp1, DSPP=dspp, RUNX2 = runx2, BMP4=bmp4, DLX1= dlx1,
                       LEF1=lef1, LRP5=lrp5, NOTCH2=notch2, HES1=hes1)
cor_matrix <- cor(cor_data)
colors <- c("red", "#66CC99", "#66CC99", "#66CC99","#66CC99","#66CC99", "#FFCC66", "#FFCC66","#FF9966", "#FF9966")
chordDiagram(cor_matrix, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1), grid.col = colors)
gap <- -1  # 设定一个适当的距离值
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] - gap, CELL_META$sector.index, niceFacing = TRUE, cex = 0.8)
}, bg.border = NA)


#BEAM分析
seurat_obj <- sub_DP
expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
cell_metadata <- seurat_obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expr_matrix))
rownames(gene_annotation) <- rownames(expr_matrix)
cds <- newCellDataSet(expr_matrix,
                      phenoData = new("AnnotatedDataFrame", data = cell_metadata),
                      featureData = new("AnnotatedDataFrame", data = gene_annotation),
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, ordering_genes$gene_id)
cds <- reduceDimension(cds, method = 'DDRTree')
cds <- orderCells(cds)
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- cds[disp.genes,]
plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 3, progenitor_method = "duplicate")
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-4)),]
plot_genes_branched_heatmap(mycds_sub_beam, branch_point = 3, num_clusters = 3, show_rownames = T)



Idents(sub_DP) <- "Cell_Type"
selected_DP1 <- c("DLX6-AS1+ DP")
sub_DP <- subset(sub_DP, idents = selected_DP1)
sub_DP <- NormalizeData(sub_DP1, normalization.method = "LogNormalize", scale.factors=10000)
sub_DP@meta.data[["sample"]]
unique_sample<- unique(sub_DP1$sample)
print(unique_sample)
sub_DP1$sample_celltype <- paste(sub_DP1@meta.data[["sample"]], sub_DP1@meta.data[["Cell_Type4"]], sep = "-")
sub_DP1@meta.data[["sample_celltype"]]
unique_sample_celltype<- unique(sub_DP1$sample_celltype)

sub_DP <- RunDEtest(srt = sub_DP, group_by = "sample", fc.threshold = 1, only.pos = FALSE)
DEG <- sub_DP@tools$DEtest_sample$AllMarkers_wilcox
DEG$label <- ifelse(DEG$p_val_adj<0.05,
                    ifelse(DEG$avg_log2FC>0, "Up", "Down"),
                    "None")
DEG$label = factor(DEG$label, levels = c("Up", "None", "Down"))
table(DEG$label)

DEG.text = NULL;
for(group1.name in unique(DEG$group1) ){
  print(group1.name)
  tmp=DEG %>% filter(group1==group1.name & p_val_adj <0.05) %>% top_n(n = 100, wt = avg_log2FC)
  DEG.text=rbind(DEG.text, tmp);
  tmp=DEG %>% filter(group1==group1.name & p_val_adj <0.05 ) %>% top_n(n = 100, wt = -avg_log2FC)
  DEG.text=rbind(DEG.text, tmp);
}
dim(DEG.text)
table(DEG.text$label, DEG.text$group1)
head(DEG.text)

p1 =  ggplot() +
  geom_jitter(data = DEG, aes(group1, y=avg_log2FC, color=label), size = 0.85,width =0.4)+
  geom_jitter(data = DEG.text, aes(x = group1, y = avg_log2FC, color = label), size = 1,width =0.4)


# 火山图展示差异基因
sub_DP <- RunDEtest(srt = sub_DP, group_by = "Cell_Type4", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = sub_DP, group_by = "Cell_Type")

sub_DP <- readRDS("~/12_new data_wntnotch_dlx6as1_paper_figure _20240701/图4.DP亚群_DEtest.rds")
DEGs <- sub_DP@tools$DEtest_Cell_Type4$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 0.5 & p_val_adj < 0.05), ]
ht <- FeatureHeatmap(
  srt = sub_DP, group.by = "Cell_Type4", features = DEGs$gene, feature_split = DEGs$group1,
  species = "Homo_sapiens", db = c("GO_BP", "KEGG"), anno_terms = TRUE,
  height = 3, width = 1
)


#DP亚群细胞丰度
freq_table <- table(sub_DP@meta.data$sample, sub_DP@meta.data$Cell_Type4)
prop_table <- prop.table(freq_table, margin = 1)
prop_table <- read_csv("prop_table.csv")#修改过细胞比例的新表
ggplot(prop_table_long, aes(x = ...1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cell_colors) +
  labs(x = "时间点", y = "比例", fill = "集群") +
  theme_minimal()

#火山图比较阴阳两群细胞
updown.markers <- FindMarkers(sub_DP, only.pos = FALSE, ident.1 = "DLX6-AS1+ DP", ident.2 = "DLX6-AS1- DP", min.pct = 0.25, logfc.threshold = 0.25)
updown.markers <- updown.markers %>% mutate(gene_name = rownames(updown.markers))
updown.markers <- updown.markers %>% mutate(NegLog10Pvalue = -log10(p_val_adj))
sign_thres = 1
updown.markers <- updown.markers %>% mutate(
  color = case_when(
    p_val_adj >= 0.05 ~ "NoChange",
    avg_log2FC > sign_thres ~ "Upregulated",
    avg_log2FC <= -sign_thres ~ "Downregulated",
    TRUE ~ "NoChange"
  )
)
volcano_plot <- ggplot(updown.markers, aes(x = avg_log2FC, y = NegLog10Pvalue, color = color)) +
  geom_point(alpha = 0.6, size = 3) +
  geom_vline(xintercept = c(-sign_thres, sign_thres), linetype = "dashed") +  # 添加垂直虚线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # 添加水平虚线
  geom_text(data = top_up_genes, aes(label = gene_name), vjust = -2, size = 3) +  # 添加基因名标签
  geom_text(data = top_down_genes, aes(label = gene_name), vjust = -2, size = 3) +  # 添加基因名标签
  
  scale_color_manual(values = c("Upregulated" = "#8B1A1A", "Downregulated" = "#104E8B", "NoChange" = "grey")) +
  theme_classic() +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  ggtitle("Volcano plot of differentially expressed genes") +
  coord_cartesian(ylim = c(0, p_val_upbound))  # 设置y轴的范围
print(volcano_plot)





