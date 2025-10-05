
#R3_Q2
scRNA_harmony <- readRDS("~/24new_human_all.rds")
CellDimPlot(srt = scRNA_harmony, 
                    group.by = c("sample"), 
                    reduction = "umap", 
                    theme_use = "theme_blank", label = F)
p_split_feature <- FeaturePlot(
  scRNA_harmony,
  features = c("KIF5C", "SOX9", "PTCH1"),
  split.by = "sample", # 核心参数
  pt.size = 0.5,
  combine = TRUE, # Seurat会自动排列
  ncol = 4 # 
)
desired_order <- c("8w", "10w", "12w", "15w", "17w", "24w", "Pulp")
sample_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2")
names(sample_colors) <- desired_order
scRNA_harmony$sample <- factor(scRNA_harmony$sample, levels = desired_order)
Idents(scRNA_harmony) <- "CellType"
p_split_vln <- VlnPlot(
  scRNA_harmony,
  features = c("S100B", "MPZ", "PLP1"),
  group.by = "CellType",
  split.by = "sample",
  cols = sample_colors, # 直接在这里指定颜色
  pt.size = 0,
  idents = "Neuronal Cell"
) + 
  theme(legend.position = "bottom")
ggsave("Vln_Neuronal Cell.png", plot = p_split_vln, width = 12, height = 7)


#R3_Q6
Idents(scRNA_harmony) <- "CellType"
all.markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25, 
                              group.by = "CellType")
top10 <- all.markers%>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
heatmap_plot <- DoHeatmap(
  scRNA_harmony, 
  features = top10$gene, 
  group.by = "CellType", # 确保热图的列是按CellType组织的
  label = TRUE # 在热图顶部显示细胞类型标签
) + NoLegend() # 移除图例，使热图更清晰

Idents(sub_DE) <- "cell_type"
DE_all.markers <- FindAllMarkers(sub_DE, only.pos = TRUE, 
                                 min.pct = 0.25, logfc.threshold = 0.25, 
                                 group.by = "cell_type")
DE_top10 <- DE_all.markers%>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

heatmap_plot <- DoHeatmap(
  sub_DE, 
  features = DE_top10$gene, 
  group.by = "cell_type", 
  label = TRUE 
) + NoLegend() 

#R3_Q8
pathway_list <- list(
  WNT = list(Ligand = wnt_ligands, Receptor = wnt_receptors),
  BMP = list(Ligand = bmp_ligands, Receptor = bmp_receptors),
  FGF = list(Ligand = fgf_ligands, Receptor = fgf_receptors),
  SHH = list(Ligand = shh_ligands, Receptor = shh_receptors),
  EDA = list(Ligand = eda_ligands, Receptor = eda_receptors),
  NOTCH = list(Ligand = notch_ligands, Receptor = notch_receptors),
  TGF = list(Ligand = tgf_ligands, Receptor = tgf_receptors)
)

gene_info_df <- do.call(rbind, lapply(names(pathway_list), function(pathway_name) {
  pathway_data <- pathway_list[[pathway_name]]
  ligands <- data.frame(Pathway = pathway_name, Gene = pathway_data$Ligand, Gene_Type = "Ligand", stringsAsFactors = FALSE)
  receptors <- data.frame(Pathway = pathway_name, Gene = pathway_data$Receptor, Gene_Type = "Receptor", stringsAsFactors = FALSE)
  rbind(ligands, receptors)
}))

de_results <- FindMarkers(sub_DEDP, 
                          ident.1 = "Epithelium Cell", 
                          ident.2 = "Dental Papilla", 
                          features = all_genes_of_interest,
                          logfc.threshold = 0,
                          min.pct = 0,
                          test.use = "wilcox",
                          assay = "RNA") # 使用标准的Wilcoxon秩和检验
de_results <- rownames_to_column(de_results, var = "Gene")

avg_expr <- AverageExpression(sub_DEDP, 
                              features = all_genes_of_interest, 
                              group.by = "CellType",
                              assay = "RNA")
avg_expr_df <- as.data.frame(avg_expr$RNA) # 假设您的assay是"RNA"
avg_expr_df <- rownames_to_column(avg_expr_df, var = "Gene")

final_table <- left_join(gene_info_df, de_results, by = "Gene")
final_table <- left_join(final_table, avg_expr_df, by = "Gene")
final_supplementary_table <- final_table %>%
  select(
    Pathway,
    Gene,
    Gene_Type,
    DE_Avg_Expr = `Epithelium Cell`, # 重命名列
    DPL_Avg_Expr = `Dental Papilla`, # 重命名列
    Log2FC = avg_log2FC,
    p_val_adj
  ) %>%
  # 按照通路和基因类型排序，使表格更美观
  arrange(Pathway, desc(Gene_Type), Gene)

Idents(scRNA_harmony) <- "CellType"
selected_DEDP <- c("Epithelium Cell", "Dental Papilla")
sub_DEDP <- subset(scRNA_harmony, idents = selected_DEDP)
sub_DEDP <- NormalizeData(sub_DEDP, normalization.method = "LogNormalize", scale.factors=10000)

# 步骤 1: 定义通路相关的基因集 (Hallmark Gene Sets)
# 使用MSigDB中广泛接受的Hallmark基因集来代表通路活性。
# 注意：这些是人类基因符号。请确保您的Seurat对象中的基因名是这个格式。
# EDA通路没有标准的Hallmark集，我们使用一个文献中公认的核心基因列表。
pathway_gene_sets <- list(
  "BMP_Pathway" = c("BMP2", "BMP3", "BMP4", "BMP5", "BMP6", "BMP7", "BMP10", "BMP15",
                    "BMPR2", 
                    "ID1", "ID2", "ID3", "ID4", "SMAD1", "SMAD6", "SMAD7", "SMAD9", "NOG", "GREM1"),
  "EDA_Pathway" = c("EDA", "EDAR", 
                    "EDARADD", "NFKBIA", "RELB", "LTB"),
  "WNT_Pathway" = c("WNT1", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8B", "WNT10A", "WNT10B", "WNT11",
                    "FZD1","FZD2","FZD3","FZD4","FZD5","FZD6","FZD7","FZD8","FZD9","FZD10","LRP5", "LRP6", 
                    "CTNNB1", "LEF1", "TCF7", "AXIN2", "NKD1"),
  "FGF_Pathway" = c("FGF2", "FGF3", "FGF8", "FGF10",
                    "FGFR1", "FGFR2", "FGFR3", "FGFR4",
                    "SPRY2", "SPRY4", "ETV4", "ETV5", "DUSP6"),
  "SHH_Pathway" = c("SHH", "IHH", "DHH", 
                    "PTCH1", "PTCH2", 
                    "SMO", "GLI1", "GLI2", "GLI3", "HHIP"),
  "NOTCH_Pathway" = c("JAG2", "JAG1", "DLL1", "DLL3", "DLL4", 
                      "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", 
                      "HES1", "HEY1", "HEYL"), 
  "TGF_Pathway" = c("TGFB1", "TGFB2", "TGFB3", "INHBA", "INHBB",
                    "TGFBR1", "TGFBR2", "ACVR1", "ACVR2A", "ACVR2B",
                    "ARID4B", "BAMBI", "BCAR3", "BMPR2", "CCN3", "CDKN1A", 
                    "CDKN2B", "CITED1", "CREBBP", "CTNNA1", "CXXC4", "EP300", "FKBP1A",
                    "FN1", "FURIN", "FUT8", "HDAC1", "HIPK2", "ID1", "ID2", "ID3", "JUN",
                    "JUNB", "KLF10", "LEFTY2", "LEMD3", "LTBP2", "MAPK1", "NCOR2", "NEDD9", 
                    "PMEPA1", "PPM1A", "PPP1R15A", "RAB31", "RHOA", "ROCK1", "RUNX1", "SERPINE1", "SKI", 
                    "SKIL", "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9", "SMURF1", 
                    "SMURF2", "SP1", "THBS1", "TRIM33", "TSC22D1", 
                    "ZFYVE16", "ZFYVE9")
)

# 步骤 2: 循环计算每个通路的模块分数
for (pathway_name in names(pathway_gene_sets)) {
  # 提取当前通路的基因列表
  genes <- pathway_gene_sets[[pathway_name]]
  
  # 确保基因存在于数据集中，避免报错
  genes_in_data <- intersect(genes, rownames(sub_DEDP@assays$RNA@data))
  
  if (length(genes_in_data) == 0) {
    warning(paste("No genes found for pathway:", pathway_name))
    next
  }
  
  # 使用 AddModuleScore 计算通路活性
  sub_DEDP <- AddModuleScore(
    object = sub_DEDP,
    features = list(genes_in_data),
    name = pathway_name, # 输出的列名将是 "pathway_name1"
    ctrl = 100, # 对照基因数量
    seed = 1 # 设置随机种子以保证结果可重复
  )
}

violin_plot <- ggplot(plot_data_long, aes(x = CellType, y = Score, fill = CellType)) +
  geom_violin(trim = FALSE, scale = "width") + # trim=FALSE确保小提琴的尾部不被切掉
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) + # 在小提琴图上叠加一个箱线图
  
  # 使用ggpubr自动添加统计学显著性检验 (Wilcoxon test)
  stat_compare_means(
    comparisons = list(c("Epithelium Cell", "Dental Papilla")),
    method = "wilcox.test",
    label = "p.signif" # 显示 *, **, ***, ****
  ) +
  
  # 使用分面，为每个通路创建一个独立的图
  facet_wrap(~ Pathway, scales = "free_y", ncol = 4) + # scales="free_y"让每个图有自己独立的Y轴范围
  
  # 自定义颜色和主题
  scale_fill_manual(values = lineage_colors) +
  
  # 美化图形
  labs(
    title = "Pathway Activity Scores in Dental Lineages",
    x = "Cell Lineage",
    y = "Pathway Activity Score"
  ) +
  theme_classic() + # 使用经典的白色背景主题
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # 标题居中加粗
    axis.text.x = element_text(angle = 45, hjust = 1),     # X轴标签旋转45度，防止重叠
    # strip.background = element_rect(fill = "grey90", color = "black"), # 分面标题的背景
    strip.text = element_text(face = "bold"), # 分面标题的文字
    legend.position = "none" # 隐藏图例，因为颜色和X轴信息重复
  )


#R3_Q9
library(Seurat)
library(dplyr)
library(AUCell)
library(msigdbr)
library(pheatmap)
library(tibble)
set.seed(123)
exprMatrix <- GetAssayData(sub_DEDP, slot = "counts") 
unique_sample <- unique(sub_DEDP$sample)
print(unique_sample) # 检查原始顺序
chronological_order <- c("8w", "10w", "12w", "15w", "17w", "24w", "Pulp") #手动排序
sub_DEDP$sample <- factor(sub_DEDP$sample, levels = chronological_order)

获取通路基因集 (使用msigdbr) 
 pathway_keywords <- c(
   "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
   "HALLMARK_NOTCH_SIGNALING",
   "HALLMARK_HEDGEHOG_SIGNALING",
   "HALLMARK_TGF_BETA_SIGNALING",
   "REACTOME_FGF_RECEPTOR_SIGNALING_PATHWAY" # 使用Reactome中的FGF通路
 )
 
 selected_pathways <- all_gene_sets %>%
   filter(gs_name %in% pathway_keywords)
 
 pathway_gene_lists <- selected_pathways %>%
   select(gs_name, gene_symbol) %>%
   group_by(gs_name) %>%
   summarise(genes = list(gene_symbol)) %>%
   deframe() # deframe() 是一个将两列tibble转换为命名列表的便捷函数

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(pathway_gene_lists, cells_rankings)
auc_matrix <- getAUC(cells_AUC)
sub_DEDP <- AddMetaData(sub_DEDP, metadata = as.data.frame(t(auc_matrix)))
heatmap_data <- sub_DEDP@meta.data %>%
  # 确保sample是因子且顺序正确
  mutate(sample = factor(sample, levels = unique(sub_DEDP$sample)[order(unique(sub_DEDP$sample))])) %>% 
  group_by(sample) %>%
  # 计算每个通路的平均值
  summarise(across(all_of(names(pathway_gene_lists)), mean)) %>%
  # 将sample列转换为行名，以满足pheatmap的输入格式
  column_to_rownames("sample") %>%
  # 转置矩阵，使通路为行，时间点为列
  t()

desired_pathway_order <- c(
  "WNT_Pathway",
  "NOTCH_Pathway",
  "SHH_Pathway",
  "EDA_Pathway",
  "BMP_Pathway",
  "TGF_Pathway",
  "FGF_Pathway"
)
heatmap_data_ordered <- heatmap_data[desired_pathway_order, ]

color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

pheatmap(
  heatmap_data_ordered,
  main = "Pathway Activity Dynamics (Custom Order)",
  fontsize_row = 11,
  fontsize_col = 11,
  
  # 关键参数：按行进行标准化（Z-score）
  # 这能凸显每个通路自身随时间的变化趋势，而不是比较通路间的绝对活性
  scale = "row", 
  
  # 关键参数：保持行列原始顺序，不进行聚类
  cluster_cols = FALSE, 
  cluster_rows = T,
  
  # 调整边框和颜色
  border_color = "grey60",
  color = color_palette,
  angle_col = 45, # 将列标签旋转45度，防止重叠
  cellwidth = 30, # 增加单元格宽度
  cellheight = 30 # 增加单元格高度
)

color_palette <- colorRampPalette(c("lightyellow", "red"))(100)
pheatmap(
  heatmap_data_ordered,
  main = "Pathway Activity Dynamics (Custom Order)",
  fontsize_row = 11,
  fontsize_col = 11,
  cluster_cols = FALSE, 
  cluster_rows = FALSE,
  display_numbers = FALSE,
  number_format = "%.3f", # 设置数值的格式，例如保留三位小数
  fontsize_number = 8, # 调整单元格内数字的大小
  # 调整边框和颜色
  border_color = "grey60",
  color = color_palette,
  angle_col = 45,
  cellwidth = 50, # 增加单元格宽度以容纳数字
  cellheight = 30)

#R3_Q17
genes_to_plot <- c("DLX6-AS1", 
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
                   "CYP1B1")
desired_order <- c("8w", "10w", "12w", "15w", "17w", "24w", "Pulp")
plot_data <- FetchData(sub_DP, vars = c(gene, "Cell_Type", "sample"))

  colnames(plot_data)[1] <- "Expression"

  my_comparisons <- list(c("DLX6-AS1+ DP", "DLX6-AS1- DP"))

  p <- ggplot(summary_data, aes(x = Cell_Type4, y = mean_expr, fill = Cell_Type4)) +
    # 绘制条形图
    # stat="identity"表示Y值直接由数据中的'mean_expr'列决定
    geom_bar(stat = "identity", color = "black", position = position_dodge()) +
    
    # 添加误差线
    geom_errorbar(
      aes(ymin = mean_expr - se, ymax = mean_expr + se),
      width = 0.2, # 误差线顶部的宽度
      position = position_dodge(0.9)
    ) +
    
    facet_wrap(~sample, nrow = 1, scales = "free_y") +
        
    labs(
      title = paste("Expression of", gene, "Across Developmental Stages"),
      x = NULL,
      y = "Mean Expression (+/- SE)"
    ) +
    theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      strip.text = element_text(face = "bold", size = 12), # 分面标题样式
      strip.background = element_blank(), # 移除分面标题的背景框
      panel.border = element_rect(colour = "black", fill=NA, size=1), # 添加面板边框
      legend.position = "none"
    )
  


