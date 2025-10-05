library(Seurat)
library(Matrix)
library(SeuratObject)

scRNA_harmony <- readRDS("~/24new_human_all.rds")
Idents(scRNA_harmony) <- "CellType"
selected_DEDP <- c("Epithelium Cell", "Dental Papilla")
selected_DP <- c("Dental Papilla")
selected_DE <- c("Epithelium Cell")
sub_DEDP <- subset(scRNA_harmony, idents = selected_DEDP)
cluster1_cells <- WhichCells(sub_DEDP, ident = "Epithelium Cell")
cluster2_cells <- WhichCells(sub_DEDP, ident = "Dental Papilla")
genes_of_interest <- c( "WNT1", 
                        "WNT3A", "WNT4",
                        "WNT5A", "WNT5B",
                        "WNT6", "WNT7A",
                        "WNT7B", "WNT8B",
                        "WNT10A", 
                        "WNT10B", "WNT11",
                        
                        "DKK1", "DKK2",
                        "DKK3", "DKK4",
                        
                        "WIF1",
                        
                        "SFRP1",  "SFRP2", 
                        "SFRP3", 
                        "SFRP4", "SFRP5", 
                        "SOSTDC1",
                        
                        "FZD1", "FZD2",
                        "FZD3",  
                        "LRP5", "LRP6",
                        
                        "JAG1", "JAG2", 
                        "DLL1", "DLL2", "DLL3",
                        
                        "NOTCH1", "NOTCH2", "NOTCH3",
                        "NOTCH4")
cluster1_expr <- FetchData(sub_DEDP, vars = genes_of_interest, cells = cluster1_cells)
cluster2_expr <- FetchData(sub_DEDP, vars = genes_of_interest, cells = cluster2_cells)
expr_data <- rbind(cluster1_expr, cluster2_expr)
expr_data_long <- melt(expr_data, id.vars = "Cluster", variable.name = "Gene", value.name = "Expression")

#点图
p <- ggplot(expr_data_long, aes(x = Gene, y = Expression, color = Cluster)) +
  geom_point(position = position_jitterdodge(), size = 1, alpha = 0.7) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  scale_color_manual(values = c("#CE5E5E", "#91BAFF")) +
  coord_cartesian(ylim = c(min(0), max(expr_data_long$Expression))) +
  labs(x = "Gene", y = "Expression", title = "Gene Expression Comparison between Two Subgroups")

#饼图
gene_cluster_summary <- aggregate(Expression ~ Gene + Cluster, data = expr_data_long, sum)
gene_cluster_summary_wide <- dcast(gene_cluster_summary, Gene ~ Cluster, value.var = "Expression")
gene_cluster_summary_wide$total <- gene_cluster_summary_wide$Cluster1 + gene_cluster_summary_wide$Cluster2
gene_cluster_summary_wide$Cluster1_percent <- gene_cluster_summary_wide$Cluster1 / gene_cluster_summary_wide$total
gene_cluster_summary_wide$Cluster2_percent <- gene_cluster_summary_wide$Cluster2 / gene_cluster_summary_wide$total
draw_pie <- function(gene_name) {
  data <- data.frame(cluster = c("Cluster1", "Cluster2"),
                     percent = c(gene_cluster_summary_wide[gene_cluster_summary_wide$Gene == gene_name, ]$Cluster1_percent,
                                 gene_cluster_summary_wide[gene_cluster_summary_wide$Gene == gene_name, ]$Cluster2_percent))
  ggplot(data, aes(x = "", y = percent, fill = cluster)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    labs(title = paste("Gene:", gene_name)) +
    scale_fill_manual(values = c("#CE5E5E", "#91BAFF")) +
    theme_void()
}

# 热图
p<- DotPlot(object = sub_DE, 
            features =c("WNT1", 
                        "WNT3A", "WNT4",
                        "WNT5A", "WNT5B",
                        "WNT6", "WNT7A",
                        "WNT7B", "WNT8B",
                        "WNT10A", 
                        "WNT10B", "WNT11",
                        
                        "SOSTDC1",
                        
                        "JAG1", "JAG2", 
                        "DLL1", "DLL2", "DLL3"),
            assay = NULL,
            cols = c("lightgrey", "blue"),
            col.min = -2.5,
            col.max = 2.5,
            dot.min = 0,
            dot.scale = 6,
            idents = NULL,
            group.by = c("sample_cluster"),
            split.by = NULL,
            cluster.idents = FALSE,
            scale = TRUE,
            scale.by = "radius",
            scale.min = NA,
            scale.max = NA)

p<- DotPlot(object = sub_DP, 
            features =c("DKK1", "DKK2", 
                        "DKK3", "DKK4",  
                        "WIF1",
                        
                        "SFRP1",  "SFRP2", 
                        "SFRP3", 
                        "SFRP4", "SFRP5", 
                        
                        "FZD1", "FZD2",
                        "FZD3",  
                        "LRP5", "LRP6",
                        
                        
                        "NOTCH1", "NOTCH2", "NOTCH3",
                        "NOTCH4"),
            assay = NULL,
            cols = c("lightgrey", "blue"),
            col.min = -2.5,
            col.max = 2.5,
            dot.min = 0,
            dot.scale = 6,
            idents = NULL,
            group.by = c("sample_cluster"),
            split.by = NULL,
            cluster.idents = FALSE,
            scale = TRUE,
            scale.by = "radius",
            scale.min = NA,
            scale.max = NA)

df<- p$data
exp_mat<-df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()
percent_mat<-df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 
row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

expression_long <- melt(exp_mat, varnames = c("Cluster", "Gene"), value.name = "Expression")
percent_long <- melt(percent_mat, varnames = c("Cluster", "Gene"), value.name = "Percent")
plot_data <- merge(expression_long, percent_long, by = c("Cluster", "Gene"))
p <- ggplot(plot_data, aes(x = Gene, y = Cluster)) +
  geom_point(aes(size = Percent, color = Expression), shape = 21, stroke = 0.5) +
  scale_size_continuous(range = c(1, 7), name = "Expression Percent") +
  scale_color_gradient(low = "#990000", high = "#FF9999", name = "Expression Level") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(title = "Dot Plot", x = "Gene", y = "Cluster")
dev.off()

#基因空间映射
tooth_molar_15w <- readRDS("~/tooth_molar_15w.rds")
tooth_molar_15w <- RunPCA(tooth_molar_15w, assay = "SCT", verbose = FALSE)
tooth_molar_15w <- FindNeighbors(tooth_molar_15w, reduction = "pca", dims = 1:30)
tooth_molar_15w <- FindClusters(tooth_molar_15w, resolution = 1.0, verbose = FALSE)
tooth_molar_15w <- RunUMAP(tooth_molar_15w, reduction = "pca", dims = 1:30)
SpatialFeaturePlot(tooth_molar_15w, features = c("WNT6", "WNT7B", "WNT10B", "JAG1"),
                   pt.size.factor = 1, min.cutoff = 0, ncol = 3, 
                   crop = F, alpha = c(0,1)) + scale_alpha()


##########################画趋势折现图###################################
counts_matrix <- GetAssayData(sub_DEDP, assay = "RNA", slot = "counts")
normalized_data_matrix <- GetAssayData(sub_DEDP, assay = "RNA", slot = "data")
clusters_of_interest <- c(                   
  "8w-Epithelium Cell",  
  
  "10w-Epithelium Cell",              
  
  "12w-Epithelium Cell",      
  
  "15w-Epithelium Cell", 
  
  "17w-Epithelium Cell",
  
  "24w-Epithelium Cell",
  
  "Pulp-Epithelium Cell"
)
cells_of_interest <- WhichCells(sub_DEDP, expression = sample_cluster %in% clusters_of_interest)
genes_of_interest <- c(
  "WNT1", 
  "WNT3A", "WNT4",
  "WNT5A", "WNT5B",
  "WNT6", "WNT7A",
  "WNT7B", "WNT8B",
  "WNT10A", 
  "WNT10B", "WNT11",
  
  "SOSTDC1",
  
  "JAG1", "JAG2", 
  "DLL1", "DLL3"
)
gene_expression_list <- list()

for (gene in genes_of_interest) {
  # 获取每个基因的数据
  data_of_interest <- FetchData(sub_DEDP, vars = gene, cells = cells_of_interest)
  
  # 添加 sample_cluster 列
  data_of_interest$sample_cluster <- sub_DEDP@meta.data[cells_of_interest, "sample_cluster", drop = FALSE]
  
  # 计算每个群体的平均表达量
  avg_expr <- data_of_interest %>%
    group_by(sample_cluster) %>%
    summarise(average_expression = mean(.data[[gene]], na.rm = TRUE)) # 确保处理缺失值
  
  # 将数据框存储到列表中，以基因名称为键
  gene_expression_list[[gene]] <- avg_expr
}

# 查看结果，例如查看某个基因的数据
print(gene_expression_list[["WNT6"]])

# 创建一个文件夹来存储输出文件
output_dir <- "图3.DE gene_expression_outputs"
dir.create(output_dir, showWarnings = FALSE)

# 导出每个基因的数据框为 CSV 文件
for (gene in names(gene_expression_list)) {
  file_name <- paste0(output_dir, "/", gene, "_expression.csv")
  write.csv(gene_expression_list[[gene]], file = file_name, row.names = FALSE)
}


