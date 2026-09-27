#load library
library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)
library(magick)
library(Cairo)
library(ComplexHeatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(EnhancedVolcano)
library(ggrepel)
library(textshaping)
library(dplyr)
library(apeglm)
library(circlize)
library(clusterProfiler)
library(pheatmap)
library(RColorBrewer)
library(stringr)

# read the count file 
counts_data0 <- read.csv("merged_counts.csv", row.names = 1)

# read sample info
colData0 <- read.csv("sample_info.tsv", sep="\t", row.names = 1)

norm <- colData0[colData0$condition == "NORM",]$sample
prim <- colData0[colData0$condition == "PRIM",]$sample
recu <- colData0[colData0$condition == "RECU",]$sample
colData <- colData0[c(norm, recu),]
counts_data <- counts_data0[,c(norm,recu)]

######### load deseq2 ########
rownames(counts_data) 
colnames(counts_data)

# make sure the rownames in colData match column names in counts_data
all(colnames(counts_data) %in% rownames(colData))
#make sure tey are in the same order 
all(colnames(counts_data) == rownames(colData))

# construct DESeq object 
dds0 <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~condition)

dds0 # 60616

# keeping genes that meet the following minimum count
keep <- rowSums(counts(dds0)) >= 100
dds <- dds0[keep,]

dds # 26715

# contrasts:
resultsNames(dds)

# set the control sample
dds$condition <- relevel(dds$condition, ref = "NORM")
dds$condition

# get normalized counts 
# Estimate size factors
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized =TRUE)

# get statistics 
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05) #setting alpha threshold = 0.05
summary(res)

# dispersion plot
plotMA(res)
plotDispEsts(dds)
vsd <- vst(dds)
plotPCA(vsd, intgroup="condition") +
  geom_point(size = 5) +
  scale_color_manual(values = c("NORM" = "black", "RECU" = "red")) +
  theme_minimal() +
  theme(
    axis.text.y   = element_text(size=12, colour = "black"),
    axis.text.x   = element_text(size=12, colour = "black"),
    axis.title.y  = element_text(size=12, colour = "black"),
    axis.title.x  = element_text(size=12, colour = "black"),
    legend.text = element_text(size=11),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth = 0.5)
  )

######### heatmap #############
 
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2])
resLFC <- as.data.frame(resLFC)

# getting gene symbol, gene name and entrez ID
resLFC$symbol <- mapIds(org.Hs.eg.db, keys = rownames(resLFC), keytype = "ENSEMBL", column = "SYMBOL")
resLFC$name <- mapIds(org.Hs.eg.db, keys=rownames(resLFC), keytype = "ENSEMBL", column = 'GENENAME', multiVals="first")
resLFC$entrez <- mapIds(org.Hs.eg.db, keys=rownames(resLFC), keytype = "ENSEMBL", column = 'ENTREZID', multiVals="first")
# removing uncharacterized genes 
dim(resLFC)
resLFC_clean <- resLFC[is.na(resLFC$symbol) == F, ]
dim(resLFC_clean) # 

# get uncharacterized genes for further inspection
unchar_resLFC <- resLFC[is.na(resLFC$symbol) == T, ]
unchar_resLFC_DE <- dplyr::filter(unchar_resLFC, (unchar_resLFC$padj < 0.05) & (abs(unchar_resLFC$log2FoldChange) > 1))
write.csv(unchar_resLFC_DE, file="uncharacterized_genes.csv")


##### Get differentially expressed genes (DEGs)
filtered <- resLFC_clean %>% 
  dplyr::filter((resLFC_clean$padj < 0.05) & (abs(resLFC_clean$log2FoldChange) > 1))

dim(filtered)
dim(filtered[filtered$log2FoldChange > 1,]) # upregulated DEGs
View(filtered[filtered$log2FoldChange > 1,])
dim(filtered[filtered$log2FoldChange < 1,]) # downregulated DEGs

##### heatmap 
rlog_out <- rlog(dds, blind=FALSE) # get normalized counts
mat <- assay(rlog_out)[rownames(filtered), rownames(colData)]
colnames(mat) <- rownames(colData)
mat.scaled <- t(apply(mat, 1, scale)) # center and scale each column (Z-score) then transpose
colnames(mat.scaled) <- colnames(mat)

# Get 50 top DE:
num_keep <- 50 # keep 25 highest and lowest log2FoldChange -> total 50
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled) - num_keep), 
                                    nrow(mat.scaled) ))

# get log2FC for each gene
l2_val <- as.matrix(filtered[rows_keep,]$log2FoldChange)
colnames(l2_val) <- "logFC"
# get baseMean for each gene
mean <- as.matrix(filtered[rows_keep,]$baseMean)
colnames(mean) <- "AveExpr"

col_logFC <- colorRamp2(c(min(l2_val), -0.8, max(l2_val)), c("mediumseagreen", "white", "red3"))
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red2"))
col_Z <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# Ensure the colData has the right levels
colData$condition <- factor(colData$condition, levels = c("NORM", "RECU"))

# the heatmap:
ha <- columnAnnotation(Condition = colData$condition,
                       col = list(Condition = c("NORM" = "grey35", "RECU" = "grey")),
                       show_annotation_name = F,
                       annotation_legend_param = list(
                         labels_gp = gpar(fontsize = 10),
                         title_gp = gpar(fontsize = 10, fontface = "bold"),
                         #at = c("co-cultured", "control"),
                         #labels = c("low", "zero", "high"),
                         #title = "Some values",
                         legend_height = unit(4, "cm"),
                         grid_width = unit(0.8, "cm")
                         #title_position = "leftbot-rot"
                       )) #+ Legend(at = colData[,]$condition, nrow = 1)

h1 <- Heatmap(mat.scaled[rows_keep,], 
              col = col_Z,
              cluster_rows = T,
              cluster_columns = T, 
              column_labels = colnames(mat.scaled), 
              name="Z-score",
              show_row_names = T,
              row_names_gp = gpar(fontsize = 8),
              show_column_names = F,
              #left_annotation = ha_rows,
              bottom_annotation = ha,
              row_title_gp = gpar(fontsize = 3.5),
              column_names_gp = gpar(fontsize = 12),
              column_title_gp = gpar(fontsize = 12),
              column_names_rot = 45,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 10),
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                # at = c(-2, -1, 0, 1, 2),
                #labels = c("low", "zero", "high"),
                #title = "Some values",
                legend_height = unit(4, "cm"),
                grid_width = unit(0.8, "cm"),
                direction = "horizontal"
                #title_position = "lefttop"
              )
) 


h2 <- Heatmap(l2_val,
              row_labels = filtered$symbol[rows_keep],
              row_names_gp = gpar(fontsize = 9),
              cluster_rows = T,
              name="Log2FC", 
              #top_annotation = ha,
              col = col_logFC,
              width = unit(10, "mm"),
              #cell_fun = function(j, i, x, y, w, h, col) { #add text to each grid
              #  grid.text(round(l2_val[i, j], 2), x, y,
              #            gp=gpar(col="grey90"))
              #},
              column_names_rot = 90,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 10),
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                at = seq(round(min(l2_val)), round(max(l2_val)), by = 2),
                #labels = c("low", "zero", "high"),
                #title = "Some values",
                # legend_height = unit(4, "cm"),
                grid_width = unit(0.8, "cm"),
                direction = "horizontal"
                #title_position = "lefttop-rot"
              )
)


h3 <- Heatmap(mean, row_labels = filtered$symbol[rows_keep],
              cluster_rows = F,
              name = "AveExpr",
              col = col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { #add text to each grid
                grid.text(round(mean[i, j], 2), x, y)
              },
              width = unit(10, "mm"),
              column_names_rot = 45,
              heatmap_legend_param = list(
                labels_gp = gpar(fontsize = 12),
                title_gp = gpar(fontsize = 12, fontface = "bold"),
                #at = seq(0, round(max(mean)), by = 5000),
                #labels = c("low", "zero", "high"),
                #title = "Some values",
                legend_height = unit(4, "cm"),
                grid_width = unit(0.8, "cm"),
                direction = "horizontal"
                #title_position = "lefttop-rot"
              ))

h1+h2+h3 # save at width 900 x height 1600
h1 + h2

######### volcano plot ########
# Check if resLFC$symbol is a character vector
resLFC$symbol <- as.character(resLFC$symbol)
options(repr.plot.width = 2, repr.plot.height = 3)
EnhancedVolcano(resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = resLFC$symbol,
                title = '',
                subtitle = '',
                caption = '', #paste0(nrow(resLFC), " genes in total"),
                titleLabSize = 16, subtitleLabSize = 16, captionLabSize = 16,
                pointSize = 6,
                labSize = 4, axisLabSize = 16,
                #xlim = c(-4, 10), ylim = c(-3, 170),
                legendLabSize = 14,
                legendIconSize = 7,
                pCutoff = 1e-10,
                cutoffLineWidth = 0.4,
                shape = 20,
                gridlines.major = F,
                gridlines.minor = T,
                #drawConnectors = T,
                #legendLabSize = 12,
                #legendIconSize = 4,
                legendPosition = "right",
                legendLabels = c("Not significant", "|Log2FC| > 1", "P-value < 0.05", "P-value < 0.05 and |Log2FC| > 1")
)

######### GSEA ########
res <- res[order(-res$stat),]
gene_list <- res$stat
names(gene_list) <- rownames(res)
gene_list

gse <- gseGO(gene_list,
              ont = "BP",
              keyType = "ENSEMBL",
              OrgDb = "org.Hs.eg.db",
              eps = 1e-300)
gse.df <- as.data.frame(gse)
fit <- gseaplot(gse, geneSetID = 1)
png("gsea_recurrent.png", res = 250, width = 2000, height = 1300)
print(fit)
dev.off()
fit

