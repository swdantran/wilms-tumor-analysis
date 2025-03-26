#BiocManager::install("ReactomePA")
#BiocManager::install("reactome.db")
install.packages("msigdbr")

library(msigdbr)
library(fgsea)
library(enrichplot)
library(scales)

### Three pathway enrichment functional databases were used:
### 1. Gene Ontology (GO) terms
### 2. Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway database
### 3. Hallmark genesets from the Molecular Signatures Database (MSigDB)



############################## GSEA GO terms ##########################

gene_list_go.prim <- resLFC_clean.prim$log2FoldChange
gene_list_go.prim
names(gene_list_go.prim) <- rownames(resLFC_clean.prim)
gene_list_go.prim <- na.omit(gene_list_go.prim)
gene_list_go.prim = sort(gene_list_go.prim, decreasing = TRUE)

gene_list_go.recu <- resLFC_clean.recu$log2FoldChange
gene_list_go.recu
names(gene_list_go.recu) <- rownames(resLFC_clean.recu)
gene_list_go.recu <- na.omit(gene_list_go.recu)
gene_list_go.recu = sort(gene_list_go.recu, decreasing = TRUE)

# biological processes
gse_bp.prim <- gseGO(gene_list_go.prim, ont = "BP", keyType = "ENSEMBL",
                OrgDb = "org.Hs.eg.db", eps = 1e-300, nPerm = 10000)
gse_bp.recu <- gseGO(gene_list_go.recu, ont = "BP", keyType = "ENSEMBL",
                     OrgDb = "org.Hs.eg.db", eps = 1e-300, nPerm = 10000)
# molecular functions
gse_mf.prim <- gseGO(gene_list_go.prim, ont = "MF", keyType = "ENSEMBL",
                OrgDb = "org.Hs.eg.db", eps = 1e-300, nPerm = 10000)
gse_mf.recu <- gseGO(gene_list_go.recu, ont = "MF", keyType = "ENSEMBL",
                     OrgDb = "org.Hs.eg.db", eps = 1e-300, nPerm = 10000)
# cellular components
gse_cc.prim <- gseGO(gene_list_go.prim, ont = "CC", keyType = "ENSEMBL",
                OrgDb = "org.Hs.eg.db", eps = 1e-300, nPerm = 10000)
gse_cc.recu <- gseGO(gene_list_go.recu, ont = "CC", keyType = "ENSEMBL",
                     OrgDb = "org.Hs.eg.db", eps = 1e-300, nPerm = 10000)
# gse_bp_simp <- simplify(gse_bp)
# View(as.data.frame(gse_bp_simp)


gse_bp_df.prim <- as.data.frame(gse_bp.prim)
gse_bp_df.prim <- gse_bp_df.prim[order(gse_bp_df.prim$NES, decreasing=TRUE), ]
gse_bp_df.prim$genes <- fetch_gene_symbols(gse_bp_df.prim$core_enrichment, "ENSEMBL")

gse_bp_df.recu <- as.data.frame(gse_bp.recu)
gse_bp_df.recu <- gse_bp_df.recu[order(gse_bp_df.recu$NES, decreasing=TRUE), ]
gse_bp_df.recu$genes <- fetch_gene_symbols(gse_bp_df.recu$core_enrichment, "ENSEMBL")

gse_mf_df.prim <- as.data.frame(gse_mf.prim)
gse_mf_df.prim <- gse_mf_df.prim[order(gse_mf_df.prim$NES, decreasing=TRUE), ]
gse_mf_df.prim$genes <- fetch_gene_symbols(gse_mf_df.prim$core_enrichment, "ENSEMBL")

gse_mf_df.recu <- as.data.frame(gse_mf.recu)
gse_mf_df.recu <- gse_mf_df.recu[order(gse_mf_df.recu$NES, decreasing=TRUE), ]
gse_mf_df.recu$genes <- fetch_gene_symbols(gse_mf_df.recu$core_enrichment, "ENSEMBL")

gse_cc_df.prim <- as.data.frame(gse_cc.prim)
gse_cc_df.prim <- gse_cc_df.prim[order(gse_cc_df.prim$NES, decreasing=TRUE), ]
gse_cc_df.prim$genes <- fetch_gene_symbols(gse_cc_df.prim$core_enrichment, "ENSEMBL")

gse_cc_df.recu <- as.data.frame(gse_cc.recu)
gse_cc_df.recu <- gse_cc_df.recu[order(gse_cc_df.recu$NES, decreasing=TRUE), ]
gse_cc_df.recu$genes <- fetch_gene_symbols(gse_cc_df.recu$core_enrichment, "ENSEMBL")


# GSEA plot
gseaplot2(gse_bp.prim, 
          geneSetID = c("GO:0000819"),
          color = "red",
          pvalue_table = TRUE)


# (OPTIONAL) select several terms for visualization
gse_bp_df_sel <- gse_bp_df[c("GO:0000819", "GO:0007059", "GO:0033046"),]
gse_bp_df_sel <- gse_bp_df_sel[order(gse_bp_df_sel$NES, decreasing=T), ]
gse_bp_df_sel$gene_count <- str_count(gse_bp_df_sel$core_enrichment, "/") + 1
gse_bp_df_sel$GeneRatio <- gse_bp_df_sel$gene_count / gse_bp_df_sel$setSize
gse_bp_df_sel$GeneRatio_str <- paste(gse_bp_df_sel$gene_count, gse_bp_df_sel$setSize, sep = " / ")
View(gse_bp_df_sel)
View(dplyr::filter(gse_bp_df_sel, gse_bp_df_sel$NES<0))


gse_mf_df_sel <- gse_mf_df[c("GO:0030527", "GO:0022824", "GO:0022835",
                             "GO:0070412", "GO:0005518", "GO:0048029"),]
gse_mf_df_sel <- gse_mf_df_sel[order(gse_mf_df_sel$NES, decreasing=T), ]
gse_mf_df_sel$gene_count <- str_count(gse_mf_df_sel$core_enrichment, "/") + 1
gse_mf_df_sel$GeneRatio <- gse_mf_df_sel$gene_count / gse_mf_df_sel$setSize
gse_mf_df_sel$GeneRatio_str <- paste(gse_mf_df_sel$gene_count, gse_mf_df_sel$setSize, sep = " / ")
View(gse_mf_df_sel)


gse_cc_df_sel <- gse_cc_df[c("GO:0062023", "GO:0030312", "GO:0031012", 
                             "GO:0005604", "GO:0031093", "GO:0031091"),]
gse_cc_df_sel <- gse_cc_df_sel[order(gse_cc_df_sel$NES, decreasing=T), ]
gse_cc_df_sel$gene_count <- str_count(gse_cc_df_sel$core_enrichment, "/") + 1
gse_cc_df_sel$GeneRatio <- gse_cc_df_sel$gene_count / gse_cc_df_sel$setSize
gse_cc_df_sel$GeneRatio_str <- paste(gse_cc_df_sel$gene_count, gse_cc_df_sel$setSize, sep = " / ")
View(gse_cc_df_sel)


####### LOLLIPOP PLOT 

ggplot(gse_bp_df_sel, aes(x = NES, y = reorder(Description, NES, decreasing = F))) +
  geom_segment(aes(xend = 0, yend = Description), color = "black") +
  geom_point(aes(size = gene_count, fill = p.adjust), shape = 21, color = "black") +
  scale_fill_gradient(low = "deepskyblue4", high = "white") + #brown3 #seagreen #deepskyblue4
  scale_size(range = c(6, 14), breaks = c(10, 40, 70, 100)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05)),
                   labels = label_wrap(150)) +
  theme_minimal() +
  theme(
    #axis.text.y = element_text(color = colors_common_gse),
    axis.line.x=element_line(color="black"),
    axis.line.y=element_line(color="black"),
    panel.grid.major = element_line(color="white"),  
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(size = 14, color = "black"),  
    axis.title = element_text(size = 14, color = "black", face = "bold"), 
    axis.title.y = element_blank(),
    legend.text = element_text(size = 14, color = "black"), 
    legend.title = element_text(size = 14, face = "bold", color = "black")  
  ) +
  labs(x = "Normalized enrichment score (NES)", y = NULL, 
       size = "Gene hit", fill = "FDR")







############################## GSEA KEGG ##############################

gene_list_kegg.prim <- resLFC_clean.prim$log2FoldChange
names(gene_list_kegg.prim) <- resLFC_clean.prim$entrez
gene_list_kegg.prim <- na.omit(gene_list_kegg.prim)
gene_list_kegg.prim = sort(gene_list_kegg.prim, decreasing = TRUE)

gene_list_kegg.recu <- resLFC_clean.recu$log2FoldChange
names(gene_list_kegg.recu) <- resLFC_clean.recu$entrez
gene_list_kegg.recu <- na.omit(gene_list_kegg.recu)
gene_list_kegg.recu = sort(gene_list_kegg.recu, decreasing = TRUE)

gse_kegg.prim <- gseKEGG(gene_list_kegg.prim,
                    organism = "hsa",
                    keyType = "ncbi-geneid",
                    nPerm = 10000,
                    eps = 1e-300,
                    minGSSize    = 1,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.1)

gse_kegg_df.prim <- gse_kegg.prim@result
View(gse_kegg_df.prim)

gse_kegg.recu <- gseKEGG(gene_list_kegg.recu,
                         organism = "hsa",
                         keyType = "ncbi-geneid",
                         nPerm = 10000,
                         eps = 1e-300,
                         minGSSize    = 1,
                         maxGSSize    = 500,
                         pvalueCutoff = 0.1)

gse_kegg_df.recu <- gse_kegg.recu@result
View(gse_kegg_df.recu)

# GSEA plot
gseaplot2(gse_kegg, 
          geneSetID = c("hsa04066"),# ,"hsa04066", "hsa04926","hsa04512","hsa04974", "hsa04820", "hsa04510", "hsa00010"),
          color = c("royalblue3"), #, "red3", "darkgoldenrod2"
          #pvalue_table = T,
          base_size = 14,
          subplots = 1:2
)

# (OPTIONAL) select several pathways for visualization
gse_kegg_df_sel <- gse_kegg_df[c("hsa04066", "hsa04820", "hsa04610", "hsa04974","hsa04926"),]
View(gse_kegg_df_sel)
gse_kegg_df_sel <- gse_kegg_df_sel[order(gse_kegg_df_sel$NES, decreasing=T), ]
gse_kegg_df_sel$gene_count <- str_count(gse_kegg_df_sel$core_enrichment, "/") + 1
gse_kegg_df_sel$geneNames <- fetch_gene_symbols(gse_kegg_df_sel$core_enrichment, "ENTREZID")
# write.csv(gse_kegg_df_sel, file="gse_kegg_df_selected.csv")
dev.off()


# CNETPLOT
cnetplot(gse_kegg,
         #categorySize="pvalue", 
         showCategory = 6,
         vertex.label.font=5,
         node_label = 'none',
         circular = F,
         layout = 'kk', #'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl'
         color.params = list(category = "grey60", #"deepskyblue4", 
                             edge = T,
                             foldChange = gene_list_kegg_heat #"red3"
         ),
         hilight.params = list(category = NULL, # select pathway name here for highlight
                               alpha_hilight = 2, 
                               alpha_no_hilight = 0.5),
         cex.params = list(gene_node = 4)
) +
  theme(
    legend.text = element_text(colour="black", size=13), # face = "bold"
    legend.title = element_text(colour="black", size=13, face = "bold"), # , vjust = 0.9
  )  + scale_color_gradient2(name='foldChange', 
                             low = '#FFC3B0', 
                             high='red3',
                             na.value = "grey80")


# HEATPLOT
heatplot(gse_kegg, #foldChange = gene_list_kegg_heat,
         showCategory = c("HIF-1 signaling pathway",
                          "Glycolysis / Gluconeogenesis")
) +
  theme(
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line = element_line(color = "black"), 
    axis.text = element_text(color = "black", size = 12),  
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.y = element_text(size = 12, vjust = 1, face = "bold"),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 10), 
    legend.title = element_text(size = 12, face = "bold")  
  )


####### LOLLIPOP PLOT 

gse_kegg_df_sel <- gse_kegg_df_sel[order(gse_kegg_df_sel$NES, decreasing=T), ]

ggplot(gse_kegg_df_sel, aes(x = NES, y = reorder(Description, NES, decreasing = F))) +
  geom_segment(aes(xend = 0, yend = Description), color = "black") +
  geom_point(aes(size = gene_count, fill = p.adjust), shape = 21, color = "black") +
  scale_fill_gradient(low = "deepskyblue3", high = "white") +
  scale_size(range = c(3, 12), breaks = c(20, 40, 60, 80)) +
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.02))) +
  theme_minimal() +
  theme(
    #axis.text.y = element_text(color = colors_common_gse),
    axis.line.x=element_line(color="black"),
    axis.line.y=element_line(color="black"),
    panel.grid.major = element_line(color="grey97"),  
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),  
    axis.text = element_text(size = 14, color = "black"),  
    axis.title = element_text(size = 14, color = "black", face = "bold"), 
    axis.title.y = element_blank(),
    legend.text = element_text(size = 14, color = "black"),  
    legend.title = element_text(size = 14, face = "bold", color = "black")  
  ) +
  labs(x = "Normalized enrichment score (NES)", y = NULL, size = "Gene hit", fill = "FDR")



######## PATHVIEW VISUALIZATION OF KEGG PATHWAYS
# Define plotting function
plot_pathway = function(pid) pathview(gene.data=gene_list_kegg, 
                                      pathway.id=pid, 
                                      species="hsa", 
                                      new.signature=FALSE)


# selected 7 up-regulated pathways for Pathview visualization
kegg_up <- c("hsa04066", "hsa04820", "hsa04974", "hsa04926", "hsa04512", "hsa00010", "hsa04510")
# selected 6 down-regulated pathways for Pathview visualization
kegg_down <- c("hsa00190", "hsa03050", "hsa03008", "hsa03030", "hsa03040", "hsa03013", "hsa03020", "hsa00280", "hsa04110")

sapply(kegg_up, function(pid) pathview(gene.data = gene_list_kegg, 
                                       pathway.id = pid, 
                                       species = "hsa",
                                       kegg.dir = "./kegg_gsea/up/raw",
                                       #kegg.native = F,
                                       #expand.node = T,
                                       low = list(gene = "blue", cpd = "green"), 
                                       mid = list(gene = "gray", cpd = "gray"), 
                                       high = list(gene = "red", cpd = "yellow"),
                                       na.col = "white"))

sapply(kegg_down, function(pid) pathview(gene.data = gene_list_kegg, 
                                         pathway.id = pid, 
                                         species = "hsa",
                                         kegg.dir = "./kegg_gsea/down/raw",
                                         #kegg.native = F,
                                         #expand.node = T,
                                         low = list(gene = "blue", cpd = "green"), 
                                         mid = list(gene = "gray", cpd = "gray"), 
                                         high = list(gene = "red", cpd = "yellow"),
                                         na.col = "white"))



############################## GSEA HALLMARKS ##############################

hallmark_gene_sets <- msigdbr(species = "human", category = "H")
human_gene_sets <- msigdbr(species = "human")
# fixing format to work with fgsea
hallmark_gene_sets_fix = split(x = hallmark_gene_sets$entrez_gene, f = hallmark_gene_sets$gs_name)
gene_list_hall <- gene_list_kegg

fgseaRes <- fgsea(pathways=hallmark_gene_sets_fix, 
                  gene_list_hall)
fgseaRes <- fgseaRes[order(fgseaRes$NES, decreasing=TRUE), ]
fgseaRes_sel <- fgseaRes[c(1:12, 45:50),]
fgseaRes_sel <- fgseaRes_sel[order(fgseaRes_sel$NES, decreasing=T), ]


###### Lolipop
ggplot(fgseaRes_sel, aes(x = NES, y = reorder(pathway, NES, decreasing = F))) +
  geom_segment(aes(xend = 0, yend = pathway), color = "black") +
  geom_point(aes(size = size, fill = padj), shape = 21, color = "black") +
  scale_fill_gradient(low="#675bb5", high="#d4d1eb") + 
  scale_size(range = c(6, 14), breaks = c(50, 100, 150, 200)) + 
  scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  theme_minimal() +
  theme(
    #axis.text.y = element_text(color = colors_sel_hm),
    axis.line.x=element_line(color="black"),
    axis.line.y=element_line(color="black"),
    panel.grid.major = element_line(color="grey97"), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"), 
    axis.text = element_text(size = 14, color = "black"), 
    axis.title = element_text(size = 14, color = "black", face = "bold"), 
    axis.title.y = element_blank(),
    legend.text = element_text(size = 14, color = "black"), 
    legend.title = element_text(size = 14, face = "bold", color = "black")  
  ) +
  labs(x = "Normalized enrichment score (NES)", y = NULL, 
       size = "Gene hit", fill = "FDR")

# Function to fetch gene symbols for Gene IDs
fetch_gene_symbols <- function(my_ids, id_type) { # id_type is either ENTREZID or ENSEMBL  
  # Split multiple NCBI Gene IDs in a cell separated by '/'  
  id_list <- unlist(strsplit(my_ids, "/"))  
  # Fetch gene symbols from the annotation database  
  gene_symbols <- mapIds(org.Hs.eg.db, keys = id_list, column = "SYMBOL", keytype = id_type)  
  # Combine gene symbols corresponding to multiple NCBI Gene IDs in a cell, separated by '/'  
  symbol_combined <- sapply(strsplit(my_ids, "/"), function(ids) paste(gene_symbols[ids], collapse = "/"))  
  return(symbol_combined)}



