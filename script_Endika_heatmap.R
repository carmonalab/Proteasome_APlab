##############################
######### Libraries ##########
##############################

#library(ellipsis)
#library(Seurat)
#library(SeuratDisk)
#library(dplyr)
#library(patchwork)
#library(BiocManager)
#library(metap)
#library(xlsx)
#library(plyr)
#library(tidyverse)
#library(ggplot2)
#library(renv)
#library(Matrix)
#library(ProjecTILs)
#library(gridExtra)
#library(GEOquery)
#library(devtools)
#library(presto)
#library(msigdbr)
#library(fgsea)
#library(dplyr)
#library(ggplot2)
#library(Scillus)

########################################################
######### Previously: Create the Seurat object #########
########################################################

# Tutorial  https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html

# Starting point: several .rds files from which I extracted only the CD8_Tex and CD8_EffectorMemory subpopulations (ovarian, liver, breast, lung_1 & lung_2)

#merged <- merge(x=lung_2, y= c(lung_1, breast, liver, ovarian))
#DefaultAssay(merged) <- "RNA"

#merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 3000, verbose = T)
#merged <- ScaleData(merged, verbose = FALSE)
#merged <- RunPCA(merged, npcs = 20, verbose = FALSE)
#merged <- RunUMAP(merged, reduction = "pca", dims = 1:20)
#merged <- FindNeighbors(merged, reduction = "pca", dims = 1:20)
#merged <- FindClusters(merged, resolution = 2.0)

#saveRDS(merged, "Pan_Tex_EffectorMemory_CD8.rds", refhook = NULL)

###############################################################################################
######### Load the Seurat object and remove lung_1 (excluded  after heatmap analysis) #########
###############################################################################################

merged <- readRDS("Tex_EffectorMemory_CD8.rds", refhook = NULL)

merged2 <- subset(merged, subset = cancer.type != "Lung_1")

#p1 <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters", pt.size = 1, label = T)
#p2 <- DimPlot(merged, reduction = "umap", group.by = "functional.cluster", pt.size = 1)
#p3 <- DimPlot(merged, reduction = "umap", group.by = "cancer.type", pt.size = 1)

#tiff('Tpex_EffectorMemory.tiff', res=600, width = 12000, height = 4000, units='px', compression = 'lzw')
#p1+p2+p3
#dev.off()

markers <- FindMarkers(merged2, group.by = "functional.cluster", ident.1 = "CD8_Tex", min.pct = 0.5)
markers <- subset(markers, subset = p_val <= 0.05)
head(markers); dim(markers)

# write.xlsx(markers, "markers_0.5.xlsx", sheetName = "min.pct 0.5", append = TRUE)

# Me quedo solo con los markers del volcano (> 0.5 y < -0.5)

markers2 <- subset(markers, subset = avg_log2FC >= 0.5)
markers3 <- subset(markers, subset = avg_log2FC <= -0.5)

markers.volcano <- rbind(markers2, markers3)
markers.list <- rownames(markers.volcano)

# Proteasome

markers.list <- c("Psma1","Psma2", "Psma3", "Psma4","Psma5","Psma6","Psma7","Psma8","Psmb1","Psmb2","Psmb3","Psmb4","Psmb5","Psmb6","Psmb7","Psmb8","Psmb9","Psmb10","Psmc1","Psmc2","Psmc3","Psmc4","Psmc5","Psmc6","Psmd1","Psmd2","Psmd3","Psmd4","Psmd5","Psmd6","Psmd7","Psmd8","Psmd9","Psmd10","Psmd11","Psmd12","Psmd13","Psmd14","Sem1","Adrm1","Psme1","Psme2","Psme3","Psme4","Psmf1")

###########################
######### Heatmap #########
###########################

# https://scillus.netlify.app/vignettes/plotting.html

pheatmap <- merged2[markers.list,] #] merged2 es contiene todos los datos, excepto de lung_1 (el que hay que eliminar)

pheatmap_Tex <- subset(pheatmap, subset = functional.cluster == "CD8_Tex")
pheatmap_Tem <- subset(pheatmap, subset = functional.cluster == "CD8_EffectorMemory")

pheatmap_Tem <- pheatmap_Tem[, sample(colnames(pheatmap_Tem), size = ncol(pheatmap_Tex), replace=F)]

pheatmap <- merge(pheatmap_Tem, pheatmap_Tex)

DefaultAssay(pheatmap) <- "RNA"

pheatmap <- NormalizeData(pheatmap, normalization.method = "LogNormalize", scale.factor = 10000)
pheatmap <- ScaleData(pheatmap, features = markers.list)

#tiff('Heatmap.tiff', res=700, width = 5000, height = 4000, units='px', compression = 'lzw')
plot_heatmap(dataset = pheatmap,
             row_font_size = 8,
             markers = markers.list,
             sort_var = c("functional.cluster","cancer.type"),
             anno_var = c("functional.cluster","cancer.type"),
             anno_colors = list(c("green3", "mediumblue"),
                                c("orange1","snow3", "mediumorchid3","springgreen2")),
             hm_limit = c(-1,0,1),
             hm_colors = c("green3", "white","mediumblue"))
#dev.off()

######################################################################################
######### GSEA using all genes from the T cells extracted from the 4 studies #########
######################################################################################

# Use all genes from the CD8A expressing T cells from the 4 studies included (withouth lung_1).

head(merged2)

markers.genes <- wilcoxauc(X = merged2, 'functional.cluster')
dplyr::count(markers.genes, group)

markers.genes <- markers.genes %>%
  dplyr::filter(pval <= 0.05) %>%
  dplyr::filter(group == "CD8_Tex") %>%
  arrange(desc(logFC)) %>%
  dplyr::select(feature, logFC)           # select only the features with p< 0.05 and order by logFC.

markers.genes$feature = toupper(markers.genes$feature)

ranks<- deframe(markers.genes)
head(ranks); dim(ranks)

barplot(sort(ranks, decreasing = T))

# Load the reference

# To do Gene set enrichment analysis, we need to have the annotated gene set first. One popular source is the MsigDB from Broad Institute.
# https://www.gsea-msigdb.org/gsea/msigdb/index.jsp

msigdbr_species()

m_df<- msigdbr(species = "Homo sapiens", category = "H")
head(m_df)

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# Run the analysis

fgseaRes<- fgsea(fgsea_sets, ranks, minSize=15, maxSize = 500, nperm= 1000)

head(fgseaRes[order(padj, -abs(NES)), ], n=10)

write.xlsx(fgseaRes[,c(1:7)], "GSEA_pathways.xlsx", sheetName = "H_4studies", append = TRUE)

# Plot the results

pathway <- "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"

tiff('HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY.tiff', res=700, width = 3000, height = 1300, units='px', compression = 'lzw')
plotEnrichment(fgsea_sets[[pathway]], ranks, ticksSize = 0.5, stats = T) +
  labs(title=pathway) +
  plotGseaTable(fgsea_sets[[pathway]], ranks, fgseaRes, gseaParam = 1) + 
  theme_bw()
dev.off()


tiff('Panel GSEA.tiff', res=700, width = 3000, height = 1300, units='px', compression = 'lzw')

dev.off()

####################################################
######### Heatmap using the gene signature #########
####################################################

gene.signature <- c("Cxcr4", "Rbpj", "Rps27a", "Psmb2", "Ldha", "Rbx1", "Psmb8", "Eno1", "Pkm", "Psmc3", "Chchd2", "Pgk1")

pheatmap2 <- merged2[gene.signature,] #] merged2 es contiene todos los datos, excepto de lung_1 (el que hay que eliminar)

pheatmap2_Tex <- subset(pheatmap2, subset = functional.cluster == "CD8_Tex")
pheatmap2_Tem <- subset(pheatmap2, subset = functional.cluster == "CD8_EffectorMemory")

pheatmap2_Tem <- pheatmap2_Tem[, sample(colnames(pheatmap2_Tem), size = ncol(pheatmap2_Tex), replace=F)]

pheatmap2 <- merge(pheatmap2_Tem, pheatmap2_Tex)

DefaultAssay(pheatmap2) <- "RNA"

pheatmap2 <- NormalizeData(pheatmap2, normalization.method = "LogNormalize", scale.factor = 10000)
pheatmap2 <- ScaleData(pheatmap2, features = gene.signature)

tiff('Heatmap_gene signature.tiff', res=700, width = 5000, height = 1500, units='px', compression = 'lzw')
plot_heatmap(dataset = pheatmap2,
             row_font_size = 8,
             markers = gene.signature,
             sort_var = c("functional.cluster","cancer.type"),
             anno_var = c("functional.cluster","cancer.type"),
             anno_colors = list(c("green3", "mediumblue"),
                                c("orange1","snow3", "mediumorchid3","springgreen2")),
             hm_limit = c(-1.5,0,1.5),
             hm_colors = c("green3", "white","mediumblue"))
dev.off()

violin <- VlnPlot(pheatmap, features = c("Sem1"), group.by = "functional.cluster")
violin$data
write.xlsx2(violin$data, "human proteasome gene expression-chaperones.xlsx", sheetName = "Psmf1", append= TRUE)
