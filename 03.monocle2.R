###Take mouse epithelial cells as an example. The chicken samples are processed in the same way.
library(Seurat)
library(cowplot)
library(dplyr)
library(monocle)
library(destiny)
library(slingshot)
library(SingleCellExperiment)
library(FateID)
library(ggplot2)
library(RColorBrewer)
source("/public/home/zhaolei/workshop/R_package/code_R/XY_code_integrate.R")
library(reticulate)
library(ReductionWrappers)
library(s2a)
library(trqwe)
EPT<-subset(all_mouse,cells=rownames(all_mouse@meta.data[all_mouse@meta.data$celltype %in% c("AT1","AT2","S","EP"),]))
EPT<-FindVariableFeatures(EPT,mean.function=ExpMean,dispersion.function=LogVMR)
EPT<-ScaleData(EPT,vars.to.regress = c('nCount_RNA'),model.use = 'linear', use.umi = FALSE)
EPT<- RunPCA(EPT, pc.genes = VariableFeatures(EPT),npc = 50)
ElbowPlot(EPT)
dev.off()
EPT<-FindNeighbors(EPT,dims = 1:19)
EPT<-FindClusters(EPT,resolution = 1) 
EPT<-RunUMAP(EPT,dims = 1:19)
EPT<-RunTSNE(EPT,dims = 1:19)
pdf("msEPT1.pdf")
DimPlot(EPT1,reduction = "umap",label = T)+theme_test()
DimPlot(EPT1,reduction = "umap",label = T,group.by="celltype",split.by = "stage",ncol=4)+NoLegend()+theme_test()
FeaturePlot(EPT1,features=c("Sftpb","Sftpc","Ascl1","Trp63","Krt5","Scgb3a2"))
FeaturePlot(EPT1,features=c("Foxj1","Sox2","Sox9","Id2","Ager","Clic5","Hopx"))
dev.off()

EPT1<-subset(EPT,idents=c(14,15,19,3,8,17,13,9,2,0,4,22,10,6,1,5))
EPT1<-subset(EPT1,cells=setdiff(rownames(EPT1@meta.data),cells))
sim <- as.SingleCellExperiment(EPT1)
colData(sim)$order <- colData(sim)$seurat_clusters
table(colData(sim)$seurat_clusters)
library(RColorBrewer)
colors <- unique(union(brewer.pal(8,'Set3'),brewer.pal(8,'Accent')))
colors <- unique(union(colors,brewer.pal(8,'Set1')))[c(1:length(unique(colData(sim)$order)))]
names(colors) <- unique(colData(sim)$order)
sce <- slingshot(sim, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP',start.clus = c("14"),end.clus =c("1","9","19"))

p3=DimPlot(object = EPT1, reduction = "umap",label=TRUE,group.by="stage") +NoLegend()+labs(title="umap")

pdf("msslingshotEP.pdf")
plot(reducedDims(sce)$UMAP[,c(1,2)], col = colors[colData(sim)$order])
lines(SlingshotDataSet(sce), lwd=2, col='black',show.constraints = TRUE)
p3
dev.off()


all_meta <- EPT1[[]][,1:9]
pseudo <- colData(sce)
pseudo <- pseudo[rownames(all_meta),]
all_meta <- cbind(all_meta,pseudo[,c("slingPseudotime_1")])
colnames(all_meta)[c(ncol(all_meta))] <- c("dis_pro_Pseudotime_1")
EPT1@meta.data <- as.data.frame(all_meta)
all_meta <- cbind(all_meta,pseudo[,c("slingPseudotime_3")])
colnames(all_meta)[c(ncol(all_meta))] <- c("dis_pro_Pseudotime_3")
EPT1@meta.data <- as.data.frame(all_meta)
all_meta <- cbind(all_meta,pseudo[,c("slingPseudotime_2")])
colnames(all_meta)[c(ncol(all_meta))] <- c("dis_pro_Pseudotime_2")
EPT1@meta.data <- as.data.frame(all_meta)

p1 = Pesudo_FeaturePlot(object = EPT1, features = c("dis_pro_Pseudotime_1"),  ncol=1,pt.size=.5,reduction="umap",label=T,cols = CustomPalette(low ="#007BBF", mid = "#FFF485",high = "#FF0000")) +NoAxes()
p2 = Pesudo_FeaturePlot(object = EPT1, features = c("dis_pro_Pseudotime_2"),  ncol=1,pt.size=.5,reduction="umap",label=T,cols = CustomPalette(low ="#007BBF", mid = "#FFF485",high = "#FF0000")) +NoAxes()
p3 = Pesudo_FeaturePlot(object = EPT1, features = c("dis_pro_Pseudotime_3"),  ncol=1,pt.size=.5,reduction="umap",label=T,cols = CustomPalette(low ="#007BBF", mid = "#FFF485",high = "#FF0000")) +NoAxes()

pdf("msEPslingshot2.pdf")
plot_grid(p1,p2,p3,nrow=2)
dev.off()

#############slingPseudotime_1

EPseuratX_all <- EPT1
EPmeta_tmp <- EPseuratX_all@meta.data
EPmeta_tmp$new_idents <- EPmeta_tmp$dis_pro_Pseudotime_1
EPmeta_tmp$new_idents <- ifelse(EPmeta_tmp$new_idents %in% NA,"score_no","score_yes")
EPseuratX_all@meta.data <- EPmeta_tmp
Idents(EPseuratX_all) <- EPseuratX_all$new_idents
EPPseudotime_1 <- subset(EPseuratX_all,idents=c("score_yes"))

EPmeta_tmp <- EPPseudotime_1@meta.data
EPmeta_tmp <- EPmeta_tmp[order(EPmeta_tmp$dis_pro_Pseudotime_1),]
EPmeta_tmp$order <- 1:nrow(EPmeta_tmp)
EPPseudotime_1$Pseudotime_all <- EPmeta_tmp[rownames(EPPseudotime_1[[]]),]$order
EPdata <- as(as.matrix(EPPseudotime_1@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = EPPseudotime_1@meta.data)
fData <- data.frame(gene_short_name = row.names(EPdata), row.names = row.names(EPdata))
fd <- new('AnnotatedDataFrame', data = fData)
EPmonocle_obj1_all <- newCellDataSet(EPdata,
                              phenoData = pd,
                              featureData = fd,
                              expressionFamily = uninormal())
                            
pData(EPmonocle_obj1_all)$Pseudotime <- pData(EPmonocle_obj1_all)$dis_pro_Pseudotime_1


EP1diff_test_res_all <- monocle::differentialGeneTest(EPmonocle_obj1_all,
              fullModelFormulaStr = "~sm.ns(dis_pro_Pseudotime_1)")


EPsig_gene_names1 <- row.names(subset(EPdiff_test_res_all1, qval < 1e-6))
bin <- Binner(cds_object=EPmonocle_obj1_all,anno_group="stage")

EPt1 <- plot_pseudotime_heatmap(EPmonocle_obj1_all[EPsig_gene_names1[1:4000]],
                num_clusters = 6,
                show_rownames = F,norm_method="vstExprs",
                add_annotation_col = bin, return_heatmap=TRUE)
pdf("msEP1pseudotime_heatmap.pdf")
EPt1
dev.off()

EPannotation_row <- data.frame(Cluster = factor(cutree(EPt1$tree_row,6)))
EPannotation_row$gene <- rownames(EPannotation_row)

low<-EPannotation_row[EPannotation_row$Cluster %in% c(1,4,6),]$gene
high<-EPannotation_row[EPannotation_row$Cluster %in% c(3,2,5),]$gene

EPt1 <- plot_pseudotime_heatmap(EPmonocle_obj1_all[EPsig_gene_names1[4001:8211]],
                num_clusters = 6,
                show_rownames = F,norm_method="vstExprs",
                cores = 1,
                add_annotation_col = bin, return_heatmap=TRUE)
pdf("msEP1pseudotime_heatmap.pdf")
EPt1
dev.off()
