####Take lungfish as an example, the rest of the samples are operated in the same way as this one
gene<-read.table("/data02/liye/project/analys_sr/29.revised_2/01.lung_Adult_gene/lungfish_HVG_gene.tsv",header=T)
gene_cell_exp <- AverageExpression(data2_Noother, group.by = 'test1', slot = 'data')
####################lungfishEP
gene_cell_exp_scale<-as.data.frame(t(apply(gene_cell_exp$RNA,1,function(x){(x-mean(x[c(1:26,28:39)]))/sd(x[c(1:26,28:39)])})))
EP.gene1<-rownames(gene_cell_exp_scale[gene_cell_exp_scale$lung_EP >=1.645,])
gene_cell_exp_scale1<-as.data.frame(t(apply(gene_cell_exp$RNA[EP.gene1,grep("EP",colnames(gene_cell_exp$RNA))],1,function(x){(x-mean(x[c(1:5,7:8)]))/sd(x[c(1:5,7:8)])})))
specific_EPgene<-rownames(gene_cell_exp_scale1[gene_cell_exp_scale1$lung_EP >=1.645,])
write(unique(gene[gene$Gene %in% specific_EPgene,]$geneID),file="/data02/liye/project/analys_sr/29.revised_2/03.lungSpecificGene/lungfish_EP_specific_gene.tsv",sep="\n")

#######################lungfishEN
gene_cell_exp_scale<-as.data.frame(t(apply(gene_cell_exp$RNA,1,function(x){(x-mean(x[c(1:25,27:39)]))/sd(x[c(1:25,27:39)])})))
EN.gene1<-rownames(gene_cell_exp_scale[gene_cell_exp_scale$lung_EN >=1.645,])
gene_cell_exp_scale1<-as.data.frame(t(apply(gene_cell_exp$RNA[EN.gene1,grep("EN",colnames(gene_cell_exp$RNA))],1,function(x){(x-mean(x[c(1:6,8:9)]))/sd(x[c(1:6,8:9)])})))
specific_ENgene<-rownames(gene_cell_exp_scale1[gene_cell_exp_scale1$lung_EN >=1.645,])
write(unique(gene[gene$Gene %in% specific_ENgene,]$geneID),file="/data02/liye/project/analys_sr/29.revised_2/03.lungSpecificGene/lungfish_EN_specific_gene.tsv",sep="\n")

#######################lungfishIM

gene_cell_exp_scale<-as.data.frame(t(apply(gene_cell_exp$RNA,1,function(x){(x-mean(x[c(1:27,29:39)]))/sd(x[c(1:27,29:39)])})))
IM.gene1<-rownames(gene_cell_exp_scale[gene_cell_exp_scale$lung_Immune >=1.645,])
gene_cell_exp_scale1<-as.data.frame(t(apply(gene_cell_exp$RNA[IM.gene1,grep("Immune",colnames(gene_cell_exp$RNA))],1,function(x){(x-mean(x[c(1:7,9:11)]))/sd(x[c(1:7,9:11)])})))
specific_IMgene<-rownames(gene_cell_exp_scale1[gene_cell_exp_scale1$lung_Immune >=1.645,])
write(unique(gene[gene$Gene %in% specific_IMgene,]$geneID),file="/data02/liye/project/analys_sr/29.revised_2/03.lungSpecificGene/lungfish_IM_specific_gene.tsv",sep="\n")

########################lungfishStromal

gene_cell_exp_scale<-as.data.frame(t(apply(gene_cell_exp$RNA,1,function(x){(x-mean(x[c(1:28,30:39)]))/sd(x[c(1:28,30:39)])})))
STR.gene1<-rownames(gene_cell_exp_scale[gene_cell_exp_scale$lung_Stromal >=1.645,])
gene_cell_exp_scale1<-as.data.frame(t(apply(gene_cell_exp$RNA[STR.gene1,grep("Stromal",colnames(gene_cell_exp$RNA))],1,function(x){(x-mean(x[c(1:7,9:11)]))/sd(x[c(1:7,9:11)])})))
specific_STRgene<-rownames(gene_cell_exp_scale1[gene_cell_exp_scale1$lung_Stromal >=1.645,])
write(unique(gene[gene$Gene %in% specific_STRgene,]$geneID),file="/data02/liye/project/analys_sr/29.revised_2/03.lungSpecificGene/lungfish_STR_specific_gene.tsv",sep="\n")

