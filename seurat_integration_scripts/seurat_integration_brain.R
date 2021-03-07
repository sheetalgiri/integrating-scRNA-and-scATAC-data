### karoline lutz - 2021.02.10
### preprocessing of atac
### rna already is preprocessed
### integrating both datasets
### coembeding for visualization


# this script follows the Seurat Integration tutorial: https://satijalab.org/seurat/archive/v3.2/atacseq_integration_vignette.html

library(Seurat)
library(ggplot2)
library(patchwork)




# read in raw data
file_atac <- "C:\\Users\\Karol\\Desktop\\singlecell\\data\\brain\\GSM4156599_brain_raw_cell_annot.atac.h5ad"
gtf_file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\gencode.vM25.annotation.gtf"
# recreated the h5ad file using anndata 0.6.22 in conda env old_scanpy
atac <- ReadH5AD(file_atac, assay = "ATAC", layers = "data")
atac 


###GeneActivityMatrix and Preprocessing of atac and matrix data

#create the raw geneactivity matrix
chromosome <- paste0('chr', c(1:19, "X", "Y"))
chromosome
activity.matrix <- Seurat::CreateGeneActivityMatrix(peak.matrix = atac[['ATAC']]@counts, annotation.file = gtf_file, 
                                            seq.levels = chromosome, upstream = 2000, verbose = TRUE)

# start building the seurat object, add the geneactivity and metadata
brain.atac <- CreateSeuratObject(counts = atac[['ATAC']], assay = "ATAC", project = "SHAREseq_ATAC")
brain.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
brain.atac@meta.data$paper.cell.type <- atac@meta.data$paper.cell.type


# filter the atac data
brain.atac$nFeature_ATAC
VlnPlot(brain.atac, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2)

# Open a pdf file
pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\brain\\brain_featureplot.pdf") 
# 2. Create a plot
VlnPlot(brain.atac, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2)
# Close the pdf file
dev.off() 


brain.atac <- subset(brain.atac, subset = nFeature_ATAC > 200 & nFeature_ATAC < 6000)
brain.atac <- subset(brain.atac, subset = nCount_ATAC < 6000)
brain.atac

# select used features
brain.atac <- FindTopFeatures(brain.atac, min.cutoff = 'q75')
VariableFeatures(brain.atac)

# filter geneactivity
DefaultAssay(brain.atac) <- "ACTIVITY"
brain.atac
VlnPlot(brain.atac, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2)
pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\brain\\brain_featureplot_activity.pdf") 
VlnPlot(brain.atac, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2)
dev.off() 


brain.atac <- FindVariableFeatures(brain.atac)
#saveRDS(brain.atac, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\preprocessed\\brain_preprocess_atac.rds")

brain.atac <- NormalizeData(brain.atac)

###Preprocessing done

brain.atac
#An object of class Seurat 
#451647 features across 3248 samples within 2 assays 
#Active assay: ACTIVITY (23606 features, 2000 variable features)
#1 other assay present: ATAC
saveRDS(brain.atac, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\preprocessed\\brain_logNorm_atac.rds")

#brain.atac = readRDS("C:\\Users\\Karol\\Desktop\\singlecell\\data\\preprocessed\\brain_logNorm_atac.rds")
#brain.atac






###continue with seurat workflow
DefaultAssay(brain.atac) <- "ACTIVITY"
brain.atac <- ScaleData(brain.atac)

brain.atac <- RunPCA(brain.atac, npcs = 50)
brain.atac <- RunUMAP(brain.atac, reduction = "pca", dims = 1:50)
p2 <- DimPlot(brain.atac, reduction = "umap", group.by = "paper.cell.type", label = TRUE) + NoLegend() + ggtitle("activitymatrix")
p2

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\brain\\brain_activity_umap.pdf") 
p2
dev.off() 


# run lsi and umap on atac data
DefaultAssay(brain.atac) <- "ATAC"
brain.atac
VariableFeatures(brain.atac) <- names(which(Matrix::rowSums(brain.atac) > 10))
brain.atac <- RunLSI(brain.atac, n = 50, scale.max = NULL)
brain.atac <- RunUMAP(brain.atac, reduction = "lsi", dims = 1:50)

#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#top10 <- head(VariableFeatures(pbmc), 10)

p1 <- DimPlot(brain.atac, reduction = "umap", group.by = "paper.cell.type", label = TRUE) + NoLegend() + ggtitle("atac")
p1
p1+p2

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\brain\\brain_atac_umap.pdf") 
p1
dev.off() 


# read in rna preprocessed data
file_rna <- "C:\\Users\\Karol\\Desktop\\singlecell\\data\\preprocessed\\brain_logNorm_rna.h5ad"
#file_rna <- "C:\\Users\\Karol\\Desktop\\singlecell\\data\\brain\\brain_rna_withemptybarcode.h5ad"
rna <- ReadH5AD(file_rna, assay = "RNA", layers = "data")
rna 
#An object of class Seurat 
#18931 features across 3260 samples within 1 assay 
#Active assay: RNA (18931 features, 0 variable features)


rna <- FindVariableFeatures(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(rna), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(rna)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)


rna <- RunPCA(rna, features = VariableFeatures(object = rna))

print(rna[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(rna, dims = 1:2, reduction = "pca")

p3 <- DimPlot(rna, reduction = "pca", group.by = "paper.cell.type", label = TRUE)
p3

p4 <- DimHeatmap(rna, dims = 1, cells = 500, balanced = TRUE)
p4

rna <- FindNeighbors(rna, dims = 1:10)
rna <- FindClusters(rna, resolution = 0.5)

rna <- RunUMAP(rna, dims = 1:10)
p5 <- DimPlot(rna, reduction = "umap", group.by = "paper.cell.type", label = TRUE) + NoLegend()
p5

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\brain\\brain_rna_umap.pdf") 
p5
dev.off() 





# compare umaps
p6 <- DimPlot(brain.atac, reduction = "umap", group.by = "paper.cell.type", label = TRUE) + NoLegend() + ggtitle("scATAC-seq-log")
p7 <- DimPlot(rna, group.by = "paper.cell.type", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
p6 + p7




# do the seurat prediction
transfer.anchors <- FindTransferAnchors(reference = rna, query = brain.atac, features = VariableFeatures(object = rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca", k.filter = NA)

celltype_ref <- as.character(rna@meta.data$paper.cell.type)
#celltype_ref

celltype.predictions <-  TransferData(anchorset = transfer.anchors, refdata = celltype_ref, 
                                      weight.reduction = brain.atac[["lsi"]])


#look at prediction scores and maybe filter out low ones
brain.atac <- AddMetaData(brain.atac, metadata = celltype.predictions)
hist(brain.atac$prediction.score.max, main = "Prediction Scores Brain Integration")
abline(v = 0.5, col = "red")

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\brain\\brain_prediction_scores.pdf") 
hist(brain.atac$prediction.score.max, main = "Prediction Scores Brain Integration", xlab ="Scores")
abline(v = 0.5, col = "red")
dev.off() 

table(brain.atac$prediction.score.max > 0.5)
#FALSE  TRUE 
#1963  1285

#maybe filter? but idk if it will be a good idea since our scores are so low anyways
brain.atac.filtered <- subset(brain.atac, subset = prediction.score.max > 0.3)
brain.atac.filtered$predicted.id <- factor(brain.atac.filtered$predicted.id, levels = levels(rna))

# compare atac and rna again
p8 <- DimPlot(brain.atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p9 <- DimPlot(rna, group.by = "paper.cell.type", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
p8 + p9

#co-embed
genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

#does not use the filtered atac data
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = brain.atac[["lsi"]])

# add the imputation to our seurat object
brain.atac[["RNA"]] <- imputation
rna$tech <- "rna"
brain.atac$tech <- "atac"
coembed <- merge(x = rna, y = brain.atac)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed@meta.data$paper.cell.type <- ifelse(!is.na(coembed$paper.cell.type), coembed$paper.cell.type, coembed$predicted.id)


p10 <- DimPlot(coembed, group.by = "tech") + ggtitle("ATAC RNA Integration (Brain)")
p11 <- DimPlot(coembed, group.by = "paper.cell.type", label = TRUE, repel = TRUE) + ggtitle("Celltypes") + NoLegend()
p10 + p11

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\brain\\brain_integration_umap.pdf", width = 10) 
p10+p11
dev.off() 



saveRDS(coembed, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\brain\\brain_final_coembed.rds")




#write out the coembed data into h5Seurat
library(SeuratDisk)
SaveH5Seurat(coembed, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\brain\\brain_coembed_nobarcode.h5Seurat")
SaveH5Seurat(brain.atac, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\preprocessed\\brain_final_atac.h5Seurat")
SaveH5Seurat(rna, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\preprocessed\\brain_final_rna.h5Seurat")


library(SeuratDisk)
#save as h5seurat
SaveH5Seurat(coembed, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\skin_coembed.h5Seurat")

saveRDS(skin.atac, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\skin_final_atac.h5Seurat")
saveRDS(rna, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\skin_final_rna.h5Seurat")


#convert coembed as h5ad
Convert("C:\\Users\\Karol\\Desktop\\singlecell\\data\\brain\\brain_coembed_nobarcode.h5Seurat", dest = "h5ad")



