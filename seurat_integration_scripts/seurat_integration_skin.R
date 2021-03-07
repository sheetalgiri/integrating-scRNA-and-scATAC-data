### karoline lutz - 2021.02.10
### SKIN DATA
### preprocessing of atac
### rna already is preprocessed
### integrating both datasets
### coembeding for visualization

### at the bottom is a slightly changed version of the CreateGeneActivityMatrix
### using the modified version or a "file too large" error will be triggered


library(Seurat)
library(ggplot2)
library(patchwork)




# read in raw data
file_atac <- "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\GSM4156597_skin.late.anagen_raw_cell_annot.atac.h5ad"
gtf_file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\gencode.vM25.annotation.gtf"
# recreated the h5ad file using anndata 0.6.22 in conda env old_scanpy
atac <- ReadH5AD(file_atac, assay = "ATAC", layers = "data")


saveRDS(atac, file = 'C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\GSM4156597_skin.late.anagen_raw_cell_annot.atac.rds')
atac <- readRDS('C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\GSM4156597_skin.late.anagen_raw_cell_annot.atac.rds')
atac

#create the raw geneactivity matrix
chromosome <- paste0('chr', c(1:19, "X", "Y"))
chromosome
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = atac[['ATAC']]@counts, annotation.file = gtf_file, 
                                            seq.levels = chromosome, upstream = 2000, verbose = TRUE)


saveRDS(activity.matrix, file = 'C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\skin_activityMatrix.rds')

# start building the seurat object, add the geneactivity and metadata
skin.atac <- CreateSeuratObject(counts = atac[['ATAC']], assay = "ATAC", project = "SHAREseq_ATAC")


skin.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
skin.atac@meta.data$paper.cell.type <- atac@meta.data$paper.cell.type




skin.atac
saveRDS(skin.atac, file = 'C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\skin_atacPlusActivity.rds')
skin.atac = readRDS('C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\skin_atacPlusActivity.rds')

skin.atac
#An object of class Seurat 
#367091 features across 34774 samples within 2 assays 
#Active assay: ATAC (344592 features, 0 variable features)
#1 other assay present: ACTIVITY



# filter the atac data
skin.atac$nFeature_ATAC
VlnPlot(skin.atac, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2)

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\skin\\skin_featureplot.pdf") 
VlnPlot(skin.atac, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2)
dev.off()

skin.atac <- subset(skin.atac, subset = nFeature_ATAC > 200 & nFeature_ATAC < 20000)
skin.atac <- subset(skin.atac, subset = nCount_ATAC < 20000)
skin.atac
#An object of class Seurat 
#367091 features across 34471 samples within 2 assays 
#Active assay: ATAC (344592 features, 0 variable features)
#1 other assay present: ACTIVITY
VlnPlot(skin.atac, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2)

# select used features
skin.atac <- FindTopFeatures(skin.atac, min.cutoff = 'q75')
skin.atac 
#An object of class Seurat 
#367091 features across 34471 samples within 2 assays 
#Active assay: ATAC (344592 features, 86204 variable features)
#1 other assay present: ACTIVITY
VariableFeatures(skin.atac)




# filter geneactivity
DefaultAssay(skin.atac) <- "ACTIVITY"
skin.atac #Active assay: ACTIVITY (22499 features, 0 variable features)
VlnPlot(skin.atac, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2)

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\skin\\skin_featureplot_activity.pdf") 
VlnPlot(skin.atac, features = c("nFeature_ATAC", "nCount_ATAC"), ncol = 2)
dev.off()

skin.atac <- FindVariableFeatures(skin.atac)
skin.atac <- NormalizeData(skin.atac)

saveRDS(skin.atac, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\skin_logNorm_atac.rds")
skin.atac
#An object of class Seurat 
#367091 features across 34471 samples within 2 assays 
#Active assay: ACTIVITY (22499 features, 2000 variable features)
#1 other assay present: ATAC


skin.atac = readRDS("C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\skin_logNorm_atac.rds")
skin.atac



###Continue with Seurat pipeline

skin.atac <- ScaleData(skin.atac)

skin.atac <- RunPCA(skin.atac, npcs = 50)
skin.atac <- RunUMAP(skin.atac, reduction = "pca", dims = 1:50)
p2 <- DimPlot(skin.atac, reduction = "umap", group.by = "paper.cell.type", label = TRUE) + NoLegend() + ggtitle("Skin GeneActivity")
p2

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\skin\\skin_activity_umap.pdf") 
p2
dev.off()


# run lsi and umap on atac data
DefaultAssay(skin.atac) <- "ATAC"
skin.atac
VariableFeatures(skin.atac) <- names(which(Matrix::rowSums(skin.atac) > 10))
skin.atac <- RunLSI(skin.atac, n = 50, scale.max = NULL)
skin.atac <- RunUMAP(skin.atac, reduction = "lsi", dims = 1:50)

#pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#top10 <- head(VariableFeatures(pbmc), 10)

p1 <- DimPlot(skin.atac, reduction = "umap", group.by = "paper.cell.type", label = TRUE) + NoLegend() + ggtitle("Skin ATAC")
p1

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\skin\\skin_atac_umap.pdf") 
p1
dev.off()







# read in rna preprocessed data
#file_rna <- "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\skin_logNorm_rna.h5ad"
file_rna <- "C:\\Users\\Karol\\Desktop\\singlecell\\server_data\\final_skin_logNorm_rna.h5ad"
rna <- ReadH5AD(file_rna, assay = "RNA", layers = "data")
rna 
#An object of class Seurat 
#21834 features across 34755 samples within 1 assay 
#Active assay: RNA (21834 features, 0 variable features)


#rna <- FindVariableFeatures(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna
#An object of class Seurat 
#21834 features across 34755 samples within 1 assay 
#Active assay: RNA (21834 features, 2000 variable features)



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
p5 <- DimPlot(rna, reduction = "umap", group.by = "paper.cell.type", label = TRUE)
p5

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\skin\\skin_rna_umap.pdf") 
p5
dev.off()






# compare umaps
p6 <- DimPlot(skin.atac, reduction = "umap", group.by = "paper.cell.type", label = TRUE) + NoLegend() + ggtitle("scATAC-seq-log (skin)")
p7 <- DimPlot(rna, group.by = "paper.cell.type", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq (skin)")
p6 + p7




# do the seurat prediction
transfer.anchors <- FindTransferAnchors(reference = rna, query = skin.atac, features = VariableFeatures(object = rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca", k.filter = NA)



#save anchors because it takes ages to calculate
saveRDS(transfer.anchors, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\transfer_anchors.rds")
#save atac and rna inbetween as well, just to make sure in case rstudio randomly closes again
saveRDS(skin.atac, file= "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\atac_tmp.rds")
saveRDS(rna, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\rna_tmp.rds")



transfer.anchors = readRDS("C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\transfer_anchors.rds")
skin.atac = readRDS("C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\atac_tmp.rds")
rna = readRDS("C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\rna_tmp.rds")




celltype_ref <- as.character(rna@meta.data$paper.cell.type)
#celltype_ref

celltype.predictions <-  TransferData(anchorset = transfer.anchors, refdata = celltype_ref, 
                                      weight.reduction = skin.atac[["lsi"]])

#look at prediction scores and maybe filter out low ones
skin.atac <- AddMetaData(skin.atac, metadata = celltype.predictions)
hist(skin.atac$prediction.score.max, main = "Prediction Scores Skin", xlab = "Scores")
abline(v = 0.5, col = "red")

pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\skin\\skin_predictionscores.pdf", width = 10) 
hist(skin.atac$prediction.score.max, main = "Prediction Scores Skin", xlab = "Scores")
abline(v = 0.5, col = "red")
dev.off()


table(skin.atac$prediction.score.max > 0.5)
#FALSE  TRUE 
#13221 21250 


#maybe filter? but idk if it will be a good idea since our scores are so low anyways
skin.atac.filtered <- subset(skin.atac, subset = prediction.score.max > 0.35)
skin.atac.filtered$predicted.id <- factor(skin.atac.filtered$predicted.id, levels = levels(rna))


# compare atac and rna again
p8 <- DimPlot(skin.atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells (skin)") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p9 <- DimPlot(rna, group.by = "paper.cell.type", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells (skin)") + 
  NoLegend()
p8 + p9

#co-embed
genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

#does not use the filtered atac data
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = skin.atac[["lsi"]])


# add the imputation to our seurat object
skin.atac[["RNA"]] <- imputation
rna$tech <- "rna"
skin.atac$tech <- "atac"
coembed <- merge(x = rna, y = skin.atac)

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed@meta.data$paper.cell.type <- ifelse(!is.na(coembed$paper.cell.type), coembed$paper.cell.type, coembed$predicted.id)


p10 <- DimPlot(coembed, group.by = "tech") + ggtitle("ATAC RNA Integration (Skin)")
p11 <- DimPlot(coembed, group.by = "paper.cell.type", label = TRUE, repel = TRUE) + ggtitle("Celltypes (Skin)") + NoLegend()
p10 + p11


pdf("C:\\Users\\Karol\\Desktop\\singlecell\\figures\\skin\\skin_integration.pdf", width = 10) 
p10+p11
dev.off()


#save coembed as rds
saveRDS(coembed, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\coembed_skin.rds")

coembed = readRDS("C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\coembed_skin.rds")



#memory.size()
#memory.size(TRUE)
#memory.limit(9999999999)

#> memory.limit()
#[1] 8103
## To increase the storage capacity
#> memory.limit(size=56000)
#[1] 56000  



library(SeuratDisk)
#save as h5seurat
SaveH5Seurat(coembed, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\skin_coembed.h5Seurat")

saveRDS(skin.atac, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\skin_final_atac.h5Seurat")
saveRDS(rna, file = "C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\skin_final_rna.h5Seurat")


#convert coembed as h5ad
Convert("C:\\Users\\Karol\\Desktop\\singlecell\\data\\skin\\preprocessed\\skin_coembed.h5Seurat", dest = "h5ad")








### used for creating the geneactivitymatrix while otherwise file too large error

library(utilities)
library(pbapply)
### NEW GENEACTIVITY FUNCTION
CreateGeneActivityMatrix <- function (peak.matrix, annotation.file, seq.levels = c(1:22, 
                                                                                   "X", "Y"), include.body = TRUE, upstream = 2000, downstream = 0, 
                                      verbose = TRUE) 
{
  #if (!PackageCheck("GenomicRanges", error = FALSE)) {
  #  stop("Please install GenomicRanges from Bioconductor.")
  #}
  #if (!PackageCheck("rtracklayer", error = FALSE)) {
  #  stop("Please install rtracklayer from Bioconductor.")
  #}
  peak.df <- rownames(x = peak.matrix)
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, 
                                                            pattern = ":", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", "start", "end")
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 
                                 0, ]) <- 1
  gtf <- rtracklayer::import(con = annotation.file)
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, 
                                     pruning.mode = "coarse")
  if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
    GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
  }
  gtf.genes <- gtf[gtf$type == "gene"]
  if (include.body) {
    gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, 
                            downstream = downstream)
  }
  else {
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, 
                                                     upstream = upstream, downstream = downstream)
  }
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, 
                                                     subject = gtf.body_prom)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 
                                    0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
  gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]
  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
  annotations <- peak.ids[, c("peak", "gene.name")]
  colnames(x = annotations) <- c("feature", "new_feature")
  # 
  #peak.matrix <- as(object = peak.matrix, Class = "matrix")
  all.features <- unique(x = annotations$new_feature)
  # I'm guessing the parallel execution requires a full matrix
  # if (nbrOfWorkers() > 1) {
  #   mysapply <- future_sapply
  # }
  # else {
  mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
  # }
  newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x) {
    features.use <- annotations[annotations$new_feature == 
                                  all.features[[x]], ]$feature
    submat <- peak.matrix[features.use, ]
    if (length(x = features.use) > 1) {
      return(Matrix::colSums(x = submat))
    }
    else {
      return(submat)
    }
  })
  newmat <- t(x = newmat)
  rownames(x = newmat) <- all.features
  colnames(x = newmat) <- colnames(x = peak.matrix)
  return(as(object = newmat, Class = "dgCMatrix"))
}
