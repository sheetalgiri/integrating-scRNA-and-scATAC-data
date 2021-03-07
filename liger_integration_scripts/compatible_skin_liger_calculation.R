library(Seurat)
library(SeuratWrappers)
library(Signac)
library(rliger)

##NOTE overriding function only to rename references to "liger" in the code to "rliger"##
RunOptimizeALS <- function(
  object,
  k,
  assay = NULL,
  split.by = 'orig.ident',
  lambda = 5,
  thresh = 1e-6,
  max.iters = 30,
  reduction.name = 'iNMF_raw',
  reduction.key = 'riNMF_',
  nrep = 1,
  H.init = NULL,
  W.init = NULL,
  V.init = NULL,
  rand.seed = 1,
  print.obj = FALSE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  if (Seurat:::IsMatrixEmpty(x = GetAssayData(object = object, slot = 'scale.data'))) {
    stop("Data is unscaled, splease scale before running", call. = FALSE)
  }
  if (is.character(x = split.by) && length(x = split.by) == 1) {
    split.by <- object[[split.by]]
  }
  split.cells <- split(x = colnames(x = object), f = split.by)
  scale.data <- lapply(
    X = split.cells,
    FUN = function(x) {
      return(t(x = GetAssayData(
        object = object,
        slot = 'scale.data',
        assay = assay
      )[, x]))
    }
  )
  # scale.data <- sapply(X = scale.data, FUN = t, simplify = FALSE)
  out <- rliger::optimizeALS(
    object = scale.data,
    k = k,
    lambda = lambda,
    thresh = thresh,
    max.iters = max.iters,
    nrep = nrep,
    H.init = H.init,
    W.init = W.init,
    V.init = V.init,
    rand.seed = rand.seed,
    print.obj = print.obj
  )
  colnames(x = out$W) <- VariableFeatures(object = object)
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = do.call(what = 'rbind', args = out$H),
    loadings = t(x = out$W),
    assay = assay,
    key = reduction.key
  )
  Tool(object = object) <- sapply(
    X = out$V,
    FUN = function(x) {
      colnames(x = x) <- VariableFeatures(object = object)
      rownames(x = x) <- colnames(x = object[[reduction.name]])
      return(t(x = x))
    },
    simplify = FALSE
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}

RunQuantileNorm <- function(
  object,
  split.by = 'orig.ident',
  reduction = 'iNMF_raw',
  reduction.name = 'iNMF',
  reduction.key = 'iNMF_',
  quantiles = 50,
  ref_dataset = NULL,
  min_cells = 20,
  knn_k = 20,
  dims.use = NULL,
  do.center = FALSE,
  max_sample = 1000,
  eps = 0.9,
  refine.knn = TRUE,
  ...
) {
  embeddings <- sapply(
    X = SplitObject(object = object, split.by = split.by),
    FUN = function(x) {
      return(Embeddings(object = x[[reduction]]))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (is.null(x = ref_dataset)) {
    num.samples <- vapply(
      X = embeddings,
      FUN = nrow,
      FUN.VALUE = integer(length = 1L)
    )
    ref_dataset <- names(x = embeddings)[which.max(x = num.samples)]
  } else if (is.numeric(x = ref_dataset)) {
    ref_dataset <- names(x = embeddings)[ref_dataset]
  }
  if (is.character(x = ref_dataset) && !ref_dataset %in% names(x = embeddings)) {
    stop("Cannot find reference dataset '", ref_dataset, "' in the split", call. = FALSE)
  }
  out <- rliger::quantile_norm(
    object = embeddings,
    quantiles = quantiles,
    ref_dataset = ref_dataset,
    min_cells = min_cells,
    knn_k = knn_k,
    dims.use = dims.use,
    do.center = do.center,
    max_sample = max_sample,
    eps = eps,
    refine.knn = refine.knn,
    ...
  )
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = out$H.norm,
    assay = DefaultAssay(object = object[[reduction]]),
    key = reduction.key
  )
  out <- as.data.frame(x = out[names(x = out) != 'H.norm'])
  object[[colnames(x = out)]] <- out
  Idents(object = object) <- 'clusters'
  object <- LogSeuratCommand(object = object)
  return(object)
}
############################################################

#load files for skin
file_activity_matrix_atac <- "/dss/dsshome1/lxc02/di82cox/data/skin_data/integration_data/skin_activityMatrix.rds"
activity.matrix <-readRDS(file_activity_matrix_atac)
preprocessed_rna <- readRDS("/dss/dsshome1/lxc02/di82cox/code/sheetal/skin_rna_preprocessed.rds")
raw_atac <- readRDS("/dss/dsshome1/lxc02/di82cox/code/sheetal/skin_atac.rds")


# atac preprocessing
atac <- CreateSeuratObject(counts = raw_atac[['ATAC']], assay = "ATAC", project = "ATAC")
atac[["RNA"]] <- CreateAssayObject(counts = activity.matrix)
atac@meta.data$paper.cell.type <- raw_atac@meta.data$paper.cell.type


# filter the atac data
atac <- subset(atac, subset = nFeature_ATAC > 200 & nFeature_ATAC < 20000)
atac <- subset(atac, subset = nCount_ATAC < 20000)
atac

# select used features
atac <- FindTopFeatures(atac, min.cutoff = 'q75')
atac 


#rna preprocessing
rna <- CreateSeuratObject(counts = preprocessed_rna[['RNA']], assay = "RNA", project = "RNA")
rna@meta.data$paper.cell.type <- preprocessed_rna@meta.data$paper.cell.type #recover cell type lost in conversion

#combine rna and atac, convert seurat to liger
combined <- merge(rna, y = atac, add.cell.ids = c("rna", "atac"), project = "scmo")

#use liger operations on seurat object using seurat-wrappers, tutorial: https://rdrr.io/github/satijalab/seurat-wrappers/f/docs/liger.Rmd 
combined <- NormalizeData(combined) #log normalize
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined, split.by = "orig.ident", do.center = FALSE)
combined <- RunOptimizeALS(combined, k = 20, lambda = 5, split.by = "orig.ident")
combined <- RunQuantileNorm(combined, split.by = 'orig.ident') #get clusters, stored in "clusters" in the objectt metadat

saveRDS(combined,"/dss/dsshome1/lxc02/di82cox/code/sheetal/compatible_skin_combined.rds")
combined <- readRDS("/dss/dsshome1/lxc02/di82cox/code/sheetal/compatible_skin_combined.rds")

#construct SNN graph
combined <- FindNeighbors(combined, reduction = 'iNMF', dims = 1:20)

#get louvain clusters, stored in "seurat_clusters" in the object metadata
combined <- FindClusters(combined, resolution = 0.7)

#Run UMAP and plot
combined <- RunUMAP(combined, dims = 1:ncol(combined[["iNMF"]]), reduction = "iNMF")
DimPlot(combined, split.by = c("orig.ident"),group.by = c("paper.cell.type") ,ncol = 1)
DimPlot(combined, group.by = c("paper.cell.type","seurat_clusters","clusters"), ncol = 3)


###########save anndata###
library(SeuratDisk)
SaveH5Seurat(combined,"compatible_skin",overwrite=TRUE)
adata <- Convert("compatible_skin.h5seurat","compatible_skin.h5ad",assay="RNA",overwrite=TRUE)
