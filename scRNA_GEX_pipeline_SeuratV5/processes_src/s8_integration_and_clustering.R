s8.integration.and.clustering <- function(s.obj, 
                           path.to.output, 
                           save.RDS.s8,
                           PROJECT, 
                           num.dim.integration,
                           num.PCA,
                           num.PC.used.in.UMAP,
                           num.PC.used.in.Clustering,
                           cluster.resolution = 0.5,
                           my_random_seed = 42,
                           umap.method = "uwot",
                           genes.to.not.run.PCA = NULL,
                           inte_pca_reduction_name = "INTE_PCA", 
                           inte_umap_reduction_name = "INTE_UMAP",
                           with.TSNE = FALSE,
                           k.filter = 200){
  
  if (is.na(k.filter) == TRUE){
    k.filter.default <- 200
  } else {
    k.filter.default <- k.filter
  }
  
  data.list <- SplitObject(s.obj, split.by = "name")
  data.list <- lapply(X = data.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})    
  
  anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:num.dim.integration, scale=F,
                                    k.filter = k.filter)## THIS IS CCA DIMENSIONS
  
  s.obj_inte <- IntegrateData(anchorset = anchors, dims = 1:num.dim.integration, k.weight = k.filter.default) ## THIS IS PCA DIMENSION
  
  ## keep the order of integration obj
  s.obj_inte <- s.obj_inte[, colnames(s.obj)]
  
  s.obj[['integrated']] <- s.obj_inte[['integrated']]
  
  s.obj@commands <- c(s.obj@commands, s.obj_inte@commands)
  
  s.obj@tools <- c(s.obj@tools, s.obj_inte@tools)
  
  DefaultAssay(s.obj) <- "integrated"
  
  s.obj <- ScaleData(s.obj, verbose = FALSE, features = row.names(s.obj))
  
  if (is.null(genes.to.not.run.PCA) == TRUE){
    s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=inte_pca_reduction_name)
  } else {
    s.obj <- FindVariableFeatures(s.obj)
    pca.genes <- setdiff(VariableFeatures(s.obj), genes.to.not.run.PCA)
    print("####################################################################")
    print(sprintf("Running PCA with %s genes, after removing %s genes", length(pca.genes), length(genes.to.not.run.PCA)))
    print("####################################################################")
    s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=inte_pca_reduction_name, 
                    features = pca.genes)
  }
  
  s.obj <- RunUMAP(s.obj, reduction = inte_pca_reduction_name, dims = 1:num.PC.used.in.UMAP, reduction.name=inte_umap_reduction_name, 
                   seed.use = my_random_seed, umap.method = umap.method)
  
  if (with.TSNE == TRUE){
    s.obj <- RunTSNE(s.obj, reduction = inte_pca_reduction_name, dims = 1:num.PC.used.in.UMAP, reduction.name="INTE_TSNE", 
                     seed.use = my_random_seed)
  }
  
  # clustering 
  s.obj <- FindNeighbors(s.obj, reduction = inte_pca_reduction_name, dims = 1:num.PC.used.in.Clustering)
  
  s.obj <- FindClusters(s.obj, resolution = cluster.resolution, random.seed = 0)
  
  if (save.RDS.s8 == TRUE){
    dir.create(file.path(path.to.output, "s8_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s8_output", 
                             paste0(PROJECT, ".output.s8.rds")))
  }
  
  return(s.obj)
}