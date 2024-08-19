s8.integration.and.clustering_V5 <- function(s.obj, 
                                             save.RDS.s8,
                                             path.to.output,
                                             use.sctransform,
                                             num.PCA,
                                             num.PC.used.in.UMAP,
                                             num.PC.used.in.Clustering,
                                             cluster.resolution,
                                             vars.to.regress = c("percent.mt"),
                                             remove.genes = NULL){
  s.obj[["RNA"]] <- split(s.obj[["RNA"]], f = s.obj$name)
  DefaultAssay(s.obj) <- "RNA"
  
  if (use.sctransform == TRUE){
    s.obj <- SCTransform(s.obj, vars.to.regress = vars.to.regress, verbose = FALSE)
    normalization.method <- "SCT" 
  } else {
    s.obj <- NormalizeData(s.obj, normalization.method = "LogNormalize") # ---> use Log Normalized
    s.obj <- FindVariableFeatures(s.obj, selection.method = "vst")
    s.obj <- ScaleData(s.obj, features = rownames(s.obj), vars.to.regress = vars.to.regress)
    normalization.method <- "LogNormalize" 
  }
  if (is.null(remove.genes) == FALSE){
    ##### remove the removed genes from the variable features, e.g TCR or BCR genes.
    s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = TRUE, reduction.name = "RNA_PCA", features = setdiff(VariableFeatures(s.obj), remove.genes))  
  } else {
    s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = TRUE, reduction.name = "RNA_PCA")      
  }

  s.obj <- RunUMAP(s.obj, dims = 1:num.PC.used.in.UMAP, reduction = "RNA_PCA", reduction.name = "umap.unintegrated")  
  print("UMAP finished")
  
  print("Start integration with CCA ...")
  s.obj <- IntegrateLayers( 
    object = s.obj,
    method = CCAIntegration,
    orig.reduction = "RNA_PCA",
    new.reduction = "integrated.cca",
    verbose = TRUE,
    normalization.method = normalization.method
  )
  
  print("Start integration with RPCA ...")
  s.obj <- IntegrateLayers(
    object = s.obj, 
    method = RPCAIntegration,
    orig.reduction = "RNA_PCA", 
    new.reduction = "integrated.rpca", 
    verbose = TRUE, 
    normalization.method = normalization.method
  )
  
  print("Start integration with Harmony ...")
  s.obj <- IntegrateLayers(
    object = s.obj, 
    method = HarmonyIntegration,
    orig.reduction = "RNA_PCA", 
    new.reduction = "harmony",
    verbose = TRUE, 
    normalization.method = normalization.method
  )
  
  print("finished all integrations.")
  all.reductions <- c(
    "integrated.cca",
    "integrated.rpca",
    "harmony"
  )
  
  for (selected.reduction in all.reductions){
    new.reduction.name <- str_replace(selected.reduction, "integrated.", "")
    s.obj <- FindNeighbors(s.obj, dims = 1:num.PC.used.in.Clustering, reduction = selected.reduction)
    s.obj <- FindClusters(s.obj, resolution = cluster.resolution, cluster.name = sprintf("%s.cluster.%s", new.reduction.name, cluster.resolution))
    s.obj <- RunUMAP(
      object = s.obj,
      reduction = selected.reduction,
      dims = 1:num.PC.used.in.UMAP, 
      reduction.name = sprintf("%s_UMAP", new.reduction.name))
  }
  
  DefaultAssay(s.obj) <- "RNA"
  s.obj <- JoinLayers(s.obj)
  DefaultAssay(s.obj) <- "SCT"
  
  if (save.RDS.s8 == TRUE){
    dir.create(file.path(path.to.output, "s8_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s8_output", 
                             paste0(PROJECT, ".output.s8.rds")))
  }
  
  return(s.obj)
}


