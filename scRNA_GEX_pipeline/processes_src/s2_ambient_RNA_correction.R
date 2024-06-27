s2_ambient_RNA_correction <- function(s.obj, 
                                      chosen.method = "decontX",
                                      path.to.output,
                                      save.RDS.s2,
                                      PROJECT){
  
  #' Function to perform the Estimation of Ambient RNA background profile in the scRNA-seq dataset
  #'
  #'
  #' @param s.obj The input Seurat object
  #' @param chosen.method The method used in the Ambient RNA estimation process. 
  #' This can be decontX or SoupX.
  #' 
  
  if (chosen.method == "decontX"){
    # Convert the Seurat object into Bioconductor::SingleCellExperiment object 
    s.obj.sce <- as.SingleCellExperiment(s.obj)
    
    # Perform decontX in each sample in this dataset
    s.obj.decontX <- decontX(s.obj.sce, batch = s.obj$name)
    
    # Add the Ambient-RNA estimation back to the original Seurat object. 
    s.obj[["decontX"]] <- CreateAssayObject(counts = s.obj.decontX@assays@data$decontXcounts)
    
    s.obj <- AddMetaData(object = s.obj, metadata = s.obj.decontX$decontX_contamination,
                         col.name = "AmbientRNA")
    
    s.obj <- AddMetaData(object = s.obj, metadata = s.obj.decontX$decontX_clusters,
                         col.name = "decontX_clusters")
    
    ambient.cluster.RNA <- list()
    ambient.contamination <- list()
    
    for (sample.umap in reducedDimNames(s.obj.decontX)){
      umap <- reducedDim(s.obj.decontX, sample.umap)
      
      plot.ambientRNA.cluster <- plotDimReduceCluster(x = s.obj$decontX_clusters,
                                                      dim1 = umap[, 1], dim2 = umap[, 2]) +
        ggtitle(sprintf("Dim. reduce cluster, Sample: %s", sample.umap))
      ambient.cluster.RNA[[paste0(sample.umap, ".plot")]] <- plot.ambientRNA.cluster 
      
    }
    for (sample.name in s.obj$name){
      ambient.contamination[[sample.name]] <- plotDecontXContamination(s.obj.decontX, batch = sample.name) +
        ggtitle(sprintf("Contamination level in each cell in Sample: %s", sample.name))
    }  
    
    # add all the ambient.cluster.RNA.plot to the Seurat object
    s.obj@misc$ambient.cluster.RNA.plot <- ambient.cluster.RNA
    s.obj@misc$ambient.contamination.plot <- ambient.contamination
    
  } else if (chosen.method == "SoupX") {
    print("not finished")
    
  } else {
    stop("User must specify either decontX or SoupX to be the Ambient RNA estimation method.")
    
  }
  if (save.RDS.s2 == TRUE){
    dir.create(file.path(path.to.output, "s2_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s2_output", 
                             paste0(PROJECT, ".output.s2.rds")))
  }
  return(s.obj)
}