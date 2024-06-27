s9.findClusterMarkers <- function(s.obj, 
                                  path.to.output, 
                                  save.RDS.s9,
                                  PROJECT, 
                                  DE.test,
                                  find.conserved.markers = FALSE){
  DefaultAssay(s.obj) <- "RNA"
  
  markers <- FindAllMarkers(object = s.obj, test.use = DE.test)
  conserved.markers <- hash()
  if (find.conserved.markers == TRUE){
    
    
    for (cluster.id in unique(s.obj$seurat_clusters)){
      conserved.markers[[cluster.id]] <- FindConservedMarkers(object = s.obj, 
                                                              ident.1 = cluster.id,
                                                              grouping.var = "stage")
    }
  } else {
    conserved.markers <- NULL
  }
  
  cluster.markers <- list(diff.markers = markers,
                          conserved.markers = conserved.markers)
  
  if (save.RDS.s9 == TRUE){
    dir.create(file.path(path.to.output, "s9_output"), showWarnings = FALSE)
    saveRDS(object = cluster.markers, 
            file = file.path(path.to.output, "s9_output", 
                             paste0(PROJECT, ".clusterMarkers.s9.rds")))
  }
  return(cluster.markers)
}