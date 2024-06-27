s10.PROGENy_pathway_analysis_between_conditions <- function( s.obj, 
                                                            path.to.output, 
                                                            progeny.params,
                                                            cluster.name, 
                                                            save.RDS.s10,
                                                            scran_or_seurat = "scran",
                                                            PROJECT, 
                                                            condition1, 
                                                            condition2){
  #### FUNCTION FORKED FROM https://github.com/CostaLab/scrna_seurat_pipeline/blob/master/data_factory.R
  Idents(s.obj) <- cluster.name
  
  s.obj <- progeny::progeny(s.obj, scale = progeny.params$scale, 
                            organism = progeny.params$species, 
                            top = progeny.params$top, 
                            perm = progeny.params$perm, 
                            return_assay = progeny.params$return_assay)
  DefaultAssay(s.obj) <-'progeny'
  s.obj <- ScaleData(s.obj, features = row.names(s.obj))
  all.pathways <- rownames(s.obj@assays$progeny)
  
  Idents(s.obj) <- "stage"
  
  if (is.factor(s.obj@meta.data[, cluster.name])){
    s.obj@meta.data[, cluster.name] <- droplevels(s.obj@meta.data[, cluster.name])
  }else{
    c_names <- unique(s.obj@meta.data[, cluster.name])
    s.obj@meta.data[, cluster.name] <- factor(s.obj@meta.data[, cluster.name],
                                              levels=sort(c_names))
  }
  
  diff.pathways <- data.frame()
  
  for (cluster.id in levels(s.obj@meta.data[[cluster.name]])){
    cells.condition1 <- row.names(subset(s.obj@meta.data, s.obj@meta.data$stage == condition1 & s.obj@meta.data[[cluster.name]] == cluster.id))
    cells.condition2 <- row.names(subset(s.obj@meta.data, s.obj@meta.data$stage == condition2 & s.obj@meta.data[[cluster.name]] == cluster.id))
    
    subset.s.obj <- subset(s.obj, cells = c(cells.condition1, cells.condition2))
    g <- as.character(subset.s.obj@meta.data$stage)
    g <- factor(g, levels=c(condition1, condition2))
    tmp <- scran::findMarkers(as.matrix(subset.s.obj@assays$progeny@data), g)[[1]] %>%
      as.data.frame() 
    tmp$cluster <- cluster.id
    r <- sapply(all.pathways, 
                function(pw) rcompanion::wilcoxonR(as.vector(subset.s.obj@assays$progeny@data[pw,]), g))
    names(r) <- to_vec(for (item in names(r)) str_replace(item, "[.]r", "")[[1]][[1]]) 
    tmp[names(r), "r"] <- r
    tmp <- tmp[names(r), ] %>% rownames_to_column("pathway")
    diff.pathways <- rbind(diff.pathways, tmp)
  }
  diff.pathways$p_val_adj <- p.adjust(diff.pathways$p.value, method = "fdr")
  diff.pathways$tag <- sapply(diff.pathways$p_val_adj, function(pval) {
    if(pval< 0.001) {
      txt <- "***"
    } else if (pval < 0.01) {
      txt <- "**"
    } else if (pval < 0.05) {
      txt <- "*"
    } else {
      txt <- "NS"
    }
    return(txt)
  })
  
  # save the differential pathway activity score analysis to the main seurat object
  s.obj@misc[[sprintf("progeny.diff.pathways_%s_vs_%s", condition1, condition2)]] <- diff.pathways
  
  if (save.RDS.s10 == TRUE){
    dir.create(file.path(path.to.output, "s10_output"), showWarnings = FALSE)
    saveRDS(object = cluster.markers, 
            file = file.path(path.to.output, "s10_output", 
                             paste0(PROJECT, ".PROGENy.s10.rds")))
  }
  return(s.obj)
}
