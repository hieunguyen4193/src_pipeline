s5.preprocess.before.cellCycle.scoring <- function(s.obj,
                                                   path.to.output,
                                                   save.RDS.s5,
                                                   PROJECT, 
                                                   use.sctransform = FALSE){
  s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern = "^mt-|^MT-")
  s.obj[["percent.ribo"]] <- PercentageFeatureSet(s.obj, pattern = "^Rpl|^Rps|^RPL|^RPS")
  
  mt.genes <- grep(pattern = "^mt-|^MT-", x = rownames(x = s.obj), value = TRUE)
  ribo.genes <- grep("^Rpl|^Rps|^RPL|^RPS",  x = rownames(x = s.obj), value = TRUE)
  
  s.obj[["percent.exclude"]]  <- PercentageFeatureSet(s.obj, features = c(mt.genes, ribo.genes))
  
  if (use.sctransform == TRUE){
    s.obj <- SCTransform(s.obj, vars.to.regress = c("percent.mt", "percent.ribo"), verbose = FALSE)
    
  } else {
    s.obj <- NormalizeData(s.obj) # ---> use Log Normalized
    s.obj <- FindVariableFeatures(s.obj, selection.method = "vst")
    s.obj <- ScaleData(s.obj, features = rownames(s.obj))
    
  }
  
  
  all.genes <- rownames(x = s.obj)
  
  s.genes <- paste0("^", cc.genes$s.genes, "$", collapse = "|")
  s.genes <- all.genes[grepl(s.genes, all.genes, ignore.case = TRUE)]
  
  g2m.genes <- paste0("^", cc.genes$g2m.genes, "$", collapse = "|")
  g2m.genes <- all.genes[grepl(g2m.genes, all.genes, ignore.case = TRUE)]
  
  # s.obj <- RunPCA(s.obj, features = c(s.genes, g2m.genes), reduction.name="RAW_PCA")
  
  if (save.RDS.s5 == TRUE){
    dir.create(file.path(path.to.output, "s5_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s5_output", 
                             paste0(PROJECT, ".output.s5.rds")))
  }
  return(s.obj)
}