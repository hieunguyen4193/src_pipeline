s3.filter <- function(s.obj,
                      PROJECT,
                      path.to.output,
                      save.RDS.s3, 
                      nFeatureRNAfloor, 
                      nFeatureRNAceiling,
                      nCountRNAfloor, 
                      nCountRNAceiling,
                      pct_mitofloor, 
                      pct_mitoceiling,
                      pct_ribofloor, 
                      pct_riboceiling,
                      ambientRNA_thres,
                      log10GenesPerUMI_thres){
  
  if (is.null(nFeatureRNAfloor) == FALSE){
    s.obj <- subset(s.obj, subset = nFeature_RNA > nFeatureRNAfloor)
  }
  
  if(is.null(nFeatureRNAceiling) == FALSE){
    s.obj <- subset(s.obj, subset = nFeature_RNA < nFeatureRNAceiling)
  }
  
  if (is.null(nCountRNAfloor) == FALSE){
    s.obj <- subset(s.obj, subset = nCount_RNA > nCountRNAfloor)
  }
  
  if(is.null(nCountRNAceiling) == FALSE){
    s.obj <- subset(s.obj, subset = nCount_RNA < nCountRNAceiling)
  }
  
  if (is.null(pct_mitofloor) == FALSE){
    s.obj <- subset(s.obj, subset = percent.mt > pct_mitofloor)
  }
  
  if (is.null(pct_mitoceiling) == FALSE){
    s.obj <- subset(s.obj, subset = percent.mt < pct_mitoceiling)
  }
  
  if (is.null(pct_ribofloor) == FALSE){
    s.obj <- subset(s.obj, subset = percent.ribo > pct_ribofloor)    
  }
  if (is.null(pct_riboceiling) == FALSE){
    s.obj <- subset(s.obj, subset = percent.ribo < pct_riboceiling)    
  }
  
  if (is.null(ambientRNA_thres) == FALSE){
    s.obj <- subset(s.obj, subset = AmbientRNA < ambientRNA_thres)    
  }
  
  if (is.null(log10GenesPerUMI_thres) == FALSE){
    s.obj <- subset(s.obj, subset = log10GenesPerUMI >= log10GenesPerUMI_thres)    
  }
  
  if (save.RDS.s3 == TRUE){
    dir.create(file.path(path.to.output, "s3_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s3_output", 
                             paste0(PROJECT, ".output.s3.rds")))
  }
  return(s.obj)
}