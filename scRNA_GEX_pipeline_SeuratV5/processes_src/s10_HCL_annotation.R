generate_scrna_HCLannotate <- function(scrna){
  ret_code = 0
  suppressPackageStartupMessages(require(scHCL))
  hcl_result <- scHCL(GetAssayData(object=scrna, slot="counts"), numbers_plot = 3)
  pattern = gsub("-",".",HCL_NAME)
  corr=hcl_result$cors_matrix[grep(pattern,rownames(hcl_result$cors_matrix),fixed=TRUE),]
  rnms <- gsub("_[a-zA-Z]+\\.$", "", rownames(corr))
  rnms <- gsub(paste0("\\.", pattern, "\\."), "", rnms)
  rnms <- gsub("_.*high", "", rnms)
  rnms <- gsub("[0-9]\\.$", "", rnms)
  rnms <- gsub("\\.[0-9]$", "", rnms)
  rnms <- gsub("[0-9]\\.$", "", rnms)
  rnms <- gsub("\\.\\.", ".", rnms)
  rnms <- gsub("\\.$", "", rnms)
  rownames(corr) <- rnms
  if(dim(corr)[1] == 0){
    logger.error("Cannot find HCL name, please check!")
    stop("exit 1")
  }
  res=rownames(corr)[max.col(t(corr))]
  scrna$HCL_annotate <- res
  Idents(object=scrna) <- "name"
  return(list(scrna, ret_code))
}
