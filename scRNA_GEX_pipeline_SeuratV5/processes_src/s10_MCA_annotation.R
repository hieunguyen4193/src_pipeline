generate_scrna_MCAannotate <- function(scrna){
  ret_code = 0
  suppressPackageStartupMessages(require(scMCA))
  mca_result <- scMCA(GetAssayData(object=scrna, slot="counts"), numbers_plot = 3)
  #pattern = paste0("(", gsub("-", "_", MCA_NAME), ")")
  pattern = paste0("\\(", MCA_NAME, "\\)")
  corr=mca_result$cors_matrix[grep(pattern,rownames(mca_result$cors_matrix)),]
  if(dim(corr)[1] == 0){
    logger.error("Cannot find MCA name, please check!")
    stop("exit 1")
  }
  res=rownames(corr)[max.col(t(corr))]
  scrna$MCA_annotate <- res
  Idents(object=scrna) <- "name"
  return(list(scrna, ret_code))
}