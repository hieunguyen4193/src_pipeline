s7.cellFactorRegressOut <- function(s.obj,
                                    path.to.output, 
                                    save.RDS.s7,
                                    PROJECT,
                                    features_to_regressOut, 
                                    regressOut_mode = "alternative"){
  if (is.null(features_to_regressOut) == FALSE){
    all.genes <- rownames(x = s.obj)
    
    s.genes <- paste0("^", cc.genes$s.genes, "$", collapse = "|")
    s.genes <- all.genes[grepl(s.genes, all.genes, ignore.case = TRUE)]
    
    g2m.genes <- paste0("^", cc.genes$g2m.genes, "$", collapse = "|")
    g2m.genes <- all.genes[grepl(g2m.genes, all.genes, ignore.case = TRUE)]
    
    ### Scale data add exclude
    s.obj$CC.Difference <- s.obj$S.Score - s.obj$G2M.Score
    if (regressOut_mode == "normal"){
      s.obj <- ScaleData(s.obj, vars.to.regress = features_to_regressOut, features = rownames(s.obj))
    } else if (regressOut_mode == "alternative"){
      s.obj <- ScaleData(s.obj, vars.to.regress = "CC.Difference", features = rownames(s.obj))
    }

    if (save.RDS.s7 == TRUE){
      dir.create(file.path(path.to.output, "s7_output"), showWarnings = FALSE)
      saveRDS(object = s.obj, 
              file = file.path(path.to.output, "s7_output", 
                               paste0(PROJECT, ".output.s7.rds")))
    }
    return(s.obj)  
  } else {
    return(s.obj)
  }
  
}
