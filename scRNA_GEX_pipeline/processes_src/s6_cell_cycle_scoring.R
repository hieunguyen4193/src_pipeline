s6.cellCycleScoring <- function(s.obj, 
                                path.to.output, 
                                save.RDS.s6,
                                PROJECT, mode = "gene_name"){
  library(org.Hs.eg.db)

  
  s.obj <- ScaleData(s.obj, features = row.names(s.obj))
  
  if (mode == "gene_name"){
    all.genes <- rownames(x = s.obj)
    s.genes <- paste0("^", cc.genes$s.genes, "$", collapse = "|")
    s.genes <- all.genes[grepl(s.genes, all.genes, ignore.case = TRUE)]
    g2m.genes <- paste0("^", cc.genes$g2m.genes, "$", collapse = "|")
    g2m.genes <- all.genes[grepl(g2m.genes, all.genes, ignore.case = TRUE)]
    
    s.obj <- CellCycleScoring(s.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    
  } else if (mode == "ensembl"){
    all.genes <- rownames(x = s.obj)
    cc.genes.ensembl <- list()
    cc.genes.ensembl$s.genes <- mapIds(org.Hs.eg.db,
                                keys=cc.genes$s.genes,
                                column="ENSEMBL",
                                keytype="SYMBOL",
                                multiVals="first")
    cc.genes.ensembl$g2m.genes <- mapIds(org.Hs.eg.db,
                                keys=cc.genes$g2m.genes,
                                column="ENSEMBL",
                                keytype="SYMBOL",
                                multiVals="first")
    
    s.obj <- CellCycleScoring(s.obj, s.features = cc.genes.ensembl$s.genes, g2m.features = cc.genes.ensembl$g2m.genes, set.ident = TRUE)
  }
  
  s.obj$G1.Score = 1 - s.obj$S.Score - s.obj$G2M.Score
  s.obj$CC.Difference <- s.obj$S.Score - s.obj$G2M.Score
  
  if("percent.mt" %ni% names(s.obj@meta.data)){ ## for existed clusters analysis
    s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, pattern = "^mt-|^MT-")
    }
  
  if("percent.ribo" %ni% names(s.obj@meta.data)){
    s.obj[["percent.ribo"]] <- PercentageFeatureSet(s.obj, pattern = "^Rpl|^Rps|^RPL|^RPS")
  }
  
  if (save.RDS.s6 == TRUE){
    dir.create(file.path(path.to.output, "s6_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s6_output", 
                             paste0(PROJECT, ".output.s6.rds")))
  }
  
  return(s.obj)
}