s4.DoubletDetection <- function(s.obj, 
                                path.to.10X.doublet.estimation,
                                save.RDS.S4,
                                path.to.output, 
                                remove_doublet = FALSE,
                                PROJECT){
  
  # estimation.10X <- read.csv("static/DoubletEstimation10X.csv")
  # model.recovered <- approxfun(x = estimation.10X$CellsRecovered, y = estimation.10X$MultipletRate, rule = 2)
  # ---> approxfun: Given a list of points, return a linearly interpolated function. 
  
  # previously we have merged all assays data into one single Seurat object, 
  # now we have to split them into separated assay and perform Doublet detection for each of them. 
  
  # path.to.10X.doublet.estimation <- "C:/main/offline/src/DoubletEstimation10X.csv"
  
  number_of_cells <- table(s.obj$name)
  
  estimation.10X <- read.csv(path.to.10X.doublet.estimation)
  
  # based on the information provided by 10X Genomics, we interpolate and approximate 
  # the Doublet proportion. The more cells you have, the higher doublet proportion you may have. 
  model.recovered <- approxfun(x = estimation.10X$CellsRecovered, 
                               y = estimation.10X$MultipletRate, rule = 2)
  
  doublet_formation_rate <- model.recovered(number_of_cells)
  
  names(doublet_formation_rate) <- names(number_of_cells)
  
  splitted.s.obj <- SplitObject(s.obj, split.by = "name")
  
  splitted.s.obj <- lapply(splitted.s.obj, function(x){
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x <- ScaleData(x)
    x <- RunPCA(x)
    x <- RunUMAP(x, dims = 1:20)
  })
  
  sweep.res.list <- lapply(splitted.s.obj, paramSweep_v3, PCs = 1:20,  sct = FALSE)
  
  sweep.stats <- lapply(sweep.res.list, summarizeSweep, GT = FALSE)
  
  bcmvn <- lapply(sweep.stats, find.pK)
  
  # a little helper function to run Doulbet detection ----------------------------
  run.Doublet.Detection <- function(s.obj, doublet_formation_rate, pN = 0.25, pK = 9.09){
    # obj.name <- names(table(s.obj$name))
    obj.name <- unique(s.obj$name)
    nExp_poi <- doublet_formation_rate[[obj.name]]
    tmp <- doubletFinder_v3(s.obj, PCs = 1:10, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    tmp[["classifications"]] <- factor(tmp[[paste("DF.classifications", pN, pK, nExp_poi, sep = "_")]][,1],
                                       levels = c("Singlet", "Doublet"))
    return(tmp)
  }
  # ------------------------------------------------------------------------------
  
  s.obj.DoubletDetection <- lapply( splitted.s.obj, 
                                    run.Doublet.Detection, 
                                    doublet_formation_rate = doublet_formation_rate)
  
  # collect all Singlet/Doublet classification from all assays
  classifications <- c()
  cells <- c()
  for(i in 1:length(s.obj.DoubletDetection)){
    cells               <- c(names(s.obj.DoubletDetection[[i]]$classifications), cells)
    classifications     <- c(as.character(s.obj.DoubletDetection[[i]]$classifications), classifications)
  }
  classifications <- factor(classifications, levels = c("Singlet", "Doublet"))
  names(classifications) <- cells
  
  # add the classification back to the original merged object s.obj
  s.obj <- AddMetaData(s.obj, classifications, col.name = "Doublet_classifications")
  
  # remove doublet
  if (remove_doublet == TRUE){
    s.obj <- subset(s.obj, Doublet_classifications == "Singlet")
  }
  
  if (save.RDS.S4 == TRUE){
    dir.create(file.path(path.to.output, "s4_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s4_output", 
                             paste0(PROJECT, ".output.s4.rds")))
  }
  
  return(s.obj)
}