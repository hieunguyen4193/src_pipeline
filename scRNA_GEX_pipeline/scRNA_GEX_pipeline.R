# __________ scRNA-seq data analysis pipeline __________
# Trong-Hieu Nguyen, trnguyen@ukaachen.de
#
#
#
#
#
#
#

`%ni%` = Negate(`%in%`)
# __________ HELPER FUNCTIONS __________

# __________ generate data table in rendered HTML file __________
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}

run_pipeline_GEX <- function(path2src,
                         path2input,
                         path.to.logfile.dir,
                         stage_lst,
                         path.to.10X.doublet.estimation,
                         MINCELLS,
                         MINGENES,
                         PROJECT,
                         remove_doublet,
                         save.RDS,
                         path.to.output,
                         rerun, 
                         DE.test,
                         num.PCA,
                         num.PC.used.in.UMAP,
                         num.PC.used.in.Clustering,
                         use.sctransform,
                         filtered.barcodes=NULL,
                         filter.thresholds,
                         input.method,
                         path.to.anno.contigs,
                         path.to.count.clonaltype,
                         cluster.resolution=0.5,
                         mode_cell_cycle_scoring = "gene_name",
                         my_random_seed = 42,
                         umap.method = "uwot",
                         ambientRNA_thres = 0.5,
                         path.to.renv = NULL,
                         features_to_regressOut=NULL,
                         regressOut_mode="alternative",
                         num.dim.integration=30,
                         remove_XY_genes,
                         with.VDJ,
                         genes.to.not.run.PCA = NULL,
                         inte_pca_reduction_name = "INTE_PCA",
                         inte_umap_reduction_name = "INTE_UMAP",
                         pca_reduction_name = NULL,
                         umap_reduction_name = NULL,
                         path.to.s3a.source = NULL,
                         path.to.h5.file = NULL,
                         path.to.h5.meta.data = NULL,
                         with.TSNE = FALSE,
                         k.filter = 200,
                         sw = NULL){
  
  # load renv.lock file
  # require(renv)
  # renv::restore(lockfile = path.to.renv)
  
  dir.create(path.to.output, showWarnings = FALSE)
  dir.create(path.to.logfile.dir, showWarnings = FALSE)
  
  # The source file for step 1 is specified in the process
  source(file.path(path2src, "s2_ambient_RNA_correction.R"))
  if (is.null(path.to.s3a.source) == FALSE){
    source(path.to.s3a.source)
  } else {
    source(file.path(path2src, "s3_first_filter.R"))    
  }
  source(file.path(path2src, "s4_Doublet_detection.R"))
  source(file.path(path2src, "s5_preprocess_before_cellCycle_scoring.R"))
  source(file.path(path2src, "s6_cell_cycle_scoring.R"))
  source(file.path(path2src, "s7_cellFactor_regressOut.R"))
  source(file.path(path2src, "s8_integration_and_clustering.R"))
  source(file.path(path2src, "s8a_cluster_without_integration.R"))
  source(file.path(path2src, "s9_findMarkers.R"))
  
  # Save SessionIinfo to human readable file
  writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))
  
  log.file <- sprintf(file.path(path.to.logfile.dir, sprintf("%s_%s.log", PROJECT, Sys.time())))
  log.file <- str_replace_all(log.file, ":", "_")
  log.file <- str_replace(log.file, " ", "_")
  
  parameter.log.file <- sprintf(file.path(path.to.logfile.dir, sprintf("%s_%s.parameter.log", PROJECT, Sys.time())))
  parameter.log.file <- str_replace_all(parameter.log.file, ":", "_")
  parameter.log.file <- str_replace(parameter.log.file, " ", "_")
  
  write(sprintf("LOG FILE GENERATED ON %s \n", Sys.time()), file=parameter.log.file, append=TRUE)
  
  if ((sw[["s8"]] == "on") & (length(stage_lst) == 1)){
    write(sprintf("There is only 1 sample in the input data, turn off INTEGRATION AND CLUSTERING process..."), file=parameter.log.file, append=TRUE)
    sw[["s8"]] <- "off"
  }
  
  write("\n \n ########## ALL INPUT PARAMETERS ########## \n", file=parameter.log.file, append=TRUE)
  write(sprintf("Using renv: path.to.renv = %s \n", path.to.renv), file=parameter.log.file, append=TRUE)
  write(sprintf("Path to src.: path2src = %s \n", path2src), file=parameter.log.file, append=TRUE)
  write(sprintf("Path to input data: path2input = %s \n", path2input), file=parameter.log.file, append=TRUE)
  write(sprintf("Path to logs dir.: path.to.logfile.dir = %s \n", path.to.logfile.dir), file=parameter.log.file, append=TRUE)
  write(sprintf("Path to 10x doublet estimation: %s \n", path.to.10X.doublet.estimation), file=parameter.log.file, append=TRUE)
  write(sprintf("MINCELLS = %s, MINGENES = %s, PROJECT_NAME = %s \n", MINCELLS, MINGENES, PROJECT), file=parameter.log.file, append=TRUE)
  # all filter thresholds
  for (crit in names(filter.thresholds)){
    write(sprintf("%s = %s\n", crit, filter.thresholds[[crit]]), file=parameter.log.file, append=TRUE)
  }
  
  write(sprintf("Path to output dir.: %s \n", path.to.output), file=parameter.log.file, append=TRUE)
  write(sprintf("Use SCTRANSFORM: %s \n", use.sctransform), file=parameter.log.file, append=TRUE)
  write(sprintf("Pct.mitochondrial: %s", filter.thresholds$pct_mitoceiling),file=parameter.log.file, append=TRUE)
  write(sprintf("Pct.AmbientRNA: %s", ambientRNA_thres),file=parameter.log.file, append=TRUE)
  write(sprintf("UMAP method: %s", umap.method),file=parameter.log.file, append=TRUE)
  write(sprintf("Differential gene expression test: %s \n", DE.test), file=parameter.log.file, append=TRUE)
  write(sprintf("LOG FILE GENERATED ON %s", Sys.time()), 
        file=log.file, 
        append=TRUE)
  
  if (sw[["s1"]] == "on"){
    write("\n \n ########## PROCESS S1: PREPROCESSING AND QUALITY CONTROL ##########", 
          file=log.file, 
          append=TRUE)
    
    # __________ Process S1: Preprocessing and Quality control __________
      if (rerun[["s1"]] == TRUE | file.exists(file.path(path.to.output, "s1_output", sprintf("%s.output.s1.rds", PROJECT)))  == FALSE){
        print("Running process S1: Preprocessing and generate QC results ...")
        
        if (input.method == "normal"){
          source(file.path(path2src, "s1_preprocessing_QC.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output,
                                     path.to.anno.contigs=path.to.anno.contigs,
                                     path.to.count.clonaltype=path.to.count.clonaltype,
                                     with.VDJ=with.VDJ)
          write("Using source file from s1_preprocessing_QC.R", 
                file=log.file, 
                append=TRUE)
          
        } else if (input.method == "normal_with_added_VDJ"){
          source(file.path(path2src, "s1_preprocessing_QC_normal_with_VDJ.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output,
                                     path.to.anno.contigs = path.to.anno.contigs,
                                     path.to.count.clonaltype = path.to.count.clonaltype,
                                     filtered.barcodes = filtered.barcodes)
          write("Using source file from s1_preprocessing_QC_normal_with_VDJ.R", 
                file=log.file, 
                append=TRUE)
          
        } else if (input.method == "filterTRAB"){
          source(file.path(path2src, "s1_preprocessing_QC_modified_filterTRAB.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output, 
                                     path.to.anno.contigs = path.to.anno.contigs,
                                     path.to.count.clonaltype = path.to.count.clonaltype,
                                     filtered.barcodes = filtered.barcodes)
          write("Using source file from s1_preprocessing_QC_modified_filterTRAB.R", 
                file=log.file, 
                append=TRUE)
          
        } else if (input.method =="filterTRAB_clonaltype"){
          source(file.path(path2src, "s1_preprocessing_QC_excluded_noClonotype_cells.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output,
                                     path.to.anno.contigs = path.to.anno.contigs,
                                     path.to.count.clonaltype = path.to.count.clonaltype,
                                     filtered.barcodes = filtered.barcodes)
          write("Using source file from s1_preprocessing_QC_excluded_noClonotype_cells.R", 
                file=log.file, 
                append=TRUE)
          write(sprintf("Using annotated contigs file: %s", path.to.anno.contigs), 
                file=log.file, 
                append=TRUE)
          write(sprintf("Using count clonaltype file: %s", path.to.count.clonaltype), 
                file=log.file, 
                append=TRUE)
        } else if (input.method == "CITESEQ"){
          source(file.path(path2src, "s1_preprocessing_QC_CITESEQ.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output,
                                     filtered.barcodes = filtered.barcodes)
          write("Using source file from s1_preprocessing_QC_CITESEQ.R", 
                file=log.file, 
                append=TRUE)
          write(sprintf("Filtered %s cells based on the result of 1st rounds", length(filtered.barcodes)), 
                file=log.file, 
                append=TRUE)
        } else if (input.method == "CITESEQ_ensembl"){
          source(file.path(path2src, "s1_preprocessing_QC_CITESEQ_ensembl.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output,
                                     filtered.barcodes = filtered.barcodes)
          write("Using source file from s1_preprocessing_QC_CITESEQ.R", 
                file=log.file, 
                append=TRUE)
          write(sprintf("Filtered %s cells based on the result of 1st rounds", length(filtered.barcodes)), 
                file=log.file, 
                append=TRUE)
        } else if (input.method == "from_txt") {
          source(file.path(path2src, "s1_preprocessing_from_txt_counts.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output)
          write("Using source file from s1_preprocessing_from_txt_counts.R", 
                file=log.file, 
                append=TRUE)
        } else if (input.method == "filterIg") {
          source(file.path(path2src, "s1_preprocessing_QC_modified_filterIg.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output, 
                                     path.to.anno.contigs = path.to.anno.contigs,
                                     path.to.count.clonaltype = path.to.count.clonaltype,
                                     filtered.barcodes = filtered.barcodes)
          write("Using source file from s1_preprocessing_QC_modified_filterIg.R", 
                file=log.file, 
                append=TRUE)
        } else if (input.method == "filterTRAB_with_hashtag_Antibody"){
          source(file.path(path2src, "s1_preprocessing_QC_modified_filterTRAB_with_hastagAntibody.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output, 
                                     path.to.anno.contigs = path.to.anno.contigs,
                                     path.to.count.clonaltype = path.to.count.clonaltype,
                                     filtered.barcodes = filtered.barcodes)
          write("Using source file from s1_preprocessing_QC_modified_filterTRAB_with_hastagAntibody.R", 
                file=log.file, 
                append=TRUE)
        } else if (input.method == "filterTCR_BCR_with_hashtag_Antibody") {
          source(file.path(path2src, "s1_preprocessing_QC_modified_filterTRAB_BCR_with_hastagAntibody.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output, 
                                     path.to.anno.contigs = path.to.anno.contigs,
                                     path.to.count.clonaltype = path.to.count.clonaltype,
                                     filtered.barcodes = filtered.barcodes)
          write("Using source file from s1_preprocessing_QC_modified_filterTRAB_BCR_with_hastagAntibody.R", 
                file=log.file, 
                append=TRUE)
        } else if (input.method == "readH5") {
          source(file.path(path2src, "s1_preprocessing_QC_readH5.R"))
          s.obj <- s1.input.raw.data(path.to.h5.file = path.to.h5.file, 
                                     path.to.h5.meta.data = path.to.h5.meta.data, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output)
          write("Using source file from s1_preprocessing_QC_readH5.R", 
                file=log.file, 
                append=TRUE)
        } else if (input.method == "remove_XY_genes"){
          source(file.path(path2src, "s1_preprocessing_QC_remove_XY_genes.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output,
                                     remove_XY_genes = remove_XY_genes)
          write("Using source file from s1_preprocessing_QC_remove_XY_genes.R", 
                file=log.file, 
                append=TRUE)
        } else if (input.method == "normal_with_hashtag_antibodies"){
          source(file.path(path2src, "s1_preprocessing_QC_with_hastagAntibody.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output, 
                                     path.to.anno.contigs = path.to.anno.contigs,
                                     path.to.count.clonaltype = path.to.count.clonaltype,
                                     filtered.barcodes = filtered.barcodes)
          write("Using source file from s1_preprocessing_QC_with_hastagAntibody.R", 
                file=log.file, 
                append=TRUE)
        } else if (input.method == "filterTRAB_rename_mm10_genes"){
          source(file.path(path2src, "s1_preprocessing_QC_modified_filterTRAB_rename_mm10_genes.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output, 
                                     path.to.anno.contigs = path.to.anno.contigs,
                                     path.to.count.clonaltype = path.to.count.clonaltype,
                                     filtered.barcodes = filtered.barcodes)
          
        } else if (input.method == "filterTRABGD_with_hashtag_Antibody"){
          source(file.path(path2src, "s1_preprocessing_QC_modified_filterTRABGD_with_hastagAntibody.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output, 
                                     path.to.anno.contigs = path.to.anno.contigs,
                                     path.to.count.clonaltype = path.to.count.clonaltype,
                                     filtered.barcodes = filtered.barcodes)
        } else if (input.method == "filterIG_with_hashtagAntibody"){
          source(file.path(path2src, "s1_preprocessing_QC_filterIG_with_hastagAntibody.R"))
          s.obj <- s1.input.raw.data(path2input = path2input, 
                                     stage_lst = stage_lst, 
                                     MINCELLS = MINCELLS, 
                                     MINGENES = MINGENES,
                                     PROJECT = PROJECT,
                                     save.RDS.s1 = save.RDS[["s1"]],
                                     path.to.output = path.to.output, 
                                     path.to.anno.contigs = path.to.anno.contigs,
                                     path.to.count.clonaltype = path.to.count.clonaltype,
                                     filtered.barcodes = filtered.barcodes)
        }

        status.message <- sprintf("New output is saved at %s", file.path(path.to.output, "s1_output", sprintf("%s.output.s1.rds", PROJECT)))
        
        } else {
        s.obj <- readRDS(file.path(path.to.output, "s1_output", sprintf("%s.output.s1.rds", PROJECT)))
        status.message <- sprintf("Seurat object was loaded from %s", file.path(path.to.output, "s1_output", sprintf("%s.output.s1.rds", PROJECT)))
      }
     write(status.message,file=log.file, append=TRUE)
  }

  # __________ Process S2: Ambient RNA correction __________
  if (sw[["s2"]] == "on"){
    write("\n \n ########## PROCESS S2: AMBIENT RNA CONTAMINATION CORRECTION ##########", 
          file=log.file, 
          append=TRUE)
    if (rerun[["s2"]] == TRUE | file.exists(file.path(path.to.output, "s2_output", sprintf("%s.output.s2.rds", PROJECT)))  == FALSE){
      print("Running process S2: Correcting ambient RNA contamination ...")
      s.obj <- s2_ambient_RNA_correction(s.obj = s.obj, 
                                         chosen.method = "decontX",
                                         path.to.output = path.to.output,
                                         save.RDS.s2 = save.RDS[["s2"]],
                                         PROJECT = PROJECT
                                         )
      
      status.message <- sprintf("New output is saved at %s", file.path(path.to.output, "s2_output", sprintf("%s.output.s2.rds", PROJECT)))
      
    } else {
      s.obj <- readRDS(file.path(path.to.output, "s2_output", sprintf("%s.output.s2.rds", PROJECT)))
      status.message <- sprintf("Seurat object was loaded from %s", file.path(path.to.output, "s2_output", sprintf("%s.output.s2.rds", PROJECT)))
    }
    write(status.message,file=log.file, append=TRUE)
  }

  # __________ Process S3: FILTER __________
  if (sw[["s3"]] == "on"){
    write("\n \n ########## PROCESS S3: FILTER ##########", 
          file=log.file, 
          append=TRUE)
    if (is.null(path.to.s3a.source) == FALSE){
      write(sprintf("Using additional filter criteria at %s", path.to.s3a.source),file=log.file, append=TRUE)
    }
    if (rerun[["s3"]] == TRUE | file.exists(file.path(path.to.output, "s3_output", sprintf("%s.output.s3.rds", PROJECT)))  == FALSE){
      s.obj <- s3.filter( s.obj,
                          PROJECT = PROJECT,
                          path.to.output,
                          save.RDS[["s3"]], 
                          nFeatureRNAfloor = filter.thresholds$nFeatureRNAfloor, 
                          nFeatureRNAceiling = filter.thresholds$nFeatureRNAceiling,
                          nCountRNAfloor = filter.thresholds$nCountRNAfloor, 
                          nCountRNAceiling = filter.thresholds$nCountRNAceiling,
                          pct_mitofloor = filter.thresholds$pct_mitofloor, 
                          pct_mitoceiling = filter.thresholds$pct_mitoceiling,
                          pct_ribofloor = filter.thresholds$pct_ribofloor, 
                          pct_riboceiling = filter.thresholds$pct_riboceiling,
                          ambientRNA_thres = filter.thresholds$ambientRNA_thres,
                          log10GenesPerUMI_thres = filter.thresholds$log10GenesPerUMI)
      status.message <- sprintf("New output is saved at %s", file.path(path.to.output, "s3_output", sprintf("%s.output.s3.rds", PROJECT)))
      
    } else {
      print("Loading existing results ... s3 ")
      s.obj <- readRDS(file.path(path.to.output, "s3_output", sprintf("%s.output.s3.rds", PROJECT)))
      
      status.message <- sprintf("Seurat object was loaded from %s", file.path(path.to.output, "s3_output", sprintf("%s.output.s3.rds", PROJECT)))
    }
    write(status.message,file=log.file, append=TRUE)
    
  } else {
    write("\n \n ########## PROCESS S3: FILTER ##########", 
          file=log.file, 
          append=TRUE)
    status.message <- "Process s3-Filter has not been run"
    print(status.message)
    write(status.message,file=log.file, append=TRUE)
  }

  # __________ Process S4: Doublet detection __________
  if (sw[["s4"]] == "on"){
    write("\n \n ########## PROCESS S4: DOUBLET DETECTION ##########", 
          file=log.file, 
          append=TRUE)
    
    if (rerun[["s4"]] == TRUE | file.exists(file.path(path.to.output, "s4_output", sprintf("%s.output.s4.rds", PROJECT)))  == FALSE){
      print("Running process S4: Doublet detection ...")
      s.obj <- s4.DoubletDetection(s.obj, 
                                   path.to.10X.doublet.estimation,
                                   save.RDS[["s4"]],
                                   path.to.output, 
                                   remove_doublet = FALSE,
                                   PROJECT = PROJECT)

      
      status.message <- sprintf("New output is saved at %s", file.path(path.to.output, "s4_output", sprintf("%s.output.s4.rds", PROJECT)))
      
    } else {
      s.obj <- readRDS(file.path(path.to.output, "s4_output", sprintf("%s.output.s4.rds", PROJECT)))
      status.message <- sprintf("Seurat object was loaded from %s", file.path(path.to.output, "s4_output", sprintf("%s.output.s4.rds", PROJECT)))
    }
    write(status.message,file=log.file, append=TRUE)
  }
  
  # __________ Process S5: Preprocessing before cell Cycle scoring __________
  if (sw[["s5"]] == "on"){
    write("\n \n ########## PROCESS S5: PREPROCESSING BEFORE CELL CYCLE SCORING ##########", 
          file=log.file, 
          append=TRUE)
    
    if (rerun[["s5"]] == TRUE | file.exists(file.path(path.to.output, "s5_output", sprintf("%s.output.s5.rds", PROJECT)))  == FALSE){
      print("Running process S5: Preprocessing before cell cycle scoring ...")
      
      s.obj <- s5.preprocess.before.cellCycle.scoring(s.obj,
                                                      path.to.output,
                                                      save.RDS[["s5"]],
                                                      PROJECT = PROJECT,
                                                      use.sctransform = use.sctransform)
        
      status.message <- sprintf("New output is saved at %s", file.path(path.to.output, "s5_output", sprintf("%s.output.s5.rds", PROJECT)))
      
    } else {
      s.obj <- readRDS(file.path(path.to.output, "s5_output", sprintf("%s.output.s5.rds", PROJECT)))
      status.message <- sprintf("Seurat object was loaded from %s", file.path(path.to.output, "s5_output", sprintf("%s.output.s5.rds", PROJECT)))
    }
    write(status.message,file=log.file, append=TRUE)
  }
  
  # __________ Process S6: Cell cycle scoring __________
  if (sw[["s6"]] == "on"){
    write("\n \n ########## PROCESS S6: CELL CYCLE SCORING ##########", 
          file=log.file, 
          append=TRUE)
    
    if (rerun[["s6"]] == TRUE | file.exists(file.path(path.to.output, "s6_output", sprintf("%s.output.s6.rds", PROJECT)))  == FALSE){
      print("Running process S6: Running cell cycle scoring ...")
      
      s.obj <- s6.cellCycleScoring(s.obj,
                                    path.to.output,
                                    save.RDS[["s6"]],
                                   PROJECT = PROJECT, mode = mode_cell_cycle_scoring)
      
      status.message <- sprintf("New output is saved at %s", file.path(path.to.output, "s6_output", sprintf("%s.output.s6.rds", PROJECT)))
      
    } else {
      print("Loading existing results ... s6")
      s.obj <- readRDS(file.path(path.to.output, "s6_output", sprintf("%s.output.s6.rds", PROJECT)))
      status.message <- sprintf("Seurat object was loaded from %s", file.path(path.to.output, "s6_output", sprintf("%s.output.s6.rds", PROJECT)))
    }
    write(status.message,file=log.file, append=TRUE)
  }
  
  # __________ Process S7: Cell cycle factor regress out __________
  if (sw[["s7"]] == "on"){
    write("\n \n ########## PROCESS S7: CELL CYCLE FACTOR REGRESS OUT ##########", 
          file=log.file, 
          append=TRUE)
    
    if (rerun[["s7"]] == TRUE | file.exists(file.path(path.to.output, "s7_output", sprintf("%s.output.s7.rds", PROJECT)))  == FALSE){
      print("Running process S7: Running cell cycle scoring ...")
      
      s.obj <- s7.cellFactorRegressOut(s.obj,
                                   path.to.output,
                                   save.RDS[["s7"]],
                                   PROJECT = PROJECT, 
                                   features_to_regressOut = features_to_regressOut,
                                   regressOut_mode = regressOut_mode)
      
      status.message <- sprintf("New output is saved at %s", file.path(path.to.output, "s7_output", sprintf("%s.output.s7.rds", PROJECT)))
      
    } else {
      s.obj <- readRDS(file.path(path.to.output, "s7_output", sprintf("%s.output.s7.rds", PROJECT)))
      status.message <- sprintf("Seurat object was loaded from %s", file.path(path.to.output, "s7_output", sprintf("%s.output.s7.rds", PROJECT)))
    }
    write(status.message,file=log.file, append=TRUE)
  }
  
  # __________ Process S8: INTEGRATION AND CLUSTERING __________
  if (sw[["s8"]] == "on"){
    write("\n \n ########## PROCESS S8: INTEGRATION AND CLUSTERING ##########", 
          file=log.file, 
          append=TRUE)
    
    if (rerun[["s8"]] == TRUE | file.exists(file.path(path.to.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))  == FALSE){
      print("Running process S8: Running integration and clustering ...")
      
      s.obj <- s8.integration.and.clustering(s.obj, 
                                             path.to.output, 
                                             save.RDS[["s8"]],
                                             PROJECT, 
                                             num.dim.integration,
                                             num.PCA,
                                             num.PC.used.in.UMAP,
                                             num.PC.used.in.Clustering,
                                             cluster.resolution,
                                             my_random_seed,
                                             umap.method,
                                             genes.to.not.run.PCA,
                                             inte_pca_reduction_name, 
                                             inte_umap_reduction_name,
                                             with.TSNE,
                                             k.filter)
      
      status.message <- sprintf("New output is saved at %s", file.path(path.to.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
      
    } else {
      s.obj <- readRDS(file.path(path.to.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
      status.message <- sprintf("Seurat object was loaded from %s", file.path(path.to.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
    }
    write(status.message, file=log.file, append=TRUE)
  } else {
    write("Process S8 has not been run!", 
          file=log.file, 
          append=TRUE)
  }
  
  # __________ Process S8: INTEGRATION AND CLUSTERING __________
  if (sw[["s8a"]] == "on"){
    write("\n \n ########## PROCESS S8a: CLUSTERING WITHOUT INTEGRATION ##########", 
          file=log.file, 
          append=TRUE)
    
    if (rerun[["s8a"]] == TRUE | file.exists(file.path(path.to.output, "s8a_output", sprintf("%s.output.s8a.rds", PROJECT)))  == FALSE){
      print("Running process S8a: Running clustering without annotation ...")
      
      s.obj <- s8a.cluster.wo.integration(s.obj, 
                                          path.to.output, 
                                          save.RDS[["s8a"]],
                                          PROJECT,
                                          num.PCA,
                                          num.PC.used.in.UMAP,
                                          num.PC.used.in.Clustering,
                                          cluster.resolution,
                                          my_random_seed,
                                          umap.method,
                                          genes.to.not.run.PCA,
                                          pca_reduction_name,
                                          umap_reduction_name)
      
      status.message <- sprintf("New output is saved at %s", file.path(path.to.output, "s8a_output", sprintf("%s.output.s8a.rds", PROJECT)))
      
    } else {
      s.obj <- readRDS(file.path(path.to.output, "s8a_output", sprintf("%s.output.s8a.rds", PROJECT)))
      status.message <- sprintf("Seurat object was loaded from %s", file.path(path.to.output, "s8a_output", sprintf("%s.output.s8a.rds", PROJECT)))
    }
    write(status.message, file=log.file, append=TRUE)
  } else {
    write("Process S8a has not been run!", 
          file=log.file, 
          append=TRUE)
  }
  
  
  # __________ Process S9: FIND CLUSTER MARKERS __________
  if (sw[["s9"]] == "on"){
    write("\n \n ########## PROCESS S9: FIND CLUSTER MARKERS ##########", 
          file=log.file, 
          append=TRUE)
    
    if (rerun[["s9"]] == TRUE | file.exists(file.path(path.to.output, "s9_output", sprintf("%s.clusterMarkers.s9.rds", PROJECT)))  == FALSE){
      print("Running process S9: Running find cluster markers ...")
      
      cluster.markers <- s9.findClusterMarkers(s.obj, 
                                                path.to.output, 
                                                save.RDS[["s9"]],
                                                PROJECT, 
                                                DE.test)
      
      status.message <- sprintf("New output is saved at %s", file.path(path.to.output, "s9_output", sprintf("%s.clusterMarkers.s9.rds", PROJECT)))
      
    } else {
      s.obj <- readRDS(file.path(path.to.output, "s9_output", sprintf("%s.clusterMarkers.s9.rds", PROJECT)))
      status.message <- sprintf("Seurat object was loaded from %s", file.path(path.to.output, "s9_output", sprintf("%s.clusterMarkers.s9.rds", PROJECT)))
    }
    write(status.message, file=log.file, append=TRUE)
  } else {
    write("Process S9-Find All markers has not been run!", 
          file=log.file, 
          append=TRUE)
  }
  
  return(s.obj)
}

  
  