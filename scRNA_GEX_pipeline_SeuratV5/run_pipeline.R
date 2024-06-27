gc()
rm(list = ls())
my_random_seed <- 42
# __________VDJ DATA ANYLYSIS PIPELINE__________
PROJECT <- "test_pbmc10k"
path.to.storage <- "/home/hieunguyen/CRC1382/src/single_cell_pipeline_23052023/test_pbmc_10k/ready-to-use"

path.to.main.output <- file.path("/home/hieunguyen/CRC1382/src/single_cell_pipeline_23052023/OUTPUT")
dir.create(path.to.main.output, showWarnings = FALSE)

source("/home/hieunguyen/CRC1382/src/single_cell_pipeline_23052023/scRNA_VDJ_pipeline/main_VDJ_pipeline.R")

path.to.VDJ.input <- file.path(path.to.storage, "VDJ")
path.to.VDJ.output <- file.path(path.to.main.output, "VDJ_output")

summarize_vdj_data(path.to.VDJ.input, 
                   path.to.VDJ.output, 
                   PROJECT, 
                   removeNA=FALSE, 
                   removeMulti=FALSE, 
                   T_or_B = "T")

# __________GEX DATA ANALYSIS PIPELINE__________
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src/single_cell_pipeline_23052023/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")

set.seed(my_random_seed)

source(file.path(path2src, "import_libraries.R"))

source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

path2input <- file.path(path.to.storage, "GEX")

stage_lst <- list()

stage_lst[["test"]] <- c( test = "pbmc10k")


MINCELLS  <- 5
MINGENES  <- 50

save.RDS <- list(s1 = TRUE,
                 s2 = TRUE,
                 s3 = TRUE,
                 s4 = TRUE,
                 s5 = TRUE,
                 s6 = TRUE,
                 s7 = FALSE,
                 s8 = TRUE,
                 s8a = TRUE,
                 s9 = TRUE)

sw <- list(s1 = "on",
           s2 = "on",
           s3 = "on",
           s4 = "on",
           s5 = "on",
           s6 = "on",
           s7 = "off",
           s8 = "off",
           s8a = "on",
           s9 = "on")

rerun <- list(s1 = FALSE, 
              s2 = FALSE,
              s3 = FALSE,
              s4 = FALSE,
              s5 = FALSE,
              s6 = FALSE,
              s7 = FALSE,
              s8 = FALSE,
              s8a = FALSE,
              s9 = FALSE)


pct_mitoceiling <- 10

input.filter.thresholds <- hash()
##### 1st round
input.filter.thresholds[["1st_round"]] <- list( nFeatureRNAfloor = NULL,
                                                nFeatureRNAceiling = NULL,
                                                nCountRNAfloor = NULL,
                                                nCountRNAceiling = NULL,
                                                pct_mitofloor = NULL,
                                                pct_mitoceiling = pct_mitoceiling,
                                                pct_ribofloor = NULL,
                                                pct_riboceiling = NULL,
                                                ambientRNA_thres = 0.5)

##### 2nd round
input.filter.thresholds[["2nd_round"]] <- list( nFeatureRNAfloor = 250,
                                                nFeatureRNAceiling = NULL,
                                                nCountRNAfloor = 2000,
                                                nCountRNAceiling = NULL,
                                                pct_mitofloor = NULL,
                                                pct_mitoceiling = pct_mitoceiling,
                                                pct_ribofloor = NULL,
                                                pct_riboceiling = NULL,
                                                ambientRNA_thres = 0.5)

##### 3rd round
input.filter.thresholds[["3rd_round"]] <- list( nFeatureRNAfloor = 2500,
                                                nFeatureRNAceiling = NULL,
                                                nCountRNAfloor = 5000,
                                                nCountRNAceiling = NULL,
                                                pct_mitofloor = NULL,
                                                pct_mitoceiling = pct_mitoceiling,
                                                pct_ribofloor = NULL,
                                                pct_riboceiling = NULL,
                                                ambientRNA_thres = 0.5)


remove_doublet <- TRUE
path.to.10X.doublet.estimation <- "/media/hieunguyen/HD0/ext_HDD/resources/DoubletEstimation10X.csv"

num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
num.dim.integration <- 25

# remove_XY_genes <- c("Xist","Jpx","Ftx","Tsx","Cnbp2")
remove_XY_genes <- NULL
filtered.barcodes <- NULL

cluster.resolution <- 1

for (analysis.round in c("1st_round")){
  filter.thresholds <- input.filter.thresholds[[analysis.round]]
  
  path.to.anno.contigs <- NULL 
  path.to.count.clonaltype <- NULL
  
  path.to.output <- file.path(path.to.main.output, analysis.round, sprintf("pct_mito_%s_%s", pct_mitoceiling, cluster.resolution))
  dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)
  
  s.obj <- run_pipeline_GEX(path2src=path2src,
                            path2input=file.path(path2input),
                            path.to.logfile.dir=file.path(path.to.output, "logs"),
                            stage_lst=stage_lst,
                            path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                            MINCELLS=MINCELLS,
                            MINGENES=MINGENES,
                            PROJECT=PROJECT,
                            remove_doublet=remove_doublet,
                            save.RDS=save.RDS,
                            path.to.output=path.to.output,
                            rerun=rerun,
                            DE.test="wilcox",
                            num.PCA=num.PCA,
                            num.PC.used.in.UMAP=num.PC.used.in.UMAP,
                            num.PC.used.in.Clustering=num.PC.used.in.Clustering,
                            use.sctransform=FALSE,
                            filtered.barcodes=filtered.barcodes,
                            filter.thresholds=filter.thresholds,
                            path.to.anno.contigs=path.to.anno.contigs,
                            path.to.count.clonaltype=path.to.count.clonaltype,
                            input.method = "filterTRAB",
                            my_random_seed = my_random_seed,
                            remove_XY_genes = remove_XY_genes,
                            num.dim.integration=num.dim.integration,
                            num.dim.cluster=num.dim.cluster,
                            cluster.resolution = cluster.resolution)
}

writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))  


# EOF