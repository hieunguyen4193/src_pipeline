
# __________ Step 1: Input data, preprocessing and QC __________
s1.input.raw.data <- function(path2input,
                              stage_lst,
                              MINCELLS,
                              MINGENES, 
                              PROJECT,
                              save.RDS.s1,
                              path.to.output,
                              filtered.barcodes = NULL,
                              path.to.anno.contigs = NULL,
                              path.to.count.clonaltype = NULL){
  #' Function to read in the raw input data. The data is assumed to have the 
  #' CellRanger output format. 
  #'
  #'
  #' @param path2input Path to the folder containing all experiments data matrices
  #' @param stage_lst A R named list specifying the stage (condition) of each experimental data
  #' @param MINCELLS Integer, the minimum number of cells required for each experiment.
  #' @param MINGENES Integer, the minimum number of genes (features) obtain in each cell in each experiment
  #' @param PROJECT String, name of the project.
  #' 
  #' What do we have in the ouptut of this process?
  #' - All QC plots are stored in the slot "misc" --> all.QC[[plot.name]]
  
  
  # Convert the traditional SEURAT object to my modified SEURAT object
  all_exprs <- Sys.glob(file.path(path2input, "*"))
  
  # assign folder names as name of the experiment data. 
  names(all_exprs) <- to_vec(for(exprs in all_exprs) basename(exprs)) 
  
  # print(paste0("Number of scRNA experiments in this dataset: ", length(all_exprs)))
  
  data.list = list() # a list containing all experiment data. 
  
  for (i in seq_along(all_exprs)){
    path_to_expr <- all_exprs[i]
    input.data <- Read10X(path_to_expr) 
    
    count.data <- input.data$`Gene Expression`
    count.adt <- input.data$`Antibody Capture`
    
    ##### remove genes 
    # nonzero <- count.data > 0
    # count_cells_per_gene <- Matrix::rowSums(nonzero)
    # keep_genes <- Matrix::rowSums(nonzero) >= MINCELLS
    # count.data <- count.data[keep_genes, ]
    
    ##### remove some TRA/TRB genes --------------------------------------------
    TR_genes_patterns <- c("Trav", "Traj", "Trac", "Trbv", "Trbd", "Trbj", "Trbc",
                           "Trgv", "Trgj", "Trgc", "Trdv", "Trdc", "Trdj") 
    
    all_genes <- row.names(count.data)
    
    genes.to.exclude <- unlist(lapply(all_genes, function(x){
      if (substr(x, 1, 4) %in% TR_genes_patterns){
        return(x)
      } else {
        return(NA)
      }
    }))
    genes.to.exclude <- subset(genes.to.exclude, is.na(genes.to.exclude) == FALSE)
    
    tra.genes <- subset(genes.to.exclude, grepl("Tra", genes.to.exclude) == TRUE)
    trb.genes <- subset(genes.to.exclude, grepl("Trb", genes.to.exclude) == TRUE)
    trg.genes <- subset(genes.to.exclude, grepl("Trg", genes.to.exclude) == TRUE)
    trd.genes <- subset(genes.to.exclude, grepl("Trd", genes.to.exclude) == TRUE)
    
    counts.tra <- count.data[which(rownames(count.data) %in% tra.genes),]
    counts.trb <- count.data[which(rownames(count.data) %in% trb.genes),]
    counts.trg <- count.data[which(rownames(count.data) %in% trg.genes),]
    counts.trd <- count.data[which(rownames(count.data) %in% trd.genes),]
    
    count.data <- count.data[-(which(rownames(count.data) %in% genes.to.exclude)),]
    ##### FINISH remove some TRA/TRB genes -------------------------------------
    
    expr.name <- names(all_exprs)[i]
    s.obj <- CreateSeuratObject(counts = count.data, 
                                min.cells = MINCELLS, 
                                min.features = MINGENES, 
                                project = PROJECT)
    
    s.obj@meta.data[, "name"] <- expr.name
    s.obj@meta.data[, "stage"] <- stage_lst[expr.name]
    
    # estimate the percentage of mapped reads to Mitochondrial and Ribosome genes
    s.obj[["percent.mt"]] <- PercentageFeatureSet(s.obj, 
                                                  pattern = "^mt-|^MT-")
    s.obj[["percent.ribo"]] <- PercentageFeatureSet(s.obj, 
                                                    pattern = "^Rpl|^Rps|^RPL|^RPS")
    
    new.count.data <- GetAssayData(s.obj)
    
    count.adt <- count.adt[, colnames(new.count.data)]
    
    adt_assay <- CreateAssayObject(counts = count.adt)
    s.obj[["ADT"]] <- adt_assay
    
    data.list[[i]] <- s.obj
    
    # clean up
    rm(s.obj)
  }
  
  s.obj <- data.list[[1]]
  
  if(length(data.list) > 1){
    # if there are more than 1 experiment, merge them all into 1 single object. 
    s.obj <- merge(x = data.list[[1]], 
                   y = unlist(data.list[2:length(data.list)]),
                   merge.data = FALSE, 
                   add.cell.ids = names(all_exprs), 
                   project = PROJECT)
  }
  s.obj[["No.Exprs"]] <- length(all_exprs)
  
  s.obj@tools[["meta_order"]] <- list(name = names(all_exprs), 
                                      stage=unique(stage_lst))
  
  # Filter barcodes, only apply in round 2 analysis
  if (is.null(filtered.barcodes) == FALSE){
    all_cells <- colnames(s.obj)
    cells.to.keep <- subset(all_cells, all_cells %in% filtered.barcodes == FALSE)
    s.obj <- subset(s.obj, cells = cells.to.keep)
  }
  
  all.QC <- list()
  
  # Number of cells obtained per experiment/sample
  all.QC$cell.counts.plot <- ggplot(s.obj@meta.data, 
                                    aes(x=name , fill=name)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    ggtitle("Number of cells in each dataset")
  
  # distribution of number of UMI in each sample
  all.QC$nCountRNA.distribution <- ggplot(s.obj@meta.data,
                                          aes(color=name, x=nCount_RNA, fill = name)) + 
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    ylab("Cell density") +
    geom_vline(xintercept = 500, color = "red") +
    xlim(0, 40000) + 
    ggtitle("Distribution of read depths in each sample")
  
  # distribution of number of features (genes detected) in each sample
  all.QC$nFeature_RNA.distribution <- ggplot(s.obj@meta.data,
                                             aes(color=name, x=nFeature_RNA, fill = name)) + 
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    ylab("Cell density") +
    geom_vline(xintercept = 1000, color = "red") +
    xlim(0, 10000) +
    ggtitle("Distribution of number of detected genes in each sample")
  
  
  # scatter plot showing the relation between cell read-depth and number of genes detected.
  ## with Mitochondrial percentage
  all.QC$nCount.vs.nFeature.MT <- ggplot(s.obj@meta.data, 
                                         aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm, formula = y ~ x) + # apply a linear regression to show the relation, if existed.
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    facet_wrap(~name) +
    ggtitle("Scatter plot: nCount_RNA vs. nFeature_RNA, cmap % Mitochondrial genes")
  
  ## with Ribosome percentage
  all.QC$nCount.vs.nFeature.Ribo <- ggplot(s.obj@meta.data, 
                                           aes(x=nCount_RNA, y=nFeature_RNA, color=percent.ribo)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm, formula = y ~ x) + # apply a linear regression to show the relation, if existed.
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    facet_wrap(~name) +
    ggtitle("Scatter plot: nCount_RNA vs. nFeature_RNA, cmap % Ribosome genes")
  
  # Complexity: 
  
  ## We can see the samples where we sequenced each cell less have a higher overall complexity, 
  # that is because we have not started saturating the sequencing for any given gene for these samples. 
  # Outlier cells in these samples might be cells that have a less complex RNA species than other cells. 
  # Sometimes we can detect contamination with low complexity cell types like red blood cells via this metric. 
  # Generally, we expect the novelty score to be above 0.80. 
  
  ## More expanations: they are looking for cells that have a low number of genes with a high number of UMI counts. 
  # This likely means that you only captured transcripts from a low number of genes, and simply sequenced transcripts 
  # from those lower number of genes over and over again. This could be because of the cell type 
  # (such as a red blood cell having little to no RNA as they mentioned), or some other strange artifact.
  
  # Compute the complexity and add it to the s.obj@meta.data
  s.obj@meta.data <- s.obj@meta.data %>% 
    mutate(log10GenesPerUMI = log10(nFeature_RNA) / log10(nCount_RNA))
  
  all.QC$complexity <- ggplot(s.obj@meta.data,
                              aes(x=log10GenesPerUMI, color = name, fill=name)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 14)) +
    geom_vline(xintercept = 0.8) +
    ggtitle("Complexity: Log10(nCount_RNA) / log10(nFeature_RNA)")
  
  # add new slot for all.QC into the existed SEURAT OBJECT. 
  s.obj@misc$all.QC <- all.QC
  
  ##############################################################################
  # __________ ADD CLONAL TYPE INFORMATION __________
  if (is.null(path.to.anno.contigs) == FALSE){
    anno.contigs <- read.csv(path.to.anno.contigs)
    
    count.clonaltype <- read.csv(path.to.count.clonaltype)
    
    anno.contigs <- anno.contigs %>% 
      rowwise %>% 
      mutate(barcode = sprintf("%s_%s", sample, tail(unlist(str_split(barcode, pattern = "_")), 1)))
    
    combined.metadata <- merge(slot(s.obj, "meta.data"), anno.contigs, 
                               by.x = 0, by.y = "barcode", all = TRUE)
    
    combined.metadata <- subset(combined.metadata, combined.metadata$Row.names %in% row.names(s.obj@meta.data))
    
    s.obj <- AddMetaData(object = s.obj, metadata = combined.metadata$CTaa,
                         col.name = "CTaa")
  }
  
  # save both: before and after filtering the genes
  if (save.RDS.s1 == TRUE){
    dir.create(file.path(path.to.output, "s1_output"), showWarnings = FALSE)
    saveRDS(object = s.obj, 
            file = file.path(path.to.output, "s1_output", 
                             paste0(PROJECT, ".output.s1.rds")))
    saveRDS(object = genes.to.exclude,
            file = file.path(path.to.output, "s1_output", 
                             paste0(PROJECT, "genes_to_be_excluded.rds")))
  }
  return(s.obj)
}

