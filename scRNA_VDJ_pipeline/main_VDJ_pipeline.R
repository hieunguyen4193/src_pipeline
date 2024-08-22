# __________ scRNA-seq VDJ data analysis pipeline __________
# Trong-Hieu Nguyen, trnguyen@ukaachen.de
#
#
#
#
#
#
#
# __________________________________________________________
rearrange_string_CTaa <- function(s){
  split_s <- str_split(s, "_")[[1]]
  
  for (i in seq(length(split_s))){
    if (grepl(";", (split_s[[i]])) == TRUE){
      selected_index <- i
      selected_string <- split_s[[i]]
    }
  }
  selected_string <- paste(sort(str_split(selected_string, ";")[[1]]), collapse = ";")
  
  if (selected_index == 1){
    output_string <- sprintf("%s_%s", selected_string, split_s[[2]])
  } else {
    output_string <- sprintf("%s_%s", split_s[[1]], selected_string)
  }
  return(output_string)
}

summarize_vdj_data <- function(path.to.input, 
                               path.to.output, 
                               PROJECT, 
                               removeNA=FALSE, 
                               removeMulti=FALSE, 
                               T_or_B,
                               threshold = 0.85){
  dir.create(path.to.output, showWarnings = FALSE)
  
  `%ni%` <- Negate(`%in%`)
  
  path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_VDJ_pipeline/src_dir"
  
  source(file.path(path.to.pipeline.src, "import_libraries.R"))
  
  all.contig.files <- Sys.glob(file.path(path.to.input, "*/filtered_contig_annotations.csv"))
  
  sample.names <- unlist(lapply(lapply(all.contig.files, dirname), basename))
  
  names(all.contig.files) <- sample.names
  
  contig_list <- lapply(all.contig.files, vroom, show_col_type = FALSE)
  
  names(contig_list) <- sample.names
  if (T_or_B == "T"){
    combined.contigs <- combineTCR(contig_list,
                                   samples = sample.names,
                                   ID = sample.names,
                                   removeNA = removeNA,
                                   removeMulti = removeMulti)
  } else if (T_or_B == "B"){
    combined.contigs <- combineBCR(contig_list,
                                   samples = sample.names,
                                   ID = sample.names,
                                   removeNA = removeNA,
                                   removeMulti = removeMulti,
                                   threshold = threshold)
  }
  
  names(combined.contigs) <- sample.names
  
  sample.clonaltype <- hash()
  
  for (sample in names(combined.contigs)){
    combined.contigs[[sample]] <- combined.contigs[[sample]]%>% 
      rowwise %>% 
      mutate(CTaa = ifelse(grepl(";", CTaa) == TRUE, rearrange_string_CTaa(CTaa), CTaa))
    
    count.clonaltype <-  combined.contigs[[sample]] %>% group_by(CTaa) %>% summarise(count = n()) %>% arrange(desc(count))
    
    sample.clonaltype[[sample]] <- count.clonaltype
    
    write.csv(combined.contigs[[sample]], file.path(path.to.output,
                                                    sprintf("annotated_contigs_clonaltype_%s.csv", sample)))
    
    write.csv(sample.clonaltype[[sample]], file.path(path.to.output,
                                                     sprintf("count_clonaltype_%s.csv", sample)))
  }
  saveRDS(object = combined.contigs, file = file.path(path.to.output, sprintf("%s_combined_contigs.rds", PROJECT)))
}

# EOF