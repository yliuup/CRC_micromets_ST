for (i in c("PAT2","PAT3","PAT4","PAT5","PAT6","PAT7","PAT8","PAT9","PAT10","PAT11")){
 
  dat_new <- get(paste0(i))
  ep_data <- as.matrix(GetAssayData(object = dat_new, slot = "data"))
  
  gene_order_file= read.table("gencode_v21_gene_pos.txt",header = FALSE, 
                              row.names = 1, sep = "\t", check.names = FALSE)
  colnames(ep_data) <- gsub(pattern = "-",replacement = "_",x = rownames(dat_new@meta.data))
  ep_data <- ep_data[rownames(ep_data) %in% rownames(gene_order_file),]
  
  write.table(ep_data,file = paste0(i,"_data_sample.matrix"),sep = "\t",quote = F)
  
  ep_sample_annotation <- data.frame(colnames(dat_new),dat_new$anno %>% as.character())
  colnames(ep_sample_annotation) <-c("cell_names","cluster")
  ep_sample_annotation$cell_names <- gsub(pattern = "-",replacement = "_",
                                          x = ep_sample_annotation$cell_names)
  write.table(ep_sample_annotation,file = paste0(i,"_sample_annotation.txt"),sep = "\t",col.names = F,row.names = F,quote=F)
  
}



library(infercnv)
i <- commandArgs(trailingOnly = T)
annotations_file=read.table(paste0(i,"_sample_annotation.txt"),
                            sep = "\t",header = FALSE, row.names = 1,
                            stringsAsFactors = FALSE, colClasses = "character")

ref_names <- unique(annotations_file$V2) 
ref_names <- ref_names[ref_names %in% ref_all]

infercnv_obj = infercnv::CreateInfercnvObject(
  raw_counts_matrix=as.matrix(read.table(paste0(i,"_data_sample.matrix"), sep = "\t",header = TRUE, row.names = 1)),
  annotations_file=read.table(paste0(i,"_sample_annotation.txt"),
                              sep = "\t",header = FALSE, row.names = 1,
                              stringsAsFactors = FALSE, colClasses = "character"),
  delim="\t", 
  gene_order_file= read.table("gencode_v21_gene_pos.txt",header = FALSE, 
                              row.names = 1, sep = "\t", check.names = FALSE),
  ref_group_names=ref_names) 


dir.create(i)
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=paste0(i), 
                             num_threads = 15,
                             cluster_by_groups=F, 
                             denoise=TRUE,
                             HMM=F)

save(infercnv_obj,file = "infercnv_obj.RData")



