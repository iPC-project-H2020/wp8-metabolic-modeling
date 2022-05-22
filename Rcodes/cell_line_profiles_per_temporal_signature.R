cell_line_profiles = function(working_dir, input_file, model_dir, dest_dir, model_name,
                   ic_score_by = NaN, 
                   aggregation_by = "median"){
  
  dest_dir = paste0(dest_dir, "/cl_aggregation_by_", aggregation_by)
  if (file.exists(dest_dir)){
    cat("Folder already exists.. \n")
  } else {
    dir.create(dest_dir, recursive= TRUE)
  }
  ## sample annotation
  sample_annotation <- read.delim(paste0(working_dir,"/scRNAseq_data/sample_annotation.txt"))

  ## gene symbol to gene GEM identifier
  model.gene.description <- read.table(file = paste0(model_dir, "/gene.tsv"), sep="\t", header = TRUE)
  gene_number <- model.gene.description$gene_number
  names(gene_number) <- model.gene.description$symbol

  aggregation_meth = get(aggregation_by)
  
  # read the input data
  dat <- read.delim(input_file)
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1]

  ## remove genes that do not have a GEM gene id mapping
  data_gene_number <- gene_number[rownames(dat)]
  sel.rm <- which(is.na(data_gene_number))
  dat <- dat[-sel.rm, ]
  # change gene id coding
  rownames(dat) <- data_gene_number[-sel.rm]
  
  rownames(sample_annotation) <- sample_annotation$CELL
  cell.not.ann <- which(!(colnames(dat) %in% sample_annotation$CELL))
  if (length(cell.not.ann)!=0){
    annDat <- sample_annotation[colnames(dat[,-cell.not.ann]), ]
    dat <- dat[, -cell.not.ann]
  } else {
    annDat <- sample_annotation[colnames(dat), ]
  }
  
  for (gr.annot in unique(annDat$ANNOT)){
    sel.gr <- which(annDat$ANNOT == gr.annot)
    cellsInGr <- annDat$CELL[sel.gr]
    datGr <- dat[, cellsInGr]
    gr.profile <- apply(datGr, 1, aggregation_meth)
    write.table(gr.profile,
                file = paste0(dest_dir,
                              "/", gr.annot, ".txt"))
    if (grepl("pdx", gr.annot)|grepl("Dox", gr.annot)){ 
      if (!file.exists(paste0(working_dir,"/scRNAseq_data/", gr.annot, ".txt"))){
        write.table(datGr, file = paste0(working_dir,"/scRNAseq_data/", gr.annot, ".txt"))
      }
    }
  }
}

