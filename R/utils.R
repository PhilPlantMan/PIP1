pkgglobalenv <- new.env(parent=emptyenv())

#TODO Write test for this
choose_directory = function(caption = 'Select data directory') {
  if (Sys.info()['sysname'] == 'Windows') {
    utils::choose.dir(caption = caption)
  } else {
    tcltk::tk_choose.dir(caption = caption)
  }
}

choose_dataset = function(){
  datasets = pantherCompatibleEnsemlDatasets()
  description = utils::select.list(datasets$description, graphics = TRUE, title = "Select species to generate for")
  dataset = datasets$dataset[datasets$description == description]
  return(dataset)
}

get_primary_ensembl_dataset = function(){
  if (is.null(pkgglobalenv$primary_ensembl_dataset)){
    pkgglobalenv$primary_ensembl_dataset = "athaliana_eg_gene"
  }
  return(pkgglobalenv$primary_ensembl_dataset)
}

get_primary_ensembl = function(){
  if (is.null(pkgglobalenv$primary_ensembl)){
    biomaRt::useMart(biomart = "plants_mart",
                     dataset = get_primary_ensembl_dataset(),
                     host = "https://plants.ensembl.org")
  }
}

get_primary_PANTHER_classification = function(){
  if (is.null(pkgglobalenv$primary_PANTHER_classification)){
    pkgglobalenv$primary_PANTHER_classification = 'subfamily'}
  return(pkgglobalenv$primary_PANTHER_classification)
}

get_primary_PANTHER_classification_column_name = function(){
  if (get_primary_PANTHER_classification() == 'subfamily'){
    PANTHER_classification_column_name = "hmmpanther_subfams"}
  if (get_primary_PANTHER_classification() == 'family'){
    PANTHER_classification_column_name = "hmmpanther"}
  return(PANTHER_classification_column_name)
}

get_PANTHER_classifications = function(){
  if (is.null(pkgglobalenv$PANTHER_classifications)){
    url = "ftp://ftp.pantherdb.org/hmm_classifications/17.0/PANTHER17.0_HMM_classifications"
    tmp <- tempfile()
    curl::curl_download(url, tmp)
    PANTHER_classifications = utils::read.table(tmp,fill = TRUE, sep = "\t",
                                                comment.char = "",quote = "",
                                                col.names = c("PANTHER_ID",
                                                              "DESCRIPTION",
                                                              "Molecular_Function",
                                                              "Biological_Process",
                                                              "Cellular_Component",
                                                              "PANTHER_Protein_Class",
                                                              "Pathways"))
    pkgglobalenv$PANTHER_classifications <- PANTHER_classifications
  }else{PANTHER_classifications = pkgglobalenv$PANTHER_classifications}
  return(PANTHER_classifications)
}

get_internal_PD_proteome_dir = function(){
  return(file.path(system.file("extdata", package = "pip1"),"PD_proteomes"))}

output_directories = function(output_dir, dataset_choice){
  parent = file.path(output_dir, paste0(Sys.Date(), "_", dataset_choice))
  counts_name = paste0(get_primary_PANTHER_classification(),"_counts")
  verified_panther_members_name = paste0("verified_", get_primary_PANTHER_classification(),"_members")
  directory_list = list(parent = parent,
                        counts = file.path(parent, counts_name),
                        verified_panther_members = file.path(parent, verified_panther_members_name))
  for (dir in directory_list){
    dir.create(dir)
  }
  pkgglobalenv$output_directories = directory_list
  return(directory_list)
}

read_spreadsheet = function(path){
  ext = tools::file_ext(path)
  if (ext == "csv"){
    df = utils::read.csv(path)
  }
  if (ext == "xlsx"){
    df = readxl::read_excel(path)
  }
  return(df)
}

stop_if_no_ensembl_gene_id_column = function(df, df_description){
  if (!"ensembl_gene_id" %in% colnames(df)){
    stop(paste("ensembl_gene_id column is required in",path))
  }
}

stop_if_no_NCBI_taxon_column = function(df, df_description){
  if (!"NCBI_taxon" %in% colnames(df)){
    stop(paste("NCBI_taxon column is required in",path))
  }
}


download_NCBI_gene_info = function(){
  tmp <- tempfile()
  curl::curl_download("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Plants/All_Plants.gene_info.gz", tmp)
  gunzipped=gzfile(tmp,'rt')
  gene_info =utils::read.csv(gunzipped, header = TRUE, sep = "\t")
  close(gunzipped)
  return(gene_info)
}

annotate_with_biomaRt = function(df, taxon_id, attributes, filter){
  dataset = pantherCompatibleEnsemlDatasets()[pantherCompatibleEnsemlDatasets()$taxon_id == taxon_id, 'dataset']
  ensembl = biomaRt::useMart(biomart = "plants_mart",
                             dataset = dataset,
                             host = "https://plants.ensembl.org")
  ensembl_annotations <- biomaRt::getBM(attributes = attributes,
                                        filters = filter,
                                        values = df[,filter],
                                        mart = ensembl)
  joined_df = dplyr::left_join(df, ensembl_annotations, by = filter)
  return(joined_df)
}

annotate_with_ensembl_gene_id_with_PANTHER = function(df, taxon_id){
  df_raw_annotated = annotate_with_biomaRt(df, taxon_id,
                                           c("ensembl_gene_id", "hmmpanther"),
                                           "ensembl_gene_id")
  df_raw_annotated = dplyr::select(df_raw_annotated,ensembl_gene_id, hmmpanther)
  contains_colon = grepl(":", df_raw_annotated$hmmpanther, fixed = TRUE)
  df.hmmpanther_subfams = df_raw_annotated[contains_colon,]
  df.hmmpanther_subfams = dplyr::rename(df.hmmpanther_subfams,
                                        hmmpanther_subfams = hmmpanther)
  df.hmmpanther = df_raw_annotated[!contains_colon,]
  df_annotated = df %>%
    dplyr::left_join(df.hmmpanther, by = "ensembl_gene_id") %>%
    dplyr::left_join(df.hmmpanther_subfams, by = "ensembl_gene_id") %>%
    dplyr::distinct()
  return(df_annotated)
}

annotate_proteome = function(proteome_path){
  df = read_spreadsheet(proteome_path)
  stop_if_no_ensembl_gene_id_column(df, path)
  stop_if_no_NCBI_taxon_column(df, path)
  taxon_id = df$NCBI_taxon[1]
  df_annotated = annotate_with_ensembl_gene_id_with_PANTHER(df, taxon_id)
}

proteomes_PANTHER_counts = function(proteome_dir, output_dir = pkgglobalenv$output_directories$counts){
  file_list = list.files(proteome_dir, pattern = ".csv|.xlsx", full.names = TRUE)
  for (file in file_list){
    message(paste("Performing counts for", file))
    annotated_proteome = annotate_proteome(file)
    file_sans_ext = tools::file_path_sans_ext(basename(file))
    counts_column = paste0(get_primary_PANTHER_classification(),"_counts")
    df.fam_counts = annotated_proteome %>% dplyr::count(!!rlang::sym(get_primary_PANTHER_classification_column_name()),
                                                        name = counts_column)
    file_path = file.path(output_dir, paste0(file_sans_ext,".csv"))
    utils::write.csv(df.fam_counts, file_path, row.names = FALSE)
    proteome_column_name = paste(file_sans_ext, counts_column, sep = ".")
    df.fam_counts = dplyr::rename(df.fam_counts, !!proteome_column_name := dplyr::all_of(counts_column))
    if (file == file_list[1]){concat_df = df.fam_counts
    }else{concat_df = dplyr::full_join(concat_df, df.fam_counts, by = get_primary_PANTHER_classification_column_name())}
  }
  concat_df = concat_df %>% dplyr::mutate_if(is.numeric,dplyr::coalesce,0)
  concat_path =  file.path(output_dir, paste0(get_primary_PANTHER_classification_column_name(),"_cross-referenced.csv"))
  utils::write.csv(concat_df, concat_path, row.names = FALSE)
  return(concat_df)
}

#TODO tidy this function
export_Euler_diagrams = function(cross_referenced_counts, output_dir = pkgglobalenv$output_directories$counts){
  df = cross_referenced_counts[,-1]
  df[df > 0] = 1
  fit <- eulerr::euler(df, shape = "ellipse")

  png(filename = file.path(output_dir,paste0(get_primary_PANTHER_classification(), '_euler.png')),height=500,width=1000)
  print(plot(fit, quantities =  list(type = c("counts")),
             fill = "transparent",
             edges = list(col = "black", lex = 2))) #,labels = rep('', dim(df.bool_cols)[2])
  dev.off()

  png(filename= file.path(output_dir,paste0(get_primary_PANTHER_classification(), '_no-label-euler.png')),height=500,width=1000)
  print(plot(fit, quantities =  list(type = c("counts")),
             fill = "transparent",
             edges = list(col = "black", lex = 2),
             labels = rep('', dim(df)[2]))) #,labels = rep('', dim(df.bool_cols)[2])
  dev.off()

  png(filename = file.path(output_dir,paste0(get_primary_PANTHER_classification(), '_venn.png')),height=500,width=700)
  print(plot(eulerr::venn(df),
             edges = list(col = "black", lex = 2),
             environment = environment())) # quantities =  list(type = c("percent", "counts"))
  dev.off()

  png(filename = file.path(output_dir,paste0(get_primary_PANTHER_classification(), '_venn no label.png')),height=500,width=700)
  print(plot(eulerr::venn(df),
             edges = list(col = "black", lex = 2),
             environment = environment(),
             labels = rep('', dim(df)[2]))) # quantities =  list(type = c("percent", "counts"))
  dev.off()
}

pull_PANTHER_members = function(df,
                                panther_column = get_primary_PANTHER_classification_column_name(),
                                ensembl = get_primary_ensembl()){
  ensembl_annotations <- biomaRt::getBM(attributes = c("ensembl_gene_id","hmmpanther"),
                                        filters = "hmmpanther",
                                        values = unique(df[,panther_column]),
                                        mart = ensembl)
  ensembl_annotations = dplyr::rename(ensembl_annotations,
                                      !!get_primary_PANTHER_classification_column_name():= hmmpanther)
  joined_df = dplyr::left_join(df, ensembl_annotations, by = get_primary_PANTHER_classification_column_name())
  return(joined_df)
}

