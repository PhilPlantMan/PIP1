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

annotate_PANTHER_mixed_taxons = function(df, df_description){
  stop_if_no_ensembl_gene_id_column(df, df_description)
  stop_if_no_NCBI_taxon_column(df, df_description)
  for (taxon_id in unique(df$NCBI_taxon)){
    stop_if_incompatable_taxon(taxon_id, df_description)
    df.subset = df[df$NCBI_taxon == taxon_id,]
    df.subset.PANTHER = annotate_with_ensembl_gene_id_with_PANTHER(df.subset, taxon_id)
    if (taxon_id == unique(df$NCBI_taxon)[1]){
      concat_df = df.subset.PANTHER
    }else{concat_df = rbind(concat_df, df.subset.PANTHER)}
  }
  return(concat_df)
}

pull_verified_PD_PANTHER_members = function(df){
  df_description = "known_PD_genes"
  annotated_df = annotate_PANTHER_mixed_taxons(df, df_description)
  annotated_df = dplyr::rename(annotated_df, verified_gene = ensembl_gene_id)
  pulled_panther_members_df = pull_PANTHER_members(annotated_df)
  filename = paste0("verified_PD_", get_primary_PANTHER_classification(),"_mems.csv")
  save_path = file.path(pkgglobalenv$output_directories$verified_panther_members, filename)
  message(save_path)
  utils::write.csv(pulled_panther_members_df, save_path, row.names = FALSE)
}
