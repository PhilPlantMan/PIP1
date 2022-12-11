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

stop_if_incompatable_taxon = function(taxon_id, df_description){
  if (!taxon_id %in% pantherCompatibleEnsemlDatasets()[,'taxon_id']){
    stop(paste("NCBI_taxon",
               taxon_id,
               "in", df_description,
               "is not in compatible datasets. See ?pantherCompatibleEnsemlDatasets for more information"))
  }
}

