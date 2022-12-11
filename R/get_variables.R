get_primary_ensembl_dataset = function(){
  if (is.null(pkgglobalenv$primary_ensembl_dataset)){
    pkgglobalenv$primary_ensembl_dataset = "athaliana_eg_gene"
  }
  return(pkgglobalenv$primary_ensembl_dataset)
}

get_primary_ensembl = function(){
  if (is.null(pkgglobalenv$primary_ensembl)){
    pkgglobalenv$primary_ensembl = biomaRt::useMart(biomart = "plants_mart",
                     dataset = get_primary_ensembl_dataset(),
                     host = "https://plants.ensembl.org")
  }
  return(pkgglobalenv$primary_ensembl)
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
                                                              "PANTHER_description",
                                                              "Molecular_function",
                                                              "Biological_process",
                                                              "Cellular_component",
                                                              "PANTHER_protein_class",
                                                              "Pathways"))
    pkgglobalenv$PANTHER_classifications <- PANTHER_classifications
  }else{PANTHER_classifications = pkgglobalenv$PANTHER_classifications}
  return(PANTHER_classifications)
}

get_NCBI_gene_info = function(){
  if (is.null(pkgglobalenv$NCBI_gene_info)){
  tmp <- tempfile()
  curl::curl_download("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Plants/All_Plants.gene_info.gz", tmp)
  gunzipped=gzfile(tmp,'rt')
  pkgglobalenv$NCBI_gene_info =utils::read.csv(gunzipped, header = TRUE, sep = "\t")
  close(gunzipped)
  }
  return(pkgglobalenv$NCBI_gene_info)
}

get_internal_PD_proteome_dir = function(){
  return(file.path(system.file("extdata", package = "pip1"),"PD_proteomes"))}

get_predictions_cache_path = function(){

  return(file.path(system.file("extdata", package = "pip1"),"predictions_cache.csv"))
}


