choose_directory = function(caption = 'Select data directory') {
  if (Sys.info()['sysname'] == 'Windows') {
    choose.dir(caption = caption)
  } else {
    tk_choose.dir(caption = caption)
  }
}

choose_dataset = function(){
  datasets = pantherCompatibleEnsemlDatasets()
  description = select.list(datasets$description, graphics = TRUE, title = "Select species to generate for")
  dataset = datasets$dataset[datasets$description == description]
  return(dataset)
}

get_PANTHER_classifications = function(){
  url = "ftp://ftp.pantherdb.org/hmm_classifications/17.0/PANTHER17.0_HMM_classifications"
  tmp <- tempfile()
  curl::curl_download(url, tmp)
  PantherClassifications = read.table(tmp,fill = TRUE, sep = "\t",
                                      comment.char = "",quote = "",
                                      col.names = c("PANTHERID",
                                                    "DESCRIPTION",
                                                    "Molecular_Function",
                                                    "Biological_Process",
                                                    "Cellular_Component",
                                                    "PANTHER_Protein_Class",
                                                    "Pathways"))
  return(PantherClassifications)
}

