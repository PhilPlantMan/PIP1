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

get_PANTHER_classifications = function(){
  url = "ftp://ftp.pantherdb.org/hmm_classifications/17.0/PANTHER17.0_HMM_classifications"
  tmp <- tempfile()
  curl::curl_download(url, tmp)
  PantherClassifications = utils::read.table(tmp,fill = TRUE, sep = "\t",
                                      comment.char = "",quote = "",
                                      col.names = c("PANTHER_ID",
                                                    "DESCRIPTION",
                                                    "Molecular_Function",
                                                    "Biological_Process",
                                                    "Cellular_Component",
                                                    "PANTHER_Protein_Class",
                                                    "Pathways"))
  return(PantherClassifications)
}

output_directories = function(output_dir, dataset_choice){
  parent = file.path(output_dir, paste0(Sys.Date(), "_", dataset_choice))
  directory_list = list(parent = parent,
                        counts = file.path(parent, "counts"))
  for (dir in directory_list){
    dir.create(dir)
    return(directory_list)
  }
}

