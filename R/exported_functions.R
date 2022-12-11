#TODO
#Cross reference will other files provided


#' run_PIP1
#'
#' Main function to generate an in silco proteome for a given species
#'
#' @param proteome_dir path of directory containing proteomes: example datasets is set as default
#' @param known_PD_genes path of .xlsx containing experimentally verified PD proteins. example datasets is set as default
#' @param output_dir path where PIP1 will provide output: dialog will appear by default
#' @param species_dataset ensemble dataset to use; dialog will appear by default
#' @param primary_PANTHER_classification string: either "subfamily" or "family"
#' @param minimum_num_proteomes integer: Minimum number of proteomes for inclusion in list A and B
#' @param PIPAC_or_PIPBD_filter string ("secreted", "secreted_gpi_or_tmhmm", "secreted_or_gpi_or_tmhmm" , "tm_or_gpi", "gpi", "tm", "targetP_not_C_or_M", "targetP_not_C_or_M_or_dash") : logic choice to group output to PIPA or PIPC (if true) or PIPB or PIPD (if false). See associated publication.
#' @param annotate_with_NCBI_geneInfo bool: if true, will download NCBI data to get gene description and synonyms. Slow
#'
#' @return Returns null
#' @export
#'
#' @examples run_PIP1()
run_PIP1 = function(proteome_dir = file.path(system.file("extdata", package = "pip1"),"PD_proteomes"),
                    known_PD_genes = read_spreadsheet(file.path(system.file("extdata", package = "pip1"),"Known_PD_genes.xlsx")),
                    output_dir = choose_directory("Choose output directory"),
                    species_dataset = choose_dataset(),
                    primary_PANTHER_classification = "subfamily",
                    minimum_num_proteomes = 2,
                    PIPAC_or_PIPBD_filter = "secreted_or_gpi_or_tmhmm",
                    annotate_with_NCBI_geneInfo = FALSE
){

  #set variables to use if used in automatic mode i.e. during testing
  if (!interactive()){
    species_dataset = "athaliana_eg_gene"
    output_dir = withr::local_tempdir()
  }
  pip_output = list()
  #set package level cache
  pkgglobalenv$primary_ensembl_dataset = species_dataset
  pkgglobalenv$primary_PANTHER_classification = primary_PANTHER_classification


  #main pipelin
  output_directories = output_directories(output_dir, species_dataset)
  pip_output$proteome_counts = proteomes_PANTHER_counts(proteome_dir)
  pip_output$pulled_PD_PANTHER_members = pull_verified_PD_PANTHER_members(known_PD_genes)
  pip_output$pulled_candidates = pull_PANTHER_members(pip_output$proteome_counts)
  pip_output$pulled_candidates_annotated = full_annotation_boil_down(pip_output$pulled_candidates)

  export_Euler_diagrams(pip_output$proteome_counts)

  #TODO final_annotations including annotate_with_NCBI_geneInfo
return(pip_output)
}



#' pantherCompatibleEnsemlDatasets
#'
#' Lists all ensembl dataset compatible with latest PANTHER database
#'
#' @return dataframe with ensembl datasets represented in PANTHER
#' @export
#'
#' @examples datasets = pantherCompatibleEnsemlDatasets
pantherCompatibleEnsemlDatasets = function(){
  if (is.null(pkgglobalenv$pantherCompatibleEnsemlDatasets)){
    ensembl = biomaRt::useMart(biomart = "plants_mart",
                               dataset = "athaliana_eg_gene",
                               host = "https://plants.ensembl.org")
    datasets <- biomaRt::listDatasets(ensembl)
    datasets$long_name = as.vector(sub("(\\w+\\s+\\w+).*", "\\1", datasets$description))
    panther_organisms = rbioapi::rba_panther_info(what = "organisms")
    panther_organisms.minimal = panther_organisms %>% dplyr::select(long_name,short_name,taxon_id)
    compatible_datasets = datasets[datasets$long_name %in% panther_organisms$long_name,]
    compatible_datasets = dplyr::left_join(compatible_datasets, panther_organisms.minimal,
                                           by = "long_name")
    pkgglobalenv$pantherCompatibleEnsemlDatasets = compatible_datasets
  }else(compatible_datasets = pkgglobalenv$pantherCompatibleEnsemlDatasets)
  return(compatible_datasets)
}

