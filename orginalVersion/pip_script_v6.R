packages <- c("plyr", "dplyr", "readxl", "stringr", "eulerr","remotes", 'tcltk', 'tools')
install.packages(setdiff(packages, rownames(installed.packages())))
if (!"ragp" %in% rownames(installed.packages())){remotes::install_github("missuse/ragp")}

if (!"biomaRt" %in% rownames(installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
  BiocManager::install("biomaRt")}
library(biomaRt)
library(plyr)
library(dplyr)
library(ragp)
library(readxl)
library(stringr)
library(eulerr)
library(tcltk)
library(tools)

choose_directory = function(caption = 'Select data directory') {
  if (Sys.info()['sysname'] == 'Windows') {
    choose.dir(caption = caption)
  } else {
    tk_choose.dir(caption = caption)
  }
}

print("Point to PIP directory")
pipeline_files_dir = choose_directory("Point to PIP directory")

print('Choose directory to save output')
wd_path = choose_directory("Choose directory to save output")
setwd(wd_path)

datasets = read.csv(file.path(pipeline_files_dir, "ensembl_panther_compat_datasets_info.csv"))


for (species in datasets$species){print(species)}

species_choice = readline(prompt="Which species? (don't include quotation marks) ")
if (!species_choice %in% datasets$species){stop("species not in list, try again")}
dataset_choice = datasets[datasets$species == species_choice, 'dataset.id']
dataset_description = datasets[datasets$species == species_choice, 'description']
minimum_num_proteomes = readline(prompt="Minimum number of proteomes for inclusion in list A and B (default = 2)? ")
minimum_num_proteomes = as.numeric(minimum_num_proteomes)
if (is.na(as.numeric(minimum_num_proteomes))){minimum_num_proteomes = 2}
fam_or_sub = readline("Would you like to pull and classify proteins based on family (f) or (default) subfamily (s)? ")
if (!fam_or_sub %in% c('f','s')){
  print("input not f or s, applying s as default")
  fam_or_sub = 's'}
if (fam_or_sub == 's'){
  sort_col = 'hmmpanther_subfams'
  sort_name = "subfamily"}
if (fam_or_sub == 'f'){
  sort_col = 'hmmpanther'
  sort_name = "family"}

to_cross_ref = list.files(file.path(pipeline_files_dir,'Cross reference with'), pattern = '.csv')
join_choice_vector = vector()
if (length(to_cross_ref)>0){
  for (file in to_cross_ref){
    print(paste("How would you like to cross reference with ",
                tools::file_path_sans_ext(file), "?", sep = ""))
    print("1: ensembl_gene_id")
    print("2: panther subfamily)(Panther annotations will be determied automatically)")
    print("3: panther family (Panther annotations will be determied automatically")
    join_choice = readline("Choose a number: ")
    if (!join_choice %in% c(1,2,3)){stop("Choose either 1, 2 or 3")}
    join_choice_vector = c(join_choice_vector,as.numeric(join_choice))
  }
}


zz=gzfile(file.path(pipeline_files_dir, "All_Plants.gene_info.gz"),'rt')
gene_info =read.csv(zz, header = TRUE, sep = "\t")
close(zz)
gene_info = gene_info[gene_info$X.tax_id == datasets[datasets$species == species_choice, 'Taxa.ID'],]
colnames(gene_info)[match('LocusTag',colnames(gene_info))] = 'ensembl_gene_id'

ensembl <- useMart(biomart = "plants_mart",
                   dataset = dataset_choice,
                   host = "https://plants.ensembl.org")
attributes = listAttributes(ensembl)


pd_proteome_filepaths = list.files(path = file.path(pipeline_files_dir, 'PD proteomes'), pattern = '.xlsx', all.files = FALSE,
                                   full.names = TRUE, recursive = FALSE,
                                   ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
verified_genes_df  = read_excel(file.path(pipeline_files_dir,'Known PD genes.xlsx'))


filter = 'secreted_or_gpi_or_tmhmm'
######################TRACKING LINE#########################

annotate_panther1 = function(df){
  panther16_classifications = read.csv(file.path(pipeline_files_dir, "PANTHER16_PD_pipeline_classifications.csv"))
  df.annotated = plyr::join(df,panther16_classifications, by = 'ensembl_gene_id')
  df.annotated = dplyr::select(df.annotated,
                               -Short.Name,
                               -Species)
  return(df.annotated)
}

annotate_gene_info = function(df){
  gene_info_condensed = gene_info[,c('Symbol', 'ensembl_gene_id', 'description', 'Synonyms')]
  df = plyr::join(df, gene_info_condensed, 'ensembl_gene_id')
}



pull_all = function(df, species_choice, sort_col = sort_col){
  panther16_classifications = read.csv(file.path(pipeline_files_dir, "PANTHER16_PD_pipeline_classifications.csv"))
  panther16_classifications = panther16_classifications[panther16_classifications$Species == species_choice,]
  panther16_classifications.matching <<- panther16_classifications[panther16_classifications[,sort_col]
                                                                   %in% df[,sort_col], ]
  df.to_return = dplyr::select(panther16_classifications.matching,
                               -Short.Name,
                               -Species)
  df.to_return <<- df.to_return[!duplicated(df.to_return$ensembl_gene_id),]
  return(df.to_return)
}

secreted_annotation <- function(df, sequence_col, peptide_id, override_NETGPI = FALSE,
                                override_phobius = FALSE,
                                override_bigGPI = FALSE){

  #get signalp predictions
  df.output <- df
  df.input <- df

  signalp <- eval(substitute(get_signalp(df.input,
                                         sequence_col,
                                         peptide_id)))

  df.output$signalp.secreted = signalp$is.signalp
  df.output$signalp.cutsite = signalp$Cmax.pos
  print("SignalP predicted")

  #get targetp predictions
  targetp <- eval(substitute(get_targetp(df.input,
                                         sequence_col,
                                         peptide_id)))
  df.output$targetp.Loc = targetp$Loc
  df.output$targetp.secreted = targetp$is.targetp
  print("TargetP predicted")

  #get phobius predictions
  if (override_phobius == TRUE){
    phobius = df.output
    phobius$is.phobius = FALSE
    phobius$cut_site = "ERROR, phobius NOT WORKING"
  }else{ phobius <- get_phobius(df.input,
                                peptide,
                                ensembl_peptide_id)}
  df.output$phobius.secreted = phobius$is.phobius
  df.output$phobius.cutsite = phobius$cut_site
  print("Phobius predicted")

  #makes consensus column
  df.output$at_least_one_prediction.secreted = (df.output$signalp.secreted +
                                                  df.output$targetp.secreted +
                                                  df.output$phobius.secreted)
  df.output$at_least_one_prediction.secreted = df.output$at_least_one_prediction.secreted >= 1

  # set.seed(5)
  # bind_cols(signalp, targetp, phobius) %>%
  #   select(is.phobius, is.targetp, is.signalp) %>%
  #   euler(shape = "ellipse", input = "disjoint") %>%
  #   plot(quantities = T)
  # knitr::include_graphics("eulerr.svg")

  # bind_cols(nsp_signalp, nsp_targetp, nsp_phobius) %>%
  #   select(is.phobius, is.targetp, is.signalp, id) %>%
  #   mutate(vote = rowSums(.[,1:3])) %>%
  #   filter(vote >= 2) %>%
  #   pull(id) -> id_nsp
  #
  # #get concencus of secretory predictions
  # at_nsp2 <- filter(at_nsp, Transcript.id %in% id_nsp)

  ##predict proline df.input
  # at_hyp <- eval(substitute(predict_hyp(data = df.input,
  #                        sequence = sequence_col,
  #                        id = peptide_id)))
  # stop()
  # df.output$hypdroxyp.AG_aa = at_hyp$AG_aa
  # df.output$hypdroxyp.total_length = at_hyp$total_length
  # df.output$hypdroxyp.longest = at_hyp$longest
  # print("hyp predicted")

  #predict Arabinogalactan identification
  scaned_at <- eval(substitute(scan_ag(data = df.input,
                                       sequence_col,
                                       peptide_id)))
  df.output$hypdroxyp.AG_aa = scaned_at$AG_aa
  df.output$hypdroxyp.total_length = scaned_at$total_length
  df.output$hypdroxyp.longest = scaned_at$longest
  print("Arabinogalactan predicted")
  #domain annotation
  # at_hmm <- get_hmm(at_nsp,
  #                   sequence,
  #                   Transcript.id)

  #gpi anchor prediction
  if (override_bigGPI == TRUE){
    at_big_gpi = df.output
    at_big_gpi$is.gpi = FALSE
    at_big_gpi$error = "ERROR, bigGPI.pred NOT WORKING"
  } else{ at_big_gpi <- eval(substitute(get_big_pi(df.input,
                                                   sequence_col,
                                                   peptide_id)))}
  df.output$bigGPI.pred = at_big_gpi$is.gpi
  print("bigGPI predicted")

  at_gpi_pred <- eval(substitute(get_pred_gpi(df.input,
                                              sequence_col,
                                              peptide_id)))
  df.output$predGPI.pred = at_gpi_pred$is.gpi
  print("predGPI predicted")

  if (override_NETGPI == FALSE){
    at_netGPI <- eval(substitute(get_netGPI(df.input,
                                            sequence_col,
                                            peptide_id)))}

  if (override_NETGPI == TRUE){
    at_netGPI = df.output
    at_netGPI$is.gpi = FALSE
    at_netGPI$omega_site = "ERROR, NETGPI NOT WORKING"
  }

  df.output$netGPI.pred = at_netGPI$is.gpi
  df.output$netGPI.omegasite = at_netGPI$omega_site
  print("netGPI predicted")

  df.output$at_least_one_prediction.gpi = (df.output$bigGPI.pred +
                                             df.output$predGPI.pred +
                                             df.output$netGPI.pred)
  df.output$at_least_one_prediction.gpi = df.output$at_least_one_prediction.gpi >= 1


  df.output$at_least_one_prediction.gpi_or_sp = (df.output$at_least_one_prediction.gpi +
                                                   df.output$at_least_one_prediction.secreted)
  df.output$at_least_one_prediction.gpi_or_sp = df.output$at_least_one_prediction.gpi_or_sp >= 1

  df.output = dplyr::select(df.output, -sequence_col)
  return(df.output)
}

predict_properties <- function(df, ensembl, peptide_id_provided = TRUE, override_NETGPI = FALSE){
  startwd = getwd()
  # setwd(find.package("customFunctions"))
  if (peptide_id_provided == FALSE){
    peptide_id<- getBM(attributes = c('ensembl_peptide_id',
                                      "ensembl_gene_id"),
                       filters = c("ensembl_gene_id"),
                       values = df$ensembl_gene_id,
                       mart = ensembl)
    df = plyr::join(df, peptide_id, by = 'ensembl_gene_id')
  }

  peptide_seqs<- getBM(attributes = c('ensembl_peptide_id',
                                      "peptide"),
                       filters = c("ensembl_peptide_id"),
                       values = df$ensembl_peptide_id,
                       mart = ensembl)

  peptide_seqs <- unique(peptide_seqs)
  # customFunc.path <- find.package("customFunctions", lib.loc = .libPaths())
  saved_annots <- read.csv(file.path(pipeline_files_dir, "predictions_cache.csv"))
  # usethis::use_data("saved_annots")
  retrieved_annots <- dplyr::filter(saved_annots, ensembl_peptide_id %in% peptide_seqs$ensembl_peptide_id)
  unretrieved_annots <- dplyr::filter(peptide_seqs, !ensembl_peptide_id %in% saved_annots$ensembl_peptide_id)
  unretrieved_annots <- unretrieved_annots[!is.na(unretrieved_annots$ensembl_peptide_id),]
  unretrieved_annots <- unretrieved_annots[unretrieved_annots$ensembl_peptide_id != "",]
  if (dim(unretrieved_annots)[1] > 0){

    new_secreted_annots <- batch_secreted_annot(unretrieved_annots, override_NETGPI = override_NETGPI)
    all_annots = rbind(new_secreted_annots, retrieved_annots)
  } else {all_annots = retrieved_annots}
  df = plyr::join(df, all_annots, by = 'ensembl_peptide_id')
  # write.csv(new_secreted_annots, "C:/Users/Phil/Google Drive/PhD/Sam and Phil Bioinformatics/R scripts/Phils stuff/Prediction storage/predictions_database.csv")
  #revert wd
  # setwd(startwd)
  return(df)}

batch_secreted_annot <- function(df, override_NETGPI = FALSE){
  numb_batches = ceiling(dim(df)[1] / 50)
  for (batch in seq(1, numb_batches)){
    index_start = ((batch - 1) * 50) +1
    if (batch == numb_batches){index_end = dim(df)[1] - ((numb_batches-1)*50) + index_start -1}
    else {index_end = batch * 50}
    df.batch = df[index_start:index_end,]
    df.batch_annotated <<- secreted_annotation(df.batch,
                                               sequence_col = 'peptide',
                                               peptide_id  = 'ensembl_peptide_id',
                                               override_NETGPI = override_NETGPI)
    head(df.batch_annotated)
    if (batch == 1){concat_df = df.batch_annotated}
    else {concat_df = rbind(concat_df, df.batch_annotated)}

    # customFunc.path <- find.package("customFunctions", lib.loc = .libPaths())
    saved_annots <- read.csv(file.path(pipeline_files_dir, "predictions_cache.csv"))
    updated_annots <- rbind(saved_annots, df.batch_annotated)

    write.csv(updated_annots, file.path(pipeline_files_dir, "predictions_cache.csv"), row.names = FALSE)

    print(paste("Batch", batch, "of", numb_batches, "completed"))
  }
  return(concat_df)
}

filter_candidates <- function(df_input, method, return_retentate = FALSE){
  if (!'at_least_one_prediction.gpi' %in% colnames(df_input)){
    df_input$at_least_one_prediction.gpi = df_input$GPI
    df_input$at_least_one_prediction.secreted = df_input$SP
    df_input$tmhmm  = df_input$TM
    remove_cols_at_end = TRUE
  }

  df <- df_input
  if (method == 'secreted'){
    df = df[df$at_least_one_prediction.secreted,]}

  if (method == 'secreted_gpi_or_tmhmm'){
    df$secreted_GPI = (df$at_least_one_prediction.secreted  + df$at_least_one_prediction.gpi)
    df$secreted_GPI = df$secreted_GPI == 2

    df$secreted_tmhmm = (df$at_least_one_prediction.secreted  + df$tmhmm)
    df$secreted_tmhmm = df$secreted_tmhmm == 2

    df$secreted_gpi_or_tmhmm = (df$secreted_GPI  + df$secreted_tmhmm)

    df$secreted_gpi_or_tmhmm  = df$secreted_gpi_or_tmhmm >= 1
    df = df[df$secreted_gpi_or_tmhmm,]}

  if (method == 'secreted_or_gpi_or_tmhmm'){
    df$at_least_one_feature = (df$at_least_one_prediction.secreted  + df$at_least_one_prediction.gpi + df$tmhmm)
    df = df[df$at_least_one_feature >= 1,]
    df = dplyr::select(df,
                       -at_least_one_feature)}


  if (method == 'tm_or_gpi'){
    df$TM_GPI <- (df$tmhmm  + df$at_least_one_prediction.gpi)
    df <- df[df$TM_GPI >= 1,]}

  if (method == 'gpi'){
    df$GPI <- df$at_least_one_prediction.gpi
    df <- df[df$GPI == 1,]}

  if (method == 'tm'){
    df$tm <- df$tmhmm
    df <- df[df$tm == 1,]}

  if (method == 'targetP_not_C_or_M'){
    df <- df[!df$targetp.Loc %in% c('C', 'M'),]}

  if (method == 'targetP_not_C_or_M_or_dash'){
    df <- df[!df$targetp.Loc %in% c('C', 'M', '_'),]}



  df <- df[!is.na(df$ensembl_gene_id),]
  if (remove_cols_at_end == TRUE){df = dplyr::select(df,
                                                     -at_least_one_prediction.gpi,
                                                     -at_least_one_prediction.secreted,
                                                     -tmhmm )}
  if (return_retentate == TRUE){
    # df_test <<-df
    # df_input_test <<- df_input
    retentate = df_input[!df_input$ensembl_gene_id %in% df$ensembl_gene_id,]
    if (remove_cols_at_end == TRUE){retentate = dplyr::select(retentate,
                                                              -at_least_one_prediction.gpi,
                                                              -at_least_one_prediction.secreted,
                                                              -tmhmm )}
    return(retentate)}else{return(df)}}


annotate_gene_ids <- function(df, ensembl){
  if (is.data.frame(df)){gene_ids = df$ensembl_gene_id}
  if (is.vector(df)){
    gene_ids = df
    df = as.data.frame(gene_ids)
    colnames(df) = 'ensembl_gene_id'
  }
  ensemb_attr<- getBM(attributes = c("ensembl_gene_id",
                                     'ensembl_peptide_id',
                                     'tmhmm'),
                      filters = c("ensembl_gene_id"),
                      values = gene_ids,
                      mart = ensembl)
  ensemb_attr$tmhmm = ensemb_attr$tmhmm == "TMhelix"
  df = plyr::join(df, ensemb_attr, 'ensembl_gene_id')
  # df = df[!duplicated(df$ensembl_gene_id),]
  return(df)
}

boil_down_predictions <- function(df, by_col = "ensembl_gene_id"){
  df[is.na(df)] <- FALSE
  boiled_down <- df %>%
    dplyr::group_by(get(by_col)) %>%
    dplyr::summarise(GPI = max(at_least_one_prediction.gpi),
                     SP = max(at_least_one_prediction.secreted),
                     TM = max(tmhmm),
                     TargetP_localisation = first(targetp.Loc))
  colnames(boiled_down)[1] = by_col
  return(boiled_down)
}



full_annotation_boil_down <- function(df, ensembl, override_NETGPI = FALSE, join_back = TRUE, reduce = TRUE){
  ensembl_attr <- annotate_gene_ids(df, ensembl)
  predictions <<- predict_properties(ensembl_attr, ensembl, override_NETGPI = override_NETGPI)
  boil_down <- boil_down_predictions(predictions)
  panther_attr <- annotate_panther1(boil_down)
  gene_info_anno <- annotate_gene_info(panther_attr)
  if (join_back == TRUE){gene_info_anno = plyr::join(df, gene_info_anno, by = 'ensembl_gene_id')}
  if(reduce == TRUE){gene_info_anno = gene_info_anno[!duplicated(gene_info_anno$ensembl_gene_id),]}
  gene_info_anno$Symbol[is.na(gene_info_anno$Symbol)] = gene_info_anno$ensembl_gene_id[is.na(gene_info_anno$Symbol)]
  return(gene_info_anno)
}

annotate_family_and_proteome_counts = function(df, known_annotated_PD_df){
  fams_or_subfams = read.csv(file.path(getwd(),"Counts cross referenced", paste(sort_col, "cross-referenced.csv") ))
  fams_or_subfams_counts = fams_or_subfams[,!str_sub(colnames(fams_or_subfams), 1,3) == "in."]
  fams_or_subfams_counts[is.na(fams_or_subfams_counts)] <- 0
  fams_or_subfams_proteome_counts = rowSums(fams_or_subfams_counts[,-1] != 0)
  sum_colname = paste(sort_name,'proteome_sum', sep = '_')
  counts_colname = paste(sort_name,'proteome_counts', sep = '_')
  fams_or_subfams_counts[,sum_colname] = rowSums(fams_or_subfams_counts[,-1])
  fams_or_subfams_counts[,counts_colname] = fams_or_subfams_proteome_counts
  df_output = plyr::join(df, fams_or_subfams_counts, by = sort_col)
  df_output$known_PD_gene = df_output$ensembl_gene_id %in% known_annotated_PD_df$ensembl_gene_id
  df_output$known_PD_subfam = df_output$hmmpanther_subfams %in% known_annotated_PD_df$hmmpanther_subfams
  df_output$known_PD_fam = df_output$hmmpanther %in% known_annotated_PD_df$hmmpanther

  panther16_classifications = read.csv(file.path(pipeline_files_dir, "PANTHER16_PD_pipeline_classifications.csv"))
  panther16_classifications = panther16_classifications[panther16_classifications$Species == species_choice,]
  all_counts = as.data.frame(table(panther16_classifications$hmmpanther))
  colnames(all_counts) = c('hmmpanther', 'family_size')
  df_output = plyr::join(df_output, all_counts, by = 'hmmpanther')

  all_counts = as.data.frame(table(panther16_classifications$hmmpanther_subfams))
  colnames(all_counts) = c('hmmpanther_subfams', 'subfamily_size')
  df_output = plyr::join(df_output, all_counts, by = 'hmmpanther_subfams')
  return(df_output)}

sort_rows = function(df, method = 'method1'){
  df$feature_count = rowSums(df[,c("GPI", "TM", "SP")])
  if ('family_proteome_counts' %in% colnames(df)){
    df = df[with(df, order(-known_PD_gene,
                           -known_PD_fam,
                           -known_PD_fam,
                           -feature_count,
                           -family_proteome_counts,
                           -family_proteome_sum)),]}
  if ('subfamily_proteome_counts' %in% colnames(df)){
    df = df[with(df, order(-known_PD_gene,
                           -known_PD_fam,
                           -known_PD_fam,
                           -feature_count,
                           -subfamily_proteome_counts,
                           -subfamily_proteome_sum)),]}
  return(df)
}

sort_cols = function(df){
  col_order = c('ensembl_gene_id',
                'Symbol',
                'Synonyms',
                'description',
                'hmmpanther',
                'family_size',
                'hmmpanther_subfams',
                'subfamily_size',
                'FAMILY_TERM',
                'SUBFAMILY_TERM',
                'GPI',
                'SP',
                'TM',
                'TargetP_localisation',
                'feature_count',
                'known_PD_gene',
                'known_PD_subfam',
                'known_PD_fam',
                'family_proteome_counts',
                'family_proteome_sum',
                'subfamily_proteome_counts',
                'subfamily_proteome_sum')
  col_orders_not_in_df = col_order[!col_order %in% colnames(df)]
  df_cols_not_in_order = colnames(df)[!colnames(df) %in% col_order]
  # for (col in col_orders_not_in_df){
  #   print(paste(col, "not in df"))
  # }
  # for (col in df_cols_not_in_order){
  #   print(paste(col, "not in df"))}
  valid__col_order = col_order[col_order %in% colnames(df)]
  valid__col_order_with_missing = c(valid__col_order,df_cols_not_in_order)
  df = df[,valid__col_order_with_missing]
  return(df)
}

#Family and family counting for each proteome

dir.create(file.path(getwd(),paste(sort_name,"counts")))

for (file in pd_proteome_filepaths){
  print(file)
  file_df = read_excel(file)
  file_suffix = gsub(".xls.*","",basename(file))
  df.with_families <- annotate_panther1(file_df)
  df.with_families <- df.with_families[!is.na(df.with_families[,sort_col]), ]
  df.fam_counts = df.with_families %>% dplyr::count(!!sym(sort_col), name = paste(sort_name,"counts",sep = "_"))
  write.csv(df.fam_counts, file.path(paste(sort_name,"counts"),
                                     paste(file_suffix,".csv", sep = "")), row.names = FALSE)
}



#Cross reference by family
cross_ref_sort_col = function(dir_path, reference_column_name = sort_col){
  file_list = list.files(path = dir_path, pattern = '.csv', all.files = FALSE,
                         full.names = FALSE, recursive = FALSE,
                         ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  for (file in file_list){
    #print(file)
    file_df = read.csv(file.path(dir_path,file))
    file_suffix = gsub(".csv*","",file)
    reference_col_pos = match(reference_column_name, colnames(file_df))
    colnames(file_df)[-reference_col_pos] = paste(file_suffix, '.',
                                                  colnames(file_df)[-reference_col_pos],
                                                  sep = '')
    #print(file_suffix)
    file_df[[paste('in.', file_suffix, sep = '')]] = TRUE
    if (file == file_list[1]){concat_df = file_df}
    else {concat_df = plyr::join(concat_df, file_df, by = reference_column_name, type = 'full')}
  }
  concat_df[is.na(concat_df)] <- ''
  concat_df = unique(concat_df)
  #summary_record['starting_proteomes', reference_column_name] <<- length(concat_df[,reference_column_name])
  write.csv(concat_df, file = file.path(getwd(),
                                        'Counts cross referenced',
                                        paste(reference_column_name,'cross-referenced.csv')),
            row.names = FALSE)
}


cross_ref_generic = function(df, join_choice_vector){
  file_list = list.files(path = file.path(pipeline_files_dir,'Cross reference with'), pattern = c('.csv'), all.files = FALSE,
                         full.names = FALSE, recursive = FALSE,
                         ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  by_cols = c('ensembl_gene_id','hmmpanther_subfams','hmmpanther')
  for (file_position in seq(file_list)){
    file_sans_ext = tools::file_path_sans_ext(file_list[file_position])
    ref_df = read.csv(file.path(pipeline_files_dir,'Cross reference with',file_list[file_position]))
    join_choice = by_cols[join_choice_vector[file_position]]
    if (!join_choice == 'ensembl_gene_id'){
      ref_df_annotated = annotate_panther1(ref_df)
      ref_df = data.frame(choice = ref_df_annotated[,join_choice])
      colnames(ref_df) = join_choice}
    ref_df[,paste(join_choice,'in',file_sans_ext, sep = ".")] = TRUE
    colnames_seq = seq(colnames(ref_df))
    last_element = colnames_seq[length(colnames_seq)]
    remaining_elements = colnames_seq[1:length(colnames_seq)-1]
    reorder_seq = c(last_element,remaining_elements)
    ref_df = ref_df[,reorder_seq]

    df = dplyr::left_join(df, ref_df, by = join_choice, suffix = c("",""))
    df[is.na(df[,paste(join_choice,'in',file_sans_ext, sep = ".")]),
       paste(join_choice,'in',file_sans_ext, sep = ".")] = FALSE
  }
  return(df)
}








dir.create(file.path(getwd(),"Counts cross referenced"))
# cross_ref_families(file.path(getwd(),"Subfamily counts"),'hmmpanther_subfams')

cross_ref_sort_col(file.path(getwd(),paste(sort_name, "counts")),sort_col)

fams_or_Subfams = read.csv(file.path(getwd(),"Counts cross referenced", paste(sort_col, "cross-referenced.csv")))
fams_or_Subfams_bool_cols = colnames(fams_or_Subfams)[str_sub(colnames(fams_or_Subfams), 1,3) == "in."]
fams_or_Subfam_bools = fams_or_Subfams[,fams_or_Subfams_bool_cols]
fams_or_Subfam_bools[is.na(fams_or_Subfam_bools)] <- FALSE

fit <- euler(fams_or_Subfam_bools, shape = "ellipse")

png(file = file.path(getwd(),"Counts cross referenced",paste(sort_name, 'euler.png')),height=500,width=1000)
print(plot(fit, quantities =  list(type = c("counts")),
           fill = "transparent",
           edges = list(col = "black", lex = 2))) #,labels = rep('', dim(df.bool_cols)[2])
dev.off()

png(file = file.path(getwd(),"Counts cross referenced",paste(sort_name, ' no-label-euler.png')),height=500,width=1000)
print(plot(fit, quantities =  list(type = c("counts")),
           fill = "transparent",
           edges = list(col = "black", lex = 2),
           labels = rep('', dim(fams_or_Subfam_bools)[2]))) #,labels = rep('', dim(df.bool_cols)[2])
dev.off()

png(file = file.path(getwd(),"Counts cross referenced",paste(sort_name, ' venn.png')),height=500,width=700)
print(plot(venn(fams_or_Subfam_bools),
           edges = list(col = "black", lex = 2),
           environment = environment())) # quantities =  list(type = c("percent", "counts"))
dev.off()

png(file = file.path(getwd(),"Counts cross referenced",paste(sort_name, ' venn no label.png')),height=500,width=700)
print(plot(venn(fams_or_Subfam_bools),
           edges = list(col = "black", lex = 2),
           environment = environment(),
           labels = rep('', dim(fams_or_Subfam_bools)[2]))) # quantities =  list(type = c("percent", "counts"))
dev.off()
# Make list of all known PD family members for species
dir.create(file.path(getwd(),paste("Verified",sort_name, "members")))


verified_genes_df.annotated = annotate_panther1(verified_genes_df)
verified.all_fams_members = pull_all(verified_genes_df.annotated,species_choice, sort_col = sort_col)
verified.all_fams_members_BD = full_annotation_boil_down(data.frame(ensembl_gene_id=verified.all_fams_members$ensembl_gene_id), ensembl)
verified.all_fams_members_BD = annotate_family_and_proteome_counts(verified.all_fams_members_BD, verified_genes_df.annotated)
verified.all_fams_members_BD = sort_rows(verified.all_fams_members_BD)
verified.all_fams_members_BD = sort_cols(verified.all_fams_members_BD)
verified.all_fams_members_BD = cross_ref_generic(verified.all_fams_members_BD, join_choice_vector)
write.csv(verified.all_fams_members_BD,
          file.path(getwd(), paste("Verified",sort_name, "members"),
                    paste(species_choice,"Verified",sort_name, "members.csv")),
          row.names = FALSE)


fams_or_Subfams$rowsum = rowSums(fams_or_Subfams[,fams_or_Subfams_bool_cols], na.rm = TRUE)

fams_or_Subfams.multiple_prots = fams_or_Subfams[fams_or_Subfams$rowsum >= minimum_num_proteomes,]


fams_or_Subfams.not_in_multiple_prots = fams_or_Subfams[fams_or_Subfams$rowsum < minimum_num_proteomes,]




genes_in_multiple_proteome_fams_or_Subfams =  pull_all(fams_or_Subfams.multiple_prots,species_choice, sort_col)
genes_not_in_multiple_proteome_fams_or_Subfams =  pull_all(fams_or_Subfams.not_in_multiple_prots,species_choice, sort_col)


genes_in_multiple_proteome_fams_or_Subfams_BD = full_annotation_boil_down(data.frame(ensembl_gene_id=genes_in_multiple_proteome_fams_or_Subfams$ensembl_gene_id), ensembl)
genes_not_in_multiple_proteome_fams_or_Subfams_BD = full_annotation_boil_down(data.frame(ensembl_gene_id=genes_not_in_multiple_proteome_fams_or_Subfams$ensembl_gene_id), ensembl)

candidates_a = filter_candidates(genes_in_multiple_proteome_fams_or_Subfams_BD, filter,return_retentate = FALSE)
candidates_b = filter_candidates(genes_in_multiple_proteome_fams_or_Subfams_BD, filter,return_retentate = TRUE)
candidates_c = filter_candidates(genes_not_in_multiple_proteome_fams_or_Subfams_BD, filter,return_retentate = FALSE)
candidates_d = filter_candidates(genes_not_in_multiple_proteome_fams_or_Subfams_BD, filter,return_retentate = TRUE)

candidates_a = annotate_family_and_proteome_counts(candidates_a, verified_genes_df.annotated)
candidates_b = annotate_family_and_proteome_counts(candidates_b, verified_genes_df.annotated)
candidates_c = annotate_family_and_proteome_counts(candidates_c,verified_genes_df.annotated)
candidates_d = annotate_family_and_proteome_counts(candidates_d,verified_genes_df.annotated)

candidates_a = sort_rows(candidates_a)
candidates_b = sort_rows(candidates_b)
candidates_c = sort_rows(candidates_c)
candidates_d = sort_rows(candidates_d)

candidates_a = sort_cols(candidates_a)
candidates_b = sort_cols(candidates_b)
candidates_c = sort_cols(candidates_c)
candidates_d = sort_cols(candidates_d)

candidates_a = cross_ref_generic(candidates_a, join_choice_vector)
candidates_b = cross_ref_generic(candidates_b, join_choice_vector)
candidates_c = cross_ref_generic(candidates_c, join_choice_vector)
candidates_d = cross_ref_generic(candidates_d, join_choice_vector)

dir.create(file.path(getwd(),"Candidate PD gene lists"))
write.csv(candidates_a, file.path("Candidate PD gene lists",
                                  paste(species_choice,"candidate_list_a.csv")), row.names = FALSE)
write.csv(candidates_b, file.path("Candidate PD gene lists",
                                  paste(species_choice,"candidate_list_b.csv")), row.names = FALSE)
write.csv(candidates_c, file.path("Candidate PD gene lists",
                                  paste(species_choice,"candidate_list_c.csv")), row.names = FALSE)
write.csv(candidates_d, file.path("Candidate PD gene lists",
                                  paste(species_choice,"candidate_list_d.csv")), row.names = FALSE)
output_textfile_name = paste("Based_on_PANTHER_",sort_name,".txt", sep = "")
output_string = paste("Candidate genes are generated and classified according to", sort_name, "identity. All genes in the genome with the identities present in the proteomes are pulled and partially classified by how many proteomes contain genes with this identity. Please see the manuscript for more details.")

fileConn<-file(file.path("Candidate PD gene lists",output_textfile_name))
writeLines(output_string, fileConn)
close(fileConn)
print('DONE')
