boil_down_predictions <- function(df, by_col = "ensembl_gene_id"){
  #df[is.na(df)] <- FALSE
  boiled_down <- df %>%
    dplyr::group_by(get(!!by_col)) %>%
    dplyr::summarise(GPI = max(at_least_one_prediction.gpi),
                     SP = max(at_least_one_prediction.secreted),
                     TM = max(tmhmm),
                     TargetP_localisation = dplyr::first(targetp.Loc)) %>%
    dplyr::rename(!!by_col := 1)
  return(boiled_down)
}

#TODO this is a legacy function that could be better integrated
annotate_gene_ids <- function(df, ensembl){
  if (is.data.frame(df)){gene_ids = df$ensembl_gene_id}
  if (is.vector(df)){
    gene_ids = df
    df = as.data.frame(gene_ids)
    colnames(df) = 'ensembl_gene_id'
  }
  ensemb_attr<- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                     'ensembl_peptide_id',
                                     'tmhmm'),
                      filters = c("ensembl_gene_id"),
                      values = gene_ids,
                      mart = ensembl)
  ensemb_attr$tmhmm = ensemb_attr$tmhmm == "TMhelix"
  df = dplyr::left_join(df, ensemb_attr, by = 'ensembl_gene_id')
  # df = df[!duplicated(df$ensembl_gene_id),]
  return(df)
}

predict_properties <- function(df, ensembl, peptide_id_provided = TRUE, override_NETGPI = FALSE){

  # setwd(find.package("customFunctions"))
  if (peptide_id_provided == FALSE){
    peptide_id<- biomaRt::getBM(attributes = c('ensembl_peptide_id',
                                      "ensembl_gene_id"),
                       filters = c("ensembl_gene_id"),
                       values = df$ensembl_gene_id,
                       mart = ensembl)
    df = dplyr::left_join(df, peptide_id, by = 'ensembl_gene_id')
  }

  peptide_seqs<- biomaRt::getBM(attributes = c('ensembl_peptide_id',
                                      "peptide"),
                       filters = c("ensembl_peptide_id"),
                       values = df$ensembl_peptide_id,
                       mart = ensembl)

  peptide_seqs <- unique(peptide_seqs)
  # customFunc.path <- find.package("customFunctions", lib.loc = .libPaths())
  saved_annots <- read_spreadsheet(get_predictions_cache_path())
  # usethis::use_data("saved_annots")
  retrieved_annots <- dplyr::filter(saved_annots, ensembl_peptide_id %in% peptide_seqs$ensembl_peptide_id)
  unretrieved_annots <- dplyr::filter(peptide_seqs, !ensembl_peptide_id %in% saved_annots$ensembl_peptide_id)
  unretrieved_annots <- unretrieved_annots[!is.na(unretrieved_annots$ensembl_peptide_id),]
  unretrieved_annots <- unretrieved_annots[unretrieved_annots$ensembl_peptide_id != "",]
  if (dim(unretrieved_annots)[1] > 0){

    new_secreted_annots <- batch_secreted_annot(unretrieved_annots, override_NETGPI = override_NETGPI)
    all_annots = rbind(new_secreted_annots, retrieved_annots)
  } else {all_annots = retrieved_annots}
  df = dplyr::left_join(df, all_annots, by = 'ensembl_peptide_id')
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
    df.batch_annotated <- secreted_annotation(df.batch,
                                               sequence_col = 'peptide',
                                               peptide_id  = 'ensembl_peptide_id',
                                               override_NETGPI = override_NETGPI)
    utils::head(df.batch_annotated)
    if (batch == 1){concat_df = df.batch_annotated}
    else {concat_df = rbind(concat_df, df.batch_annotated)}

    saved_annots <- utils::read.csv(get_predictions_cache_path())
    updated_annots <- rbind(saved_annots, df.batch_annotated)

    utils::write.csv(updated_annots, get_predictions_cache_path(), row.names = FALSE)

    print(paste("Batch", batch, "of", numb_batches, "completed"))
  }
  return(concat_df)
}

secreted_annotation <- function(df, sequence_col, peptide_id, override_NETGPI = FALSE,
                                override_phobius = FALSE,
                                override_bigGPI = FALSE){

  #get signalp predictions
  df.output <- df
  df.input <- df

  signalp <- eval(substitute(ragp::get_signalp(df.input,
                                         sequence_col,
                                         peptide_id)))

  df.output$signalp.secreted = signalp$is.signalp
  df.output$signalp.cutsite = signalp$Cmax.pos
  print("SignalP predicted")

  #get targetp predictions
  targetp <- eval(substitute(ragp::get_targetp(df.input,
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
  }else{ phobius <- ragp::get_phobius(df.input,
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


  scaned_at <- eval(substitute(ragp::scan_ag(data = df.input,
                                       sequence_col,
                                       peptide_id)))
  df.output$hypdroxyp.AG_aa = scaned_at$AG_aa
  df.output$hypdroxyp.total_length = scaned_at$total_length
  df.output$hypdroxyp.longest = scaned_at$longest
  print("Arabinogalactan predicted")
  #domain annotation


  #gpi anchor prediction
  if (override_bigGPI == TRUE){
    at_big_gpi = df.output
    at_big_gpi$is.gpi = FALSE
    at_big_gpi$error = "ERROR, bigGPI.pred NOT WORKING"
  } else{ at_big_gpi <- eval(substitute(ragp::get_big_pi(df.input,
                                                   sequence_col,
                                                   peptide_id)))}
  df.output$bigGPI.pred = at_big_gpi$is.gpi
  print("bigGPI predicted")

  at_gpi_pred <- eval(substitute(ragp::get_pred_gpi(df.input,
                                              sequence_col,
                                              peptide_id)))
  df.output$predGPI.pred = at_gpi_pred$is.gpi
  print("predGPI predicted")

  if (override_NETGPI == FALSE){
    at_netGPI <- eval(substitute(ragp::get_netGPI(df.input,
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
