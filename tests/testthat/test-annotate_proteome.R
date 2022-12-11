test_that("PD proteomes are annotated with family/subfamily correctly", {
  internal_PD_proteome_dir = get_internal_PD_proteome_dir()
  Brault_proteome = file.path(internal_PD_proteome_dir,"Brault-2018-PD-proteome.xlsx")
  Brault_annotated = annotate_proteome(Brault_proteome)
  hmmpanther_subfam_entry = Brault_annotated[Brault_annotated$ensembl_gene_id == "AT1G51570", "hmmpanther_subfams"]
  expect_equal(hmmpanther_subfam_entry$hmmpanther_subfams, "PTHR45707:SF52")
})
