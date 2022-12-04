test_that("panther classification returns values", {

  results = get_PANTHER_classifications()

  expect_equal(colnames(results) , c("PANTHERID",
                                     "DESCRIPTION",
                                     "Molecular_Function",
                                     "Biological_Process",
                                     "Cellular_Component",
                                     "PANTHER_Protein_Class",
                                     "Pathways") )
  expect_gt(nrow(results), 1000)
})
