test_that("pantherCompatibleEnsemlDatasets all have hmmpanther attribute", {

  datasets = pantherCompatibleEnsemlDatasets()
  results = vector()

  for(dataset in datasets$dataset){
    ensembl = biomaRt::useMart(biomart = "plants_mart",
                      dataset = dataset,
                      host = "https://plants.ensembl.org")
    attributes = biomaRt::listAttributes(ensembl)
    results = c(results, "hmmpanther" %in% attributes$name)
    print(paste(dataset,"- hmmpanther present = ",  "hmmpanther" %in% attributes$name))
  }

  expect_equal(length(results[results == TRUE]), length(datasets$dataset))
})
