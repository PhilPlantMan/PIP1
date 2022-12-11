test_that("Output directory paths are correct", {

  test_output_dir = withr::local_tempdir()
  test_dataset = "athaliana_eg_gene"
  test_directories = output_directories(test_output_dir, test_dataset)


  expect_equal(test_directories$parent, file.path(test_output_dir, paste0(Sys.Date(), "_", test_dataset)))
})
