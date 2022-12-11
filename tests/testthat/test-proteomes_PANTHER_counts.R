test_that("Counts are performed on internal proteomes", {

  test_output_dir = withr::local_tempdir()
  test_dataset = "athaliana_eg_gene"
  test_directories = output_directories(test_output_dir, test_dataset)

  internal_PD_proteome_dir = file.path(system.file("extdata", package = "pip1"),"PD_proteomes")
  proteome_counts = proteomes_PANTHER_counts(internal_PD_proteome_dir, output_dir = test_directories$counts)

  single_proteome_count_path = file.path(test_directories$counts,"Brault-2018-PD-proteome.csv")
  expect_equal(file.exists(single_proteome_count_path), TRUE)

  cross_referenced_path = file.path(test_directories$counts,"hmmpanther_subfams_cross-referenced.csv")
  expect_equal(file.exists(cross_referenced_path), TRUE)

  export_Euler_diagrams(proteome_counts,output_dir = test_directories$counts)
  first_euler_path = file.path(test_directories$counts,"subfamily_venn.png")
  expect_equal(file.exists(first_euler_path), TRUE)

  pulled_all_subfam_members = pull_PANTHER_members(proteome_counts)
  expect_gt(nrow(pulled_all_subfam_members), nrow(proteome_counts))
  expect_equal(ncol(pulled_all_subfam_members), ncol(proteome_counts) + 1)
})

