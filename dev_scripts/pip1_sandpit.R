library(devtools)

known_PD_genes = read_spreadsheet(file.path(system.file("extdata", package = "pip1"),"Known_PD_genes.xlsx"))

pip_output = run_PIP1()
compatible_datasets = pantherCompatibleEnsemlDatasets()

pulled_PD_PATHERS= pip_output$pulled_PD_PANTHER_members
pulled_candidates = pip_output$pulled_candidates

full_annotation_boil_down_candidates = pip_output$pulled_candidates_annotated
