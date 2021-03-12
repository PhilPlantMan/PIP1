# PIP (Plasmodesmata *in silico* Proteome)
`PIP` is a R-based tool for generating and classifying candidate genes encoding plasmodesmata (PD) proteins based on proteomics and predictions of structural features. See [publication] for more information  
Compatible species with this pipeline are those listed in both Ensemble Plants and PANTHER16 databases. See Supplemental table 2 of the associated publication for a list of compatible plant species  
### Prerequisites
This pipeline has only been tested in Windows 10. There might be issues with other operating systems.  
This pipeline was written with R version 4.0.2.  
All dependencies should install automatically.  
If you are unfamiliar with R, please install [RStudio](https://rstudio.com/products/rstudio/download/).  
## Getting started
Download and unzip PIP directory.zip
#### Inside PIP directory are:  
- PD proteomes - folder with PD experimental proteomes. User editable  
- All_Plants.gene_info.gz - data required for pipeline function  
- ensembl_panther_compat_datasets_info.csv - data required for pipeline function  
- Known PD genes.xlsx   - user editable spreadsheet with known PD genes in Arabidopsis thaliana. Genes encoding PD proteins from other species might be added to this list.  
- PANTHER16_PD_pipeline_classifications.csv - data required for pipeline function  
- pip_script.R - R script to perform PIP pipeline  
- predictions_cache.csv - data required for pipeline function  
#### Running the pipeline
1.	Add or remove any proteomes you desire in the PD proteomes folder. Compatible gene IDs only (see ensembl_panther_compat_datasets_info.csv inside the PIP directory for compatible reference genome versions)  
2.	Add or amend any genes in ‘Known PD genes.xlsx’.  Compatible gene IDs only (as above)  
3.	Load the pip_script.R into RStudio (or from the terminal)  
4.	If running in RStudio, run the script in Echo mode (Windows: Ctrl+Shift+Enter)  
5.	Accept or respond to any prompts to install packages  
6.	A popup will appear asking you where the PIP directory is saved on your computer  
7.	A popup will appear asking you where you would like to save the output.  
8.	The console will list compatible species. You will be asked in the console to select a specie.
9.	The console will ask you ‘Minimum number of proteomes for inclusion in list A and B (default = 2)?’. See Figure 1 in the manuscript for more details.  
10.	The script will now run. This may take several hours. The slowest stage is prediction of protein features (SP, GPI and TM, see the associated publication for more details). Predictions are cached; this saves time if PIP is re-ran.   
#### PIP output
On completion, in the folder that you set as the output destination (step 5 above) you will have the following folders:  
- Candidate PD gene lists - Folder with the 4 lists of PD candidates (A-D) (See figure 1 of the associated publication)  
- Subfamily counts - Folder with the number of appearances of each PANTHER subfamily present in each PD proteome  
- Counts cross referenced - A concatenated version of above with Venn and Euler diagrams comparing PANTHER subfamilies represented in each proteome  
- Verified subfamily members - Orthologues (based on PANTHER subfamily membership) of known PD-localising genes (listed in ‘Known PD genes.xlsx’) in the species selected.   
## Troubleshooting
The most common errors are associated with connection problems (e.g. any error that mentions ‘curl’). This is often caused by problems with your internet connection or with servers that host the 3rd party prediction tools. Try again later and at different times of the day. Predictions are stored in case this happens partially through a run.
If problems are encountered early in running PIP, it is likely that your R and package versions do not align with the scripts. Working session information below may help you diagnose your problem.

[Session info]  
R version 4.0.2 (2020-06-22)  
Platform: x86_64-w64-mingw32/x64 (64-bit)  
Running under: Windows 10 x64 (build 19042)  
Attached packages:  
eulerr_6.1.0, stringr_1.4.0, readxl_1.3.1, ragp_0.3.2.0002 dplyr_1.0.5, plyr_1.8.6, biomaRt_2.44.4   
loaded via a namespace (and not attached):  
Rcpp_1.0.5, lattice_0.20-41, prettyunits_1.1.1, assertthat_0.2.1, utf8_1.1.4, mime_0.10, BiocFileCache_1.12.1, R6_2.5.0, cellranger_1.1.0, stats4_4.0.2, RSQLite_2.2.3, httr_1.4.2, pillar_1.5.1         rlang_0.4.10, progress_1.2.2, curl_4.3, rstudioapi_0.13, data.table_1.14.0, blob_1.2.1, S4Vectors_0.26.1, Matrix_1.2-18, polyclip_1.10-0, bit_4.0.4, compiler_4.0.2, polylabelr_0.2.0, pkgconfig_2.0.3, askpass_1.1,BiocGenerics_0.34.0, openssl_1.4.3, tidyselect_1.1.0, tibble_3.0.3, IRanges_2.22.2, XML_3.99-0.5, fansi_0.4.2, crayon_1.4.1, dbplyr_2.1.0,withr_2.4.1,MASS_7.3-53.1, rappdirs_0.3.3, grid_4.0.2, lifecycle_1.0.0, DBI_1.1.1, magrittr_2.0.1, cli_2.3.1, stringi_1.5.3, cachem_1.0.4, remotes_2.2.0, xml2_1.3.2, seqinr_4.2-5, ellipsis_0.3.1, generics_0.1.0, vctrs_0.3.6, xgboost_1.3.2.1, tools_4.0.2, ade4_1.7-16, bit64_4.0.5, Biobase_2.48.0, glue_1.4.2, purrr_0.3.4, hms_1.0.0, parallel_4.0.2, fastmap_1.1.0, AnnotationDbi_1.50.3, BiocManager_1.30.10  memoise_2.0.0       


## Built With

* [biomartr](https://github.com/ropensci/biomartr) - Peptide sequence retrieval
* [ragp](https://github.com/missuse/ragp) - Protein feature prediction

## Bug reports
If you encounter problems with `PIP`, open an [issue](https://github.com/PhilPlantMan/PIP/issues)

## Publication

**********
