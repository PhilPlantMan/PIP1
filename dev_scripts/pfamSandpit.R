packages <- c("plyr", "dplyr", "readxl", "stringr", "eulerr","remotes", 'tcltk', 'tools')
install.packages(setdiff(packages, rownames(installed.packages())))
if (!"ragp" %in% rownames(installed.packages())){remotes::install_github("missuse/ragp")}

if (!"biomaRt" %in% rownames(installed.packages())){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
  BiocManager::install("biomaRt")}
if (!"PFAM.db" %in% rownames(installed.packages())){
  BiocManager::install("PFAM.db")}
library(biomaRt)
library(plyr)
library(dplyr)
library(ragp)
library(readxl)
library(stringr)
library(eulerr)
library(tcltk)
library(tools)
library(PFAM.db)


setwd("F:/Cleaned whole PD pipeline/pip_directory_v6/221204 thesis corrections/pfam annotation sandpit")

athalPipA = read.csv("Arabidopsis thaliana candidate_list_a.csv")

ensemblMoss <- useMart(biomart = "plants_mart",
                   dataset = "ppatens_eg_gene",
                   host = "https://plants.ensembl.org")
ensemblArab <- biomaRt::useMart(biomart = "plants_mart",
                   dataset = "athaliana_eg_gene",
                   host = "https://plants.ensembl.org")

phytozome_v13 <- useMart(biomart = "phytozome_mart",
                         dataset = "phytozome",
                         host = "https://phytozome-next.jgi.doe.gov")


datasets <- listDatasets(phytozome_v13)
attributesPhyto = listAttributes(phytozome_v13)
attributesArab = biomaRt::listAttributes(ensemblArab)
ArabFilters = biomaRt::listFilters(ensemblArab)
attributesMoss = listAttributes(ensemblMoss)

pfam<- getBM(attributes = c('ensembl_gene_id',
                                  "pfam"),
                   filters = c("ensembl_gene_id"),
                   values = athalPipA$ensembl_gene_id,
                   mart = ensemblArab)

AC2DE <- as.list(PFAMDE)
pfamID = as.list(PFAMID)
pfamLoad = as.list(PFAMLOAD)
ewnsembls = listEnsembl()

phtytoFilters = listFilters(phytozome_v13)
phtyo<- getBM(attributes = c('organism_name',"organism_tla"
                            ),
             #filters = c("ensembl_gene_id"),
             #values = athalPipA$ensembl_gene_id,
             mart = phytozome_v13)

arab_genes = getBM(attributes = c("gene_name1","panther_id","panther_desc"),
      filters = "organism_id",
      values = "447 ",
      mart = phytozome_v13)


ensemblePantherARab<- getBM(attributes = c('ensembl_gene_id',
                            "hmmpanther"),
             filters = c("ensembl_gene_id"),
             values = athalPipA$ensembl_gene_id,
             mart = ensemblArab)

ensemblePantherMoss<- getBM(attributes = c('ensembl_gene_id',
                                           "hmmpanther"),
                            #filters = c("ensembl_gene_id"),
                            #values = athalPipA$ensembl_gene_id,
                            mart = ensemblMoss)
##################Get panther descriptions########################
library(curl)
url = "ftp://ftp.pantherdb.org/hmm_classifications/17.0/PANTHER17.0_HMM_classifications"
tmp <- tempfile()
curl_download(url, tmp)
PantherClassifications = read.table(tmp,fill = TRUE, sep = "\t",
                                    comment.char = "",quote = "",
                                    col.names = c("PANTHERID",
                                                  "DESCRIPTION",
                                                  "Molecular_Function",
                                                  "Biological_Process",
                                                  "Cellular_Component",
                                                  "PANTHER_Protein_Class",
                                                  "Pathways"))
linearProteinClass.string = paste(PantherClassifications$PANTHER_Protein_Class,";")
linearProteinClass.split = strsplit(linearProteinClass.string,";")
uniqueProteinClassifications = unique(linearProteinClass.split)
