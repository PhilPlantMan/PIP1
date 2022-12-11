library(biomaRt)
library(rbioapi)
library(stringr)
library(dplyr)
library(stringr)

#arbidposis athaliana_eg_gene
#ppatens ppatens_eg_gene


userListElementSelection  = function(elementList, instructionString = ""){
  print(instructionString)
  for (i in seq(elementList)){
    print(paste0(i, ": ", elementList[i]))
  }
  selectionIndex= readline(paste0(instructionString," ","\n[Select list element by number]"))
  selectionIndex = as.integer(selectionIndex)
  if (selectionIndex > length(elementList)){
    print("Number is larger that length of list, please try again")
    userListElementSelection(elementList, instructionString)
  }
  return(elementList[selectionIndex])
}


ensembl <- useMart(biomart = "plants_mart",
                   dataset = "ppatens_eg_gene",
                   host = "https://plants.ensembl.org")
attributes = listAttributes(ensembl)
datasets <- listDatasets(ensembl)
datasets$Species = as.vector(sub("(\\w+\\s+\\w+).*", "\\1", datasets$description))

peptide_id<- getBM(attributes = c('ensembl_peptide_id',
                                  "ensembl_gene_id"),
                   #filters = c("ensembl_gene_id"),
                   #values = df$ensembl_gene_id,
                   mart = ensembl)


## 1 We get the available annotation datasets in PANTHER (we need to select one of them to submit an enrichment request)
annots <- rba_panther_info(what = "datasets")
organisms = rba_panther_info(what = "organisms")
family_list = rba_panther_info(what = "families")
family_pages = family_list$pages_count

for (page in seq(family_pages+1)){
  page_numeric = as.numeric(page)
  family_temp = rba_panther_info(what = "families", families_page = page_numeric)$familiy
  if (page == 1){familys = family_temp
  }else{familys = rbind(familys, family_temp)}
}

species_tree = rba_panther_info(what = "species_tree")

#ensemble and panther overlapping species
ensemble_species = as.vector(sub("(\\w+\\s+\\w+).*", "\\1", datasets$description))
inBoth = organisms[organisms$long_name %in% ensemble_species,"long_name"]

#get Taxon from species names
species = "Physcomitrella patens"
taxonId = as.numeric(organisms[organisms$long_name == species, "taxon_id"])

#get annotations from panther
panther_classifications_df = data.frame()
panther_annotations = rba_panther_mapping(tail(peptide_id$ensembl_gene_id),taxonId)
genes = panther_annotations$mapped_genes$gene
for (i in seq(length(genes))){
  panther_accession = genes[[i]]$accession
  panther_accession_split = unlist(stringr::str_split(panther_accession, "\\|"))
  Short.Name = panther_accession_split[1]
  ensembl_accession = panther_accession_split[2]
  ensembl_accession = unlist(stringr::str_split(ensembl_accession, "="))[2]
  panther_classifications_df[i, 'Short.Name'] = Short.Name
  panther_classifications_df[i, 'panther_accession'] = panther_accession
  panther_classifications_df[i, 'ensembl_accession'] = ensembl_accession
}




alreadinpip1 = c("Amborella trichopoda","Arabidopsis thaliana","Brachypodium distachyon","Brassica napus",
                 "Capsicum annuum","Chlamydomonas reinhardtii","Cucumis sativus","Glycine max","Helianthus annuus",
                 "Manihot esculenta","Medicago truncatula","Physcomitrella patens","Populus trichocarpa",
                 "Prunus persica","Selaginella moellendorffii","Setaria italica","Solanum lycopersicum",
                 "Solanum tuberosum","Sorghum bicolor","Triticum aestivum","Vitis vinifera","Zea mays")
notInPip1 = inBoth[!inBoth %in% alreadinpip1 ]
notInNewList = alreadinpip1[!alreadinpip1 %in% inBoth ]
notInPip1df = data.frame(notInPIP1 = notInPip1)



library(curl)
url = "ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/"
h = new_handle(dirlistonly=TRUE)
con = curl(url, "r", h)
tbl = read.table(con, stringsAsFactors=TRUE, fill=TRUE)
close(con)
head(tbl)
latestVersion = unlist(regmatches(as.vector(tbl$V1[1]), regexpr("_", as.vector(tbl$V1[1])), invert = TRUE))[1]


tmp <- tempfile()
curl_download("ftp://ftp.pantherdb.org/sequence_classifications/current_release/species", tmp)
panther_species = read.table(tmp,fill = TRUE, col.names = c("FTPname", "short_name", "Identifier", "taxon_id"))
panther_species = left_join(panther_species, organisms, by = "short_name")

panther_species$long_name[panther_species$long_name=="Physcomitrella patens"] = "Physcomitrium patens"
organismsInEnsemble = panther_species[panther_species$long_name %in% ensemble_species,]


zz=gzfile("D:/221001 gene info/All_Plants.gene_info.gz",'rt')
gene_info =read.csv(zz, header = TRUE, sep = "\t")
close(zz)
gene_info.narrow = gene_info %>%
  select(X.tax_id, GeneID, Symbol, LocusTag, Synonyms)


for (i in seq(organismsInEnsemble$FTPname)){
  plant = organismsInEnsemble$FTPname[i]
  species = organismsInEnsemble$long_name[i]
  shortName = organismsInEnsemble$short_name[i]
  taxonId = organismsInEnsemble$taxon_id.x[i]
  plant_url = paste0(url, latestVersion,"_",plant)

  gene_info.taxon = gene_info.narrow[gene_info.narrow$X.tax_id == taxonId,]

  tmp <- tempfile()
  print(paste0("Downloading ",plant, " (", latestVersion,"_",plant,")"))
  curl_download(plant_url, tmp, quiet = F)
  plant_df_temp = read.table(tmp, sep = "\t", fill = TRUE, col.names  = c("panther_accession", "uniprotsptrembl", "id2", "hmmpanther_subfams",
                                                                          "FAMILY_TERM", "SUBFAMILY_TERM","molecular_function",
                                                                          "biological_process","cellular_components",
                                                                          "protein_class","pathway"),
                             quote = "")
  plant_df_temp = dplyr::select(plant_df_temp, -molecular_function,-biological_process,-cellular_components,
                                -protein_class,-pathway)
  plant_df_temp$Spcies = species
  ensemblDataset = datasets[datasets$Species==species,"dataset"]


  if (length(ensemblDataset)>1){
    ensemblDataset = userListElementSelection(ensemblDataset, paste0("Multiple possible datasets found. Which ensembl dataset does the Panther dataset ",
                                                                     shortName, " correspond to?"))}

  ensembl <- useMart(biomart = "plants_mart",
                     dataset = ensemblDataset,
                     host = "https://plants.ensembl.org")
  attributes = listAttributes(ensembl)

  attributesPertenant = c("ensembl_gene_id",
                          "uniprotsptrembl","uniprotswissprot", "external_gene_name", "entrezgene_id")
  presentAttributes = attributesPertenant[attributesPertenant %in% attributes$name]
  uniprotsptrembl<- getBM(attributes = presentAttributes,
                          mart = ensembl)
  trembl = uniprotsptrembl[,c("ensembl_gene_id", "uniprotsptrembl")]
  swiss = uniprotsptrembl[,c("ensembl_gene_id", "uniprotswissprot")]


  plant_df_tempTempbl = left_join(plant_df_temp, trembl, by = "uniprotsptrembl")  %>%
    select(ensembl_gene_id, panther_accession) %>%
    rename(ensembl_gene_id.trembl = ensembl_gene_id)

  plant_df_tempSwiss = plant_df_temp %>%
    rename(uniprotswissprot = uniprotsptrembl) %>%
    left_join(swiss, by = "uniprotswissprot") %>%
    select(ensembl_gene_id, panther_accession) %>%
    rename(ensembl_gene_id.swiss = ensembl_gene_id)


  if ("external_gene_name" %in% presentAttributes){
    external_id = uniprotsptrembl %>%
      select(ensembl_gene_id,external_gene_name) %>%
      mutate(external_gene_name = toupper(external_gene_name)) %>%
      distinct(external_gene_name, .keep_all = TRUE)

    plant_df_tempExternal = plant_df_temp %>%
      rename(external_gene_name = id2) %>%
      mutate(external_gene_name = toupper(external_gene_name)) %>%
      left_join(external_id, by = "external_gene_name") %>%
      select(ensembl_gene_id, panther_accession) %>%
      rename(ensembl_gene_id.external = ensembl_gene_id)

  }else{presentAttributes = plant_df_temp
  presentAttributes$ensembl_gene_id.external = ""
  }

  plant_df_temp = left_join(plant_df_temp, plant_df_tempTempbl, by = 'panther_accession') %>%
    left_join(plant_df_tempSwiss, by = 'panther_accession') %>%
    left_join(plant_df_tempExternal, by = 'panther_accession')

  plant_df_temp = plant_df_temp %>%
    rename(Symbol = id2) %>%
    left_join(gene_info.taxon, by = "Symbol") %>%
    rename(entrezgene_id = GeneID)


  if (entrezgene_id %in% presentAttributes){
  entrezEnsembl<- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"),
                          filters = c("entrezgene_id"),
                          values = plant_df_temp$entrezgene_id,
                          mart = ensembl)

  entrezEnsembl = entrezEnsembl %>% rename(ensembl_gene_id.entrez = ensembl_gene_id)
  plant_df_temp = left_join(plant_df_temp,entrezEnsembl, by = "entrezgene_id" )
  } else{plant_df_temp$entrezgene_id = ""
  }

  if (plant == organismsInEnsemble$name[1]){
    plant_df = plant_df_temp
  }else{plant_df = rbind(plant_df, plant_df_temp)}
}




plant_df$ensembl_gene_id = plant_df$ensembl_gene_id.trembl
plant_df$ensembl_gene_id[is.na(plant_df$ensembl_gene_id.trembl)] = plant_df$ensembl_gene_id.swiss[is.na(plant_df$ensembl_gene_id.trembl)]
plant_df.na = plant_df[is.na(plant_df$ensembl_gene_id),]

plant_df.na = plant_df.na %>% tidyr::separate(panther_accession, c(NA, "codes", NA), sep = "\\|",remove = F)
plant_df.na = plant_df.na %>% tidyr::separate(codes, c("left", "right"), sep = "=",remove = F,extra = "merge")

plant_df.na.Gene_OrderedLocusName = plant_df.na %>%
  select(-ensembl_gene_id)  %>%
  filter(left == "Gene_OrderedLocusName") %>%
  mutate(right = toupper(right)) %>%
  rename(ensembl_gene_id.Gene_OrderedLocusName = right) %>%
  select(panther_accession, ensembl_gene_id.Gene_OrderedLocusName)
plant_df = left_join(plant_df, plant_df.na.Gene_OrderedLocusName, by = "panther_accession")
plant_df$ensembl_gene_id[is.na(plant_df$ensembl_gene_id)] = plant_df$ensembl_gene_id.Gene_OrderedLocusName[is.na(plant_df$ensembl_gene_id)]


plant_df.na.Gene_ORFName = plant_df.na %>%
  select(-ensembl_gene_id)  %>%
  filter(left == "Gene_ORFName") %>%
  rename(ensembl_gene_id.Gene_ORFName = right) %>%
  select(panther_accession, ensembl_gene_id.Gene_ORFName)
plant_df = left_join(plant_df, plant_df.na.Gene_ORFName, by = "panther_accession")
plant_df$ensembl_gene_id[is.na(plant_df$ensembl_gene_id)] = plant_df$ensembl_gene_id.Gene_ORFName[is.na(plant_df$ensembl_gene_id)]


plant_df.na.ensembl = plant_df.na %>%
  select(-ensembl_gene_id)  %>%
  filter(left == "EnsemblGenome") %>%
  rename(ensembl_gene_id.EnsemblGenome = right) %>%
  select(panther_accession, ensembl_gene_id.EnsemblGenome)
plant_df = left_join(plant_df, plant_df.na.ensembl, by = "panther_accession")
plant_df$ensembl_gene_id[is.na(plant_df$ensembl_gene_id)] = plant_df$ensembl_gene_id.EnsemblGenome[is.na(plant_df$ensembl_gene_id)]

plant_df$ensembl_gene_id[is.na(plant_df$ensembl_gene_id)] = plant_df$ensembl_gene_id.external[is.na(plant_df$ensembl_gene_id)]

plant_df$ensembl_gene_id[is.na(plant_df$ensembl_gene_id)] = plant_df$ensembl_gene_id.entrez[is.na(plant_df$ensembl_gene_id)]

plant_df.na.remaining = plant_df[is.na(plant_df$ensembl_gene_id),]


hmmpanther = str_split_fixed(plant_df$hmmpanther_subfams, ":", 2)
plant_df[,'hmmpanther'] = hmmpanther[,1]
write.csv(plant_df, "D:/Cleaned whole PD pipeline/pip_directory_v6/PANTHER_PD_pipeline_classifications.csv", row.names = F)

plant_df = read.csv("D:/Cleaned whole PD pipeline/pip_directory_v6/PANTHER_PD_pipeline_classifications.csv")


stop()
#______________________________________________________________________________________

ensembl <- useMart(biomart = "plants_mart",
                   dataset = "jregia_eg_gene",
                   host = "https://plants.ensembl.org")

attributes = listAttributes(ensembl)

uniprotsptrembl<- getBM(attributes = c("ensembl_gene_id",
                                       "entrezgene_accession", "strand", "entrezgene_id"),
                        #filters = c("ensembl_gene_id"),
                        #values = df$ensembl_gene_id,
                        mart = ensembl)
test = uniprotsptrembl[uniprotsptrembl$uniprotsptrembl != "",]
colnames(plant_df_temp)[colnames(plant_df_temp)=="uniprotsptrembl"] = "uniprotswissprot"
plant_df_temp = left_join(plant_df_temp, uniprotsptrembl, by = "uniprotswissprot")


tmp <- tempfile()
curl::curl_download("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Plants/All_Plants.gene_info.gz", tmp)
zz=gzfile(tmp,'rt')
gene_info =read.csv(zz, header = TRUE, sep = "\t")
close(zz)
gene_test = gene_info[gene_info$X.tax_id == 51240,]
# ensemblUniprot <- useMart("unimart", dataset="uniprot")
# attributesEnsembl = listAttributes(ensemblUniprot)
# datasetsEnsembl <- listDatasets(ensembl)


tmp <- tempfile()
curl::curl_download("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz", tmp)
zz=gzfile(tmp,'rt')
gene2ensembl =read.csv(zz, header = TRUE, sep = "\t")
close(zz)

pantherCompatibleEnsemlDatasets = pantherCompatibleEnsemlDatasets()
gene_info.clean = gene_info[gene_info$X.tax_id %in% pantherCompatibleEnsemlDatasets$taxon_id,]
gene2ensembl.clean = gene2ensembl[gene2ensembl$X.tax_id =="3702",]
gene2ensembl.clean = gene2ensembl[gene2ensembl$X.tax_id %in% pantherCompatibleEnsemlDatasets$taxon_id,]
gene_info_ensemble = left_join(gene_info.clean, gene2ensembl.clean, by = "X.tax_id")

gene_test$Symbol = toupper(gene_test$Symbol)
gene_info_na_join = plant_df.na.remaining %>%
  rename(Symbol = id2) %>%
  mutate(Symbol = toupper(Symbol))  %>%
  left_join(gene_info, by = "Symbol")


zz=gzfile("D:/221001 gene info/gene2ensembl.gz",'rt')
gene2ensembl =read.csv(zz, header = TRUE, sep = "\t")
close(zz)
