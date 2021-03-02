# PD candidate generator
R pipeline for generating and classifying plasmodesmata (PD) candidates based on proteomics and feature prediction.  
The pipeline pulls all genes for a given species that are in a subfamily present in any of the proteomes included as an input.  
Genes are classified by how many proteomes the subfamily appears in and whether the gene has features that are overrepresented in known PD genes  
Please see ******* for more information
![schematic of PD pipeline](https://github.com/PhilPlantMan/PD-candidate-generator/blob/main/Figure%201.png)

## Getting Started

### Prerequisites

This pipeline has only been tested in Windows 10. Issues with other operating systems are likely  
This pipeline was written with R version 4.0.2.  
All other dependencies should install automatically  
If you are unfamiliar with R, please install [RStudio](https://rstudio.com/products/rstudio/download/)

### Installation and execution

* Download and unzip "directory.zip"
* Follow the instructions provided in the README within

## Built With

* [biomartr](https://github.com/ropensci/biomartr) - Peptide sequence retrieval
* [ragp](https://github.com/missuse/ragp) - Protein feature prediction

## Publication

**********
