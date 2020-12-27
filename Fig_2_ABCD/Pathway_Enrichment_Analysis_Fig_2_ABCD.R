####################################################################
## This code is to reproduce the data and images of Figures 2A-C  ##
## Author: Wasim Aftab                                            ##
####################################################################

## Clear variables
rm(list = ls())

## Clear screen
cat('\014')

## Installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
list.of.packages <- c("clusterProfiler", "org.Mm.eg.db", "ReactomePA")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  BiocManager::install(new.packages)

## Installing CRAN packages
list.of.packages <- c("readxl", "writexl", "rstudioapi")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

## Loading packages
library(ReactomePA)
library(clusterProfiler)
library(org.Mm.eg.db)
library(readxl)
library(rstudioapi)

## Chdir to source dir
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

## Helper function to compute pathway enrichment
compute_GO_Pathway <- function(Uniprot_list, bg_uniprot) {
  eg <-
    clusterProfiler::bitr(
      Uniprot_list,
      fromType = "UNIPROT",
      toType = "ENTREZID",
      OrgDb = "org.Mm.eg.db"
    )

  ## Convert background Uniprots to ENTREZ IDs
  eg_bg <-
    clusterProfiler::bitr(
      bg_uniprot$Uniprot ,
      fromType = "UNIPROT",
      toType = "ENTREZID",
      OrgDb = "org.Mm.eg.db"
    )
  x <-
    ReactomePA::enrichPathway(
      gene = eg$ENTREZID,
      organism = "mouse",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe = eg_bg$ENTREZID,
      readable = TRUE,
      minGSSize = 1
    )
  return(x)
}


## WT Exclusive Pathway (Fig_2_ABC)
temp <- read_xlsx('Protein_List_WT_AROM+_SL.xlsx', 2)
Uniprot_list <- as.matrix(temp$WT)
idx <- union(which(is.na(Uniprot_list)), 1)
Uniprot_list_WT_Exclusive <- Uniprot_list[-idx]
bg_uniprot <-
  as.data.frame(read.table("Mouse_Proteome_Uniprot.txt", header = T))
x <- compute_GO_Pathway(Uniprot_list_WT_Exclusive, bg_uniprot)
idx <- which(x@result$p.adjust < 0.05)
if (length(idx) != 0){
  print(ReactomePA::emapplot(x, showCategory = 30, color = "p.adjust"))
  writexl::write_xlsx(
    x@result[idx,],
    '../Cytoscape_Data_Generation/Sample_Data_for_pathway_analysis/Fig2ABC_Pathway_Enrichment.xlsx',
    col_names = TRUE
  )
}


## AROM+ Exclusive Pathway (Fig_2D)
temp <- read_xlsx('Protein_List_WT_AROM+_SL.xlsx', 2)
Uniprot_list <- as.matrix(temp$`AROM+`)
idx <- union(which(is.na(Uniprot_list)), 1)
Uniprot_list_AROM_Exclusive <- Uniprot_list[-idx]
bg_uniprot <-
  as.data.frame(read.table("Mouse_Proteome_Uniprot.txt", header = T))
x <- compute_GO_Pathway(Uniprot_list_AROM_Exclusive, bg_uniprot)
idx <- which(x@result$p.adjust < 0.05)
if (length(idx) != 0){
  print(ReactomePA::emapplot(x, showCategory = 30, color = "p.adjust"))
  writexl::write_xlsx(
    x@result[idx,],
    '../Cytoscape_Data_Generation/Sample_Data_for_pathway_analysis/Fig2D_Pathway_Enrichment.xlsx',
    col_names = TRUE
  )
}
