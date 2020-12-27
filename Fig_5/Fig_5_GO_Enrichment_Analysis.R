####################################################################
## This code is to reproduce the data and images of Figures 5B-D  ##
## Author: Wasim Aftab                                            ##
####################################################################

## Clear variables
rm(list = ls())

## Clear screen
cat('\014')

## Installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
list.of.packages <- c("clusterProfiler", "org.Mm.eg.db")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  BiocManager::install(new.packages)

## Installing CRAN packages
list.of.packages <- c("readxl", "rstudioapi")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

## Loading packages
library(clusterProfiler)
library(org.Mm.eg.db)
library(readxl)
library(rstudioapi)

## Chdir to source dir
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

## Read background Uniprots
bg_uniprot <-
  as.data.frame(read.table("Data/Mouse_Proteome_Uniprot.txt",
                           header = TRUE))

## Convert background Uniprots to ENTREZ IDs
eg_bg <-
  clusterProfiler::bitr(
    bg_uniprot$Uniprot ,
    fromType = "UNIPROT",
    toType = "ENTREZID",
    OrgDb = "org.Mm.eg.db"
  )

## Helper function to compute GO enrichment
compute_GO_enrichment <- function(eg, eg_bg) {
  ego <-
    clusterProfiler::enrichGO(
      gene = eg$ENTREZID,
      OrgDb = "org.Mm.eg.db",
      universe = eg_bg$ENTREZID,
      ont = "BP",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      readable = TRUE,
      minGSSize = 1,
      maxGSSize = 500
    )
  return(ego)
}


## Fig 5B
T <-
  readxl::read_xls('Data/Top_MLP_hits_Red_Green_Tub_WT_Dscrmn_Mass_LIMMA_FC_less_0_Pvalue_DC.xls')
eg <-
  clusterProfiler::bitr(T$Uniprot,
                        fromType = "UNIPROT",
                        toType = "ENTREZID",
                        OrgDb = "org.Mm.eg.db")

ego <- compute_GO_enrichment(eg, eg_bg)

if (min(ego@result$p.adjust) < 0.05) {
  p <- enrichplot::cnetplot(ego, showCategory = 20)
  plot(p)
  idx <- which(ego@result$p.adjust < 0.05)
  writexl::write_xlsx(
    ego@result[idx,],
    '../Cytoscape_Data_Generation/Sample_Data_for_GO_analysis/Fig5B_GO_Enrichment.xlsx',
    col_names = TRUE
  )
}


## Fig 5C
T <-
  readxl::read_xls(
    'Data/Top_MLP_hits_Brown_Green_Tub_WT_Dscrmn_Mass_LIMMA_FC_less_0_Pvalue_DC.xls'
  )
eg <-
  clusterProfiler::bitr(T$Uniprot,
                        fromType = "UNIPROT",
                        toType = "ENTREZID",
                        OrgDb = "org.Mm.eg.db")

ego <- compute_GO_enrichment(eg, eg_bg)

if (min(ego@result$p.adjust) < 0.05) {
  p <- enrichplot::cnetplot(ego, showCategory = 20)
  plot(p)
  idx <- which(ego@result$p.adjust < 0.05)
  writexl::write_xlsx(
    ego@result[idx,],
    '../Cytoscape_Data_Generation/Sample_Data_for_GO_analysis/Fig5C_GO_Enrichment.xlsx',
    col_names = TRUE
  )
}

## Fig 5D
T <-
  readxl::read_xls(
    'Data/Top_MLP_hits_Tub3N_Tub2_WT_Dscrmn_Mass_LIMMA_FC_less_0_Pvalue_DC.xls'
  )
eg <-
  clusterProfiler::bitr(T$Uniprot,
                        fromType = "UNIPROT",
                        toType = "ENTREZID",
                        OrgDb = "org.Mm.eg.db")

ego <- compute_GO_enrichment(eg, eg_bg)

if (min(ego@result$p.adjust) < 0.05) {
  p <- enrichplot::cnetplot(ego, showCategory = 20)
  plot(p)
  idx <- which(ego@result$p.adjust < 0.05)
  writexl::write_xlsx(
    ego@result[idx,],
    '../Cytoscape_Data_Generation/Sample_Data_for_GO_analysis/Fig5D_GO_Enrichment.xlsx',
    col_names = TRUE
  )
}

