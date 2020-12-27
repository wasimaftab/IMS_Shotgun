################################################################################
## This code is to reproduce the data and image of Figure  1B   
## Since the missing values are imputed randomly (from normal distribution), 
## One can notice minute changes in numbers associated with fold change and p-values 
## across multiple runs of the code. However, the major patterns in the data do not change. 
## Author: Wasim Aftab
################################################################################
cat('\014')
rm(list = ls())

## Installing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
list.of.packages <- c("limma", "qvalue")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

## Installing CRAN packages
list.of.packages <- c("dplyr", "stringr", "MASS", "matlab", "plotly", "htmlwidgets", "rstudioapi", "webshot", "matrixStats")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if (is.null(webshot:::find_phantom())){
  webshot::install_phantomjs()
}

library(dplyr)
library(stringr)
library(MASS)
library(matlab)
library(plotly)
library(limma)
library(qvalue)
library(htmlwidgets)
library(rstudioapi)

##Chdir to source dir
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))
cur_dir <- getwd()

source("limma_helper_functions_Fig1B.R")

##Load the proteingroups file
myFilePath <- "proteinGroups_180416.txt"
proteingroups <-
  as.data.frame(read.table(myFilePath, header = TRUE, sep = "\t"))
###Kill the code if proteingroups does not contain crap columns#######################################
temp <-
  select(
    proteingroups,
    matches("(Reverse|Potential.contaminant|Only.identified.by.site)")
  )
if (!nrow(temp) * ncol(temp)) {
  stop("File error, It does not contain crap...enter another file with crap")
}
######################################################################################################
##Display data to faciliate choice of treatment and control
temp <- select(proteingroups, matches("(ibaq|lfq)"))
print(names(temp))

# #remove "+" identifeid rows from proteingroups##################################################
idx <- NULL
# temp <- as.matrix(select(proteingroups, matches("(Only.identified.by.site|Reverse|Potential.contaminant)")))
temp <- select(proteingroups, matches("(Only.identified.by.site|Reverse|Potential.contaminant)"))
for (i in 1:ncol(temp)){
  index <- which(unlist(!is.na(match(temp[,i], "+"))))
  idx <- union(idx, index)
}
proteingroups <- proteingroups[-idx, ] # removing indexed rows
# ################################################################################################
# 
#Extrat Uniprot and gene symbols
Uniprot <- character(length = nrow(proteingroups))
Symbol <- character(length = nrow(proteingroups))
ProteinNames <- proteingroups$Protein.names

for (i in 1:nrow(proteingroups)) {
  temp <- as.character(proteingroups$Fasta.headers[i])
  splits <- unlist(strsplit(temp, '\\|'))
  Uniprot[i] <- splits[2]
  splits <- unlist(str_match(splits[3], "GN=(.*?) PE="))
  Symbol[i] <- splits[2]
}

#Extract required data for Limma
treatment <-
  readline('Enter treatment name(case insensitive) as it appeared in the iBAQ/LFQ column= ')
control <-
  readline('Enter control name(case insensitive) as it appeared in the iBAQ/LFQ column= ')
ibaq <- readinteger_binary('Enter 1 for iBAQ or 0 for LFQ= ')

if (ibaq) {
  temp <-
    select(proteingroups, matches(paste('^.*', "ibaq", '.*$', sep = '')))
  treatment_reps <- data_sanity_check(temp, 'treatment', treatment)
  control_reps <- select(temp, matches(control))
  control_reps <- data_sanity_check(temp, 'control', control)
  data <-
    cbind(treatment_reps,
          control_reps,
          select(proteingroups, matches("^id$")),
          Uniprot,
          Symbol,
          ProteinNames)
} else {
  temp <-
    select(proteingroups, matches(paste('^.*', "lfq", '.*$', sep = '')))
  treatment_reps <- select(temp, matches(treatment))
  treatment_reps <- data_sanity_check(temp, 'treatment', treatment)
  control_reps <- select(temp, matches(control))
  control_reps <- data_sanity_check(temp, 'control', control)
  data <-
    cbind(treatment_reps,
          control_reps,
          select(proteingroups, matches("^id$")),
          Uniprot,
          Symbol,
          ProteinNames)
}

print(names(data))
rep_treats <-
  readinteger("Enter the number of treatment replicates=")
rep_conts <- readinteger("Enter the number of control replicates=")
FC_Cutoff <- readfloat("Enter the fold change cut off=")

## Find out Blank rows, i.e. proteins with all zeros in treatment and in control, see followig example
# iBAQ.Mrpl40_1 iBAQ.Mrpl40_2 iBAQ.Mrpl40_3 iBAQ.Kgd4_1 iBAQ.Kgd4_2 iBAQ.Kgd4_3  id Uniprot Symbol
#-------------------------------------------------------------------------------------------------
#       0             0             0           0           0           0        84  Q02888  INA17
temp <-
  as.matrix(rowSums(apply(data[, 1:(rep_treats + rep_conts)], 2, as.numeric)))
idx <- which(temp == 0)
if (length(idx)) {
  data <- data[-idx,] # removing blank rows
}

# ## Find out outliers, i.e. proteins with all zeros in treatment and all/some values in control, see followig example
# ## Uniprot  Symbol  treat_1   treat_2   treat_3   contrl_1    contrl_2    contrl_3
# ##-----------------------------------------------------------------------------------
# ## P25554   SGF29	    0	        0	        0        2810900	     0	        0
# temp <- as.matrix(rowSums(data[,1:rep_treats]))
temp <-
  as.matrix(rowSums(apply(data[, 1:rep_treats], 2, as.numeric)))
idx <- which(temp == 0)
if (length(idx)) {
  outliers <- data[idx,]
  filename_outliers <- "Exclusive_proteins_WT"
    # paste("Outliers_treatment_", treatment, "_", control, sep = "")
  data <- data[-idx, ] # removing indexed rows
}

# ## Find out outliers, i.e. proteins with all zeros in control and all/some values in treatment, see followig example
# iBAQ.Mrpl40_1 iBAQ.Mrpl40_2 iBAQ.Mrpl40_3 iBAQ.Kgd4_1 iBAQ.Kgd4_2 iBAQ.Kgd4_3  id   Uniprot   Symbol
##-----------------------------------------------------------------------------------------------------
#     662810        505600        559130        0           0           0        79   P38845    CRP1
# temp <-
#   as.matrix(rowSums(data[, (rep_treats + 1):(rep_conts + rep_treats)]))
temp <-
  as.matrix(rowSums(apply(data[, (rep_treats + 1):(rep_conts + rep_treats)], 2, as.numeric)))
idx <- which(temp == 0)
if (length(idx)) {
  outliers_control <- data[idx,]
  filename_outliers_control <- "Exclusive_proteins_AROM"
    # paste("Outliers_control_", treatment, "_", control, sep = "")
  data <- data[-idx, ] # removing indexed rows
}

#Impute data
data_limma <- log2(as.matrix(data[c(1:(rep_treats + rep_conts))]))
data_limma[is.infinite(data_limma)] <- NA
nan_idx <- which(is.na(data_limma))

fit <- fitdistr(c(na.exclude(data_limma)), "normal")
mu <- as.double(fit$estimate[1])
sigma <- as.double(fit$estimate[2])
sigma_cutoff <- 6
new_width_cutoff <- 0.3
downshift <- 1.8
width <- sigma_cutoff * sigma
new_width <- width * new_width_cutoff
new_sigma <- new_width / sigma_cutoff
new_mean <- mu - downshift * sigma
imputed_vals_my = rnorm(length(nan_idx), new_mean, new_sigma)
data_limma[nan_idx] <- imputed_vals_my

##Limma main code
design <-
  model.matrix( ~ factor(c(rep(2, rep_treats), rep(1, rep_conts))))
colnames(design) <- c("Intercept", "Diff")
res.eb <- eb.fit(data_limma, design, data$Symbol)
Sig_FC_idx <-
  union(which(res.eb$logFC < (-FC_Cutoff)), which(res.eb$logFC > FC_Cutoff))
Sig_Pval_mod_idx <- which(res.eb$p.mod < 0.05)
Sig_Pval_ord_idx <- which(res.eb$p.ord < 0.05)
Sig_mod_idx <-  intersect(Sig_FC_idx, Sig_Pval_mod_idx)
Sig_ord_idx <-  intersect(Sig_FC_idx, Sig_Pval_ord_idx)
categ_Ord <- rep(c("Not Significant"), times = length(data$Symbol))
categ_Mod <- categ_Ord
categ_Mod[Sig_mod_idx] <- "Significant"
categ_Ord[Sig_ord_idx] <- "Significant"

dat <-
  cbind(
    res.eb,
    categ_Ord,
    categ_Mod,
    NegLogPvalMod = (-log10(res.eb$p.mod)),
    NegLogPvalOrd = (-log10(res.eb$p.ord))
  )

dat <- select(dat, -matches("^gene$"))
# dat <- select(dat, -matches("^Uniprot$"))
##Save the data file
final_data <-
  cbind(data$Uniprot,
        data$ProteinNames,
        data$Symbol,
        data[,1:(rep_treats+rep_conts)],
        dat)

colnames(final_data)[1] <- c("Uniprot")
colnames(final_data)[2] <- c("ProteinNames")
colnames(final_data)[3] <- c("Symbol")
filename_final_data <-
  readline('Enter a filename for final data= ')

##Create plotly object and save plot as html
filename_mod <-
  readline('Enter a filename for limma plot= ')
# filename_ord <-
#   readline('Enter a filename for ordinary t-test plot= ')

# display_plotly_figs(final_data, FC_Cutoff, filename_mod, filename_ord)
display_plotly_figs(final_data, FC_Cutoff, filename_mod)


## Write final data
write.table(
  final_data,
  paste(filename_final_data, '.tsv', sep = ''),
  sep = '\t',
  row.names = FALSE,
  col.names = TRUE
)

## Write outliers in treatment
write.table(
  outliers,
  paste(filename_outliers, '.tsv', sep = ''),
  sep = '\t',
  row.names = FALSE,
  col.names = TRUE
)

## Write outliers in control
write.table(
  outliers_control,
  paste(filename_outliers_control, '.tsv', sep = ''),
  sep = '\t',
  row.names = FALSE,
  col.names = TRUE
)

setwd(cur_dir)