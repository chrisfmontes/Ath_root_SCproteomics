# SCP (Single Cell Proteomics) package training script

library(scp)
library(ggplot2)
library(magrittr)
library(dplyr)
library(impute)
library(scater)
library(sva)

#Load gene annotation file
AthGeneSymbols <- read.csv("Ath_symbol_file.csv", header = TRUE, stringsAsFactors = FALSE)
#Load MQ data
mySCdata <- read.table(file = "../DART_processed_files/output/ev_updated.txt", header = TRUE, sep = "\t", comment.char = "") %>%
  mutate(UID = paste0("UID", id))

#Load (generate) experiment metadata
raw.files <- unique(mySCdata$Raw.file)
nruns <- length(unique(mySCdata$Raw.file))
SClabellingmap <- read.csv("../TMTlabelMap.csv", header = TRUE)%>%select(4,5,6)
mySCmetadata <- data.frame(Raw.file = rep(raw.files, each = 18),
                           Channel = rep(paste0("Reporter.intensity.", 1:18), times = nruns),
                           PlexNumber = rep(paste0("Plex_", 1:nruns), each = 18),
                           TMTlabel = rep(1:18, times = nruns)) %>%
  left_join(read.csv("../TMTlabelMap.csv", header = TRUE)%>%select(4,5,6),
            by = c("PlexNumber","TMTlabel")) %>%
  select(Raw.file, Channel, SampleType) %>%
  mutate(SampleType = case_when(Channel == "Reporter.intensity.1" ~ "Carrier", TRUE ~ as.character(SampleType)),
         SampleType = case_when(Channel == "Reporter.intensity.2" ~ "Reference", TRUE ~ as.character(SampleType)),
         SampleType = case_when(Channel == "Reporter.intensity.3" | Channel == "Reporter.intensity.4" ~ "Unused", TRUE ~ as.character(SampleType)))


#Create the SingleCellEperiment object
scp <- readSCP(featureData = mySCdata,
               colData = mySCmetadata,
               channelCol = "Channel",
               batchCol = "Raw.file",
               removeEmptyCols = TRUE)

mean(dims(scp)[1,])
sum(dims(scp[, scp$SampleType %in% c("cortex", "endodermis"), ])[2,])

# Clean missing data, change zeroes to NA
scp <- zeroIsNA(scp, i = 1:nruns)
scp 

#Filter-out contaminants, decoy, and low confidence (high PIF)
scp_bckp1 <- scp
#scp <- scp_bckp1
scp <- filterFeatures(scp,
                      ~ Reverse != "+" &
                        Potential.contaminant != "+" &
                        !is.na(PIF) & PIF > 0.8)
# scp <- filterFeatures(scp,
#                       ~ Reverse != "+" &
#                         Potential.contaminant != "+" &
#                         !is.na(PIF))
scp

#get average number of PSM per cell
mean(dims(scp)[1,])
sum(dims(scp[, scp$SampleType %in% c("cortex", "endodermis"), ])[2,])

#Filter by number of detected features
dims(scp)
keepAssay <- dims(scp)[1, ] > 150
scp <- scp[, , keepAssay]

#Check our dataset
scp
mean(dims(scp)[1,])
sum(dims(scp[, scp$SampleType %in% c("cortex", "endodermis"), ])[2,])

#Filter by SCP metrics, sample-to-carrier ratio (SCR)
table(colData(scp)[,"SampleType"]) # to see different channel types

scp <- computeSCR(scp,
                  i = 1:length(scp),
                  colvar = "SampleType",
                  carrierPattern = "Carrier",
                  samplePattern = "cortex|endodermis",
                  sampleFUN = "mean",
                  rowDataName = "MeanSCR")

#we can plot the average SCR distribution
rbindRowData(scp, i = 1:length(scp)) %>%
  data.frame %>%
  ggplot(aes(x = MeanSCR)) +
  geom_histogram() +
  geom_vline(xintercept = c(1/60, 1/135, 0.16,0.3,1),
             lty = c(4, 2, 1,1,1)) +
  scale_x_log10()
#we filter those PSM with mean SCR > 0.3 (the expected mean SCR was 1/135 = 0.007, dashed line on the graph)
#but it looks like our cells are in a 1:60 = 0.016 SCR (mis-quantification or mass per cell wrong?). We'll remove
#the data with too high SCR to avoid using misquantified/misprocessed wells
scp_bckp2 <-  scp
#scp <- scp_bckp2
scp <- filterFeatures(scp,
                      ~ !is.na(MeanSCR) &
                        MeanSCR < 0.3)
                        #MeanSCR < 0.16)
                        #MeanSCR < 1)
#Let's see how is our dataset now
scp
mean(dims(scp)[1,])
sum(dims(scp[, scp$SampleType %in% c("cortex", "endodermis"), ])[2,])

#plot the filtered SCR distribution
rbindRowData(scp, i = 1:length(scp)) %>%
  data.frame %>%
  ggplot(aes(x = MeanSCR)) +
  geom_histogram() +
  geom_vline(xintercept = c(1/60, 1/135, 0.3),
             lty = c(4, 2, 1)) +
  scale_x_log10()
# Control for FDR -> convert PEP to qvalue. It seems that PEP is too conservative
scp <- pep2qvalue(scp,
                  i = 1:length(scp),
                  PEP = "dart_PEP",
                  rowDataName = "qvalue_PSMs") # we can do this on PSM
scp <- pep2qvalue(scp,
                  i = 1:length(scp),
                  PEP = "dart_PEP",
                  groupBy = "Leading.razor.protein",
                  rowDataName = "qvalue_proteins") #Also, we can do PEP to qvalue on proteins (or peptides)
scp_bckp3 <-  scp
#scp <- scp_bckp3
scp <- filterFeatures(scp,
                      ~ qvalue_proteins < 0.01) # we filter PSM, to control for protein FDR at 1%
# scp <- filterFeatures(scp,
#                       ~ qvalue_PSMs < 0.01)
mean(dims(scp)[1,])
sum(dims(scp[, scp$SampleType %in% c("cortex", "endodermis"), ])[2,])
scp
# 6. Process PSM data
#We normalize SC relative reporter ion intensity (RII) to reference channel
scp <- divideByReference(scp,
                         i = 1:length(scp),
                         colvar = "SampleType",
                         samplePattern = ".",
                         refPattern = "Reference")

# 7. Aggregate PSM into peptide data
scp <- aggregateFeaturesOverAssays(scp,
                                   i = 1:length(scp),
                                   fcol = "Modified.sequence", #Column to use for peptide ID 
                                   name = paste0("peptides_", names(scp)),    #name to create for the aggregated dataset
                                   fun = matrixStats::colMedians, na.rm = TRUE) #function to use when aggregating PSM intensity data
scp
mean(dims(scp)[1,55:108])
sum(dims(scp[, scp$SampleType %in% c("cortex", "endodermis"), ])[2,1:54])
scp_bckp4 <- scp
#create a scp object to count proteins per cell before RII filtering
scp_pre_RII <- aggregateFeaturesOverAssays(scp,
                                  i = 1:54,
                                  fcol = "Leading.razor.protein", #Column to use for peptide ID
                                  name = paste0("proteins_", names(scp)[1:54]),    #name to create for the aggregated dataset
                                  fun = matrixStats::colMedians, na.rm = TRUE) #function to use when aggregating PSM intensity data
scp_pre_RII
mean(dims(scp_pre_RII)[1,109:162])

# 8. Join SCoPE2 sets in one assay
scp <- joinAssays(scp,
                  i = c(55:length(scp)),
                  name = "peptides") # We join the aggregated datasets into one dataset named "peptides"

#Join pre RII filter proteins for summary
scp_pre_RII <- joinAssays(scp_pre_RII,
                  i = 109:162,
                  name = "proteins")
mean(dims(scp)[1,55:108])

#plot(scp)

# 9. Filter single-cells by median RII and CV per cell
# We subset the data, extracting only SC data (samples and empty wells, but not empty channels, carrier, or reference)
scp <- scp[, scp$SampleType %in% c("cortex", "endodermis"), ]
table(colData(scp)[,"SampleType"])
mean(dims(scp)[1,55:108])
dims(scp)[1,109]
scp

# Filter by median RII
medians <- colMedians(assay(scp[["peptides"]]), na.rm = TRUE)
scp$MedianRI <- medians

# visualize median RII per sample type
colData(scp) %>%
  data.frame %>%
  ggplot() +
  aes(x = MedianRI, 
      y = SampleType,
      fill = SampleType) +
  geom_violin()+
  geom_vline(xintercept = c(0.165, 0.005, 10)) +
#  geom_boxplot() +
  scale_x_log10() #we can assess if there's SC data with median RII similar to empty wells, those need remove
scp_bckp5 <- scp #make a backup of scp object
#scp <- scp_bckp5
scp
#We will analyze both populations of the bi-modal distribution separately, and all together
#scp_upper <- scp[, !is.na(scp$MedianRI) & scp$MedianRI > 0.1 & scp$MedianRI < 9 ]
#scp_upper <- scp[, !is.na(scp$MedianRI) & scp$MedianRI > 0.1]
#scp_lower <- scp[, !is.na(scp$MedianRI) & scp$MedianRI < 0.21 & scp$MedianRI > 0.017]
#scp_lower <- scp[, !is.na(scp$MedianRI) & scp$MedianRI < 0.21]
#scp <- scp[, !is.na(scp$MedianRI) & scp$MedianRI > 0.005 & scp$MedianRI < 10 ]

table(colData(scp)[,"SampleType"])
mean(dims(scp)[1,55:108])
mean(dims(scp)[1,1:54])
dims(scp)[1,109]
scp

table(colData(scp_upper)[,"SampleType"])
mean(dims(scp_upper)[1,55:108])
mean(dims(scp_upper)[1,1:54])
dims(scp_upper)[1,109]
scp_upper

table(colData(scp_lower)[,"SampleType"])
mean(dims(scp_lower)[1,55:108])
mean(dims(scp_lower)[1,1:54])
dims(scp_lower)[1,109]
scp_lower

#Now we plot the filtered subsets
colData(scp) %>%
  data.frame %>%
  ggplot() +
  aes(x = MedianRI, 
      y = SampleType,
      fill = SampleType) +
  geom_violin()+
#  geom_boxplot() +
  scale_x_log10()

colData(scp_upper) %>%
  data.frame %>%
  ggplot() +
  aes(x = MedianRI, 
      y = SampleType,
      fill = SampleType) +
  geom_violin()+
  #  geom_boxplot() +
  scale_x_log10()

colData(scp_lower) %>%
  data.frame %>%
  ggplot() +
  aes(x = MedianRI, 
      y = SampleType,
      fill = SampleType) +
  geom_violin()+
  #  geom_boxplot() +
  scale_x_log10()
scp_bckp6 <- scp
#scp <- scp_bckp6
#join features into proteins for summary
scp_post_RII <- aggregateFeaturesOverAssays(scp,
                                           i = 1:54,
                                           fcol = "Leading.razor.protein", #Column to use for peptide ID
                                           name = paste0("proteins_", names(scp)[1:54]),    #name to create for the aggregated dataset
                                           fun = matrixStats::colMedians, na.rm = TRUE) #function to use when aggregating PSM intensity data
scp_post_RII
mean(dims(scp_post_RII)[1,110:163])

scp_post_RII <- joinAssays(scp_post_RII,
                          i = 110:163,
                          name = "proteins")

scp_post_RII
#Filter by median CV per cell
scp <- medianCVperCell(scp,
                       i = 1:54,
                       groupBy = "Leading.razor.protein",
                       nobs = 4, #only calculate CV if there are x or more peptides per protein
                       norm = "div.median",
                       na.rm = TRUE,
                       colDataName = "MedianCV") #calculates each protein CV across its peptides, then cell median CV
scp_upper <- medianCVperCell(scp_upper,
                       i = 1:54,
                       groupBy = "Leading.razor.protein",
                       nobs = 4, #only calculate CV if there are 5 or more peptides per protein
                       norm = "div.median",
                       na.rm = TRUE,
                       colDataName = "MedianCV") #calculates each protein CV across its peptides, then cell median CV
scp_lower <- medianCVperCell(scp_lower,
                       i = 1:54,
                       groupBy = "Leading.razor.protein",
                       nobs = 4, #only calculate CV if there are 5 or more peptides per protein
                       norm = "div.median",
                       na.rm = TRUE,
                       colDataName = "MedianCV") #calculates each protein CV across its peptides, then cell median CV

#Let's plot the CV distribution
getWithColData(scp, "peptides") %>%
  colData %>%
  data.frame %>%
  ggplot(aes(x = MedianCV,
             y = SampleType,
             fill = SampleType)) +
  geom_violin()+
  geom_boxplot() +
  geom_vline(xintercept = c(0.65, 0.82))

getWithColData(scp_upper, "peptides") %>%
  colData %>%
  data.frame %>%
  ggplot(aes(x = MedianCV,
             y = SampleType,
             fill = SampleType)) +
  geom_violin()+
  geom_boxplot() +
  geom_vline(xintercept = c(0.65, 0.87))

getWithColData(scp_lower, "peptides") %>%
  colData %>%
  data.frame %>%
  ggplot(aes(x = MedianCV,
             y = SampleType,
             fill = SampleType)) +
  geom_violin()+
  geom_boxplot() +
  geom_vline(xintercept = c(0.65, 0.8))
scp_bckp7 <- scp
scp_up_bckp1 <- scp_upper
scp_lo_bckp1 <- scp_lower
scp <- scp_bckp7

######## if you want to test the best CV cutoff on your dataset try the following function#######
source("loopCV.R")
loop_CVs(scp_object = scp,first = 0.85,last = 0.65,step = -0.01)
#################################################################################################

####### If you know the best cutoff, define it as j in the following code ######################
j <- 0.68
data <- scp
#data <- data[, !is.na(data$MedianCV) & data$MedianCV < 0.68, ] #keep only cells with median CV <0.68
data <- data[, !is.na(data$MedianCV) & data$MedianCV < j, ] #keep only cells with median CV < j
data

#how many cells and peptides do we have after filtering
data <- data[, data$SampleType != "Blank", ] #remove the info for empty wells

# 10. Process peptide data

#Normalization of peptide RII

## Divide columns by median (sample-loading normalization)
data <- sweep(data, 
              i = "peptides",
              MARGIN = 2,
              FUN = "/",
              STATS = colMedians(assay(data[["peptides"]]), na.rm = TRUE),
              name = "peptides_norm_col")
## Divide rows by mean (multiplex normalization)
data <- sweep(data,
              i = "peptides_norm_col",
              MARGIN = 1,
              FUN = "/",
              STATS = rowMeans(assay(data[["peptides_norm_col"]]),  na.rm = TRUE),
              name = "peptides_norm")

#Remove peptides with high missing rate
## Replace NaN with NA
peptides_norm <- as.data.frame(assay(data, "peptides_norm"))
assay(data[["peptides_norm"]])[is.nan(assay(data[["peptides_norm"]]))] <- NA
data <- filterNA(data,
                 i = "peptides_norm",
                 pNA = 0.99) #remove peptides with 99% missing data
peptides_norm <- as.data.frame(assay(data, "peptides_norm"))
#log2 data transformation
data <- logTransform(data,
                     base = 2,
                     i = "peptides_norm",
                     name = "peptides_log")

# 11. Aggregate data into proteins
data <- aggregateFeatures(data,
                          i = "peptides_log",
                          name = "proteins",
                          fcol = "Leading.razor.protein",
                          fun = matrixStats::colMedians, na.rm = TRUE)
data
## Replace NaN with NA
proteins <- as.data.frame(assay(data, "proteins"))
assay(data[["proteins"]])[is.nan(assay(data[["proteins"]]))] <- NA
# data <- aggregateFeaturesOverAssays(data,
#                          i = "peptides_log",
#                          name = paste0("POSTproteins_cell_",c(1:75)),
#                          fcol = "Leading.razor.protein",
#                          fun = matrixStats::colMedians, na.rm = TRUE)
# 12. Processing the protein data

# Normalization of protein RII

## Center columns with median (SLN)
data <- sweep(data, i = "proteins",
              MARGIN = 2,
              FUN = "-",
              STATS = colMedians(assay(data[["proteins"]]),
                                 na.rm = TRUE),
              name = "proteins_norm_col")
## Center rows with mean (between plex normalization)
data <- sweep(data, i = "proteins_norm_col",
              MARGIN = 1,
              FUN = "-",
              STATS = rowMeans(assay(data[["proteins_norm_col"]]),
                               na.rm = TRUE),
              name = "proteins_norm")

#Imputation using KNN algorithm
data[["proteins_norm"]] %>%
  assay %>%
  is.na %>%
  mean #this will reveal % of missing values on dataset
#scp_bckp8 <- data

######## Testing removing proteins with >98% missing values #################
data <- filterNA(data,
                 i = "proteins_norm",
                 pNA = 0.98) #remove proteins with 98% missing data
###########################################################################
#data <- scp_bckp8
data[["proteins_norm"]] %>%
  assay %>%
  is.na %>%
  mean #this will reveal % of missing values on dataset

data <- impute(data,
               i = "proteins_norm",
               name = "proteins_imptd",
               method = "knn",
               k = 3, rowmax = 1, colmax= 1, #will impute missing values using 3 neighbors
               maxp = Inf, rng.seed = 1234)

# prot_norm <- getWithColData(data,"proteins_norm")
# library(impute)
# assay(prot_norm) <- impute.knn(data = t(assay(prot_norm)),
#                                  k = 3,
#                                  rowmax = 1,
#                                  colmax = 1,
#                                  maxp = Inf,
#                                  rng.seed = 1234) %>%
#   .[[1]] %>%
#   t()
# data <- addAssay(data,
#                  y = prot_norm,
#                  name = "proteins_imptd") #apply the batch nomr factors on dataset

data[["proteins_imptd"]] %>%
  assay %>%
  is.na %>%
  mean

#Batch correction using ComBat function from "sva" package
sce <- getWithColData(data, "proteins_imptd") #extract the processed data into a new object

batch <- colData(sce)$Raw.file #extract the different run identification
model <- model.matrix(~ SampleType, data = colData(sce)) #create a fitting model (I'm not sure how this works)
assay(sce) <- ComBat(dat = assay(sce),
                     batch = batch,
                     mod = model) #generate batch normalization factors
data <- addAssay(data,
                 y = sce,
                 name = "proteins_batchC") #apply the batch nomr factors on dataset

data <- addAssayLinkOneToOne(data, 
                             from = "proteins_imptd",
                             to = "proteins_batchC") #connect the corrected dataset to it's parent dataset
#plot(data)

# 13. Dimension reduction, PCA and UMAP using R package "scater"

#PCA
data[["proteins_batchC"]] <- runPCA(data[["proteins_batchC"]],
                                    ncomponents = 5,
                                    ntop = Inf,
                                    scale = TRUE,
                                    exprs_values = 1,
                                    name = "PCA") #compute PCA
myPCA <- plotReducedDim(data[["proteins_batchC"]],
                        dimred = "PCA",
                        colour_by = "SampleType",
                        point_alpha = 1) #plot PCA

pdf(paste0("SCPCA_CV",j,"_.pdf"), width = 6, height = 3.5)
print(myPCA)
dev.off()

#UMAP
data[["proteins_batchC"]] <- runUMAP(data[["proteins_batchC"]],
                                     ncomponents = 2,
                                     ntop = Inf,
                                     scale = TRUE,
                                     exprs_values = 1,
                                     n_neighbors = 3,
                                     dimred = "PCA", #use information from previously generated PCA
                                     name = "UMAP") #calculate UMAP
myUMAP <- plotReducedDim(data[["proteins_batchC"]],
                         dimred = "UMAP",
                         colour_by = "SampleType",
                         point_alpha = 1) #plot UMAP

pdf(paste0("SCUMAP_CV",j,"_.pdf"), width = 6, height = 3.5)
print(myUMAP)
dev.off()
data
# Save the results before exporting the data
scp_0.68 <- data
scp_0.7 <- data
scp_0.71 <- data
scp_0.72 <- data
scp_0.73 <- data
scp_1.0 <- data
scp_upper <- data
scp_lower <- data

#################################################################################
##########################Now,let's export the results ##########################
library(stringr)
library(tidyr)
library(tibble)

export_data <- function(myData = data, outname = scp_out, DE_test = ttest) {
  myData <- data
  mySCres <- as.data.frame(assay(myData, "proteins_batchC")) %>%
    rownames_to_column(var = "AGI")
  
  myCells <- data.frame(Raw.file = str_extract(colnames(mySCres), pattern = "11.*QE2-[0-9]{2}"),
                        Channel = str_extract(colnames(mySCres), pattern = "Reporter.intensity.[0-9]{1,2}")) %>%
    left_join(mySCmetadata, by = c("Raw.file","Channel")) %>%
    na.omit
  
  #colnames(mySCres) <- c("AGI", myCells$SampleType)
  mySCres %>%
    `colnames<-`(c("AGI", myCells$SampleType)) %>%
    pivot_longer(!AGI, names_to = "cell_type", values_to = "log2int") %>%
    mutate(AGI = as.factor(AGI), cell_type = as.factor(cell_type)) -> mySCres_long
  write.table(mySCres, "SCP0.68_MeanSCRbelow0.3__rem0.98missing_results_log2Int.csv", row.names = FALSE, sep = ",", quote = FALSE)
  
  mySCres <- data.frame(assay(myData, "proteins_batchC")) %>%
    rownames_to_column(var = "AGI") %>%
    `colnames<-`(c("AGI", myCells$SampleType)) %>%
    tibble(.name_repair = "unique") %>%
    select(AGI, matches("cortex"), matches("endodermis")) -> mySCres_wide
  write.table(mySCres_wide, "SCP0.68_MeanSCRbelow0.3__rem0.98missing_results_log2Int_wider.csv", row.names = FALSE, sep = ",", quote = FALSE)  
  
  ################################################################################
  ########################## Let's do some significance test  #####################
  # first, two-sample t-test
  library(rstatix)
  stat_ttest <- mySCres_long %>%
    group_by(AGI) %>%
    t_test(log2int ~ cell_type,var.equal = FALSE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance(cutpoints = c(0,0.05,0.1,1), symbols = c("**","*","ns")) %>%
    inner_join(mySCres_wide, by = "AGI") %>%
    rename(Protein = AGI) %>%
    mutate(AGI = str_sub(Protein,start = 1, end = 9),
           cortex_AVG = rowMeans(across(starts_with("cortex"))),
           endo_AVG = rowMeans(across(starts_with("endo"))),
           log2FC = c(cortex_AVG-endo_AVG)) %>%
    left_join(AthGeneSymbols%>%select(AGI,symbol), by = "AGI") %>%
    select(1,AGI,symbol,3:6,9:11,log2FC,matches("AVG"))
  stat_ttest
  write.table(stat_ttest, "SCP0.68_MeanSCRbelow0.3__rem0.98missing_results_ttest_BH_welch.csv", row.names = FALSE, sep = ",", quote = FALSE)
  
  #Let's now try PoissonSeq
  library(PoissonSeq)
  pdata <- mySCres_wide %>%
    column_to_rownames(var = "AGI") %>%
    2^.
  y <- c(rep(1,times = length(mySCres_wide%>%select(matches("cortex")))),
         rep(2, times = length(mySCres_wide%>%select(matches("endo")))))
  pseq<-PS.Main(dat=list(n=pdata,
                         y=y,
                         type="twoclass",
                         pair=FALSE,
                         gname=row.names(pdata)),
                para=list(ct.sum=0,ct.mean=0))#,div=15))
  pseq |>
    select(2,4,5) |>
    rename(UID = gname) |>
    inner_join(mySCres_wide, by = c("UID" = "AGI")) %>%
    mutate(GeneID = substr(.[]$UID,start = 1, stop = 9),
           log2FC = rowMeans(across(starts_with("cortex")))-rowMeans(across(starts_with("endo"))),
           adj_pval = p.adjust(p=pval, method = "BH"),
           .after = UID) -> myPseqResults #%>%
  #left_join(y = read.delim(file = paste0(TMT_NEAT_DIR,"/gene_annotations/ath.csv"),header = TRUE,sep = ",",stringsAsFactors = FALSE),
  #          by = c("GeneID" = "AGI"))
  
  write.table(myPseqResults, "SCP_lower_results_Pseq_BH.csv", row.names = FALSE, sep = ",", quote = FALSE)
}





