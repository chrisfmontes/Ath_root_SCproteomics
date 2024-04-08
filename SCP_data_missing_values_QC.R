library(dplyr)
library(magrittr)
library(scpdata)
library("scp")
library(wesanderson)



# loading the data
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

Specht_data <- specht2019v3() %>%
  removeAssay("peptides") %>%
  removeAssay("proteins")
assign(x = "SCP_V3", value = Specht_data)

#Create the SingleCellEperiment object
scp <- readSCP(featureData = mySCdata,
               colData = mySCmetadata,
               channelCol = "Channel",
               batchCol = "Raw.file",
               removeEmptyCols = TRUE)

scp <- selectRowData(scp, c(
  "Sequence", "Leading.razor.protein", "Reverse", 
  "Potential.contaminant", "PEP"
))
SCP_V3 <- selectRowData(SCP_V3, c(
  "Sequence", "Leading.razor.protein", "Reverse", 
  "Potential.contaminant", "PEP"
))

# basic filters on the data
## 1.
scp <- filterFeatures(scp,
                      ~ Reverse != "+" &
                        Potential.contaminant != "+" &
                        PEP < 0.01)
SCP_V3 <- filterFeatures(SCP_V3,
                      ~ Reverse != "+" &
                        Potential.contaminant != "+" &
                        PEP < 0.01)
## 2.
scp <- zeroIsNA(scp, i = 1:nruns)
SCP_V3 <- zeroIsNA(SCP_V3, i = 1:177)

## 3.
table(colData(scp)[,"SampleType"]) # to see different channel types
table(colData(SCP_V3)[,"SampleType"]) # to see different channel types

scp <- subsetByColData( 
  scp, scp$SampleType %in% c("endodermis", "cortex")
)
SCP_V3 <- subsetByColData( 
  SCP_V3, SCP_V3$SampleType %in% c("Monocyte", "Macrophage")
)

## 4.
scp <- filterNA(scp, i = names(scp), pNA = 0.9999) 
scp <- dropEmptyAssays(scp) 
SCP_V3 <- filterNA(SCP_V3, i = names(SCP_V3), pNA = 0.9999) 
SCP_V3 <- dropEmptyAssays(SCP_V3) 

## 5.
scp <- aggregateFeatures( 
  scp, i = names(scp), name = paste0("peptides_", names(scp)),
  fcol = "Sequence", fun = colMedians
)
SCP_V3 <- aggregateFeatures( 
  SCP_V3, i = names(SCP_V3), name = paste0("peptides_", names(SCP_V3)),
  fcol = "Sequence", fun = colMedians
)
## 6.
scp <- joinAssays(
  scp, i = grep("^peptides_", names(scp)), name = "peptides"
)
SCP_V3 <- joinAssays(
  SCP_V3, i = grep("^peptides_", names(SCP_V3)), name = "peptides"
)
# generate a logic data matrix (empty = FALSE, has value = TRUE) to assess sensitivity and missingness
scp_peps <- getWithColData(scp, "peptides")
assay(scp_peps) <- ifelse(is.na(assay(scp_peps)), FALSE, TRUE)
scp <- addAssay(scp, scp_peps, "peps_logic")

V3_peps <- getWithColData(SCP_V3, "peptides")
assay(V3_peps) <- ifelse(is.na(assay(V3_peps)), FALSE, TRUE)
SCP_V3 <- addAssay(SCP_V3, V3_peps, "peps_logic")

# report missing values per cell-type and calculate sensitivity
#reportMissingValues(scp, "peptides", by = scp$SampleType)
reportMissingValues(scp, "peps_logic", by = scp$SampleType)
reportMissingValues(SCP_V3, "peps_logic", by = SCP_V3$SampleType)

scpMissingValues <- reportMissingValues(scp, "peps_logic", by = scp$SampleType)
V3_MissingValues <- reportMissingValues(SCP_V3, "peps_logic", by = SCP_V3$SampleType)

# Jackard index by cell-type (features detected consistently across cells)
ji_scp <- jaccardIndex(scp, "peps_logic", by = scp$SampleType)
ji_V3 <- jaccardIndex(SCP_V3, "peps_logic", by = SCP_V3$SampleType)

#plot the results
library("ggplot2")
Jacc <- ggplot(ji_scp) +
  aes(x = jaccard,
      y = by,
      fill = by) +
  ylab(label = "Cell type")+
  xlab(label = "Jaccard index")+
  labs(fill = "Cell type")+
  theme_light(base_size = 20) +
  geom_violin() +
  scale_fill_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete"))+
  #geom_histogram() +
  #facet_grid(~ by)
  coord_flip()

pdf(file = "Jaccard_index.pdf", width =12, height =8, compress = FALSE)
print(Jacc)
dev.off()

#Cumulative sensitivity curve
csc_scp <- cumulativeSensitivityCurve(scp, "peps_logic", by = scp$SampleType,
                                  batch = scp$Set, niters = 10, 
                                  nsteps = 50)
csc_V3 <- cumulativeSensitivityCurve(SCP_V3, "peps_logic", by = SCP_V3$SampleType,
                                      batch = SCP_V3$Set, niters = 10, 
                                      nsteps = 50)
#Plot the results
(plCSC_scp <- ggplot(csc_scp) +
    aes(x = SampleSize, y = Sensitivity, colour = by) +
    geom_point(size = 3,
               stroke = 0.5))+
  scale_color_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete"))

(plCSC_V3 <- ggplot(csc_V3) +
    aes(x = SampleSize, y = Sensitivity, colour = by) +
    geom_point(size = 3,
               stroke = 0.5))+
scale_color_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete"))

# Add asymptotic regression curve to better predict sensitivity
predCSC_scp <- predictSensitivity(csc_scp, nSample = 1:378)
predCSC_V3 <- predictSensitivity(csc_V3, nSample = 1:378)

cumSens <- plCSC_scp + 
  theme_light(base_size = 20) +
  xlab(label = "Sample size") +
  ylab(label = "Protein detection sensitivity\n(in number of proteins)")+
  labs(color = "Cell type") +
  geom_line(data = predCSC_scp) +
  geom_hline(yintercept = 5322.5,
             lty = 1) + 
  geom_hline(yintercept = 1098.5025,
             lty = 4)+
  scale_y_continuous(breaks = c(seq(from = 0, to = 5000, by =1000))) +
  scale_color_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete"))

pdf(file = "Cumulative_Sensitivity.pdf", width = 12, height = 8, compress = FALSE)
print(cumSens)
dev.off()

plCSC_V3 + 
  theme_light(base_size = 20) +
  geom_line(data = predCSC_V3) +
  geom_hline(yintercept = 5161.5,
             lty = 1) + 
  geom_hline(yintercept = 943.48,
             lty = 4)+
  scale_color_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete"))

# Predict sensitivity if we had infinite cells
predictSensitivity(csc_scp, nSamples = Inf)
predictSensitivity(csc_V3, nSamples = Inf)

#Plot total sensitivity
bothMissingValues <- rbind(scpMissingValues,V3_MissingValues) |>
  mutate(Dataset = factor(c("Montes", "Montes", "Specht", "Specht")))

scpMissingValues <- mutate(scpMissingValues, CellType = factor(rownames(scpMissingValues)))

TotSens <- ggplot(data = scpMissingValues,
       aes(y = TotalSensitivity,
           x = rownames(scpMissingValues),
           fill = rownames(scpMissingValues))) +
  xlab(label = "Cell type")+
#  ylab(label = "Total sensitivity\n(number of proteins)")+
  labs(fill = "Cell type")+
#  ylim(0,6000) +
  theme_light(base_size = 20) +
  geom_bar(stat = "identity")+
  scale_y_continuous(name = "Total sensitivity\n(number of proteins)",
                     breaks = seq(from = 0, to = 6000, by = 1000),
                     limits = c(0,6000))+
  scale_fill_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete"))

pdf(file = "Total_sensitivity.pdf", width = 12, height = 8, compress = FALSE)
print(TotSens)
dev.off()

ggplot(data = bothMissingValues,
       aes(y = TotalSensitivity,
           x = rownames(bothMissingValues),
           fill = rownames(bothMissingValues))) +
  ylim(0,6000) +
  theme_light() +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete"))

#Plot local sensitivity

LocSens <- ggplot(data = scpMissingValues,
       aes(y = LocalSensitivityMean,
           x = Completeness*100))+
  geom_point(aes(size = NumberCells, color = CellType)) +
  theme_light(base_size = 20)+
  scale_size_continuous(limits = c(1,700), breaks = c(100, 300, 700)) +
  xlim(0,25)+
  ylim(0,2000)+
  geom_errorbar(aes(
    ymin = LocalSensitivityMean - LocalSensitivitySd,
    ymax = LocalSensitivityMean + LocalSensitivitySd,
    color = CellType))+
  geom_point(data = data.frame(Completeness = c(mean(scpMissingValues$Completeness)),
                               LocalSensitivityMean = c(mean(scpMissingValues$LocalSensitivityMean))),
             size =10,
             shape = 5,
             color = "#5B1A18")+
  xlab(label = "Completeness")+
  ylab(label = "Local sensitivity\n(number of proteins)")+
  labs(size = "Number of cells", color = "Cell type")+
  scale_color_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete"))

pdf(file = "Local_sensitivity.pdf", width = 12, height = 8, compress = FALSE)
print(LocSens)
dev.off()

ggplot(data = bothMissingValues,
       aes(y = LocalSensitivityMean,
           x = Completeness*100,
           color = Dataset))+
  geom_point(aes(size = NumberCells)) +
  theme_classic()+
  scale_size_continuous(limits = c(1,1500), breaks = c(300, 1000, 1500)) +
  xlim(0,25)+
  ylim(0,2000)+
  geom_errorbar(aes(
    ymin = LocalSensitivityMean - LocalSensitivitySd,
    ymax = LocalSensitivityMean + LocalSensitivitySd))+
  geom_point(data = data.frame(Completeness = c(mean(bothMissingValues$Completeness[1:2]),
                                                mean(bothMissingValues$Completeness[3:4])),
                               LocalSensitivityMean = c(mean(bothMissingValues$LocalSensitivityMean[1:2]),
                                                          mean(bothMissingValues$LocalSensitivityMean[3:4])),
                               Dataset = factor(c("Montes","Specht"))),
    size =10,
    shape = 5
  )+
  scale_color_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete"))
