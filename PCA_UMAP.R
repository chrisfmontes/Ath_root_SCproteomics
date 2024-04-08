library(tibble)
library(dplyr)
library(scater)
library(wesanderson)

## Let's make a subset of only top N most diverse proteins
#read our DE results table

read.csv("DE_results/SCP0.68_MeanSCRbelow0.3__rem0.98missing_results_ttest_BH_welch.csv",
         header = TRUE, stringsAsFactors = FALSE) |>
  slice_max(order_by = abs(log2FC), n = 775) |>
  select(1) -> diverse_proteins

#Subset the scp dataset
mostDiverse <- data[unlist(diverse_proteins),,]
mostDiverse

#PCA
mostDiverse[["proteins_batchC"]] <- runPCA(mostDiverse[["proteins_batchC"]],
                                           ncomponents = 5,
                                           ntop = Inf,
                                           scale = TRUE,
                                           exprs_values = 1,
                                           name = "PCA") #compute PCA
myPCA <- plotReducedDim(mostDiverse[["proteins_batchC"]],
                        dimred = "PCA",
                        colour_by = "SampleType",
                        point_alpha = 0.9,
                        point_size = 5) #plot PCA
plotPCA <- myPCA +
  theme_classic(base_size = 20) +
  scale_color_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete")) +
  labs(color = "Cell type")

pdf(file = "PCA_top25percMostDiverse.pdf", width = 12, height = 8, compress = FALSE)
print(plotPCA)
dev.off()

#UMAP
mostDiverse[["proteins_batchC"]] <- runUMAP(mostDiverse[["proteins_batchC"]],
                                            ncomponents = 2,
                                            ntop = Inf,
                                            scale = TRUE,
                                            exprs_values = 1,
                                            n_neighbors = 3,
                                            dimred = "PCA", #use information from previously generated PCA
                                            name = "UMAP") #calculate UMAP
myUMAP <- plotReducedDim(mostDiverse[["proteins_batchC"]],
                         dimred = "UMAP",
                         colour_by = "SampleType",
                         point_alpha = 0.9,
                         point_size = 5) #plot UMAP
plotUMAP <- myUMAP +
  theme_classic(base_size = 20) +
  scale_color_manual(values = wes_palette(name = "GrandBudapest1",type = "discrete")) +
  labs(color = "Cell type")

pdf("UMAP_top25percMostDiverse.pdf", width = 12, height = 8, compress = FALSE)
print(plotUMAP)
dev.off()
