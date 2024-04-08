library(tibble)
library(dplyr)
library(viridis)
library(ComplexHeatmap)
library(dendsort)

myDE <- read.csv("SCP0.68_MeanSCRbelow0.3__rem0.98missing_results_ttest_BH_welch.csv",
                 header = TRUE,
                 stringsAsFactors = FALSE) |>
  filter(p.adj < 0.1) |>
  select(Protein)

mySC_LogDEMatrix <- mygnames |>
  inner_join(myDE, by = c("AGI" = "Protein")) |>
  distinct(AGI, .keep_all = TRUE) |>
  column_to_rownames(var = "AGI") |>
  as.matrix()

myNormLogSCDEMatrix <- t(scale(x = t(mySC_LogDEMatrix)))

hclog <- hclust(as.dist(1-cor(t(mySC_logMatrix), method = "pearson")))

dendDEl <- dendsort(hcDEl)

hm <- Heatmap(myNormLogSCDEMatrix,
              name = "Abundance\n(Z-score)",
              col = inferno(3),
              border = T,
              row_title = "Proteins",
              column_title = "Cell type",
              column_title_side = "bottom",
              show_column_dend = T,
              column_dend_height = unit(2, "cm"),
              column_names_gp = gpar(fontsize = 10),
              row_dend_width = unit(6, "cm"),
              show_row_names = F,
              cluster_rows = dendDEl)

pdf(file = "HM_DEzNorm.pdf",width = 13, height = 10,compress = FALSE)
draw(hm)
dev.off()