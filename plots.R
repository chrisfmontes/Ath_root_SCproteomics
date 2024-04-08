library(ggplot2)
library(ggrepel)

myData <- read.csv("SCP0.68_MeanSCRbelow0.3__rem0.98missing_results_ttest_BH_welch.csv",
                   header = TRUE,
                   stringsAsFactors = FALSE)
myData$DE <- "NS"
myData$DE[myData$p.adj<0.1 & myData$log2FC < 0] <- "Endodermis"
myData$DE[myData$p.adj<0.1 & myData$log2FC > 0] <- "Cortex"
myData$DE <- factor(myData$DE, levels = c("Cortex", "Endodermis", "NS"))

MAplot <- ggplot(myData,
                 aes(x = (cortex_AVG+endo_AVG)*0.5,
                     y = log2FC,
                     label = Selected))+
  geom_point(color=alpha('black', 0.3), shape=21, size=2, aes(fill=factor(DE))) +
  scale_fill_manual(values=alpha(c('#FD6467','#56B4E9','#FDFD96'), 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = alpha("black", 0.7), linewidth = 1)+
  geom_text_repel(
    size = 3,
    max.overlaps = 100,
    box.padding = 0.35,
    point.padding = 1,
    min.segment.length = 0)+
  xlab(label = bquote("A (Average"~log[2]~"Intensity)"))+
  ylab(label = bquote("M ("*log[2]~"Fold-change)"))+
  labs(fill = element_blank())+
  theme_classic()+
  theme(axis.title.x = element_text(size =20),
        axis.title.y = element_text(size =20),
        legend.text = element_text(size=15),
        axis.text.x= element_text(size=15),
        axis.text.y= element_text(size=15))
print(MAplot)

pdf(file = "MA_plot.pdf", width = 10, height = 6,compress = FALSE)
print(MAplot)
dev.off()
