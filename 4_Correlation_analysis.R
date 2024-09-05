library("ChIPseeker")
library("org.Hs.eg.db")
library("RMariaDB")

#AnnotationHub
library("AnnotationHub")

#txdb <- makeTxDbFromEnsembl("Homo_sapiens",release=109)

AnnotationHub <- AnnotationHub()
query(AnnotationHub, c("Gencode", "TxDb", "Homo sapiens", "hg38")) 
txdb <- AnnotationHub[['AH75143']] #gencode v24

#enhancer.bed
setwd("/newdisk/jinyuan/r/Chipseeker/enhancer")
enhancer_region <- readPeakFile("encode_ELS_hg38.bed")


enhancer_anno <- annotatePeak(enhancer_region, 
                         tssRegion=c(-2000, 2000),
                         TxDb=txdb, 
                         #genomicAnnotationPriority = c("Promoter", "Intergenic"), #c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic")
                         addFlankGeneInfo=TRUE, 
                         flankDistance=3000, 
                         ignoreOverlap=FALSE, 
                         sameStrand=FALSE)

anno_df <- as.data.frame(enhancer_anno)

write.csv(anno_df, file = "enhancer_annotation.csv", row.names = FALSE)


cols <- c("#6baed6", "#ffc107", "#74c476", "#e31a1c", "#fdbb84", "#a6bddb", "#f03b20", "#a1d99b", "#ffff33", "#fd8d3c", "#e6550d", "#bcbddc")

par(mfrow=c(1,2))
plotAnnoPie(enhancer_anno, piecol=cols, legendpos="topright", main="enhancer Region Annotation", fontface="bold", line=-1)

#----------------------------------------------------------------------------------------

tss <- promoters(txdb, upstream=2000, downstream=2000, use.names=TRUE)
tss_df <- as.data.frame(tss)
write.csv(tss_df, file = "GENCODE_tss_annotation.csv", row.names = FALSE)


library(ggplot2)
library(ggthemes)
setwd("/newdisk/jinyuan/r/linear_regression/promoter")
num_df  <- read.csv("RNAvsATAC.csv" , header = T, row.names = 1)
#RNA-ATAC-8h vs 0h------------------------------------------------------------------------
correlation <- cor.test(num_df$log2_0.8.RNA, num_df$log2_0.8.ATAC, method = "pearson")
p = ggplot(num_df, aes(x=log2_0.8.RNA, y=log2_0.8.ATAC)) +
  geom_point(color = "blue") + 
  geom_smooth(method = lm, formula = y ~ x, color = "red") +
  annotate(geom = "text", x = max(num_df$log2_0.8.RNA)*0.9, y = max(num_df$log2_0.8.ATAC)*0.9, 
           label = paste0("R =", round(correlation$estimate, 2), "; p =", round(correlation$p.value, 4)), 
           hjust = 1, vjust = 1)+
  labs(title = "Correlation of expression and ATAC (8h vs 0h)", x = "expression: 8h vs 0h, log2 ratio", y = "ATAC: 8h vs 0h, log2 ratio") +
  theme(axis.title.x = element_text(family = "timesnewroman", size = 14),
        axis.title.y = element_text(family = "timesnewroman", size = 14),
        plot.title = element_text(family = "timesnewroman", size = 18, face = "bold"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )
ggsave("RNA-ATAC-8hvs0h.jpeg", plot = p, dpi = 300)
#RNA-ATAC-16h vs 0h------------------------------------------------------------------------
correlation <- cor.test(num_df$log2_0.16.RNA, num_df$log2_0.16.ATAC, method = "pearson")
p = ggplot(num_df, aes(x=log2_0.16.RNA, y=log2_0.16.ATAC)) +
  geom_point(color = "blue") +
  geom_smooth(method = lm, formula = y ~ x, color = "red") +
  annotate(geom = "text", x = max(num_df$log2_0.16.RNA)*0.9, y = max(num_df$log2_0.16.ATAC)*0.9, 
           label = paste0("R =", round(correlation$estimate, 2), "; p =", round(correlation$p.value, 4)), 
           hjust = 1, vjust = 1)+
  labs(title = "Correlation of expression and ATAC (16h vs 0h)", x = "expression: 16h vs 0h, log2 ratio", y = "ATAC: 16h vs 0h, log2 ratio") +
  theme(axis.title.x = element_text(family = "timesnewroman", size = 14),
        axis.title.y = element_text(family = "timesnewroman", size = 14),
        plot.title = element_text(family = "timesnewroman", size = 18, face = "bold"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )
ggsave("RNA-ATAC-16hvs0h.jpeg", plot = p, dpi = 300)

#RNA-ChIP-8h vs 0h------------------------------------------------------------------------
num_df  <- read.csv("RNAvsChIP.csv" , header = T, row.names = 1)
correlation <- cor.test(num_df$log2_0.8.RNA, num_df$log2_0.8.ChIP, method = "pearson")
p = ggplot(num_df, aes(x=log2_0.8.RNA, y=log2_0.8.ChIP)) +
  geom_point(color = "blue") +
  geom_smooth(method = lm, formula = y ~ x, color = "red") +
  annotate(geom = "text", x = max(num_df$log2_0.8.RNA)*0.9, y = max(num_df$log2_0.8.ChIP)*0.9, 
           label = paste0("R =", round(correlation$estimate, 2), "; p =", round(correlation$p.value, 4)), 
           hjust = 1, vjust = 1)+
  labs(title = "Correlation of expression and ChIP (8h vs 0h)", x = "expression: 8h vs 0h, log2 ratio", y = "ChIP: 8h vs 0h, log2 ratio") +
  theme(axis.title.x = element_text(family = "timesnewroman", size = 14),
        axis.title.y = element_text(family = "timesnewroman", size = 14),
        plot.title = element_text(family = "timesnewroman", size = 18, face = "bold"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )
ggsave("RNA-ChIP-8hvs0h.jpeg", plot = p, dpi = 300)
#RNA-ChIP-16h vs 0h------------------------------------------------------------------------
correlation <- cor.test(num_df$log2_0.16.RNA, num_df$log2_0.16.ChIP, method = "pearson")
p = ggplot(num_df, aes(x=log2_0.16.RNA, y=log2_0.16.ChIP)) +
  geom_point(color = "blue") + 
  geom_smooth(method = lm, formula = y ~ x, color = "red") +
  annotate(geom = "text", x = max(num_df$log2_0.16.RNA)*0.9, y = max(num_df$log2_0.16.ChIP)*0.9, 
           label = paste0("R =", round(correlation$estimate, 2), "; p =", round(correlation$p.value, 4)), 
           hjust = 1, vjust = 1)+
  labs(title = "Correlation of expression and ChIP (16h vs 0h)", x = "expression: 16h vs 0h, log2 ratio", y = "ChIP: 16h vs 0h, log2 ratio") + 
  theme(axis.title.x = element_text(family = "timesnewroman", size = 14),
        axis.title.y = element_text(family = "timesnewroman", size = 14),
        plot.title = element_text(family = "timesnewroman", size = 18, face = "bold"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )
ggsave("RNA-ChIP-16hvs0h.jpeg", plot = p, dpi = 300)

#ChIP-ATAC-8h vs 0h------------------------------------------------------------------------
num_df  <- read.csv("ChIPvsATAC.csv" , header = T, row.names = 1)
correlation <- cor.test(num_df$log2_0.8.ChIP, num_df$log2_0.8.ATAC, method = "pearson")
p = ggplot(num_df, aes(x=log2_0.8.ChIP, y=log2_0.8.ATAC)) +
  geom_point(color = "blue") +
  geom_smooth(method = lm, formula = y ~ x, color = "red") +
  annotate(geom = "text", x = max(num_df$log2_0.8.ChIP)*0.9, y = max(num_df$log2_0.8.ATAC)*0.9, 
           label = paste0("R =", round(correlation$estimate, 2), "; p =", round(correlation$p.value, 4)), 
           hjust = 1, vjust = 1)+
  labs(title = "Correlation of ChIP and ATAC (8h vs 0h)", x = "ChIP: 8h vs 0h, log2 ratio", y = "ATAC: 8h vs 0h, log2 ratio") +
  theme(axis.title.x = element_text(family = "timesnewroman", size = 14),
        axis.title.y = element_text(family = "timesnewroman", size = 14),
        plot.title = element_text(family = "timesnewroman", size = 18, face = "bold"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1) 
  )
ggsave("ChIP-ATAC-8hvs0h.jpeg", plot = p, dpi = 300)
#ChIP-ATAC-16h vs 0h------------------------------------------------------------------------
correlation <- cor.test(num_df$log2_0.16.ChIP, num_df$log2_0.16.ATAC, method = "pearson")
p = ggplot(num_df, aes(x=log2_0.16.ChIP, y=log2_0.16.ATAC)) +
  geom_point(color = "blue") + 
  geom_smooth(method = lm, formula = y ~ x, color = "red") +
  annotate(geom = "text", x = max(num_df$log2_0.16.ChIP)*0.9, y = max(num_df$log2_0.16.ATAC)*0.9, 
           label = paste0("R =", round(correlation$estimate, 2), "; p =", round(correlation$p.value, 4)), 
           hjust = 1, vjust = 1)+
  labs(title = "Correlation of ChIP and ATAC (16h vs 0h)", x = "ChIP: 16h vs 0h, log2 ratio", y = "ATAC: 16h vs 0h, log2 ratio") +
  theme(axis.title.x = element_text(family = "timesnewroman", size = 14),
        axis.title.y = element_text(family = "timesnewroman", size = 14),
        plot.title = element_text(family = "timesnewroman", size = 18, face = "bold"), 
        panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
ggsave("ChIP-ATAC-16hvs0h.jpeg", plot = p, dpi = 300)