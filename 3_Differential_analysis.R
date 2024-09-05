##DESeq2
BiocManager::install(c("R.utils","DESeq2","tidyverse","biomaRt","curl"))
library("R.utils")
library("DESeq2")
library("tidyverse")
setwd("/newdisk/jinyuan/r/RNA_seq/deseq2_edgeR")
#0-8h---------------------------------------------------------------------------------
mycounts <- read.csv("enhancer_reads.csv",header=TRUE,row.names=1,check.names = FALSE)

names(mycounts) <- c("Group_0h-1","Group_0h-2","Group_0h-3","Group_8h-1","Group_8h-2","Group_8h-3","Group_16h-1","Group_16h-2","Group_16h-3") 
head(mycounts)

mycounts_new <- mycounts[,c(1:3,4:6)]
head(mycounts_new)
##factor()
condition <- factor(c(rep("Group_0h",3),rep("Group_8h",3)), levels = c("Group_0h","Group_8h"))

colData <- data.frame(row.names=colnames(mycounts_new), condition)

mycounts_new[is.na(mycounts_new)] <- 0

dds <- DESeqDataSetFromMatrix(countData = mycounts_new, colData = colData, design= ~ condition)
dds <- DESeq(dds)

res = results(dds, contrast=c("condition", "Group_8h", "Group_0h")) 
write.csv(res,file="enhancer-0h-8h_results.csv")


##pvalue
res = res[order(res$pvalue),]
##table()
table(res$padj<0.05)
##subset()
diff_gene_deseq2 <-subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
##dim()
dim(diff_gene_deseq2)
library('biomaRt')
library("curl") 
##biomaRt
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(diff_gene_deseq2)
##getBM()
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
ensembl_gene_id<-rownames(diff_gene_deseq2)
## cbind()
diff_gene_deseq2<-cbind(ensembl_gene_id,diff_gene_deseq2)
colnames(diff_gene_deseq2)[1]<-c("ensembl_gene_id")
##merge()
diff_name<-merge(diff_gene_deseq2,mms_symbols,by="ensembl_gene_id")
write.csv(diff_name,file="0h-8h-id.csv")
#8-16h--------------------------------------------------------------------------------------------------

mycounts_new <- mycounts[,c(4:6,7:9)]
head(mycounts_new)
##factor()
condition <- factor(c(rep("Group_8h",3),rep("Group_16h",3)), levels = c("Group_8h","Group_16h"))

colData <- data.frame(row.names=colnames(mycounts_new), condition)

mycounts_new[is.na(mycounts_new)] <- 0

dds <- DESeqDataSetFromMatrix(countData = mycounts_new, colData = colData, design= ~ condition)
dds <- DESeq(dds)

res = results(dds, contrast=c("condition", "Group_16h", "Group_8h")) 
write.csv(res,file="enhancer-8h-16h_results.csv")
#0-16h---------------------------------------------------------------------------------------------------

mycounts_new <- mycounts[,c(1:3,7:9)]
head(mycounts_new)

condition <- factor(c(rep("Group_0h",3),rep("Group_16h",3)), levels = c("Group_0h","Group_16h"))

colData <- data.frame(row.names=colnames(mycounts_new), condition)

mycounts_new[is.na(mycounts_new)] <- 0

dds <- DESeqDataSetFromMatrix(countData = mycounts_new, colData = colData, design= ~ condition)
dds <- DESeq(dds) 

res = results(dds, contrast=c("condition", "Group_16h", "Group_0h")) 
write.csv(res,file="enhancer-0h-16h_results.csv")


library("edgeR")
setwd("/newdisk/jinyuan/r/RNA_seq/deseq2_edgeR")
#0-8h---------------------------------------------------------------------
rawdata <- read.csv("ATAC_reads.csv",header=TRUE,row.names=1,check.names = FALSE) 
head(rawdata) 

names(rawdata) <- c("Group_0h","Group_8h","Group_16h") 

group <- factor(c("Group_0h","Group_8h")) 

y <- DGEList(counts=rawdata[,c(1,2)],genes=rownames(rawdata),group = group) 
###TMM
y<-calcNormFactors(y) 
y$samples 

#bcv <- 0.1 
#bcv <- 0.2 
bcv <- 0.4 
et <- exactTest(y, dispersion=bcv^2) 
topTags(et) 
summary(de <- decideTestsDGE(et)) 

DE <- et$table 
write.csv(DE, "ATAC-0-8h-result.csv")

#8-16h---------------------------------------------------------------------------

group <- factor(c("Group_8h","Group_16h")) 

y <- DGEList(counts=rawdata[,c(2,3)],genes=rownames(rawdata),group = group) 
###TMM
y<-calcNormFactors(y) 
y$samples 

#bcv <- 0.1 
#bcv <- 0.2 
bcv <- 0.4 
et <- exactTest(y, dispersion=bcv^2) 
topTags(et) 
summary(de <- decideTestsDGE(et)) 

DE <- et$table 
write.csv(DE, "ATAC-8-16h-result.csv")

#0-16h---------------------------------------------------------------------------

group <- factor(c("Group_0h","Group_16h")) 

y <- DGEList(counts=rawdata[,c(1,3)],genes=rownames(rawdata),group = group) 
###TMM
y<-calcNormFactors(y) 
y$samples 

#bcv <- 0.1 
#bcv <- 0.2 
bcv <- 0.4 
et <- exactTest(y, dispersion=bcv^2) 
topTags(et) 
summary(de <- decideTestsDGE(et)) 

DE <- et$table 
write.csv(DE, "ATAC-0-16h-result.csv")


setwd("/newdisk/jinyuan/r/padj")
pvalue = read.csv("pvalue.csv",header=TRUE,row.names=NULL,check.names = FALSE)
names(pvalue) <- c("C0-8","C8-16","C0-16","A0-8","A8-16","A0-16") 
df = pvalue
for ( i in c(1:6))
{
  padj = p.adjust(c(pvalue[ ,i]),method ="BH")
  df = data.frame(df,padj)
}
write.csv(df,"padj.csv")
#padj = p.adjust(c(pvalue[ ,1]),method ="bonferroni")


#pheatmap
setwd("/newdisk/jinyuan/r/RNA_seq/pheatmap")

df <- read.csv("alway_updown_chayi_739_190.csv" , header = T, row.names = 1)

pheatmap(df[,1:3], scale = "row", 
             color = colorRampPalette(c("#3951A2", "#D9D6D6", "#A80326"))(100), 
             cluster_cols = F, cluster_rows = TRUE, clustering_method = "ward.D2", 
             cutree_rows = 2 ,#cutree_cols = 1 ,
             #legend_breaks = c(-2, 0, 2), 
             legend_labels = c("low", "median", "high"), 
             border = FALSE, show_rownames = FALSE, angle_col = 45, 
             #main = "Heatmap of eRNA expression TPM values",
             filename = "eRNA_heatmap.jpeg", width = 10,height = 8)


row_cluster <- cutree(p$tree_row, k=1)
newOrder <- df[p$tree_row$order,]
newOrder[,ncol(newOrder)+1]=row_cluster[match(rownames(newOrder), names(row_cluster))]
colnames(newOrder)[ncol(newOrder)]="Cluster"
write.csv(newOrder, "cluster_order.csv")


clustered_df <- read.csv("RNA_upandchayi.csv" , header = T, row.names = 1)

z_score_RNA <- t(apply(clustered_df[,1:3], 1, function(x) (x-mean(x))/sd(x)))
z_score_Chip <- t(apply(clustered_df[,4:6], 1, function(x) (x-mean(x))/sd(x)))
z_score_ATAC <- t(apply(clustered_df[,7:9], 1, function(x) (x-mean(x))/sd(x)))
z_scored_df = cbind(z_score_RNA, z_score_Chip, z_score_ATAC)


pheatmap(z_scored_df[,c(1:9)], scale = "none", #kmeans_k = 4,
             color = colorRampPalette(c("#3951A2", "#D9D6D6", "#A80326"))(100), 
             cluster_cols = F, cluster_rows = T, clustering_method = "ward.D2", 
             cutree_rows = 1 ,cutree_cols = 3 ,
             #legend_breaks = c(-2, 0, 2), 
             legend_labels = c("low", "median", "high"), 
             border = FALSE, show_rownames = FALSE, angle_col = 45, 
             #main = "Heatmap of eRNA expression TPM values",
             filename = "heatmap.jpeg", width = 10,height = 8)

######################################
#ClusterProfiler Go KEGG
######################################


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("readxl")
library("ggprism")
library("ReactomePA")
library("DOSE")
library("dplyr")


setwd(dir = "/newdisk/jinyuan/r/RNA_seq/ClusterProfiler")

EnsembleID <- read_excel("gene_id.xlsx",col_names = TRUE)


#'ENTREZID','ENSEMBL','REFSEQ'
geneList=bitr(EnsembleID$Gene_ID,
              fromType="ENSEMBL",
              toType=c("ENTREZID"),
              OrgDb="org.Hs.eg.db") 

#Go
ego2 <- enrichGO(gene          = geneList$ENTREZID,
                 keyType       = "ENTREZID",
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",       #可选"CC"、"BP"、"MF"、"ALL"
                 maxGSSize     = 500,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05 )
#KEGG
kk2 <- enrichKEGG(gene    = geneList$ENTREZID,
                  organism    = "hsa",
                  maxGSSize    = 500,
                  keyType      = "kegg",
                  pvalueCutoff = 0.05,
                  qvalueCutoff  = 0.05)

theme=theme_prism(base_size = 17,base_fontface = "plain",base_rect_size = 15)

barplot(ego2, showCategory=15, horiz=FALSE)+theme

ggsave("goBP_15.png", device = 'png',units = "cm",width = 30,height = 12,dpi = 300)

barplot(kk2, showCategory=10)+theme

#mutate(kk2, qscore = -log(p.adjust, base=10)) %>% 
#  barplot(x="qscore")

dotplot(kk2, showCategory=10) +theme
dotplot(ego2, showCategory=10) +theme


write.csv(x = geneList, file = "birt_id.csv")
write.csv(x = kk2,file = "kegg.csv")
write.csv(x = ego2,file="goBP.csv")



###############
######GSEA#####
###############
setwd("/newdisk/jinyuan/r/RNA_seq/ClusterProfiler/GSEA")

EnsembleID = read_excel("929_always_0-16h.xlsx",col_names = TRUE)
head(EnsembleID)


#'ENTREZID','ENSEMBL','REFSEQ'
geneList=bitr(EnsembleID$Gene_ID,
              fromType="ENSEMBL",
              toType=c("ENTREZID"),
              OrgDb="org.Hs.eg.db") 
head(geneList)


colnames(EnsembleID)[1]="ENSEMBL"
ID_bitr <- EnsembleID %>% 
  inner_join(geneList,by="ENSEMBL")
head(ID_bitr)

EnsembleID_sort <- ID_bitr %>% 
  arrange(desc(logFC))


GSEA = EnsembleID_sort[,c(3,2)]
GSEA <- na.omit(GSEA)
GSEA = as.data.frame(GSEA)
head(GSEA)
geneList = GSEA$logFC 
names(geneList) <- GSEA$ENTREZID 

# GSEA_GO
gsea <- gseGO(geneList, 
              keyType = "ENTREZID",
              ont = "BP",
              OrgDb = "org.Hs.eg.db",
              minGSSize = 10,
              maxGSSize = 500,
              eps = 1e-10, 
              pAdjustMethod = "BH",
              verbose = TRUE,
              seed = FALSE,
              by = "fgsea",
              pvalueCutoff = 1)

write.csv(x = gsea, file = "GSEA_GOBP.csv")
write.csv(x = ID_bitr, file = "ID_bitr.csv")

# visualize
p6<-gseaplot2(gsea,56,
              title = gsea$Description[56],
              color = "#FA5860",
              base_size = 12, 
              rel_heights = c(1.5, 0.5, 1), 
              subplots = 1:3,
              pvalue_table = F, 
              ES_geom = "line" 
)
p6
ggsave("GSEA.png",units = "cm",height = 12,width = 20)




KEGG_database="hsa"
gsea <- gseKEGG(geneList, organism = KEGG_database, pvalueCutoff = 0.05, minGSSize = 10,maxGSSize = 500,eps = 0, pAdjustMethod = "BH")

kegg_gmt <- read.gmt('...txt')
gsea <- GSEA(geneList,  
             TERM2GENE = kegg_gmt, 
             pvalueCutoff = 1,  
             pAdjustMethod = 'BH')


#####################
#####volcano_plot######
#####################
library(readxl)
library(ggplot2)
library(ggprism)


setwd("/newdisk/jinyuan/r/RNA_seq/Volcano_plot/gene")
dataset <- read_excel("all_gene_0-16h_gai.xlsx")


dataset$information <- factor(x=dataset$information,levels = c("stable","up_always","down_always"))


head(dataset)

ggplot(
  dataset, aes(x = logFC, y = Padj_log10, colour=information)) +
  geom_point(alpha=0.8, size=2.5) +
  scale_color_manual(values=c("#d2dae2","#e41a1c","#377eb8"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_prism(base_size = 15, base_line_size = 0.8,base_fontface = "plain")+
  labs(x="log2(fold change)",
       y="-log10 (p-adjust)")+
  theme(plot.title = element_text(), 
        legend.position="top", 
        legend.title = element_blank())+
  scale_y_continuous(expand = c(0, 0),limits = (c(0,300)))
  

ggsave("volcano_plot0-16h_gai.png",device = "png",units = "cm",width = 18,height = 15)


