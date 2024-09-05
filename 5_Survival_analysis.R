##enhancer_target_gene
setwd("/newdisk/jinyuan/r/RNA_seq/target_gene")
eRNA = read.csv("eRNA_RPM.csv",header=1)
gene = read.csv("gene_TPMdayu1.csv",header=1)
data = data.frame()
nrow(eRNA)
nrow(gene)

for ( gene_number in c(1:12860)) 
{
  for (eRNA_number in c(1:5316))
  {
    gene_chr = gene[gene_number,3]
    gene_start = gene[gene_number,4]
    gene_end = gene[gene_number,5]
    gene_strand = gene[gene_number,6]
    gene_TPM = gene[gene_number,7:15]
    
    eRNA_chr = eRNA[eRNA_number,2]
    eRNA_mid = eRNA[eRNA_number,5]
    eRNA_RPM = eRNA[eRNA_number,6:14]

    if(gene_strand == "-"){gene_TSS <- gene_end}else{gene_TSS <- gene_start}
    ab = abs(eRNA_mid - gene_TSS)
    if (eRNA_chr == gene_chr)
    {
      if (ab <= 1000000)
      {
        ##spearman
        cor = cor.test(t(eRNA_RPM),t(gene_TPM),method = "spearman")
        rho = cor$estimate
        pvalue = cor$p.value
        if (complete.cases(pvalue))
          ##&& 0.3 <= abs(rho) && pvalue <0.05)
        {
          mer = cbind(gene[gene_number,1:5],eRNA[eRNA_number,1:5],rho,pvalue)
          data = rbind(data,mer)
        }
      }
    }
  }
}

padj = p.adjust(c(data[,12]),method ="BH")
df = data.frame(data,padj)
write.csv(df,"eRNA_target_gene.csv")


library("survival")
library("survminer")
library("UCSCXenaTools")
library("dplyr")
setwd("/newdisk/jinyuan/r/KM_survival")

##XenaData$XXXX
as.data.frame(tail(sort(table(XenaData$XenaHostNames)))) 


LIHC_cohort = XenaData %>% 
  filter( XenaHostNames == "tcgaHub") %>% # select TCGA Hub
  XenaScan("TCGA Liver Cancer")   # select LIHC cohort
LIHC_cohort


cli_query = LIHC_cohort %>% 
  filter(DataSubtype == "phenotype") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>% 
  XenaDownload()
cli = XenaPrepare(cli_query)

cli = cli$LIHC_survival.txt
head(cli)


ge = LIHC_cohort %>% 
  filter(DataSubtype == "gene expression RNAseq", Label == "IlluminaHiSeq")

ge
#---------------------------
gene_name = "CDKN1A"

gene_expression = fetch_dense_values(host = ge$XenaHosts, 
                          dataset = ge$XenaDatasets, 
                          identifiers = gene_name,
                          use_probeMap = TRUE) %>% 
  .[1, ]
head(gene_expression)


merged_data = tibble(sample = names(gene_expression),
                     expression = as.numeric(gene_expression)) %>% 
  left_join(cli, by = "sample") %>% 
  filter(substr(sample, 14, 15) == "01") %>%  # Keep only 'Primary Tumor'
  select(sample, expression, OS.time, OS) %>% 
  rename(time = OS.time, 
         event = OS)
head(merged_data)


#surv_cutpoint
res.cut <- surv_cutpoint(merged_data, time = "time", event = "event",
                         variables = "expression")
summary(res.cut)

plot(res.cut, "expression", palette = "npg")


res.cat <- surv_categorize(res.cut)
head(res.cat)

fit <- survfit(Surv(time, event) ~expression, data = res.cat)

res.cox = coxph(Surv(time, event) ~expression, data = res.cat)
summary(res.cox)

#
p = ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata", 
           #linetype = "strata", 
           surv.median.line = "hv", 
           ggtheme = theme_bw(), 
           palette = c("#D62828", "#2E9FDF"))


filename = paste(gene_name, ".png", sep = "")

ggsave(filename,dpi = 600)