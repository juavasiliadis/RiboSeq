#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("edgeR")
library(magrittr)
library(dplyr)
library(edgeR)
install.packages("ggplot2")
library(ggplot2)

df = read.csv('/home/ulya/RiboSeq/HSE_RiboSeq_HT/01. RiboSeq_RNASeq_HCC_counts.tsv', sep = '\t')
df = replace(df, is.na(df), 0)
RNA = df %>% select(ends_with("RNA"))
RNA = cbind(df$geneSymbol, RNA)
RPF = df %>% select(ends_with("RPF"))
RPF = cbind(df$geneSymbol, RPF)

# Диф. экспрессия для данных RNA
y = DGEList(counts = RNA[,2:21], genes = RNA[,1])
patient = factor(c('LC001', 'LC001', 'LC033', 'LC033', 'LC034', 'LC034', 'LC501', 'LC501', 'LC502', 'LC502', 'LC505', 'LC505', 'LC506', 'LC506', 'LC507', 'LC507', 'LC508','LC508', 'LC509', 'LC509'))
tissue = factor(c('N', 'T', 'N', 'T', 'N', 'T', 'N', 'T', 'N', 'T','N', 'T', 'N', 'T', 'N', 'T', 'N', 'T', 'N', 'T'))
data.frame(Sample=colnames(y), patient, tissue)
design = model.matrix(~ patient + tissue)
rownames(design) = colnames(y)
design
y = estimateDisp(y, design, robust=TRUE)
y$common.dispersion
fit = glmFit(y, design)
lrt = glmLRT(fit)
topTags(lrt)
colnames(design)
o = order(lrt$table$PValue)
cpm(y)[o[1:15],]
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")

# Диф. экспрессия для данных RPF 
y1 = DGEList(counts = RPF[,2:21], genes = RPF[,1])
patient = factor(c('LC001', 'LC001', 'LC033', 'LC033', 'LC034', 'LC034', 'LC501', 'LC501', 'LC502', 'LC502', 'LC505', 'LC505', 'LC506', 'LC506', 'LC507', 'LC507', 'LC508','LC508', 'LC509', 'LC509'))
tissue = factor(c('N', 'T', 'N', 'T', 'N', 'T', 'N', 'T', 'N', 'T','N', 'T', 'N', 'T', 'N', 'T', 'N', 'T', 'N', 'T'))
data.frame(Sample=colnames(y1), patient, tissue)
design1 = model.matrix(~ patient + tissue)
rownames(design) = colnames(y1)
design1
y1 = estimateDisp(y1, design, robust=TRUE)
y1$common.dispersion
fit1 = glmFit(y1, design1)
lrt1 = glmLRT(fit1)
topTags(lrt1)
colnames(design1)
o = order(lrt1$table$PValue)
cpm(y1)[o[1:15],]
summary(decideTests(lrt1))
plotMD(lrt1)
abline(h=c(-1, 1), col="blue")

# Volcano plots
p = ggplot(data = df, aes(x = lrt$table$logFC, y = -log10(lrt$table$PValue))) + geom_point()
p1 = ggplot(data = df, aes(x = lrt1$table$logFC, y = -log10(lrt1$table$PValue))) + geom_point()
p + labs(x = 'Fold-change', y = '-log10(p-value)')
p1 + labs(x = 'Fold-change', y = '-log10(p-value)')
# DE results do not match

# Выберем диф. экспрессируемые гены со значимым значением p-value < 0.05`
sign_genes_RNA = list()
sign_genes_RNA = lrt$genes$genes[lrt$table$PValue < 0.05]
sign_genes_RPF = list()
sign_genes_RPF = lrt1$genes$genes[lrt1$table$PValue < 0.05]

library(VennDiagram)
# Построим диаграмму Венна
venn.diagram(
  x = list(sign_genes_RNA, sign_genes_RPF),
  category.names = c("RNA" , "RPF"),
  filename = 'venn_diagramm.png',
  output=TRUE
)

# В результате у RNA и RPF идентифицировано пересечение по 4044 диф. экспрессируемым генам. 
