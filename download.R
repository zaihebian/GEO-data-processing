8月30日,本文件下载的是用getGEO下载实验数据，而直接下载实验数据的压缩包进行读取的数据,均为RNAseq
library("FactoMineR")
library("factoextra") 
library(GEOquery)
library(limma)
library(sva)
options(stringsAsFactors = F)
workingdir = "C:/Users/Dell/Documents/R/GEO"
setwd(workingdir)
memory.limit(400000)


# ----------------------------------------------------------------------------------
# GSEseries: GSE94660   个样本,  个受试者
# ----------------------------------------------------------------------------------
# 只处理expr
expr <- read.table(file = 'GSE94660_RPKM_normalized.txt', sep = "\t", header = T, stringsAsFactors = F,fill = T)


# 只处理pdata
alldata = getGEO('GSE94660')
pdata0 = pData(alldata[[1]])
pdata = pdata0[,c("title","source_name_ch1")]
colnames(pdata) = c("patient","group")
group_list = pdata[,2]

# 将两者对应起来，顺序以pdata为准
group_list[group_list == "Non-neoplastic liver tissue"] = "HBV"
group_list[group_list == "Tumor tissue"] ="HCC"
save(expr,pdata,group_list,file = "GSE94660_expr_pdata_group.RData")
# ----------------------------------------------------------------------------------
# END                 END                 END                   END
# ----------------------------------------------------------------------------------




# ----------------------------------------------------------------------------------
# GSEseries: GSE130823   个样本,  个受试者
# ----------------------------------------------------------------------------------
# 只处理expr
alldata = getGEO('GSE130823')

expr <- exprs(alldata$GSE130823_series_matrix.txt.gz)
probes = rownames(expr)
library(biomaRt)
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
annot <- getBM(
  attributes = c('agilent_wholegenome',
                 "hgnc_symbol",
                 'ensembl_gene_id',
                 'entrezgene_id',
                 'external_gene_name'),
  filters = 'agilent_wholegenome',
  values = probes,
  mart = ensembl)

annot = annot[annot$hgnc_symbol != '',]
annot = annot[!duplicated(annot$agilent_wholegenome),]
annot = annot[!duplicated(annot$hgnc_symbol),]


expr = expr[annot$agilent_wholegenome,]
rownames(expr) = annot$hgnc_symbol


# 只处理pdata

pdata0 = pData(alldata[[1]])
pdata = pdata0[,c("title","source_name_ch1")]
colnames(pdata) = c("patient","group")
group_list = pdata[,2]

# 将两者对应起来，顺序以pdata为准
group_list[group_list %in% c("GastricHighGradeIntraepithelialNeoplasia","GastricLowGradeIntraepithelialNeoplasia")] = "Neoplasia"
group_list[group_list == "IntestinalGastricCancer"] ="GastricCancer"
save(expr,pdata,group_list,file = "GSE130823_expr_pdata_group.RData")
load(file = "GSE130823_expr_pdata_group.RData")
# ----------------------------------------------------------------------------------
# END                 END                 END                   END
# ----------------------------------------------------------------------------------




# ----------------------------------------------------------------------------------
# GSEseries: GSE55696   个样本,  个受试者
# ----------------------------------------------------------------------------------
# 只处理expr
alldata = getGEO('GSE55696')

expr <- exprs(alldata$GSE55696_series_matrix.txt.gz)
probes = rownames(expr)
library(biomaRt)
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
annot <- getBM(
  attributes = c('agilent_wholegenome',
                 "hgnc_symbol",
                 'ensembl_gene_id',
                 'entrezgene_id',
                 'external_gene_name'),
  filters = 'agilent_wholegenome',
  values = probes,
  mart = ensembl)

annot = annot[annot$hgnc_symbol != '',]
annot = annot[!duplicated(annot$agilent_wholegenome),]
annot = annot[!duplicated(annot$hgnc_symbol),]


expr = expr[annot$agilent_wholegenome,]
rownames(expr) = annot$hgnc_symbol


# 只处理pdata

pdata0 = pData(alldata[[1]])
pdata = pdata0[,c("title","source_name_ch1")]
colnames(pdata) = c("patient","group")
group_list = pdata[,2]

# 将两者对应起来，顺序以pdata为准
group_list[group_list == "inflammation"] ="CG"
save(expr,pdata,group_list,file = "GSE55696_expr_pdata_group.RData")

# ----------------------------------------------------------------------------------
# GSEseries: GSE144424   个样本,  个受试者
# ----------------------------------------------------------------------------------
# 只处理expr
expr <- read.table(file = 'GSE144424_Counts_RNA_MCW_NEB.txt', sep = "\t", header = T, stringsAsFactors = F,fill = T)


# 只处理pdata
alldata = getGEO('GSE144424')
pdata0 = pData(alldata[[1]])

exprsamples = colnames(expr)[7:90]
exprsamples = gsub("H","",exprsamples)
pdataSamples = pdata0$title
pdataSamples = gsub("_RNA-seq","",pdataSamples)

pdata0 = pdata0[match(exprsamples,pdataSamples),]

pdata = pdata0[,c("title","condition:ch1")]
colnames(pdata) = c("patient","group")
group_list = pdata[,2]


save(expr,pdata,group_list,file = "GSE144424_expr_pdata_group.RData")
# ----------------------------------------------------------------------------------
# END                 END                 END                   END
# ----------------------------------------------------------------------------------

all_TCGA_genes = rownames(all.tpm.final)
save(all_TCGA_genes,file = "all_TCGA_genes.RData")

aaa =  expr[match(all_TCGA_genes,rownames(expr)),]
rownames(aaa) =all_TCGA_genes
aaa[is.na(aaa)] = 0


