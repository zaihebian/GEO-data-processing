a = read.table('GSE5281_series_matrix.txt.gz',
               sep = '\t',quote = '',fill = T,
               comment.char = '!',header = T)
#################################################################
#_______________________________________________________________#
#                          下载数据                             #
#_______________________________________________________________#
#################################################################
downGSE <- function(studyID = 'GSE1009',destdir = '.')
{
  library(GEOquery)
  eSet <- getGEO(studyID,destdir= destdir,getGPL = F,AnnotGPL = F)
  exprSet = exprs(eSet[[1]])
  pdata = pData(eSet[[1]])
  write.csv(exprSet,paste0(studyID,"_exprSet.csv"))#或者用save
# write.csv(pdata,paste0(studyID,"_metadata.csv"))
  return(eSet)
}
set1 <- c('GSE132903','GSE118553','GSE5281','GSE37264','GSE36980')
set2 <- c('GSE20168','GSE68719')
set3 <- c('GSE33000')
allset <- c(set1,set2,set3)
for(i in allset)
{
  downGSE(i)
}
#################################################################
#_______________________________________________________________#
#                          处理数据                             #
#_______________________________________________________________#
#################################################################
control_names <- c("diagnosis: ND","Temporal_Cortex control ","MTGcontrol","CNT (exon level)","non-AD_TC",
                   "disease state: control" ,"\tC","disease status: non-demented")
case_names <- c("diagnosis: AD","Temporal_Cortex AD ","MTGaffected","AD (exon level)","AD_TC",
                "disease state: Parkinson's disease","\tP","disease status: Huntington's disease")
# 变量设置区
k <- 4
studyID <- allset[k]
i <- as.character(k)
control_name <- control_names[k]
case_name <- case_names[k]
setwd("C:/Users/dell/Documents/R/Zhouzhike/rawdata/AD")
# 1)将数据读进来
datatmp <- read.csv(paste(studyID,"_exprSet.csv",sep = ''))
datatmp[1:4,1:4]
dim(datatmp)
rownames(datatmp) <- datatmp[,1]
datatmp <- datatmp[-1]
# 2）选取需要的样本（即选取列）,并按先控制，后病例的顺序来存放
type_index <- read.csv(paste(studyID,".csv",sep = ''))
table(type_index$type) #查看分类如何区分
unique(type_index$type)
datatmp_control <- datatmp[,type_index$type== control_name]
datatmp_case <- datatmp[,type_index$type== case_name]
datatmp <- cbind(datatmp_control,datatmp_case)
#### 此处需要手动修改 ############
# 3)探针号码转ID,illuminaHumanv4可以替换为各种其他报名，查看类使用ls(包名)
# 1&2==== illuminaHumanv4
# 3 ==== hgu133plus2
# 4 ===== 
# 5 ===hugene10sttranscriptcluster
# 6 ===hgu133a
library(huex10sttranscriptcluster.db)    #### 此处需要手动修改 ############
#ids = toTable(huex10sttranscriptclusterSYMBOL)#20000多个
library(illuminaHumanv4.db)    
#ids = toTable(illuminaHumanv4SYMBOL)#20000多个
library(hgu133plus2.db)
ids = toTable(hgu133plus2SYMBOL)
library(hugene10sttranscriptcluster.db)
ids = toTable(hugene10sttranscriptclusterSYMBOL)
library(hgu133a.db)
ids = toTable(hgu133aSYMBOL)
#length(unique(ids$symbol))
#tail(sort(table(ids$symbol)))
#table(sort(table(ids$symbol)))
#plot(table(sort(table(ids$symbol))))

# 4)第一次筛选,选出id能匹配到的探针行
datatmp <- datatmp[rownames(datatmp)%in% ids$probe_id,]
ids <- ids[match(rownames(datatmp),ids$probe_id),]
# 5)第二次筛选,选出重复基因最大表达量的探针
tmp <- by(datatmp,ids$symbol,function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
datatmp = datatmp[rownames(datatmp)%in%probes,]
# 6）找出symbol_name
symbol_name = as.data.frame(ids[match(rownames(datatmp),ids$probe_id),2])
names(symbol_name) = 'symbol_name'
rownames(datatmp) <- symbol_name[,1]
datatmp[1:4,1:4]
# 7) 保存数据
var1_name <- paste('symbol_name_',i,sep = '')
var2_name <- paste('data_',i,sep='')
assign(var1_name,symbol_name)
assign(var2_name,datatmp)
#######手动执行这一步
save(data_4,symbol_name_4,file = paste('data_',i,'.RData',sep=''))
write.csv(data_k,file = paste('processed_',allset[k],'.csv',sep = ''))
#load(file = paste('data_',i,'.RData',sep=''))
#到这里得到去除批次效应前并排序的datatmp全数据
#数据融合
for(j in c(1,2,3))
{
  load(file = paste('data_',as.character(j),'.RData',sep=''))
}
data_3 <- log2(data_3+1)
ave <- colMeans(data_5)
plot(ave)
data_1t = cbind(symbol_name_1,data_1)
data_2t = cbind(symbol_name_2,data_2)
data_3t = cbind(symbol_name_3,data_3)
data_4t = cbind(symbol_name_4,data_4)
data_5t = cbind(symbol_name_5,data_5)
# 融合两组数据
mergeRawData <- merge(data_1t,data_2t,by = 'symbol_name')
mergeRawData <- merge(mergeRawData,data_3t,by = 'symbol_name')
mergeRawData <- merge(mergeRawData,data_4t,by = 'symbol_name')
mergeRawData <- merge(mergeRawData,data_5t,by = 'symbol_name')
rownames(mergeRawData) <- mergeRawData[,1]
mergeRawData <- mergeRawData[-1]
dim(mergeRawData)
allmeanv <- colMeans(mergeRawData)
plot(allmeanv)
write.csv(mergeRawData,file = 'all_AD_rawData.csv')
#################################################################
#_______________________________________________________________#
#                       去除批次效应                            #
#_______________________________________________________________#
#################################################################
library(limma)
batch <- c(rep('PD1',29),rep('PD2',73))
batch <- c(rep('AD1',195),rep('AD2',69),rep('AD3',28),rep('AD4',16),rep('AD5',29))
#去掉离群值
#mergeRawCounts <- mergeRawCounts[,-which(colnames(mergeRawCounts)=='GSM119649')]
mergeData <- removeBatchEffect(mergeRawData,batch)
allmeanv <- colMeans(mergeData)
plot(allmeanv)
write.csv(mergeData,file = 'all_PD_Data.csv')
#################################################################
#_______________________________________________________________#
#                       准备差异分析                            #
#_______________________________________________________________#
#################################################################
group1 <- c(rep('control',98),rep('case',97))
group2 <- c(rep('control',24),rep('case',45))
group3 <- c(rep('control',12),rep('case',16))
group4 <- c(rep('control',8),rep('case',8))
group5 <- c(rep('control',19),rep('case',10))
group6 <- c(rep('control',15),rep('case',14))
group7 <- c(rep('control',44),rep('case',29))
group8 <- c(rep('control',157),rep('case',157))
group_list <- c(group1,group2,group3,group4,group5)
group_list <- c(group6,group7)
group_list <- group8
mergeData <- data_8
### 查看数据（慎用，运行耗时）
par(cex = 0.7)
n.sample=ncol(mergeData)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(mergeData, col = cols,main="expression value",las=2)
#做design矩阵
design <- model.matrix(~0+relevel(factor(group_list),'control'))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(mergeData)
#做contrast矩阵
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix 
#################################################################
#_______________________________________________________________#
#                       进行差异分析                            #
#_______________________________________________________________#
#################################################################
##step1
fit <- lmFit(mergeData,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)
write.csv(nrDEG,file = paste('result_PD.csv',sep = ''))
