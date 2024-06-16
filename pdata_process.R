###########################################################
#AD part 
setwd("C:/Users/dell/Documents/R/Zhouzhike/")
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
set1 <- c('GSE132903','GSE118553','GSE5281','GSE37264','GSE36980')
set2 <- c('GSE20168','GSE68719')
set3 <- c('GSE33000')
# 1
studyID ='GSE132903'

eSet <- getGEO(studyID,destdir= '.',getGPL = F,AnnotGPL = F)
pdata = pData(eSet[[1]])
pdata = pdata[,c(2,34,35,36)]
colnames(pdata) = c('GSM','AD','Age','Gender')

AD_tmp = pdata$AD
AD_tmp[AD_tmp == 'AD']=1
AD_tmp[AD_tmp == 'ND']=0
pdata$AD = AD_tmp

gender_tmp = pdata$Gender
gender_tmp[gender_tmp == 'male']=1
gender_tmp[gender_tmp == 'female']=0
pdata$Gender = gender_tmp

pdata1 = pdata
save(pdata1,file = 'pdata1.RData')
#2 
studyID ='GSE118553'

eSet <- getGEO(studyID,destdir= '.',getGPL = F,AnnotGPL = F)
pdata = pData(eSet[[1]])
pdata = pdata[,c(2,36,35,37)]
colnames(pdata) = c('GSM','AD','Age','Gender')


type_index <- read.csv(paste(studyID,".csv",sep = ''))
pdata_control <- pdata[type_index$type== "Temporal_Cortex control ",]
pdata_case <- pdata[type_index$type== "Temporal_Cortex AD ",]
pdata <- rbind(pdata_control,pdata_case)

AD_tmp = pdata$AD
AD_tmp[AD_tmp == 'AD']=1
AD_tmp[AD_tmp == 'control']=0
pdata$AD = AD_tmp

gender_tmp = pdata$Gender
gender_tmp[gender_tmp == 'MALE']=1
gender_tmp[gender_tmp == 'FEMALE']=0
pdata$Gender = gender_tmp

pdata2 = pdata
save(pdata2,file = 'pdata2.RData')

#3
studyID ='GSE5281'

eSet <- getGEO(studyID,destdir= '.',getGPL = F,AnnotGPL = F)
pdata = pData(eSet[[1]])

type_index <- read.csv(paste(studyID,".csv",sep = ''))
pdata_control <- pdata[type_index$type== "MTGcontrol",]
pdata_case <- pdata[type_index$type== "MTGaffected",]
pdata <- rbind(pdata_control,pdata_case)

pdata = pdata[,c(2,30,21,19)]
colnames(pdata) = c('GSM','AD','Age','Gender')

AD_tmp = pdata$AD
AD_tmp[AD_tmp == 'cortical layer III neurons']=1
AD_tmp[AD_tmp == 'AD control']=0
pdata$AD = AD_tmp

gender_tmp = pdata$Gender
gender_tmp[gender_tmp == 'Sex: male聽'|gender_tmp =='sex: male聽']=1
gender_tmp[gender_tmp == 'sex: female聽'|gender_tmp =='Sex: female聽']=0
pdata$Gender = gender_tmp

AD_tmp = pdata$Age
pp = strsplit(AD_tmp,split=' ')
pdata$Age = unlist(lapply(pp,function(i)return(i[2])))
pdata3 = pdata
save(pdata3,file = 'pdata3.RData')

#4
studyID ='GSE37264'
eSet <- getGEO(studyID,destdir= '.',getGPL = F,AnnotGPL = F)
pdata = pData(eSet[[1]])

pdata = pdata[,c(2,40,39,42)]
colnames(pdata) = c('GSM','AD','Age','Gender')

AD_tmp = pdata$AD
AD_tmp[AD_tmp == "Alzheimer's disease (AD)"]=1
AD_tmp[AD_tmp == "control"]=0
pdata$AD = AD_tmp

gender_tmp = pdata$Gender
gender_tmp[gender_tmp == 'male']=1
gender_tmp[gender_tmp == 'female']=0
pdata$Gender = gender_tmp

pdata4 = pdata
save(pdata4,file = 'pdata4.RData')

5#
studyID ='GSE36980'

eSet <- getGEO(studyID,destdir= '.',getGPL = F,AnnotGPL = F)
pdata = pData(eSet[[1]])
pdata = pdata[,c(2,8,39,40)]
colnames(pdata) = c('GSM','AD','Age','Gender')

pdata_control <- pdata[pdata$AD== "Temporal cortex of non-AD brain",]
pdata_case <- pdata[pdata$AD== "Temporal cortex of AD brain",]
pdata <- rbind(pdata_control,pdata_case)

AD_tmp = pdata$AD
AD_tmp[AD_tmp == 'Temporal cortex of AD brain']=1
AD_tmp[AD_tmp == 'Temporal cortex of non-AD brain']=0
pdata$AD = AD_tmp

gender_tmp = pdata$Gender
gender_tmp[gender_tmp == 'M']=1
gender_tmp[gender_tmp == 'F']=0
pdata$Gender = gender_tmp

pdata5 = pdata
save(pdata5,file = 'pdata5.RData')
########################合并
pdata_AD = rbind(pdata1,pdata2,pdata3,pdata4,pdata5)
write.csv(pdata_AD,'pdata_AD.csv')
##############################################################
# PD part
6#
studyID ='GSE20168'

eSet <- getGEO(studyID,destdir= '.',getGPL = F,AnnotGPL = F)
pdata = pData(eSet[[1]])
pdata = pdata[,c(2,39,37,40)]
colnames(pdata) = c('GSM','AD','Age','Gender')

AD_tmp = pdata$AD
AD_tmp[AD_tmp == "Parkinson's disease"]=1
AD_tmp[AD_tmp == 'control']=0
pdata$AD = AD_tmp

gender_tmp = pdata$Gender
gender_tmp[gender_tmp == 'male']=1
gender_tmp[gender_tmp == 'female']=0
pdata$Gender = gender_tmp

pdata6 = pdata
save(pdata6,file = 'pdata6.RData')
7#
studyID ='GSE68719'

eSet <- getGEO(studyID,destdir= '.',getGPL = F,AnnotGPL = F)
pdata = pData(eSet[[1]])
pdata = pdata[,c(2,1,59,60)]
colnames(pdata) = c('GSM','AD','Age','Gender')

AD_tmp <- pdata %>% select(AD)%>%
  separate(AD,c("keep","drop"),sep = '_')%>% 
  select(-drop)
pdata$AD = AD_tmp$keep

AD_tmp = pdata$AD
AD_tmp[AD_tmp == "P"]=1
AD_tmp[AD_tmp == 'C']=0
pdata$AD = AD_tmp

gender_tmp = pdata$Gender
gender_tmp[gender_tmp == 'male']=1
gender_tmp[gender_tmp == 'female']=0
pdata$Gender = gender_tmp

pdata7 = pdata
save(pdata7,file = 'pdata7.RData')
########################合并
pdata_PD = rbind(pdata6,pdata7)
write.csv(pdata_PD,'pdata_PD.csv')

##############################################################
# HD part
8#
studyID ='GSE33000'

eSet <- getGEO(studyID,destdir= '.',getGPL = F,AnnotGPL = F)
pdata = pData(eSet[[1]])
pdata = pdata[,c(2,49,48,50)]
colnames(pdata) = c('GSM','AD','Age','Gender')

pdata_control <- pdata[pdata$AD== "non-demented",]
pdata_case <- pdata[pdata$AD== "Huntington's disease",]
pdata <- rbind(pdata_control,pdata_case)

Age_tmp <- pdata %>% select(Age)%>%
  separate(Age,c("keep","drop"),sep = ' ')%>% 
  select(-drop)
pdata$Age = Age_tmp$keep

AD_tmp = pdata$AD
AD_tmp[AD_tmp == "Huntington's disease"]=1
AD_tmp[AD_tmp == "non-demented"]=0
pdata$AD = AD_tmp

gender_tmp = pdata$Gender
gender_tmp[gender_tmp == 'male']=1
gender_tmp[gender_tmp == 'female']=0
pdata$Gender = gender_tmp

pdata8 = pdata
save(pdata8,file = 'pdata8.RData')
pdata_HD = pdata8
write.csv(pdata_HD,'pdata_HD.csv')