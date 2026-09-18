#PRADA Genomics Quality Evaluation, used for WCPG 2026 poster


library(data.table)
library(ggplot2)
library(ggrepel)


#projectFolderPath<-"/scratch/prj/sgdp_nanopore/Projects/prada_jz" #remote
projectFolderPath<-"/Users/jakz/Documents/work_rstudio/pgxrex" #local

##read sample results
dSample<-fread(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis","samples.tsv"))

dSample[is.na(evaluationRatio),evaluationRatio:=1] #fallback

dSample.filtered<-dSample[evaluationRatio>=1,] #filter


#remove duplicate samples
dSample.n<-nrow(dSample.filtered)
dSample.filtered$MERGEID<-1:dSample.n

dSample.filtered<-dSample.filtered[order(-evaluationRatio,MERGEID),]
dSample.unique<-dSample.filtered[, .(MERGEID = head(MERGEID,1)), by = c("analysis","mostCredibleReference")]
dSample.unique<-dSample.filtered[dSample.unique, on=c(MERGEID=c("MERGEID"))]

#dSequencing <- fread(file.path(projectFolderPath,"work","pradaApp","unified-plots","sampleMetaTot.tsv")) #we don't have to read in this separately now

#dSample.unique[dSequencing, on=c("analysis","barcode"), c("depth_on","depth_off"):=list(i.sdepth_q050_bed, i.sdepth_q050_nobed)]
dSample.unique[,label:=paste0(analysis,".",mostCredibleReference)]

colnamesAllGeneCalls <- c("GC_ABCG2","GC_CACNA1S","GC_CFTR","GC_CYP2B6","GC_CYP2C19","GC_CYP2C9","GC_CYP2D6","GC_CYP3A4","GC_CYP3A5","GC_CYP4F2","GC_DPYD","GC_G6PD","GC_IFNL3","GC_NUDT15","GC_RYR1","GC_SLCO1B1","GC_TPMT","GC_UGT1A1","GC_VKORC1")
colnamesAllGeneCalls.antidepressant <- c("GC_CYP2B6","GC_CYP2C19","GC_CYP2D6")

dSample.unique[,c("nCalls")] <-rowSums(dSample.unique[,..colnamesAllGeneCalls],na.rm = T)
dSample.unique[,c("nCalls.dep")] <-rowSums(dSample.unique[,..colnamesAllGeneCalls.antidepressant],na.rm = T)
dSample.unique[,sample:=mostCredibleReference]


#Samtools stats

#genomeLengthSettingBp <- 3088269832

#pilot10
folderpathProject<-file.path(projectFolderPath,"work","pradaApp","pilot10")
cProjectSamtoolsStatsMeta<-c()
cProjectSamtoolsStatsReadLength<-c()
con = file(file.path(folderpathProject,paste0("PradaA.haplotagged.bam.ontarget.stats.sn.txt")),open="r")
cProjectSamtoolsStatsMeta["PradaA"] <- list(readLines(con = con,encoding = "UTF-8", warn = F))
close(con)
iYield <- grep("bases mapped \\(cigar\\):",cProjectSamtoolsStatsMeta["PradaA"][[1]])[[1]]
matches<-gregexpr("(\\d+)",cProjectSamtoolsStatsMeta["PradaA"][[1]][iYield])
cYield<-regmatches(cProjectSamtoolsStatsMeta["PradaA"][[1]][iYield], matches)[[1]]
dSample.unique[analysis=="p10" & barcode=="PradaA",yield.ontarget:=eval(cYield)]

cProjectSamtoolsStatsReadLength["PradaA"] <- list(fread(file.path(folderpathProject,paste0("PradaA.haplotagged.bam.ontarget.stats.rl.txt")),header = F, data.table = T))
cLengthData<-(cProjectSamtoolsStatsReadLength["PradaA"][[1]])[,totLength:=V1*V2]
cLengthData<-cLengthData[order(-V1),][,cumLength:=cumsum(totLength)]
totalAssemblyLength<-max(cLengthData$cumLength)
cN50<-cLengthData[cumLength>eval(totalAssemblyLength)/2,][1,V1]
dSample.unique[analysis=="p10" & barcode=="PradaA",n50.ontarget:=eval(cN50)]


