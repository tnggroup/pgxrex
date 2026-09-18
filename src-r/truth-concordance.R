#Truth concordance analysis with unified plotting, 13/07/2026

#remotes::install_github("tnggroup/pgxrex")
#remotes::install_github("tnggroup/pgxrex",ref = 'jz_dev')
library(pgxrex)
library(data.table)
library(vcfR)

projectFolderPath<-"/scratch/prj/sgdp_nanopore/Projects/prada_jz"
#projectFolderPath<-"/Users/jakz/Documents/work_rstudio/pgxrex" #local

dAnalysis <- fread(file.path(projectFolderPath,"data","pradaApp","prada_sequencing_run_database.analysis.tsv.txt")) #these are database files containing the collected folder configurations for all runs/analyses.
dAnalysis[,code:=analysis_id] #harmonise with pgxrex naming
filenameSampleTemplate<-file.path(projectFolderPath,"data","pradaApp","prada_sequencing_run_database.sample.tsv.txt")
filenameSamplePreviousResults<-file.path(projectFolderPath,"work","pradaApp","per-sample-analysis","samples.tsv")

#use existing result if possible
usingPreviousResults<-F
if(file.exists(filenameSamplePreviousResults)) {
    catl("Using existing sample results")
    dSample<-fread(filenameSamplePreviousResults)
    usingPreviousResults<-T
  } else {
    catl("Initiating sample data from template")
    dSample <- fread(file.path(filenameSampleTemplate))
  }

#reference H (HG002)
dSampleVCF<-read.vcfR(file.path(projectFolderPath,"work","pgx","downsampled-bam-runs","truth-set","HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"), verbose = F)
dSampleVCF.H<-as.data.frame(getFIX(dSampleVCF))
dSampleVCF.H<-dSampleVCF.H[dSampleVCF.H$CHROM=="chr1",] #subset to chr1
dSamplePGX.H<-fread(file.path(projectFolderPath,"work","pradaApp","downsampled-bam-runs","pgxCallsAggCustom_multi-hg002-100_HG002_new_hg19_merged_hg38.tsv"))
setDT(dSampleVCF.H)
setkeyv(dSampleVCF.H,c("CHROM","POS","REF","ALT"))

#reference A
dSampleVCF<-read.vcfR(file.path(projectFolderPath,dAnalysis[analysis_id=="p10-wgs-A",c("pathAnalysisOutput")],"output","Prada_WGS_Sample_A",paste0("Prada_WGS_Sample_A",".filtered.vcf.gz")), verbose = F)
dSampleVCF.A<-as.data.frame(getFIX(dSampleVCF))
dSampleVCF.A<-dSampleVCF.A[dSampleVCF.A$CHROM=="chr1",] #subset to chr1
dSamplePGX.A<-fread(file.path(projectFolderPath,"work","pradaApp","pilot10","pgxCallsAggCustom_p10-wgs-A_Prada_WGS_Sample_A.tsv"))
setDT(dSampleVCF.A)
setkeyv(dSampleVCF.A,c("CHROM","POS","REF","ALT"))

#reference B
dSampleVCF<-read.vcfR(file.path(projectFolderPath,dAnalysis[analysis_id=="p10-wgs-B",c("pathAnalysisOutput")],"output","Prada_WGS_Sample_B",paste0("Prada_WGS_Sample_B",".filtered.vcf.gz")), verbose = F)
dSampleVCF.B<-as.data.frame(getFIX(dSampleVCF))
dSampleVCF.B<-dSampleVCF.B[dSampleVCF.B$CHROM=="chr1",] #subset to chr1
dSamplePGX.B<-fread(file.path(projectFolderPath,"work","pradaApp","pilot10","pgxCallsAggCustom_p10-wgs-B_Prada_WGS_Sample_B.tsv"))
setDT(dSampleVCF.B)
setkeyv(dSampleVCF.B,c("CHROM","POS","REF","ALT"))

#reference C
dSampleVCF<-read.vcfR(file.path(projectFolderPath,dAnalysis[analysis_id=="p10-wgs-C",c("pathAnalysisOutput")],"output","Prada_WGS_Sample_C",paste0("Prada_WGS_Sample_C",".filtered.vcf.gz")), verbose = F)
dSampleVCF.C<-as.data.frame(getFIX(dSampleVCF))
dSampleVCF.C<-dSampleVCF.C[dSampleVCF.C$CHROM=="chr1",] #subset to chr1
dSamplePGX.C<-fread(file.path(projectFolderPath,"work","pradaApp","pilot10","pgxCallsAggCustom_p10-wgs-C_Prada_WGS_Sample_C.tsv"))
setDT(dSampleVCF.C)
setkeyv(dSampleVCF.C,c("CHROM","POS","REF","ALT"))

print("References read")

#we could re-run the analysis and data collection step here if needed
pgxrexObj<-PgxrexClass()
pgxrexObj$folderpathWork<-file.path(projectFolderPath,"work","pradaApp","pilot13") #the folder where this is executed
pgxrexObj$applicationCoverageRegions <- fread(file = file.path(projectFolderPath,"data","roughApplicationCoverageRegionsAsOfPilot3.tsv"), na.strings = c(".",
                                                                               NA, "NA", ""), encoding = "UTF-8", check.names = T,
                                              fill = T, blank.lines.skip = T, data.table = F, nThread = 6,
                                              showProgress = F)


origWD<-getwd()

#sync database and object representations
pgxrexObj$analysisMeta<-as.data.frame(dAnalysis)
rownames(pgxrexObj$analysisMeta)<-pgxrexObj$analysisMeta$code
pgxrexObj$sampleMeta<-as.data.frame(dSample)
rownames(pgxrexObj$sampleMeta)<-paste0(pgxrexObj$sampleMeta$analysis,"_",pgxrexObj$sampleMeta$barcode)

if(is.null(pgxrexObj$analysisSettingsList)) pgxrexObj$analysisSettingsList<-list() #in case these are null due to not reading in the raw data
if(is.null(pgxrexObj$sampleSettingsList)) pgxrexObj$sampleSettingsList<-list()

#read in previous results
for(iAnalysis in 1:nrow(pgxrexObj$analysisMeta)){
  #iAnalysis<-2
  #cAnalysisID<-dAnalysis[iAnalysis,c("analysis_id")] #not used
  cPilotID<-pgxrexObj$analysisMeta[iAnalysis,c("pilot_id")]
  if(is.null(cPilotID)) next
  if(nchar(cPilotID)<1) next
  #file.exists(file.path(projectFolderPath,"work","pradaApp",cPilotID,"analysisMeta.tsv"))

  #test
  # targetFolderpath <- file.path(projectFolderPath,"work","pradaApp",cPilotID)
  # cFilepath<-file.path.coalesce(targetFolderpath,"analysisMeta.tsv")
  # analysisMeta.toAdd <- fread(file = cFilepath, na.strings = c(".",
  #                                                              NA, "NA", ""), encoding = "UTF-8", check.names = T,
  #                             fill = T, blank.lines.skip = T, data.table = F, nThread = 6,
  #                             showProgress = F)
  #
  # analysisMeta<-pgxrexObj$analysisMeta
  # analysisMeta[analysisMeta.toAdd$code,colnames(analysisMeta.toAdd)] <- analysisMeta.toAdd


  pgxrexObj$readPrintData(file.path(projectFolderPath,"work","pradaApp",cPilotID))
}

#complementary analyses - per analysis - This is time-consuming - skip if existing results
#produces a temporary version of sample results
for(iAnalysis in 1:nrow(pgxrexObj$analysisMeta)){

  #iAnalysis<-2
  cAnalysisID<-pgxrexObj$analysisMeta[iAnalysis,c("code")]
  catl(paste0("Complementary analysis for ",cAnalysisID))
  if(is.null(cAnalysisID)) next
  if(nchar(cAnalysisID)<1) next

  cSamples<-pgxrexObj$sampleMeta[pgxrexObj$sampleMeta$analysis==cAnalysisID,]

  if(any(names(cSamples)=="sdepth_q050.CYP2B6")){
    if(all(is.finite(
      unlist(cSamples[,c(
      "sdepth_q050.CYP2B6",
      "sdepth_q050.CYP2C19",
      "sdepth_q050.CYP2D6"
    )])
    ))) next
  }

  pgxrexObj$collectAnalysisDepthData(cAnalysisID)
  pgxrexObj$computeDepthDataStatistics(filePathBed = file.path(projectFolderPath,"data/bed/pgx.grch38.5k.0p7percent.bed"))


  for(iSample in 1:nrow(cSamples)){
    #iSample<-1

    cUniqueSampleId<-rownames(cSamples)[iSample]
    #depth statistics if available (computed again above)
    if(!is.null(pgxrexObj$sampleSettingsList[[cUniqueSampleId]]$sequencingStatsOriginalRegionsTable)){
      cStats<-pgxrexObj$sampleSettingsList[[cUniqueSampleId]]$sequencingStatsOriginalRegionsTable
      pgxrexObj$sampleMeta[cUniqueSampleId,c(
        paste0("sdepth_q050.CYP2B6"),
        paste0("sdepth_q050.CYP2C19"),
        paste0("sdepth_q050.CYP2D6")
      )]<-list(
        cStats[cStats$label_region=="CYP2B6",c("sdepth_q050")],
        cStats[cStats$label_region=="CYP2C19",c("sdepth_q050")],
        cStats[cStats$label_region=="CYP2D6",c("sdepth_q050")]
      )

      catl("Extra depth information calculated for ",cUniqueSampleId)

    }

    #cleanup
    pgxrexObj$sampleSettingsList[[cUniqueSampleId]]$sequencingDepthTable<-NULL
    pgxrexObj$sampleSettingsList[[cUniqueSampleId]]$sequencingDepthRegionsTable<-NULL
    pgxrexObj$sampleSettingsList[[cUniqueSampleId]]$sequencingDepthRegionsTableCustom<-NULL
    pgxrexObj$sampleSettingsList[[cUniqueSampleId]]$sequencingStatsOriginalRegionsTable<-NULL

  }


  #pgxrexObj$sampleSettingsList[["p2-gtube_barcode01"]]$sequencingDepthRegionsTable

  shru::writeFile(pgxrexObj$sampleMeta,file=file.path(projectFolderPath,"work","pradaApp","per-sample-analysis","samples.tsv"),nThreads = 5)

}



#this only depends on the file/database metadata and does not index existing barcodes again
for(iSample in 1:nrow(dSample)){
  #iSample<-69
  cSampleID<-pgxrexObj$sampleMeta[iSample,c("barcode")]
  cUniqueSampleId<-rownames(pgxrexObj$sampleMeta)[iSample]
  cAnalysisID<-pgxrexObj$sampleMeta[iSample,c("analysis")]
  cPilotID<-pgxrexObj$analysisMeta[cAnalysisID,c("pilot_id")]
  cPathAnalysisOutput<-pgxrexObj$analysisMeta[cAnalysisID,c("pathAnalysisOutput")]

  catl(paste0("Sample ",cSampleID))

  # dSampleVCF.result<-shru::supermunge(filePaths = file.path(projectFolderPath,cPathAnalysisOutput,"output",cSampleID,paste0(cSampleID,".filtered.vcf.gz")),traitNames = cSampleID,writeOutput = F)
  # dSampleVCF<-dSampleVCF.result$last

  knownReference<-NA

  folderPathSample <- file.path(projectFolderPath,cPathAnalysisOutput,"output",cSampleID)
  if(!file.exists(folderPathSample)) warning("Sample folder does not exist. May be due to a folder naming issue.")
  filePathVCF <- file.path(projectFolderPath,cPathAnalysisOutput,"output",cSampleID,paste0(cSampleID,".filtered.vcf.gz"))
  if(file.exists(filePathVCF)){
    dSampleVCF<-read.vcfR(filePathVCF, verbose = F)
    dSampleVCF.fix<-as.data.frame(getFIX(dSampleVCF))
    dSampleVCF.fix<-dSampleVCF.fix[dSampleVCF.fix$CHROM=="chr1",] #subset to chr1
    setDT(dSampleVCF.fix)
    setkeyv(dSampleVCF.fix,c("CHROM","POS","REF","ALT"))

  } else {
    #assume hg002
    dSampleVCF.fix <- dSampleVCF.H
    knownReference<-"H"
  }


  if(any(colnames(pgxrexObj$sampleMeta)=="mostCredibleReference")){
    if(!is.na(pgxrexObj$sampleMeta[iSample,c("mostCredibleReference")])) next #only process empty rows
  }


  mostCredibleReference<-"NA"
  ratioOfMostCredibleReference<-0
  mConcordantREF.credible<-NA
  mConcordantALT.credible<-NA
  mDiscordantREF.credible<-NA
  mDiscordantALT.credible<-NA

  if(is.na(knownReference)){
    cat(paste0(", establishing credible reference"))
    for(cComparison in list("H","A","B","C")){
      #cComparison<-"H"
      dSampleVCF.REF<-dSampleVCF.H
      if(cComparison=="H") dSampleVCF.REF<-dSampleVCF.H
      if(cComparison=="A") dSampleVCF.REF<-dSampleVCF.A
      if(cComparison=="B") dSampleVCF.REF<-dSampleVCF.B
      if(cComparison=="C") dSampleVCF.REF<-dSampleVCF.C

      colnames(dSampleVCF.REF)<-paste0(colnames(dSampleVCF.REF),".REF")

      dSampleVCF.fix.merged <- dSampleVCF.fix[dSampleVCF.REF, on=.(CHROM=CHROM.REF,POS=POS.REF)]

      mREF<-nrow(dSampleVCF.REF)
      mVCF<-nrow(dSampleVCF.fix)
      mConcordantREF<-nrow(dSampleVCF.fix.merged[REF==REF.REF,])
      mConcordantALT<-nrow(dSampleVCF.fix.merged[ALT==ALT.REF,])
      mDiscordantREF<-nrow(dSampleVCF.fix.merged[REF!=REF.REF,])
      mDiscordantALT<-nrow(dSampleVCF.fix.merged[ALT!=ALT.REF,])

      evaluationRatio<-(mConcordantREF+1000*mConcordantALT-mDiscordantREF-1000*mDiscordantALT)/mVCF
      if(is.finite(evaluationRatio) && evaluationRatio>ratioOfMostCredibleReference){
        mostCredibleReference<-cComparison
        ratioOfMostCredibleReference<-evaluationRatio
        mConcordantREF.credible<-mConcordantREF
        mConcordantALT.credible<-mConcordantALT
        mDiscordantREF.credible<-mDiscordantREF
        mDiscordantALT.credible<-mDiscordantALT
      }

      pgxrexObj$sampleMeta[iSample,c(
        #paste0("mREF.",cComparison),
        paste0("mVCF.",cComparison),
        paste0("mConcordantREF.",cComparison),
        paste0("mConcordantALT.",cComparison),
        paste0("mDiscordantREF.",cComparison),
        paste0("mDiscordantALT.",cComparison)
      )]<-list(mVCF,mConcordantREF,mConcordantALT,mDiscordantREF,mDiscordantALT)

    }
  } else {
    mostCredibleReference<-knownReference
  }

  cat(paste0(": ",mostCredibleReference))

  pgxrexObj$sampleMeta[iSample,c("mostCredibleReference")]<-list(
                                   mostCredibleReference)

  if(is.na(knownReference)){
    pgxrexObj$sampleMeta[iSample,c(
                    "evaluationRatio",
                    "mVCF",
                    "mConcordantREF.credible",
                    "mConcordantALT.credible",
                    "mDiscordantREF.credible",
                    "mDiscordantALT.credible")]<-list(
                      ratioOfMostCredibleReference,
                      mVCF,
                      mConcordantREF.credible,
                      mConcordantALT.credible,
                      mDiscordantREF.credible,
                      mDiscordantALT.credible)
  }

  filePathPGX<-file.path(projectFolderPath,"work","pradaApp", cPilotID ,paste0("pgxCallsAggCustom_",cAnalysisID,"_",cSampleID,".tsv")) #the pilot folder has to have the same name as the pilot ID.
  if(file.exists(filePathPGX)){
    dSamplePGX<-fread(filePathPGX)
    dSamplePGX.REF<-dSamplePGX.H
    if(mostCredibleReference=="H") dSamplePGX.REF<-dSamplePGX.H
    if(mostCredibleReference=="A") dSamplePGX.REF<-dSamplePGX.A
    if(mostCredibleReference=="B") dSamplePGX.REF<-dSamplePGX.B
    if(mostCredibleReference=="C") dSamplePGX.REF<-dSamplePGX.C

    colnames(dSamplePGX.REF)<-paste0(colnames(dSamplePGX.REF),".REF")
    dSamplePGX.merged <- dSamplePGX[dSamplePGX.REF, on=.(gene=gene.REF,chromosome_top=chromosome_top.REF)]
    allGeneCalls<-dSamplePGX.merged[,.(gene)]
    matchingGeneCalls<-dSamplePGX.merged[diplotype_name_top==diplotype_name_top.REF,.(gene)]
    matchingCond<-unlist(allGeneCalls) %in% unlist(matchingGeneCalls)

    allGeneCalls.columns<-paste0("GC_",unlist(allGeneCalls))
    pgxrexObj$sampleMeta[iSample,allGeneCalls.columns]<-0
    pgxrexObj$sampleMeta[iSample,allGeneCalls.columns[matchingCond]]<-1

    cat(paste0(" 🧬"))

  }

  shru::writeFile(pgxrexObj$sampleMeta,file=file.path(projectFolderPath,"work","pradaApp","per-sample-analysis","samples.tsv"),nThreads = 5)

}



#plots

##read back samples
dSample<-fread(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis","samples.tsv"))

##add back evaluationRatio - this is added in above now
# dSample[mostCredibleReference=="H",evaluationRatio:=(mConcordantREF.H+1000*mConcordantALT.H-mDiscordantREF.H-1000*mDiscordantALT.H)/mVCF.H]
# dSample[mostCredibleReference=="A",evaluationRatio:=(mConcordantREF.A+1000*mConcordantALT.A-mDiscordantREF.A-1000*mDiscordantALT.A)/mVCF.A]
# dSample[mostCredibleReference=="B",evaluationRatio:=(mConcordantREF.B+1000*mConcordantALT.B-mDiscordantREF.B-1000*mDiscordantALT.B)/mVCF.B]
# dSample[mostCredibleReference=="C",evaluationRatio:=(mConcordantREF.C+1000*mConcordantALT.C-mDiscordantREF.C-1000*mDiscordantALT.C)/mVCF.C]
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

#
# #data format of sdepth vars
# dSample.unique$sdepth_q050.CYP2B6<-as.numeric(dSample.unique$sdepth_q050.CYP2B6)
# dSample.unique$sdepth_q050.CYP2C19<-as.numeric(dSample.unique$sdepth_q050.CYP2C19)
# dSample.unique$sdepth_q050.CYP2D6<-as.numeric(dSample.unique$sdepth_q050.CYP2D6)
#
# #cumulative sums
# dSample.unique<-dSample.unique[order(sdepth_q050.CYP2B6,MERGEID),]
# #dSample.unique$GC_CYP2B6.cum<-cumsum(dSample.unique$GC_CYP2B6)
# #dSample.unique$sdepth_q050.CYP2B6<-as.factor(dSample.unique$sdepth_q050.CYP2B6)
# dSample.unique<-dSample.unique[order(sdepth_q050.CYP2C19,MERGEID),]
# #dSample.unique$GC_CYP2C19.cum<-cumsum(dSample.unique$GC_CYP2C19)
# #dSample.unique$sdepth_q050.CYP2C19<-as.factor(dSample.unique$sdepth_q050.CYP2C19)
# dSample.unique<-dSample.unique[order(sdepth_q050.CYP2D6,MERGEID),]
# #dSample.unique$GC_CYP2D6.cum<-cumsum(dSample.unique$GC_CYP2D6)
# #dSample.unique$sdepth_q050.CYP2D6<-as.factor(dSample.unique$sdepth_q050.CYP2D6)

library(ggplot2)
library(ggrepel)

ggplot(dSample.unique, aes(x=sdepth_q050_bed, y=nCalls, color=sample, label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsAll.png")))

ggplot(dSample.unique, aes(x=sdepth_q050_bed, y=nCalls.dep, color=sample, label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsAntidepressant.png")))

#custom CYP2B6

ggplot(dSample.unique, aes(x=sdepth_q050.CYP2B6, y=GC_CYP2B6, color=sample, label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsCYP2B6.png")))

dSample.unique<-dSample.unique[order(sdepth_q050.CYP2B6,MERGEID),]
ggplot(dSample.unique, aes(x= sdepth_q050.CYP2B6, y=cumsum(GC_CYP2B6), label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsCYP2B6.cum.png")))

dSample.unique<-dSample.unique[order(-sdepth_q050.CYP2B6,MERGEID),]
ggplot(dSample.unique, aes(x= cumsum(abs(GC_CYP2B6-1)), y=cumsum(GC_CYP2B6), label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsCYP2B6.roc.png")))

#custom CYP2C19
ggplot(dSample.unique, aes(x=sdepth_q050.CYP2C19, y=GC_CYP2C19, color=sample, label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsCYP2C19.png")))

dSample.unique<-dSample.unique[order(sdepth_q050.CYP2C19,MERGEID),]
ggplot(dSample.unique, aes(x=sdepth_q050.CYP2C19, y=cumsum(GC_CYP2C19), label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsCYP2C19.cum.png")))

dSample.unique<-dSample.unique[order(-sdepth_q050.CYP2C19,MERGEID),]
ggplot(dSample.unique, aes(x= cumsum(abs(GC_CYP2C19-1)), y=cumsum(GC_CYP2C19), label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsCYP2C19.roc.png")))


#custom CYP2D6
ggplot(dSample.unique, aes(x=sdepth_q050.CYP2D6, y=GC_CYP2D6, color=sample, label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsCYP2D6.png")))

dSample.unique<-dSample.unique[order(sdepth_q050.CYP2D6,MERGEID),]
ggplot(dSample.unique, aes(x=sdepth_q050.CYP2D6, y=cumsum(GC_CYP2D6), label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsCYP2D6.cum.png")))

dSample.unique<-dSample.unique[order(-sdepth_q050.CYP2D6,MERGEID),]
ggplot(dSample.unique, aes(x= cumsum(abs(GC_CYP2D6-1)), y=cumsum(GC_CYP2D6), label=label)) +
  geom_point() +
  geom_line() +
  geom_text_repel(size = 2) +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCallsCYP2D6.roc.png")))
#
# for(iGene in 1:length(colnamesAllGeneCalls)){
#   #iGene<-1
#   cCol<-colnamesAllGeneCalls[iGene]
#   ggplot(dSample.unique, aes(x=sdepth_q050_bed, y=UQ(as.name(cCol)), color=sample, label=label)) +
#     geom_point() +
#     geom_line() +
#     geom_text_repel(size = 2) +
#     theme_light()
#   ggsave(file.path(projectFolderPath,"work","pradaApp","per-sample-analysis",paste0("trueCalls_",cCol,".png")))
# }
