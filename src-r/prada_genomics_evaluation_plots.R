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

#genetic scores from pgxrex_pgx_score_simulation
dScore<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.individual_drug.txt")))


#Score density plots

#select score
dScore[,score:=wscore]
#dScore[,score:=score+1] #shift score +1

setorder(dScore, score, weight, id)
dScore[,order:=.I]

lDrug<-c("amitriptyline","citalopram","escitalopram","fluvoxamine","paroxetine","sertraline","venlafaxine","vortioxetine")


ggplot(data = dScore, aes(x=score)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  labs(y = "Density",
       x = "Score"
       #title = ""
  ) +
  theme_light() +
  geom_density()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_density.png")), units = "mm", width = 300, height = 150, scale = 0.40)

dplyr::tally(
  dplyr::group_by(
    dScore,
    score
    #ancestry
  )
)

ggplot(data = dScore, aes(x=weight)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  labs(y = "Density",
       x = "Weight"
       #title = ""
  ) +
  theme_light() +
  geom_density()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("weight_density.png")), units = "mm", width = 300, height = 150, scale = 0.40)

#per drug
for(iDrug in 1:length(lDrug)){
  #iDrug<-1
  sDrug<-lDrug[iDrug]

  ggplot(data = dScore[drug_name==sDrug,], aes(x=score)) + #[drug_name=="sertraline",]
    #geom_density(position="stack")
    #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
    labs(y = "Density",
         x = "Score",
         #color="Ancestry"
         #title = ""
    ) +
    theme_light() +
    geom_density()
  ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_density.",sDrug,".png")), units = "mm", width = 300, height = 150, scale = 0.40)

}


ggplot(data = dScore, aes(x=score, group=ancestry, color=ancestry, fill = ancestry)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  labs(y = "Density (relative)",
       x = "Score"
       #title = ""
  ) +
  theme_light() +
  geom_density(position="fill") +
  scale_fill_discrete(name='Ancestry',labels=c("Sub-Saharan African", "African American/Afro-Caribbean","Latino","Central/South Asian","East Asian","European","Near Eastern")) +
  scale_color_discrete(name='Ancestry',labels=c("Sub-Saharan African", "African American/Afro-Caribbean","Latino","Central/South Asian","East Asian","European","Near Eastern"))
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("ancestry_score_density.png")), units = "mm", width = 300, height = 150, scale = 0.40)


#per drug
lDrug<-c("amitriptyline","citalopram","escitalopram","fluvoxamine","paroxetine","sertraline","venlafaxine","vortioxetine")
for(iDrug in 1:length(lDrug)){
  #iDrug<-1
  sDrug<-lDrug[iDrug]

  ggplot(data = dScore[drug_name==sDrug,], aes(x=score, group=ancestry, color=ancestry, fill = ancestry)) + #[drug_name=="sertraline",]
    #geom_density(position="stack")
    #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
    labs(y = "Density (relative)",
         x = "Score",
         #color="Ancestry"
         #title = ""
    ) +
    theme_light() +
    geom_density(position="fill") +
    scale_fill_discrete(name='Ancestry',labels=c("Sub-Saharan African", "African American/Afro-Caribbean","Latino","Central/South Asian","East Asian","European","Near Eastern")) +
    scale_color_discrete(name='Ancestry',labels=c("Sub-Saharan African", "African American/Afro-Caribbean","Latino","Central/South Asian","East Asian","European","Near Eastern"))
  ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("ancestry_score_density.",sDrug,".png")), units = "mm", width = 300, height = 150, scale = 0.40)

}

# analysis of relative drug priority difference across scores
dScore<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.individual_drug.txt")))
dPheScore<-shru::readFile(file.path(projectFolderPath,"data","petrushka",paste0("final_predicted_score_petrushka_arm_(to KCL).csv")))
colnames(dPheScore)<-tolower(colnames(dPheScore))
lPheDrug<-tolower(colnames(dPheScore)[2:16])
lAllDrug<-unique(c(lPheDrug,lDrug))
dtAllDrug<-data.table(drug_name=lAllDrug)

dPheScore.long<-reshape(dPheScore,
                   direction = "long",
                   idvar = c("pid"),
                   varying = lPheDrug,
                   v.names = "score",
                   timevar = "drug_name",
                   times = lPheDrug
)

dPheScore.long[,score:=as.numeric(score)]

#make ranks - relative phenotypic score
lId<-unique(dPheScore.long$pid)
for(iId in 1:length(lId)){
  #iId<-1
  cId<-lId[iId]
  cScoreD<-dPheScore.long[pid==eval(cId),]
  cScoreD.score<-cScoreD[,
                         .(
                           id=head(.SD, 1)$pid,
                           ethnicity=head(.SD, 1)$ethnicity
                         )
                         ,by=c('score')]

  setorder(cScoreD.score, -score, na.last = T)
  cScoreD.score[,rank:=.I]

  cScoreD[cScoreD.score, on=c("score"),rank:=i.rank]

  dPheScore.long[cScoreD, on=.(pid=pid,drug_name=drug_name), rank:=i.rank]
}

#harmonise/map ethnicity/ancestry between the datasets
unique(dPheScore.long$ethnicity)
unique(dScore$ancestry)
dPheScore.long[ethnicity=="White/Caucasian",mapped_ancestry:="eur"]
dPheScore.long[ethnicity=="Other",mapped_ancestry:="cassas"]
dPheScore.long[ethnicity=="Asian",mapped_ancestry:="eas"]
dPheScore.long[ethnicity=="African/Caribbean",mapped_ancestry:="amcarafr"]


#add ancestry specific score means
dScore[,score.ancestry.mean:=mean(wscore,na.rm=T),by=c("ancestry")]

#add missing drugs
lSimId<-unique(dScore$id)
for(iSimId in 1:length(lSimId)){
  #iSimId<-1
  cSimId<-lSimId[iSimId]

  cScoreD<-dScore[id==eval(cSimId),]
  cScoreD.shared<-cScoreD[1,]
  cAllDrug<-dtAllDrug
  cAllDrug<-cAllDrug[cScoreD, on=c("drug_name"), match:=1][is.na(match),][,c("id","ancestry","score.ancestry.mean"):=list(eval(cSimId),eval(cScoreD.shared$ancestry),eval(cScoreD.shared$score.ancestry.mean))][,.(drug_name,id,ancestry,score.ancestry.mean)]
  dScore<-rbindlist(list(dScore,cAllDrug),use.names = T,fill = T)
}

#assign default values for missing scores
dScore[is.na(wscore),wscore:=score.ancestry.mean]
#dScore[is.na(wscore),wscore:=4]
dScore[is.na(wscore2),wscore2:=1]

dScore[,wscore.max:=max(wscore,na.rm=T), by=c('id')]
dScore[,wscore2.max:=max(wscore2,na.rm=T), by=c('id')]
dScore[!is.na(weight),wscore.max.gonly:=max(wscore,na.rm=T), by=c('id')]
dScore[!is.na(weight),wscore2.max.gonly:=max(wscore2,na.rm=T), by=c('id')]

#make ranks - wscore
lSimId<-unique(dScore$id)
for(iSimId in 1:length(lSimId)){
  #iSimId<-1
  cSimId<-lSimId[iSimId]
  cScoreD<-dScore[id==eval(cSimId),]
  cScoreD.score<-cScoreD[,
    .(
    id=head(.SD, 1)$id,
    ancestry=head(.SD, 1)$ancestry,
    wscore
    )
    ,by=c('wscore')]

  setorder(cScoreD.score, -wscore, na.last = T)
  cScoreD.score[,rank:=.I]

  cScoreD[cScoreD.score, on=c("wscore"),rank:=i.rank]

  dScore[cScoreD, on=.(id=id,drug_name=drug_name), rank:=i.rank]
}

#score differences as in the
dScore[,wscore.diff:=wscore-wscore.max]

#output long version of the diffscores and ranks
setorder(dScore, ancestry, id, drug_name, na.last = T)
shru::writeFile(dScore,file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.individual_drug.scoreranks.txt")))


#ancestral sumstats - this is the only thing we have to match the phenotypic scores with
dScore<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.individual_drug.scoreranks.txt")))
dScore.agg<-dScore[,.(
  #ancestry,
  #drug_name,
  rec_count = .N,
  avg_wscore=mean(wscore,na.rm=T),
  avg_rank=mean(rank,na.rm=T),
  max_wscore=max(wscore,na.rm=T),
  max_rank=max(rank,na.rm=T),
  min_wscore=min(wscore,na.rm=T),
  min_rank=min(rank,na.rm=T)
), by=c('ancestry','drug_name')]


#stats of the phenotypic scores in relation to the genetics based simulated sumstats
dPheScore.long[dScore.agg, on=.(mapped_ancestry=ancestry,drug_name=drug_name),c("avg_genetic_rank","rank.diff"):=list(i.avg_rank,rank-i.avg_rank)]
setorder(dPheScore.long, mapped_ancestry, pid, drug_name, na.last = T)
shru::writeFile(dPheScore.long,file.path(projectFolderPath,"data","petrushka",paste0("final_predicted_score_petrushka_arm.scorediff.txt")))

dPheScore.long.agg<-dPheScore.long[,.(
  #drug_name,
  rec_count = .N,
  avg_rank=mean(rank,na.rm=T),
  avg_genetic_rank=mean(avg_genetic_rank,na.rm=T),
  avg_rankdiff=mean(rank.diff,na.rm=T)
  ), by=c('mapped_ancestry','drug_name')]
setorder(dPheScore.long.agg, mapped_ancestry, drug_name, na.last = T)
shru::writeFile(dPheScore.long.agg,file.path(projectFolderPath,"data","petrushka",paste0("final_predicted_score_petrushka_arm.scorediff.agg.txt")))

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


