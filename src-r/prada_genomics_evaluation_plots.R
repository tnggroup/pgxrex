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
  #geom_bar()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_density.png")), units = "mm", width = 300, height = 150, scale = 0.40)

ggplot(data = dScore, aes(x=score)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  labs(y = "Count",
       x = "Score"
       #title = ""
  ) +
  theme_light() +
  #geom_density()
  geom_bar()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_bar.png")), units = "mm", width = 300, height = 150, scale = 0.40)

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


#per drug and ancestry - faceted
ggplot(data = dScore, aes(x=score)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Density",
       x = "Score",
       #color="Ancestry"
       #title = ""
  ) +
  facet_grid(drug_name ~ ancestry, scales = "free_y") +
  theme_light() +
  #geom_bar()
  geom_density()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_density.facet.png")), units = "mm", width = 300, height = 150, scale = 1.25)

ggplot(data = dScore, aes(x=score)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Count",
       x = "Score",
       #color="Ancestry"
       #title = ""
  ) +
  facet_grid(drug_name ~ ancestry) + #scales = "free_y"
  theme_light() +
  geom_bar()
  #geom_density()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_bar.facet.png")), units = "mm", width = 300, height = 150, scale = 1.25)

#per drug - faceted
ggplot(data = dScore, aes(x=score)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Density",
       x = "Score",
       #color="Ancestry"
       title = "Per drug"
  ) +
  facet_grid(rows = "drug_name", scales = "free_y") +
  theme_light() +
  #geom_bar()
  geom_density()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_density.facet.drug.png")), units = "mm", width = 200, height = 300, scale = 0.75)

ggplot(data = dScore, aes(x=score)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Count",
       x = "Score",
       #color="Ancestry"
       title = "Per drug"
  ) +
  facet_grid(rows = "drug_name") + #scales = "free_y"
  theme_light() +
  geom_bar()
  #geom_density()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_bar.facet.drug.png")), units = "mm", width = 200, height = 300, scale = 0.75)

#per ancestry - faceted
ggplot(data = dScore, aes(x=score)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Density",
       x = "Score",
       #color="Ancestry"
       title = "Per ancestry"
  ) +
  facet_grid(rows = "ancestry", scales = "free_y") +
  theme_light() +
  #geom_bar()
  geom_density()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_density.facet.ancestry.png")), units = "mm", width = 200, height = 300, scale = 0.75)

ggplot(data = dScore, aes(x=score)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Count",
       x = "Score",
       #color="Ancestry"
       title = "Per ancestry"
  ) +
  facet_grid(rows = "ancestry") + #scales = "free_y"
  theme_light() +
  geom_bar()
#geom_density()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_bar.facet.ancestry.png")), units = "mm", width = 200, height = 300, scale = 0.75)


#
# #old per drug
# for(iDrug in 1:length(lDrug)){
#   #iDrug<-1
#   sDrug<-lDrug[iDrug]
#
#   ggplot(data = dScore[drug_name==sDrug,], aes(x=score)) + #[drug_name=="sertraline",]
#     #geom_density(position="stack")
#     #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
#     labs(y = "Density",
#          x = "Score",
#          #color="Ancestry"
#          #title = ""
#     ) +
#     theme_light() +
#     geom_density()
#   ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("score_density.",sDrug,".png")), units = "mm", width = 300, height = 150, scale = 0.40)
#
# }


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

#  *****analysis of relative drug priority difference across scores******
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

#  ***** END analysis of relative drug priority difference across scores******

#  *****simulation of genetics in phenotypic participants******

dScore<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.individual_drug.txt")))
dScore.individual<-dScore[,.(
  ancestry=head(.SD, 1)$ancestry
)
,by=c('id')]


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



#add missing drugs
lSimId<-unique(dScore$id)
for(iSimId in 1:length(lSimId)){
  #iSimId<-1
  cSimId<-lSimId[iSimId]

  cScoreD<-dScore[id==eval(cSimId),]
  cScoreD.shared<-cScoreD[1,]
  cAllDrug<-dtAllDrug
  cAllDrug<-cAllDrug[cScoreD, on=c("drug_name"), match:=1][is.na(match),][,c("id","ancestry"):=list(eval(cSimId),eval(cScoreD.shared$ancestry))][,.(drug_name,id)]
  dScore<-rbindlist(list(dScore,cAllDrug),use.names = T,fill = T)
}

#assign default values for missing scores
dScore[is.na(wscore),wscore:=4] #i.e. no effect of the genetic score
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

#score differences as in the petrushka file
dScore[,wscore.diff:=wscore-wscore.max]


simulated.petrushka<-as.data.table(matrix(NA,0,0))

#simulate/assign/map to simulated gentotypes
lId<-unique(dPheScore.long$pid)
setkeyv(dPheScore.long,cols=c("pid"))
for(isim in 1:100){
  #isim<-1
  sample.eur<-data.table(pid=unique(dPheScore.long[mapped_ancestry=="eur",]$pid), gid=sample(x = dScore.individual[ancestry=="eur",]$id,size = length(unique(dPheScore.long[mapped_ancestry=="eur",]$pid))))
  sample.cassas<-data.table(pid=unique(dPheScore.long[mapped_ancestry=="cassas",]$pid), gid=sample(x = dScore.individual[ancestry=="cassas",]$id,size = length(unique(dPheScore.long[mapped_ancestry=="cassas",]$pid))))
  sample.eas<-data.table(pid=unique(dPheScore.long[mapped_ancestry=="eas",]$pid), gid=sample(x = dScore.individual[ancestry=="eas",]$id,size = length(unique(dPheScore.long[mapped_ancestry=="eas",]$pid))))
  sample.amcarafr<-data.table(pid=unique(dPheScore.long[mapped_ancestry=="amcarafr",]$pid), gid=sample(x = dScore.individual[ancestry=="amcarafr",]$id,size = length(unique(dPheScore.long[mapped_ancestry=="amcarafr",]$pid))))
  simulated.petrushka<-rbindlist(
  list(
    simulated.petrushka,
    rbindlist(list(sample.eur,sample.cassas,sample.eas,sample.amcarafr))[,sim:=eval(isim)]
  ), use.names = T,fill = T)
}
setkeyv(simulated.petrushka,cols=c("sim","pid","gid"))
simulated.petrushka2<-simulated.petrushka[dPheScore.long, on=.(pid=pid), allow.cartesian=TRUE]

setkeyv(simulated.petrushka2,cols=c("sim","pid","gid","drug_name"))
setkeyv(dScore,cols=c("id","drug_name"))
simulated.petrushka2[dScore[,.(drug_name,id,ancestry,wscore,wscore.max,wscore.diff,rank)], on=.(gid=id,drug_name=drug_name),c("gid","gscore","gscore.max","gscore.diff","grank"):=list(i.id,i.wscore,i.wscore.max,i.wscore.diff,i.rank)]
setorder(simulated.petrushka2,sim,pid,gid,mapped_ancestry,drug_name)

simulated.petrushka2[,rank.diff:=rank-grank]


#stats of the phenotypic scores in relation to the random genetics based simulated scores

shru::writeFile(simulated.petrushka2,file.path(projectFolderPath,"data","petrushka",paste0("final_predicted_score_petrushka_arm.scorediff.txt")))

#read back
simulated.petrushka2<-shru::readFile(file.path(projectFolderPath,"data","petrushka",paste0("final_predicted_score_petrushka_arm.scorediff.txt")))
setorder(simulated.petrushka2,pid,gid,mapped_ancestry,drug_name)

#counts to get percentages
simulated.petrushka2[,count_ancestry_drug := .N, by=c('mapped_ancestry','drug_name')]

pscore_mean<-mean(simulated.petrushka2$score,na.rm = T)
pscore_sd<-sd(simulated.petrushka2$score,na.rm = T)
gscore_mean<-mean(simulated.petrushka2$gscore,na.rm = T)
gscore_sd<-sd(simulated.petrushka2$gscore,na.rm = T)

simulated.petrushka2[,cscore0p5:=((score-eval(pscore_mean))/eval(pscore_sd)+0.5*(gscore-eval(gscore_mean))/eval(gscore_sd))/1.5]
simulated.petrushka2[,cscore1:=((score-eval(pscore_mean))/eval(pscore_sd)+1*(gscore-eval(gscore_mean))/eval(gscore_sd))/2]
simulated.petrushka2[,cscore1p5:=((score-eval(pscore_mean))/eval(pscore_sd)+1.5*(gscore-eval(gscore_mean))/eval(gscore_sd))/2.5]

#make ranks - combined scores
setorder(simulated.petrushka2,sim,pid,gid,-cscore0p5,na.last = T)
simulated.petrushka2[,rank_cscore0p5:=seq_len(.N), by=c("sim","pid","gid")][is.na(cscore0p5),rank_cscore0p5:=NA_integer_]
setorder(simulated.petrushka2,sim,pid,gid,-cscore1,na.last = T)
simulated.petrushka2[,rank_cscore1:=seq_len(.N), by=c("sim","pid","gid")][is.na(cscore1),rank_cscore1:=NA_integer_]
setorder(simulated.petrushka2,sim,pid,gid,-cscore1p5,na.last = T)
simulated.petrushka2[,rank_cscore1p5:=seq_len(.N), by=c("sim","pid","gid")][is.na(cscore1p5),rank_cscore1p5:=NA_integer_]

#rank.diff.limit.lt<-2

#rank groups
#simulated.petrushka2[,match_fuzzy1:=abs(rank.diff)<2]
simulated.petrushka2[,mismatch_fuzzy1_gsuperior:=abs(rank.diff)>=2 & rank.diff<1] #genetics is ranked higher
simulated.petrushka2[,mismatch_fuzzy1_ginferior:=abs(rank.diff)>=2 & rank.diff>-1] #genetics is ranked lower
simulated.petrushka2[,match_fuzzy1_group:="Match"][mismatch_fuzzy1_gsuperior==TRUE,match_fuzzy1_group:="Genetics ranks higher"][mismatch_fuzzy1_ginferior==TRUE,match_fuzzy1_group:="Genetics ranks lower"]

simulated.petrushka2[,rank.diff.cscore0p5:=rank-rank_cscore0p5]
simulated.petrushka2[,rank.diff.cscore0p5_gsuperior:=abs(rank.diff.cscore0p5)>=2 & rank.diff.cscore0p5<1] #genetics is ranked higher
simulated.petrushka2[,rank.diff.cscore0p5_ginferior:=abs(rank.diff.cscore0p5)>=2 & rank.diff.cscore0p5>-1] #genetics is ranked lower
simulated.petrushka2[,rank.diff.cscore0p5_group:="Match"][rank.diff.cscore0p5_gsuperior==TRUE,rank.diff.cscore0p5_group:="Composite ranks higher"][rank.diff.cscore0p5_ginferior==TRUE,rank.diff.cscore0p5_group:="Composite ranks lower"]

simulated.petrushka2[,rank.diff.cscore1:=rank-rank_cscore1]
simulated.petrushka2[,rank.diff.cscore1_gsuperior:=abs(rank.diff.cscore1)>=2 & rank.diff.cscore1<1] #genetics is ranked higher
simulated.petrushka2[,rank.diff.cscore1_ginferior:=abs(rank.diff.cscore1)>=2 & rank.diff.cscore1>-1] #genetics is ranked lower
simulated.petrushka2[,rank.diff.cscore1_group:="Match"][rank.diff.cscore1_gsuperior==TRUE,rank.diff.cscore1_group:="Composite ranks higher"][rank.diff.cscore1_ginferior==TRUE,rank.diff.cscore1_group:="Composite ranks lower"]

simulated.petrushka2[,rank.diff.cscore1p5:=rank-rank_cscore1p5]
simulated.petrushka2[,rank.diff.cscore1p5_gsuperior:=abs(rank.diff.cscore1p5)>=2 & rank.diff.cscore1p5<1] #genetics is ranked higher
simulated.petrushka2[,rank.diff.cscore1p5_ginferior:=abs(rank.diff.cscore1p5)>=2 & rank.diff.cscore1p5>-1] #genetics is ranked lower
simulated.petrushka2[,rank.diff.cscore1p5_group:="Match"][rank.diff.cscore1p5_gsuperior==TRUE,rank.diff.cscore1p5_group:="Composite ranks higher"][rank.diff.cscore1p5_ginferior==TRUE,rank.diff.cscore1p5_group:="Composite ranks lower"]


shru::writeFile(simulated.petrushka2,file.path(projectFolderPath,"data","petrushka",paste0("final_predicted_score_petrushka_arm.composite.txt")))


# #test
# quantile(simulated.petrushka2[drug_name=="mirtazapine",]$gscore, na.rm = T)
# quantile(simulated.petrushka2$gscore, na.rm = T)

simulated.petrushka2[,count_ancestry_drug := .N, by=c('mapped_ancestry','drug_name')]

simulated.petrushka2.agg<-simulated.petrushka2[,.(
  #drug_name,
  rec_count = .N,
  #avg_rank=mean(rank,na.rm=T),
  #avg_genetic_rank=mean(grank,na.rm=T),
  avg_rankdiff = mean(rank.diff,na.rm=T),
  avg_rankdiff_abs = mean(abs(rank.diff),na.rm=T),
  match_count = sum(rank.diff==0),
  match_count_fuzzy1 = sum(abs(rank.diff)<2),
  match_count_fuzzy2 = sum(abs(rank.diff)<3)
), by=c('mapped_ancestry','sim','pid','gid')]
setorder(simulated.petrushka2.agg, mapped_ancestry, pid, gid, na.last = T)
#shru::writeFile(simulated.petrushka2.agg,file.path(projectFolderPath,"data","petrushka",paste0("final_predicted_score_petrushka_arm.scorediff.agg.txt")))


ggplot(data = simulated.petrushka2, aes(x=rank.diff)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Density",
       x = "Delta Rank",
       #color="Ancestry"
       #title = ""
  ) +
  facet_grid(drug_name ~ mapped_ancestry, scales = "free_y") + #scales = "free_y"
  theme_light() +
  #geom_bar()
  geom_density() #aes(y = after_stat(count))
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_density_rankdiff.facet.png")), units = "mm", width = 300, height = 150, scale = 1.25)


ggplot(data = simulated.petrushka2, aes(x=abs(rank.diff))) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Density",
       x = "Delta Rank",
       #color="Ancestry"
       #title = ""
  ) +
  facet_grid(drug_name ~ mapped_ancestry, scales = "free_y") + #scales = "free_y"
  theme_light() +
  #geom_bar()
geom_density()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_density_rankdiff_abs.facet.png")), units = "mm", width = 300, height = 150, scale = 1.25)


ggplot(data = simulated.petrushka2, aes(x=score,y=gscore)) + #[drug_name=="sertraline",]
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Genetic score",
       x = "Phenotypic score",
       #color="Ancestry"
       #title = ""
  ) +
  facet_grid(drug_name ~ mapped_ancestry) + #scales = "free_y"
  theme_light() +
  #geom_bar()
  geom_point()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_density_scores.facet.png")), units = "mm", width = 300, height = 300, scale = 1.25)



setting.ranks.lt<-4

drugsToIncludeInPlots<-unique(simulated.petrushka2$drug_name)
#drugsToIncludeInPlots<-lDrug

#rank changes bars - genetict vs phenotypic
label_data <- simulated.petrushka2[
  drug_name %in% drugsToIncludeInPlots & (rank<eval(setting.ranks.lt) | grank<eval(setting.ranks.lt)),
  .(n = .N),
  by = .(mapped_ancestry, drug_name, match_fuzzy1_group)
]
label_data[, n_ancestry_drug_name := sum(n),
           by = .(mapped_ancestry, drug_name)]

ggplot(data = simulated.petrushka2[
  drug_name %in% drugsToIncludeInPlots & (rank<eval(setting.ranks.lt) | grank<eval(setting.ranks.lt)),], aes(x=drug_name,fill=match_fuzzy1_group)) + #[drug_name=="sertraline",]
  scale_fill_manual(values = c(
    "Genetics ranks higher" = "gold1",
    "Genetics ranks lower" = "coral1",
    "Match" = "green2"
  )) +
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Count",
       x = "Drug",
       #color="Ancestry"
       #title = ""
  ) +
  facet_grid(rows = "mapped_ancestry", scales = "free_y") + #scales = "free_y"
  theme_light() +
  #geom_bar()
  geom_bar()+
  geom_text(
    data = label_data,
    aes(
      y = n,
      label = scales::percent(n / n_ancestry_drug_name,accuracy = 1)
    ),
    position = position_stack(vjust = 0.5)
  )
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_bar_changes.facet.png")), units = "mm", width = 300, height = 300, scale = 1.00)


#rank changes bars - combined 0p5 vs phenotypic
label_data <- simulated.petrushka2[
  drug_name %in% drugsToIncludeInPlots & (rank<eval(setting.ranks.lt) | rank_cscore0p5<eval(setting.ranks.lt)),
  .(n = .N),
  by = .(mapped_ancestry, drug_name, rank.diff.cscore0p5_group)
]
label_data[, n_ancestry_drug_name := sum(n),
           by = .(mapped_ancestry, drug_name)]

ggplot(data = simulated.petrushka2[
  drug_name %in% drugsToIncludeInPlots & (rank<eval(setting.ranks.lt) | rank_cscore0p5<eval(setting.ranks.lt)),], aes(x=drug_name,fill=rank.diff.cscore0p5_group)) + #[drug_name=="sertraline",]
  scale_fill_manual(values = c(
    "Composite ranks higher" = "gold1",
    "Composite ranks lower" = "coral1",
    "Match" = "green2"
  )) +
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Count",
       x = "Drug",
       #color="Ancestry"
       #title = ""
  ) +
  facet_grid(rows = "mapped_ancestry", scales = "free_y") + #scales = "free_y"
  theme_light() +
  #geom_bar()
  geom_bar()+
  geom_text(
    data = label_data,
    aes(
      y = n,
      label = scales::percent(n / n_ancestry_drug_name,accuracy = 1)
    ),
    position = position_stack(vjust = 0.5)
  )
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_bar_changes.combo0p5facet.png")), units = "mm", width = 300, height = 300, scale = 1.00)


#rank changes bars - combined 1 vs phenotypic
label_data <- simulated.petrushka2[
  drug_name %in% drugsToIncludeInPlots & (rank<eval(setting.ranks.lt) | rank_cscore1<eval(setting.ranks.lt)),
  .(n = .N),
  by = .(mapped_ancestry, drug_name, rank.diff.cscore1_group)
]
label_data[, n_ancestry_drug_name := sum(n),
           by = .(mapped_ancestry, drug_name)]

ggplot(data = simulated.petrushka2[
  drug_name %in% drugsToIncludeInPlots & (rank<eval(setting.ranks.lt) | rank_cscore1<eval(setting.ranks.lt)),], aes(x=drug_name,fill=rank.diff.cscore1_group)) + #[drug_name=="sertraline",]
  scale_fill_manual(values = c(
    "Composite ranks higher" = "gold1",
    "Composite ranks lower" = "coral1",
    "Match" = "green2"
  )) +
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Count",
       x = "Drug",
       #color="Ancestry"
       #title = ""
  ) +
  facet_grid(rows = "mapped_ancestry", scales = "free_y") + #scales = "free_y"
  theme_light() +
  #geom_bar()
  geom_bar()+
  geom_text(
    data = label_data,
    aes(
      y = n,
      label = scales::percent(n / n_ancestry_drug_name,accuracy = 1)
    ),
    position = position_stack(vjust = 0.5)
  )
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_bar_changes.combo1facet.png")), units = "mm", width = 300, height = 300, scale = 1.00)


#rank changes bars - combined 1p5 vs phenotypic
label_data <- simulated.petrushka2[
  drug_name %in% drugsToIncludeInPlots & (rank<eval(setting.ranks.lt) | rank_cscore1p5<eval(setting.ranks.lt)),
  .(n = .N),
  by = .(mapped_ancestry, drug_name, rank.diff.cscore1p5_group)
]
label_data[, n_ancestry_drug_name := sum(n),
           by = .(mapped_ancestry, drug_name)]

ggplot(data = simulated.petrushka2[
  drug_name %in% drugsToIncludeInPlots & (rank<eval(setting.ranks.lt) | rank_cscore1p5<eval(setting.ranks.lt)),], aes(x=drug_name,fill=rank.diff.cscore1p5_group)) + #[drug_name=="sertraline",]
  scale_fill_manual(values = c(
    "Composite ranks higher" = "gold1",
    "Composite ranks lower" = "coral1",
    "Match" = "green2"
  )) +
  #geom_density(position="stack")
  #scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
  labs(y = "Count",
       x = "Drug",
       #color="Ancestry"
       #title = ""
  ) +
  facet_grid(rows = "mapped_ancestry", scales = "free_y") + #scales = "free_y"
  theme_light() +
  #geom_bar()
  geom_bar()+
  geom_text(
    data = label_data,
    aes(
      y = n,
      label = scales::percent(n / n_ancestry_drug_name,accuracy = 1)
    ),
    position = position_stack(vjust = 0.5)
  )
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_bar_changes.combo1p5facet.png")), units = "mm", width = 300, height = 300, scale = 1.00)




plotdat<-tidyr::pivot_longer(simulated.petrushka2[drug_name %in% lDrug,],
  cols = c(cscore0p5,cscore1,cscore1p5),
  names_to = "Score",
  values_to = "Value"
)
plotdat<-as.data.table(as.data.frame(plotdat))


ggplot(plotdat,
       aes(x = score,
           y = Value,
           colour = Score)) +
  geom_count() +
  scale_size_area(max_size = 12) +
  labs(
    x = "PETRUSHKA score",
    y = "Combined score"
  ) +
  theme_light()
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_sumscore_count.png")), units = "mm", width = 300, height = 300, scale = 1.00)


ggplot(plotdat,
       aes(x = factor(score),
           y = Value,
           colour = Score)) +
  geom_count() +
  scale_size_area(max_size = 12) +
  facet_wrap(drug_name ~ mapped_ancestry) +
  theme_light() +
  labs(
    x = "PETRUSHKA score",
    y = "Combined score",
    size = "Count"
  )
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_sumscore_count.facet.png")), units = "mm", width = 300, height = 300, scale = 1.00)

ggplot(plotdat,
       aes(x = factor(score),
           y = Value,
           colour = Score)) +
  geom_count() +
  scale_size_area(max_size = 12) +
  facet_wrap(~ drug_name) +
  theme_light() +
  labs(
    x = "PETRUSHKA score",
    y = "Combined score",
    size = "Count"
  )
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_sumscore_count.facet.drug.png")), units = "mm", width = 300, height = 300, scale = 1.00)

ggplot(plotdat,
       aes(x = factor(score),
           y = Value,
           colour = Score)) +
  geom_count() +
  scale_size_area(max_size = 12) +
  facet_wrap(~ mapped_ancestry) +
  theme_light() +
  labs(
    x = "PETRUSHKA score",
    y = "Combined score",
    size = "Count"
  )
ggsave(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("petrushka_sim_sumscore_count.facet.ancestry.png")), units = "mm", width = 300, height = 300, scale = 1.00)




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


