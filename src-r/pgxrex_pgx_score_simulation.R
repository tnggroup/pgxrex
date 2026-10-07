#Simulation of PGx genotypes with CPIC recommendations from pgxrex, 21/09/2026

#remotes::install_github("tnggroup/pgxrex")
#remotes::install_github("tnggroup/pgxrex",ref = 'jz_dev')
library(pgxrex)
library(data.table)


projectFolderPath<-"/scratch/prj/sgdp_nanopore/Projects/prada_jz"
#projectFolderPath<-"/Users/jakz/Documents/work_rstudio/pgxrex" #local

simpleReadTextFile<-function(path){
  con = file(path,open="r")
  res<-readLines(con = con,encoding = "UTF-8", warn = F)
  close(con)
  return(res)
}

readCpicPgxAndSimulateByAncestry<-function(ancestry="eur",nSim=2000){
  qString <- paste0("SELECT * FROM prada.harmonised_cpic_pgx_ancestry_",ancestry)
  q <- DBI::dbSendQuery(pgxrexObj$pgxrexApplicationDAO$connection,qString)
  res<-DBI::dbFetch(q)
  DBI::dbClearResult(q)

  setDT(res)
  dFreq<-res[genesymbol=="CYP2B6" | genesymbol=="CYP2C19" | genesymbol=="CYP2D6",]
  #View(dFreq)

  #Deal with missing data.
  #quantile(dFreq$consensus_allele_frequency,na.rm = T)
  #IMPUTE WITH 0 FREQUENCY
  dFreq[is.na(consensus_allele_frequency),consensus_allele_frequency:=0.0]
  #View(dFreq[genesymbol=="CYP2B6",])

  dFreq[,consensus_allele_frequency:=consensus_allele_frequency+1e-9] #add in baseline risk for everything to avoid 0 - will affect weighting equally

  #simulate n participants

  lGenes<-unique(dFreq$genesymbol)
  dSim<-as.data.frame(matrix(NA,0,0))

  for(iSimulatedIndividual in 1:nSim){
    #iSimulatedIndividual<-1
    dSim[iSimulatedIndividual,c("n")]<-c(iSimulatedIndividual)
    for(iGene in 1:length(lGenes)){
      #iGene<-1
      dGene<-dFreq[genesymbol==eval(lGenes[iGene]),]
      sDiplotype<-sample(x = dGene$diplotype,size = 1,prob = dGene$consensus_allele_frequency)
      dSim[iSimulatedIndividual,c(paste0("diplotype.",lGenes[iGene]))]<-c(sDiplotype)
    }

  }

  setDT(dSim)
  dSim
}


pgxrexObj<-PgxrexClass()
pgxrexObj$connectPgxrexDatabase(usernameToUse="tng_prada_system", dbnameToUse="prada_central_dev",passwordToUse=simpleReadTextFile("/users/k2481717/secure/tng_prada_system"))
#pgxrexObj$connectPgxrexDatabase(usernameToUse="tng_prada_system", dbnameToUse="prada_central_dev") #local

if(!file.exists(file.path(projectFolderPath,"data","pradaApp","simulated.eur.txt"))){
  dSim<-readCpicPgxAndSimulateByAncestry("eur")
  shru::writeFile(dSim,file.path(projectFolderPath,"data","pradaApp","simulated.eur.txt"))
}

if(!file.exists(file.path(projectFolderPath,"data","pradaApp","simulated.eas.txt"))){
  dSim<-readCpicPgxAndSimulateByAncestry("eas")
  shru::writeFile(dSim,file.path(projectFolderPath,"data","pradaApp","simulated.eas.txt"))
}

if(!file.exists(file.path(projectFolderPath,"data","pradaApp","simulated.cassas.txt"))){
  dSim<-readCpicPgxAndSimulateByAncestry("cassas")
  shru::writeFile(dSim,file.path(projectFolderPath,"data","pradaApp","simulated.cassas.txt"))
}

if(!file.exists(file.path(projectFolderPath,"data","pradaApp","simulated.afr.txt"))){
  dSim<-readCpicPgxAndSimulateByAncestry("afr")
  shru::writeFile(dSim,file.path(projectFolderPath,"data","pradaApp","simulated.afr.txt"))
}

if(!file.exists(file.path(projectFolderPath,"data","pradaApp","simulated.amlathis.txt"))){
  dSim<-readCpicPgxAndSimulateByAncestry("amlathis")
  shru::writeFile(dSim,file.path(projectFolderPath,"data","pradaApp","simulated.amlathis.txt"))
}

if(!file.exists(file.path(projectFolderPath,"data","pradaApp","simulated.amcarafr.txt"))){
  dSim<-readCpicPgxAndSimulateByAncestry("amcarafr")
  shru::writeFile(dSim,file.path(projectFolderPath,"data","pradaApp","simulated.amcarafr.txt"))
}

if(!file.exists(file.path(projectFolderPath,"data","pradaApp","simulated.gme.txt"))){
  dSim<-readCpicPgxAndSimulateByAncestry("gme")
  shru::writeFile(dSim,file.path(projectFolderPath,"data","pradaApp","simulated.gme.txt"))
}

dSim.eur<-shru::readFile(file.path(projectFolderPath,"data","pradaApp","simulated.eur.txt"))
dSim.eas<-shru::readFile(file.path(projectFolderPath,"data","pradaApp","simulated.eas.txt"))
dSim.cassas<-shru::readFile(file.path(projectFolderPath,"data","pradaApp","simulated.cassas.txt"))
dSim.afr<-shru::readFile(file.path(projectFolderPath,"data","pradaApp","simulated.afr.txt"))
dSim.amlathis<-shru::readFile(file.path(projectFolderPath,"data","pradaApp","simulated.amlathis.txt"))
dSim.amcarafr<-shru::readFile(file.path(projectFolderPath,"data","pradaApp","simulated.amcarafr.txt"))
dSim.gme<-shru::readFile(file.path(projectFolderPath,"data","pradaApp","simulated.gme.txt"))

#run pgxrex interpretation and scoring

#set ancestry dataset!!!!
dSim<-dSim.eur
labelAncestry<-"eur"

dSim<-as.data.frame(dSim)
#append to existing results
if(file.exists(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.",labelAncestry,".txt")))){
  dRes<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.",labelAncestry,".txt")))
} else dRes<-as.data.frame(matrix(NA,0,0))

lastIndex<-max(unique(dRes$i))

for(iSimulatedIndividual in 1:nrow(dSim)){
  #iSimulatedIndividual<-1
  if(iSimulatedIndividual<=lastIndex) next
  qString <- "DROP TABLE IF EXISTS t_gene_diplotype_input"
  q <- DBI::dbSendStatement(pgxrexObj$pgxrexApplicationDAO$connection,qString)
  DBI::dbClearResult(q)
  qString <- paste0("CREATE TEMP TABLE IF NOT EXISTS t_gene_diplotype_input AS
  SELECT 'CYP2D6' AS gene, '",dSim[iSimulatedIndividual,]$diplotype.CYP2D6,"' AS diplotype
	UNION ALL
	SELECT 'CYP2B6','",dSim[iSimulatedIndividual,]$diplotype.CYP2B6,"'
	UNION ALL
	SELECT 'CYP2C19', '",dSim[iSimulatedIndividual,]$diplotype.CYP2C19,"'")
  q <- DBI::dbSendQuery(pgxrexObj$pgxrexApplicationDAO$connection,qString)
  res<-DBI::dbFetch(q)
  DBI::dbClearResult(q)
  qString <- paste0("SELECT * FROM prada.get_application_recommendation()")
  q <- DBI::dbSendQuery(pgxrexObj$pgxrexApplicationDAO$connection,qString)
  res<-DBI::dbFetch(q)
  DBI::dbClearResult(q)
  res$i<-dSim[iSimulatedIndividual,]$n
  dRes<-rbindlist(list(dRes,res),use.names = T, fill=T)

  if(iSimulatedIndividual %% 10 == 0){
    shru::writeFile(dRes,file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.",labelAncestry,".txt")))
  }
  if(iSimulatedIndividual %% 10 == 0) catl("Progress: ",iSimulatedIndividual)
}
shru::writeFile(dRes,file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.",labelAncestry,".txt")))


#postprocess app recommendations

dRes.eur<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.eur.txt")))
dRes.eur[,ancestry:="eur"]
dRes.eas<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.eas.txt")))
dRes.eas[,ancestry:="eas"]
dRes.cassas<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.cassas.txt")))
dRes.cassas[,ancestry:="cassas"]
dRes.afr<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.afr.txt")))
dRes.afr[,ancestry:="afr"]
dRes.amlathis<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.amlathis.txt")))
dRes.amlathis[,ancestry:="amlathis"]
dRes.amcarafr<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.amcarafr.txt")))
dRes.amcarafr[,ancestry:="amcarafr"]
dRes.gme<-shru::readFile(file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.gme.txt")))
dRes.gme[,ancestry:="gme"]


dRes<-rbindlist(list(dRes.eur,dRes.eas,dRes.cassas,dRes.afr,dRes.amlathis,dRes.amcarafr,dRes.gme),use.names = T,fill = T)
dRes[,id:=paste0("",ancestry,i)]

lSimId<-unique(dRes$id)

#HERE!!! there is an error
dRes.individual_drug<-data.frame(matrix(NA,0,0))
for(iSimulatedIndividual in 1:length(lSimId)){
  #iSimulatedIndividual<-1
  idSimulatedIndividual<-lSimId[iSimulatedIndividual]
  cRes<-dRes[id==eval(idSimulatedIndividual),]
  #View(cRes)

  cRes.aggregate<-cRes[,.(
    id=head(.SD, 1)$id,
    ancestry=head(.SD, 1)$ancestry,
    prada_rec_count = .N,
    prada_start_dose=min(prada_start_dose,na.rm=T),
    prada_target_dose=min(prada_target_dose,na.rm=T),
    prada_titration_speed=min(prada_titration_speed,na.rm=T),
    prada_switch1_drug=max(prada_switch1_drug,na.rm=T),
    prada_switch1_gene=max(prada_switch1_gene,na.rm=T),
    prada_switch2_drug=max(prada_switch2_drug,na.rm=T),
    prada_switch2_gene=max(prada_switch2_gene,na.rm=T),
    prada_tdm=max(prada_tdm,na.rm=T),
    prada_avg_cpiclevel_num=mean(prada_cpiclevel_num,na.rm=T)
    ), by=c('drug_name')]

  cRes.aggregate[,score:=1.0*prada_start_dose+1.0*prada_target_dose+1.0*(prada_titration_speed/2)-1.0*prada_switch1_drug-1.0*prada_switch1_gene-0.5*prada_switch2_drug-0.5*prada_switch2_gene-1.0*prada_tdm][,weight:=(prada_avg_cpiclevel_num/8)] #no sequencing depth information

  cRes.aggregate[,id:=eval(idSimulatedIndividual)]

  dRes.individual_drug<-rbindlist(list(dRes.individual_drug,cRes.aggregate),use.names = T,fill = T)

}

shru::writeFile(dRes.individual_drug,file.path(projectFolderPath,"work","pradaApp","simulation",paste0("simulatedPgx.individual_drug.txt")))


#see the prada_genomics_evaluation_plots.R file for plots of these


#
# dRes.individual_drug[!is.finite(score),score:=NA_real_]
#
# #stats
# dRes.individual_drug$score.std<-scale(dRes.individual_drug[,c("score")],center = T,scale = T)
#
# quantile(dRes.individual_drug$score, na.rm = T)
# mean(dRes.individual_drug$score, na.rm = T)
# sd(dRes.individual_drug$score,na.rm = T)
# quantile(dRes.individual_drug$score.std, na.rm = T)
#
# #petrushka stats
# dPetrushka<-shru::readFile(file.path(projectFolderPath,"data","pradaApp","final_predicted_score_petrushka_arm.csv"))
# dPetrushka.datacolumnnames<-colnames(dPetrushka)[2:(ncol(dPetrushka)-1)]
# dPetrushka.d<-dPetrushka[,..dPetrushka.datacolumnnames]
# dPetrushka.d.allScores<-rowSums(dPetrushka.d,na.rm = T)
# dPetrushka.d.allScores<-data.table(score=dPetrushka.d.allScores)
# dPetrushka.d.allScores$score.std<-scale(dPetrushka.d.allScores[,c("score")],center = T,scale = T)
#
# quantile(dPetrushka.d.allScores$score, na.rm = T)
# mean(dPetrushka.d.allScores$score, na.rm = T)
# sd(dPetrushka.d.allScores$score,na.rm = T)
# quantile(dPetrushka.d.allScores$score.std, na.rm = T)
