
library(pgxrex)
library(data.table)

settingProjectFolder<-"/Users/jakz/Documents/work_rstudio/pgxrex"

pradaO<-PgxrexClass()
pradaO$connectPgxrexDatabase(usernameToUse="tng_prada_system", dbnameToUse="prada_central_dev")
pradaO$computeGenomeCoverage( writeToThisBedPath = file.path(settingProjectFolder,"data","bed","pgx.grch38.5k.0p7percent.bed"),nPrioritisedCnv=0, nPrioritisedSnp=0, nPrioritisedTotal = 5000)
pradaO$computeGenomeCoverage( writeToThisBedPath = file.path(settingProjectFolder,"data","bed","pgx_cnv.grch38.5k.2p1percent.bed"),nPrioritisedSnp=0, nPrioritisedTotal = 5000)
#pradaO$computeGenomeCoverage( writeToThisBedPath = "pgx_cnv_mddeur.grch38.5k.1p3percent.bed",nPrioritisedTotal = 5000) #1e6 CNV weighting
pradaO$computeGenomeCoverage( writeToThisBedPath = file.path(settingProjectFolder,"data","bed","pgx_cnv_mddeur.grch38.5k.2p6percent.bed"),nPrioritisedTotal = 5000)
pradaO$computeGenomeCoverage( writeToThisBedPath = file.path(settingProjectFolder,"data","bed","pgx_cnv_mddeur.grch38.25k.4p5percent.bed"),nPrioritisedTotal = 25000)
pradaO$computeGenomeCoverage( writeToThisBedPath = file.path(settingProjectFolder,"data","bed","pgx_cnv_mddeur.grch38.50k.7p3percent.bed"),nPrioritisedTotal = 50000)
#pradaO$computeGenomeCoverage( writeToThisBedPath = file.path(settingProjectFolder,"data","bed","pgx_cnv_mddeur.grch38.75k.10p0percent.bed"),nPrioritisedTotal = 75000)
#pradaO$computeGenomeCoverage( writeToThisBedPath = file.path(settingProjectFolder,"data","bed","pgx_cnv_mddeur.grch38.100k.12p8percent.bed"),nPrioritisedTotal = 100000)


pradaO$computeGenomeCoverage( writeToThisBedPath = file.path(settingProjectFolder,"data","bed","pgx.grch38.padded_100k.bed"),nPrioritisedCnv=0, nPrioritisedSnp=0, nPrioritisedTotal = 1000, writePaddedStrands = T, paddingGeneBp = 100000)

pradaO$computeGenomeCoverage( writeToThisBedPath = file.path(settingProjectFolder,"data","bed","pgx.grch38.padded_100k.ss.bed"),nPrioritisedCnv=0, nPrioritisedSnp=0, nPrioritisedTotal = 1000, writePaddedStrands = F, paddingGeneBp = 100000)
originalRegions<-pradaO$applicationCoverageRegionsFiltered[,c("chr","bp1","id","bp2","chr_name")]
dChromosome<-pradaO$pgxrexApplicationDAO$selectChromosomeInformation()
originalRegions[dChromosome,on=c(chr='number'),c('chr_size'):=list(i.sizebp)]
newRegions<-as.data.frame(matrix(NA,0,0))

for(iChromosome in 1:nrow(dChromosome)){
  #iChromosome<-1
  originalRegions.chr<-originalRegions[chr==iChromosome,]
  originalRegions.chr<-originalRegions.chr[order(bp1,bp2,id),]
  if(nrow(originalRegions.chr)<1) next
  smallestBp1<-originalRegions.chr[1,bp1]
  smallestId<-originalRegions.chr[1,id]
  largestBp2<-originalRegions.chr[nrow(originalRegions.chr),bp2]
  largestChrSize<-originalRegions.chr[nrow(originalRegions.chr),chr_size]
  largestId<-originalRegions.chr[nrow(originalRegions.chr),id]
  if(is.finite(smallestBp1)) newRegions[nrow(newRegions)+1,c("chr","bp1","bp2","id")]<-c(iChromosome,1,smallestBp1,paste0(smallestId,"_LEFT"))
  if(is.finite(largestBp2) && is.finite(largestChrSize) && largestChrSize>largestBp2) newRegions[nrow(newRegions)+1,c("chr","bp1","bp2","id")]<-c(iChromosome,largestBp2+1,largestChrSize,paste0(largestId,"_RIGHT"))

  newBp1<-NA_integer_
  if(nrow(originalRegions.chr)>0){
    for(iRegion in 1:nrow(originalRegions.chr)){
      #iRegion<-1
      #iRegion<-2
      newBp2<-originalRegions.chr[iRegion,bp1]-1
      newId<-originalRegions.chr[iRegion,id]

      if(!is.na(newBp1) && newBp1 < newBp2) newRegions[nrow(newRegions)+1,c("chr","bp1","bp2","id")]<-c(iChromosome,newBp1,newBp2,paste0(newId,"_LEFT"))

      newBp1<-originalRegions.chr[iRegion,bp2]

    }
  }
}

setDT(newRegions)
newRegions[,chr:=as.integer(chr)]
newRegions[,bp1:=as.integer(bp1)]
newRegions[,bp2:=as.integer(bp2)]
newRegions<-newRegions[order(chr,bp1,bp2,id),]
newRegions[dChromosome,on=c(chr='number'),c('chr_name'):=list(i.name)]

bedDf<-as.data.frame(newRegions)
bedDf$id<-gsub("[^a-zA-z1-9_-]","_",bedDf$id) #create safe IDs
setDT(bedDf)
data.table::setorderv(bedDf,
                      cols = c("chr","bp1","id","bp2"),
                      order =c(1,1,1,1)
)
bedDf<-bedDf[,.(chr_name,bp1,bp2,id)][,bp1:=bp1-1] #includes fix of bp1 indexing for BED

fwrite(bedDf,file = file.path(settingProjectFolder,"data","bed","pgx.grch38.padded_100k.ss.inverted.bed") ,append = F,sep = "\t",encoding = "UTF-8",col.names = F)



#defaults for testing
# writeToThisBedPath=NULL
# writePaddedStrands=FALSE
# paddingGeneBp=10000
# paddingVariantCnvBp=10000
# paddingVariantSnpBp=10000
# paddingGroupFinal=10000
# nPrioritisedGene=300
# nPrioritisedCnv=100
# nPrioritisedSnp=200000
# nPrioritisedTotal=25000
# wGene=1e20
# wVariantCnv=1e7
# wVariantSnp=1
# useLabelIdsForGenes=TRUE
# verbose=TRUE

#test
# writeToThisBedPath = NULL
# nPrioritisedCnv=0
# nPrioritisedSnp=0
# nPrioritisedTotal = 1000
# writePaddedStrands = F
# paddingGeneBp = 100000
# pgxrexApplicationDAO<-pgxrexO$pgxrexApplicationDAO

#View(pradaO$applicationCoverageRegionsFiltered[pradaO$applicationCoverageRegionsFiltered$label=="CYP2D6",])
# #check
# bedDf1<-pradaO$applicationCoverageRegionsFiltered
# data.table::setorderv(bedDf1,
#                       cols = c("chr","bp1","id","bp2"),
#                       order =c(1,1,1,1)
# )
# View(bedDf1)
#
# bedDf<-pradaO$applicationCoverageRegionsFilteredPaddedStrands
# data.table::setorderv(bedDf,
#                       cols = c("chr","bp1","id","strand","bp2"),
#                       order =c(1,1,1,1,1)
# )
# View(bedDf)

#PGX:                                   The coverage of the current selection is  20296834 bp or  0.006572235
#PGX + CNV:                             The coverage of the current selection is  66340186 bp or  0.02148134
#nPrioritisedTotal = 5000 (no CNV's) :  The coverage of the current selection is  38946202 bp or  0.01261101
#nPrioritisedTotal = 5000 :             The coverage of the current selection is  79617339 bp or  0.02578056
#nPrioritisedTotal = 25000 :            The coverage of the current selection is  139605404 bp or  0.04520505
#nPrioritisedTotal = 50000 :            The coverage of the current selection is  224963114 bp or  0.07284438
#nPrioritisedTotal = 75000 :            The coverage of the current selection is  309284802 bp or  0.1001482
#nPrioritisedTotal = 100000 :           The coverage of the current selection is  394041240 bp or  0.1275929
