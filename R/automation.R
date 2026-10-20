

#main function to run automated procedures
PgxrexClass$methods(
  automation=function(
    command=NULL
    ){

    if(!is.null(command)){
      contextDatabaseList["command"]<<-command
    }

    shared()

    command<-fcoalesce(contextDatabaseList$command,"about") #default to print and exit
    if(command=="about"){
      catl("PGX Rex - A software to work with ONT long read sequencing and pharmacogenomics. Please specify the command to run.")
    } else if(command=="automatic_routine"){

      #check for analysis to run

      analyses <- ifelse(
        nchar(folderpathWork)>0,
        list.files(path = folderpathWork,pattern = ".+\\.pgxrex\\.analysis\\.conf\\.txt$"),
        list.files(pattern = ".+\\.pgxrex\\.analysis\\.conf\\.txt$")
      )

      if(length(analyses)>0){

        analysisConfigurationFilepath<<-ifelse(
          nchar(folderpathWork)>0,
          file.path(folderpathWork,analyses[1]),
          analyses[1])
        analysis()

      } else {
        #no analyses
        catl("There are no analyses to run. A template analysis config file has been created.")
        defaultAnalysiOptions <- list(
          #`analysisId`="myAnalysis", #this is not put in the template as default
          `folderpathSequencing`="",
          `folderpathWFPGX`=""
          #`folderpathMosdepth`="" #if no value - assume no mosdepth has been run and run automatically
        )

        cTemplateFilePath<-ifelse(nchar(folderpathWork)>0,
                                  file.path(folderpathWork,"TEMPLATE.pgxrex.analysis.conf.txt")
                                  ,
                                  "TEMPLATE.pgxrex.analysis.conf.txt"
                                  )
        if(!file.exists(cTemplateFilePath)){
          dfOpts<-data.table(names=names(defaultAnalysiOptions), values=defaultAnalysiOptions)
          fwrite(dfOpts,cTemplateFilePath,sep = "=",row.names = FALSE,col.names = FALSE)
          catl("A template analysis config file has been created.")
        }

      }


    } else {
      catl("Misspecified command. You specified the command: ",command)
    }

  }
)

#Default command line parsing
#not automated - needs to be run by any command line application
PgxrexClass$methods(
  parseCommandLine=function(assignToContextDatabase=TRUE){

    #command line options
    optionParser <- OptionParser()

    optionParser <- add_option(
      object = optionParser,
      opt_str = c("-c", "--command"),
      type = "character",
      default = "automatic_routine",
      help = "Command to run."
    )

    optionParser <- add_option(
      object = optionParser,
      opt_str = c("-s", "--settingsFilePath"),
      type = "character",
      default = "settings.conf.txt",
      help = "Path to settings-file."
    )

    optionParser <- add_option(
      object = optionParser,
      opt_str = c("-p", "--protectedFilePath"),
      type = "character",
      default = "protected.conf.txt",
      help = "Path to protected settings-file (passwords etc.)"
    )

    parsedCLOptions <- parse_args(optionParser)

    if(assignToContextDatabase){
      contextDatabaseList[names(parsedCLOptions)]<<-parsedCLOptions
    }

    return(parsedCLOptions)

  }
)

#Default configuration file parsing
#This deals both with the standard configuration file and the protected file
#not automated - needs to be run by any command line application
PgxrexClass$methods(
  parseConfigurationFile=function(assignToContextDatabase=TRUE){

    parsedCFOptions <- list(
      `pgxrexId`="mySite",
      `pgxrexDbHost`="localhost",
      `pgxrexDbUsername`="tng_prada_system",
      `pgxrexDbName`="prada_local",
      `pgxrexDbPort`=65432,
      `redcapBaseUrl`="https://redcap.medsci.ox.ac.uk/api/",
      `ncores`=6,
      `folderpathWork`="" #default to current folder
    ) #defaults


    if(file.exists(contextDatabaseList$settingsFilePath)){
      parsedCFOptions <- pgxrex::readMetadata(filePath = contextDatabaseList$settingsFilePath)
    } else {
      if(!file.exists("settings.conf.txt")){
        dfOpts<-data.table(names=names(parsedCFOptions), values=parsedCFOptions)
        fwrite(dfOpts,"settings.conf.txt",sep = "=",row.names = FALSE,col.names = FALSE)
      }
    }

    if(assignToContextDatabase){
      contextDatabaseList[names(parsedCFOptions)]<<-parsedCFOptions
    }


    parsedPCFOptions <- list(
      `pgxrexDbPassword`="XXX_replace_with_password_XXX",
      `redcapApiToken`="XXX_replace_with_API_token_XXX"
    ) #defaults

    if(assignToContextDatabase){
      contextDatabaseList[names(parsedPCFOptions)]<<-parsedPCFOptions
    }

  }
)

#set up before running automation procedure(s)
PgxrexClass$methods(
  shared=function(){
    #contextDatabaseList<-pgxrexObj$contextDatabaseList
    catl("***PGX Rex automated routines***")
    if(any(names(contextDatabaseList)=="ncores")) {
      nThread <<- as.integer(contextDatabaseList["ncores"])
    }
    catl0("ncores=",nThread)

    if(any(names(contextDatabaseList)=="folderpathWork")) {
      if(nchar(contextDatabaseList["folderpathWork"])>0) folderpathWork <<-  normalizePath(file.path(contextDatabaseList["folderpathWork"]),mustWork = T)
    }
    catl0("folderpathWork=",folderpathWork)

    #connect database
    connectPgxrexDatabase(
      hostToUse=contextDatabaseList$pgxrexDbHost,
      usernameToUse=contextDatabaseList$pgxrexDbUsername,
      passwordToUse=contextDatabaseList$pgxrexDbPassword,
      dbnameToUse=contextDatabaseList$pgxrexDbName,
      portToUse=as.integer(contextDatabaseList$pgxrexDbPort),
      askForPassword=FALSE
    )
    catl("Database connected")

    redcapDAO<<-PgxrexRedcapUtilityClass(base_url=contextDatabaseList$redcapBaseUrl, apiToken=contextDatabaseList$redcapApiToken)
    catl("RedCap configured")


  }
)

PgxrexClass$methods(
  analysis=function(){
    #contextDatabaseList<-pgxrexObj$contextDatabaseList
    catl("Running analysis using settings from:",analysisConfigurationFilepath)

    parsedAnalysisOptions <- pgxrex::readMetadata(filePath = analysisConfigurationFilepath)


    #determine analysis id/code/name
    analysisId<<-paste0(contextDatabaseList$pgxrexId,"_unknown")
    if(any(names(parsedAnalysisOptions)=="analysisId")){
      analysisId<<-parsedAnalysisOptions[["analysisId"]]
    }
    catl("Analysis is:",analysisId)


    #general analysis data setup - requires the database
    computeGenomeCoverage( #default settings
      nPrioritisedGene=300,
      nPrioritisedCnv=0,
      nPrioritisedSnp=0,
      verbose = TRUE
    )
    catl("Genomic regions established.")

    #setup analysis working directory
    analysisFolderpathWork<<-file.path.coalesce(folderpathWork,"analysis",analysisId) #default
    dir.create(analysisFolderpathWork,recursive = TRUE)
    setwd(analysisFolderpathWork)



    #run analysis data collection - not tested!
    addAnalysisSetting(
      settingLabel = analysisId,
      folderPathAnalysisSequencingRaw = parsedAnalysisOptions["folderpathSequencing"],
      folderPathAnalysisOutputRaw = parsedAnalysisOptions["folderpathWFPGX"],
      folderPathDepthAnalysisOutputRaw = file.path(analysisFolderpathWork,"mosdepth")
      )
    collectAnalysisCallData(analysisId)
    catl("Sequencing and basecall data collected.")
    if(testFlag.offline){
      readPrintData(analysisFolderpathWork)
    } else {
      #run mosdepth for each sequenced individual in the analysis
      sampleMeta.analysis<-sampleMeta[sampleMeta$analysis==analysisId,]
      if(nrow(sampleMeta.analysis)>0){

        analysisFolderpathWork.mosdepth<-file.path.coalesce(analysisFolderpathWork,"mosdepth")
        dir.create(analysisFolderpathWork.mosdepth,recursive = TRUE)
        setwd(analysisFolderpathWork.mosdepth)

        for(iSample in 1:nrow(sampleMeta.analysis)){
          #iSample<-1
          cBarcode<-sampleMeta.analysis[iSample,c("barcode")]
          if(file.exists(paste0(cBarcode,".per-base.bed.gz"))) next
          wrapper.mosdepth(
            label = cBarcode,bamFilePath = file.path(parsedAnalysisOptions["folderpathWFPGX"],"output",cBarcode,paste0(cBarcode,".haplotagged.bam")),
            threads = nThread,
            mosdepthPath = ""
          )
        }

        setwd(analysisFolderpathWork)

      }

      collectAnalysisDepthData(analysisId)
      catl("Sequencing depth data collected.")
      computeDepthDataStatistics() #do we need the bed-file here?
      catl("Sequencing depth data statistics computed.")
      computeCallStatistics()
      catl("Sequencing and basecall data statistics computed.")

    }

    setwd(folderpathLaunch)


    #reporting
    #TODO - when our table is set up


    #phenoconversion
    browser()
    recs<-redcapDAO$exportRecords(vForms=list("concomitant_medication_log"))
    #recs<-redcapDAO$exportRecords()
    #resp_raw(recs)
    lRecs <- resp_body_json(recs)
    #length(lRecs)
    #lRecs[1]
    #lRecs[2]
    dtData <- data.table::rbindlist(lRecs)


    ##read phenoconversion data
    if(testFlag.offline){

    } else {
      #TODO
    }


    catl("block")
  }
)
