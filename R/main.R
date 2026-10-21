
#install package with
#devtools::install_github("tnggroup/prada")
#devtools::install_github("tnggroup/prada",ref = 'jz_dev')

#reference class, also called R5
# use browser() and source to debug in RStudio
PgxrexClass <- setRefClass("Pgxrex",
                                           fields = list(
                                             contextDatabaseList = "ANY",
                                             nThread = "numeric",
                                             folderpathLaunch="ANY",
                                             folderpathWork = "ANY",
                                             pgxrexApplicationDAO = "ANY",
                                             redcapDAO = "ANY",


                                             paddingPRSAnchorBp = "numeric",
                                             paddingGeneBp = "numeric",
                                             applicationCoverageRegions = 'ANY',
                                             applicationCoverageRegionsFiltered = 'ANY',
                                             applicationCoverageRegionsFilteredPaddedStrands = 'ANY',

                                             #analysis settings
                                             analysisConfigurationFilepath = 'ANY',
                                             analysisId = 'ANY',
                                             analysisFolderpathWork = 'ANY',
                                             analysisSettingsList = 'ANY',
                                             sampleSettingsList = 'ANY',
                                             analysisMeta = 'ANY',
                                             sampleMeta = 'ANY',


                                             testFlag.offline = 'ANY',
                                             testFlag.prepared.analysis = 'ANY'
                                           ),
                                           methods = list
                                           (
                                             #this is the constructor as per convention
                                             initialize=function(
                                                                  settingsFilePath=NULL,
                                                                  protectedFilePath=NULL
                                                                 )
                                             {

                                               #defaults
                                               contextDatabaseList<<-list()

                                               if(!is.null(settingsFilePath)){
                                                 contextDatabaseList["settingsFilePath"]<<-settingsFilePath
                                               } else {
                                                 contextDatabaseList["settingsFilePath"]<<-"settings.conf.txt"
                                               }

                                               if(!is.null(protectedFilePath)){
                                                 contextDatabaseList["protectedFilePath"]<<-protectedFilePath
                                               } else {
                                                 contextDatabaseList["protectedFilePath"]<<-"protected.conf.txt"
                                               }

                                               nThread <<- 6
                                               folderpathLaunch <<- getwd()
                                               folderpathWork <<- file.path("")
                                               pgxrexApplicationDAO <<- NULL
                                               redcapDAO <<- NULL


                                               paddingPRSAnchorBp<<-10000
                                               paddingGeneBp<<-10000

                                               applicationCoverageRegions <<- NULL
                                               applicationCoverageRegionsFiltered <<- NULL
                                               applicationCoverageRegionsFilteredPaddedStrands <<- NULL

                                               analysisConfigurationFilepath<<-NULL
                                               analysisId<<-NULL
                                               analysisFolderpathWork<<-NULL
                                               analysisSettingsList<<-list()
                                               sampleSettingsList<<-list()

                                               analysisMeta<<-as.data.frame(matrix(data = NA,nrow = 0,ncol = 0))
                                               sampleMeta<<-as.data.frame(matrix(data = NA,nrow = 0,ncol = 0))


                                               testFlag.offline<<-FALSE
                                               testFlag.prepared.analysis<<-FALSE
                                             }
                                           )
)

# we can add more methods after creating the ref class (but not more fields!)

#this is standardised and hard-coded - replace the dao with another for a custom connection
#connectPgxrexDatabase("tng_prada_system")
PgxrexClass$methods(
  connectPgxrexDatabase=function(hostToUse=NULL,usernameToUse=NULL,passwordToUse=NULL,dbnameToUse=NULL,portToUse=NULL, askForPassword=TRUE){
    if(is.null(passwordToUse) && askForPassword) passwordToUse <- rstudioapi::askForPassword(prompt = c("Enter database password for user: ",usernameToUse))
    if(is.null(hostToUse)) hostToUse <- pgxrexCentralDBDefaultHost
    if(is.null(dbnameToUse)) dbnameToUse <- pgxrexCentralDBDefaultDbName
    if(is.null(usernameToUse)) usernameToUse <- pgxrexCentralDBDefaultUsername
    if(is.null(portToUse)) portToUse <- pgxrexCentralDBDefaultPort
    pgxrexApplicationDAO <<- pgxrex::PgxrexPgDatabaseUtilityClass(host=hostToUse, dbname=dbnameToUse, user=usernameToUse, port=portToUse, password= passwordToUse)
  }
)



