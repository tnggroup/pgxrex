


PgxrexPgDatabaseUtilityClass <- setRefClass("PgxrexPgDatabaseUtility",
                                             fields = list(
                                               dbname = "character",
                                               user = "character",
                                               host = "character",
                                               port = "numeric",
                                               statusCode = "numeric",
                                               cacheDb = "ANY",
                                               savedCachePath = "character",
                                               connection = "ANY"
                                             ),
                                             methods = list
                                             (
                                               #this is the constructor as per convention
                                               initialize=function(dbname, host, user, port, password, doCache=T, askForPassword=F)
                                               {

                                                 dbname <<- dbname
                                                 host <<- host
                                                 user <<- user
                                                 port <<- port
                                                 if(askForPassword){
                                                   connection <<- dbConnect(RPostgres::Postgres(),
                                                                            dbname = dbname,
                                                                            host = host,
                                                                            port = port,
                                                                            user = user,
                                                                            password = rstudioapi::askForPassword(prompt = "Enter database password for specified user."))
                                                 } else {
                                                   connection <<- dbConnect(RPostgres::Postgres(),
                                                             dbname = dbname,
                                                             host = host,
                                                             port = port,
                                                             user = user,
                                                             password = password)
                                                 }
                                                 savedCachePath <<- "./pgDatabaseUtilityCache.Rds"

                                                 if(doCache==T) cache()
                                               }
                                             )
)

# we can add more methods after creating the ref class (but not more fields!)

# cache frequently and rarely changing data from the server
PgxrexPgDatabaseUtilityClass$methods(
  cache=function(refresh=F, write=T, cacheSaveFilePath=NULL){
    if(is.null(cacheSaveFilePath)){
      cacheSaveFilePath <- savedCachePath
    }

    if(file.exists(cacheSaveFilePath) && refresh == F){
      cacheDb <<- readRDS(file=cacheSaveFilePath)
    } else {
      cacheDb <<- c()

      #store cache data here

      if(write) {
        saveRDS(cacheDb,file = cacheSaveFilePath)
        savedCachePath <<- cacheSaveFilePath
      }

    }
  }
)

#this is not injection safe
PgxrexPgDatabaseUtilityClass$methods(
  setSearchPath=function(searchPath){
    #searchPath<-"dat_cohort"
    c <- paste0("SET search_path TO ",paste(searchPath, sep = ","))
    #cat(c)
    q <- dbSendQuery(connection, c)
  }
)

PgxrexPgDatabaseUtilityClass$methods(
  setSearchPath.standard=function(){
    setSearchPath("\"$user\", public")
  }
)

PgxrexPgDatabaseUtilityClass$methods(
  getPGTempSchema=function(cohort){
    q <- dbSendQuery(connection,
                     "SELECT nspname FROM pg_namespace WHERE oid  =  pg_my_temp_schema()")
    res<-dbFetch(q)
    dbClearResult(q)
    if(length(res)>0) return(res[[1]]) else return(NA_character_)
  }
)

# PgxrexPgDatabaseUtilityClass$methods(
#   getPgxrexChromosome=function(){
#     q <- dbSendQuery(connection,
#                      "SELECT * FROM prada.chromosome")
#     res<-dbFetch(q)
#     dbClearResult(q)
#     if(length(res)>0) return(res[[1]]) else return(NA_character_)
#   }
# )

PgxrexPgDatabaseUtilityClass$methods(
  selectHarmonisedCombinedPgxGene=function(
  ){
    q <- dbSendQuery(connection,
                     "SELECT * FROM prada.harmonised_combined_pgx_gene")
    res<-dbFetch(q)
    dbClearResult(q)
    return(res)
  }
)


PgxrexPgDatabaseUtilityClass$methods(
  selectApplicationCoverageRegions=function(
    paddingGeneBp=10000,
    paddingVariantCnvBp=10000,
    paddingVariantSnpBp=5000,
    nPrioritisedGene=300,
    nPrioritisedCnv=100,
    nPrioritisedSnp=200000,
    nPrioritisedTotal=25000,
    wGene=1e20,
    wVariantCnv=1e7,
    wVariantSnp=1
  ){

    # paddingGeneBp integer DEFAULT 10000,
    # paddingVariantCnvBp integer DEFAULT 10000,
    # paddingVariantSnpBp integer DEFAULT 5000,
    # nPrioritisedGene integer DEFAULT 300,
    # nPrioritisedCnv integer DEFAULT 100,
    # nPrioritisedSnp integer DEFAULT 200000,
    # nPrioritisedTotal integer DEFAULT 25000,
    # wGene double precision DEFAULT 1e20,
    # wVariantCnv double precision DEFAULT 1e7,
    # wVariantSnp double precision DEFAULT 1

    qString <- "SELECT * FROM prada.get_coverage_regions(
                        paddingGeneBp=>$1,
                        paddingVariantCnvBp=>$2,
                        paddingVariantSnpBp=>$3,
                        nPrioritisedGene=>$4,
                        nPrioritisedCnv=>$5,
                        nPrioritisedSnp=>$6,
                        nPrioritisedTotal=>$7,
                        wGene=>$8,
                        wVariantCnv=>$9,
                        wVariantSnp=>$10
                     )"
    lQArguments<-list(
      paddingGeneBp,
      paddingVariantCnvBp,
      paddingVariantSnpBp,
      nPrioritisedGene,
      nPrioritisedCnv,
      nPrioritisedSnp,
      nPrioritisedTotal,
      wGene,
      wVariantCnv,
      wVariantSnp
    )
    q <- dbSendQuery(connection,qString,lQArguments)
    res<-dbFetch(q)
    dbClearResult(q)

    q <- dbSendQuery(connection,
                     "SELECT * FROM t_coverage_region")
    res<-dbFetch(q)
    dbClearResult(q)
    return(res)
  }
)


PgxrexPgDatabaseUtilityClass$methods(
  selectFilteredApplicationCoverageRegions=function(){
    qString <- "SELECT * FROM prada.filter_coverage_regions()"
    q <- dbSendQuery(connection,qString)
    res<-dbFetch(q)
    dbClearResult(q)

    q <- dbSendQuery(connection,"SELECT * FROM t_coverage_region_filtered")
    res<-dbFetch(q)
    dbClearResult(q)
    return(res)
  }
)

PgxrexPgDatabaseUtilityClass$methods(
  selectChromosomeInformation=function(){
    qString <- "SELECT * FROM prada.chromosome"
    q <- dbSendQuery(connection,qString)
    res<-dbFetch(q)
    dbClearResult(q)
    return(res)
  }
)


#should be used with care as it passes the data by value rather than reference
PgxrexPgDatabaseUtilityClass$methods(
  importDataAsTable=function(schema_name, table_name, df, temporary = T, replace=F){

    if(is.null(schema_name) | temporary==T){
      consensus_name<-table_name
    } else {
      consensus_name<-Id(schema = schema_name, table = table_name)
    }

    if(dbExistsTable(conn = connection, name = consensus_name) & replace==T) dbRemoveTable(conn = connection, name = consensus_name, temporary= temporary)

    if(!is.null(df)){
      if(!dbExistsTable(conn = connection, name = consensus_name)) dbCreateTable(conn = connection, name = consensus_name, fields = df, temporary = temporary)

      dbAppendTable(conn = connection, name = consensus_name, value = df)
    }

  }
)
