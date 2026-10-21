#PGX Rex TEST

library(pgxrex)
#library(data.table)
#remotes::install_github("tnggroup/pgxrex",ref = 'jz_dev')

#setwd(utils::getSrcDirectory(function(){})[1]) #force wd in case of previous changes
#source("/users/k2481717/project/pgxrex/src-r/pgxrex-test.R")


PgxrexClass$trace("automation") #DEBUG

pgxrexObj<-PgxrexClass()
pgxrexObj$parseConfigurationFile()

#pgxrexObj$contextDatabaseList

pgxrexObj$testFlag.offline<-FALSE #TEST
pgxrexObj$testFlag.prepared.analysis<-FALSE #TEST

#sbrowser()

#pgxrexObj$trace("automation")
pgxrexObj$automation(command = "automatic_routine")

