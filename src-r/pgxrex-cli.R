#PGX Rex CLI

library(pgxrex)
#library(data.table)

pgxrexObj<-PgxrexClass()
pgxrexObj$parseCommandLine()
pgxrexObj$parseConfigurationFile()
pgxrexObj$automation()





