
#base_url<-"https://redcap.medsci.ox.ac.uk/api/"
PgxrexRedcapUtilityClass <- setRefClass("PgxrexRedcapUtility",
                                            fields = list(
                                              base_url = "character",
                                              apiToken = "character"

                                            ),
                                            methods = list
                                            (
                                              #this is the constructor as per convention
                                              initialize=function(base_url, apiToken, askForapiToken=F)
                                              {
                                                base_url <<- base_url
                                                apiToken <<- apiToken
                                                if(askForapiToken) apiToken <<- rstudioapi::askForPassword(prompt = "Enter API token for specified user.")





                                              }
                                            )
)


# Export Records
PgxrexRedcapUtilityClass$methods(
  exportRecords=function(vRecords=NULL, vForms=NULL, returnParsed=TRUE, returnWide=TRUE, filterRecord=NULL){

    #browser()
    recordsArg<-ifelse(is.null(vRecords),
                       '',
                       paste(vRecords,sep = ",",collapse = ",")
    )
    formsArg<- ifelse(is.null(vForms),
                      '',
                      paste(vForms,sep = ",",collapse = ",")
    )

    req <- request(base_url)

    req <- req |>
      req_headers(Accept = "application/json") |>
      req_body_form(
        token=apiToken,
        content='record',
        format='json',
        type='eav', #'flat'
        records = recordsArg,
        forms = formsArg
      )

    resp <- req |>
      req_perform()

    if(returnParsed){
      lRecs <- resp_body_json(resp)
      dtData <- rbindlist(lRecs)

      if(!is.null(filterRecord)){
        dtData <- dtData[record %in% filterRecord]
      }

      if(returnWide){

        if(any(colnames(dtData)=="redcap_repeat_instrument") && any(colnames(dtData)=="redcap_repeat_instance")){
          dtData <- dcast(dtData,  record + redcap_event_name + redcap_repeat_instrument + redcap_repeat_instance ~ field_name, value.var = "value")
        } else {
          dtData <- dcast(dtData,  record + redcap_event_name ~ field_name, value.var = "value")
        }
      }
      return(dtData)
    } else {
      return(resp)
    }

    #https://rguides.dev/guides/httr2-and-apis/#json-responses

  }
)

#test of export
#testFormExport<-pgxrexObj$redcapDAO$exportRecords(vForms = c("concomitant_medication_log","diet_lifestyle_exposures"))


# Import Records
PgxrexRedcapUtilityClass$methods(
  importRecords=function(dataToEnter.longVector, returnParsed=TRUE){

    dataToEnter.json<-jsonlite::toJSON(dataToEnter.longVector,auto_unbox = T)

    req <- request(base_url)

    req <- req |>
      req_headers(Accept = "application/json") |>
      req_body_form(
        token=apiToken,
        content='record',
        action='import',
        format='json',
        type='eav', #'flat'
        overwriteBehavior = 'normal',
        forceAutoNumber = 'false',
        data = dataToEnter.json,
        returnContent = 'ids',
        returnFormat = 'json'
      )

    resp <- req |>
      req_error(is_error = \(resp) FALSE) |>
      req_perform()

    if(returnParsed){
      lRecs <- resp_body_json(resp)
      lIds <- unlist(lRecs)
      return(lIds)
    } else {
      return(resp)
    }

  }
)


PgxrexRedcapUtilityClass$methods(
  convertToLongRecordVectorForREDCapImport=function(dfDataToEnter, recordColumnName="record", redcap_event_name=NULL){

    dfDataToEnter<-as.data.frame(dfDataToEnter)

    # dfDataToEnter<-as.data.frame(matrix(NA,0,0))
    # dfDataToEnter[1,c("record","meta_id_analysis","analysis_label","meta_id","participant_id")]<-c(1,1,"test1",1,"sample1_test1")
    # dfDataToEnter[2,c("record","meta_id_analysis","analysis_label","meta_id","participant_id")]<-c(2,1,"test1",1,"sample2_test1")
    #recordColumnName<-"record"
    #redcap_event_name<-"genetics_pgx_arm_1"

    dataToEnter<-c()
    if(!is.null(redcap_event_name)){
      for(iRecord in 1:nrow(dfDataToEnter)){
        for(iVariable in 1:ncol(dfDataToEnter)){
          if(colnames(dfDataToEnter)[iVariable]==recordColumnName) next
          dataToEnter[[length(dataToEnter)+1]]<-list(
            record=dfDataToEnter[iRecord,recordColumnName],
            field_name=colnames(dfDataToEnter)[iVariable],
            value=dfDataToEnter[iRecord,iVariable],
            redcap_event_name=redcap_event_name
          )
        }
      }
    } else {
      for(iRecord in 1:nrow(dfDataToEnter)){
        for(iVariable in 1:ncol(dfDataToEnter)){
          if(colnames(dfDataToEnter)[iVariable]==recordColumnName) next
          dataToEnter[[length(dataToEnter)+1]]<-list(
            record=dfDataToEnter[iRecord,recordColumnName],
            field_name=colnames(dfDataToEnter)[iVariable],
            value=dfDataToEnter[iRecord,iVariable]
            #redcap_event_name=redcap_event_name
          )
        }
      }
    }

    return(dataToEnter)

  }
)


#
#
# #test code for reference
# library(jsonlite)
# library(httr2)
#
# dataToEnter<-c()
# dataToEnter[[1]]<-list(
#   record=1,
#   field_name="meta_id_analysis",
#   value="1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[2]]<-list(
#   record=1,
#   field_name="analysis_label",
#   value="test1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[3]]<-list(
#   record=1,
#   field_name="meta_id",
#   value="1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[4]]<-list(
#   record=1,
#   field_name="participant_id",
#   value="sample1_test1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[5]]<-list(
#   record=1,
#   field_name="meta_id_gene_cyp2b6",
#   value="1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[6]]<-list(
#   record=2,
#   field_name="meta_id_analysis",
#   value="1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[7]]<-list(
#   record=2,
#   field_name="analysis_label",
#   value="test1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[8]]<-list(
#   record=2,
#   field_name="meta_id",
#   value="2",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[9]]<-list(
#   record=2,
#   field_name="participant_id",
#   value="sample2_test1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[10]]<-list(
#   record=2,
#   field_name="meta_id_gene_cyp2b6",
#   value="2",
#   redcap_event_name="genetics_pgx_arm_1"
# )
#
# pgxrexObj$redcapDAO$importRecords(dataToEnter.longVector = dataToEnter)
#
#
# # dataToEnter.json<-jsonlite::toJSON(dataToEnter,auto_unbox = T)
# #
# # req <- request(pgxrexObj$redcapDAO$base_url)
# #
# # req <- req |>
# #   req_headers(Accept = "application/json") |>
# #   req_body_form(
# #     token=pgxrexObj$redcapDAO$apiToken,
# #     content='record',
# #     action='import',
# #     format='json',
# #     type='eav', #'flat'
# #     overwriteBehavior = 'normal',
# #     forceAutoNumber = 'false',
# #     data = dataToEnter.json,
# #     returnContent = 'ids',
# #     returnFormat = 'json'
# #   )
# #
# # resp <- req |>
# #   req_error(is_error = \(resp) FALSE) |>
# #   req_perform()
# #
# # lRecs <- resp_body_json(resp)
# # lIds <- unlist(lRecs)
