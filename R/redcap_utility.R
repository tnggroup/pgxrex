
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
  exportRecords=function(vRecords=NULL, vForms=NULL, returnParsed=TRUE){

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
      return(dtData)
    } else {
      return(resp)
    }

    #https://rguides.dev/guides/httr2-and-apis/#json-responses

  }
)


# Import Records
PgxrexRedcapUtilityClass$methods(
  importRecords=function(dataToEnter.json.string, returnParsed=TRUE){

    dataToEnter.json<-jsonlite::toJSON(dataToEnter.json.string,auto_unbox = T)

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


#test code for reference
# library(jsonlite)
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
#
# dataToEnter[[5]]<-list(
#   record=2,
#   field_name="meta_id_analysis",
#   value="1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[6]]<-list(
#   record=2,
#   field_name="analysis_label",
#   value="test1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[7]]<-list(
#   record=2,
#   field_name="meta_id",
#   value="2",
#   redcap_event_name="genetics_pgx_arm_1"
# )
# dataToEnter[[8]]<-list(
#   record=2,
#   field_name="participant_id",
#   value="sample2_test1",
#   redcap_event_name="genetics_pgx_arm_1"
# )
#
#
# dataToEnter.json<-jsonlite::toJSON(dataToEnter,auto_unbox = T)
#
# req <- request(pgxrexObj$redcapDAO$base_url)
#
# req <- req |>
#   req_headers(Accept = "application/json") |>
#   req_body_form(
#     token=pgxrexObj$redcapDAO$apiToken,
#     content='record',
#     action='import',
#     format='json',
#     type='eav', #'flat'
#     overwriteBehavior = 'normal',
#     forceAutoNumber = 'false',
#     data = dataToEnter.json,
#     returnContent = 'ids',
#     returnFormat = 'json'
#   )
#
# resp <- req |>
#   req_error(is_error = \(resp) FALSE) |>
#   req_perform()
#
# lRecs <- resp_body_json(resp)
# lIds <- unlist(lRecs)
