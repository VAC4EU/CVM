###################################################################
# ASSIGN PARAMETERS DESCRIBING THE DATA MODEL OF THE INPUT FILES
###################################################################

# assign -ConcePTION_CDM_tables-: it is a 2-level list describing the ConcePTION CDM tables, and will enter the function as the first parameter. the first level is the data domain (in the example: 'Diagnosis' and 'Medicines') and the second level is the list of tables that has a column pertaining to that data domain 
datasources_prescriptions <- c('CPRD',"PHARMO")
thisdatasource_has_prescriptions <- ifelse(thisdatasource %in% datasources_prescriptions,TRUE,FALSE)


ConcePTION_CDM_tables <- vector(mode="list")

files<-sub('\\.csv$', '', list.files(dirinput))
for (i in 1:length(files)) {
  if (str_detect(files[i],"^VACCINES"))  ConcePTION_CDM_tables[["VaccineATC"]][[(length(ConcePTION_CDM_tables[["VaccineATC"]]) + 1)]]<-files[i]
}
files<-sub('\\.csv$', '', list.files(dirinput))
for (i in 1:length(files)) {
  if (str_detect(files[i],"^EVENTS")) {
    ConcePTION_CDM_tables[["Diagnosis"]][[(length(ConcePTION_CDM_tables[["Diagnosis"]]) + 1)]]<-files[i]
    #TODO check
    ConcePTION_CDM_tables[["Diagnosis_free_text"]][[(length(ConcePTION_CDM_tables[["Diagnosis_free_text"]]) + 1)]]<-files[i]}
  else{if (str_detect(files[i],"^MEDICINES")) ConcePTION_CDM_tables[["Medicines"]][[(length(ConcePTION_CDM_tables[["Medicines"]]) + 1)]]<-files[i] }
}

alldomain<-names(ConcePTION_CDM_tables)

ConcePTION_CDM_EAV_tables <- vector(mode="list")
EAV_table<-c()
for (i in 1:length(files)) {
  if (str_detect(files[i],"^SURVEY_OB")) {
    ConcePTION_CDM_EAV_tables[["Diagnosis"]][[(length(ConcePTION_CDM_EAV_tables[["Diagnosis"]]) + 1)]]<-list(list(files[i], "so_source_table", "so_source_column"))
    ConcePTION_CDM_EAV_tables[["Diagnosis_free_text"]][[(length(ConcePTION_CDM_EAV_tables[["Diagnosis_free_text"]]) + 1)]]<-list(list(files[i], "so_source_table", "so_source_column"))
    EAV_table<-append(EAV_table,files[i])
  }
  else{
    if (str_detect(files[i],"^MEDICAL_OB")){
      ConcePTION_CDM_EAV_tables[["Diagnosis"]][[(length(ConcePTION_CDM_EAV_tables[["Diagnosis"]]) + 1)]]<-list(list(files[i], "mo_source_table", "mo_source_column"))
      ConcePTION_CDM_EAV_tables[["Diagnosis_free_text"]][[(length(ConcePTION_CDM_EAV_tables[["Diagnosis_free_text"]]) + 1)]]<-list(list(files[i], "mo_source_table", "mo_source_column"))
      EAV_table<-append(EAV_table,files[i])
    }
  }
}

if(length(ConcePTION_CDM_EAV_tables)!=0) {
  for (t in  names(ConcePTION_CDM_EAV_tables)) {
    ConcePTION_CDM_EAV_tables_retrieve = ConcePTION_CDM_EAV_tables [[t]]
  }
}else { ConcePTION_CDM_EAV_tables_retrieve <- vector(mode="list")
}

ConcePTION_CDM_codvar <- vector(mode="list")
person_id <- vector(mode="list")
date<- vector(mode="list")
ConcePTION_CDM_datevar<-vector(mode="list")

for (ds in ConcePTION_CDM_tables[["VaccineATC"]]) {
  ConcePTION_CDM_codvar[["VaccineATC"]][[ds]] <- "vx_atc"
  person_id[["VaccineATC"]][[ds]] <- "person_id"
  date[["VaccineATC"]][[ds]] <- "vx_admin_date"
  ConcePTION_CDM_datevar[["VaccineATC"]][[ds]] <- "vx_admin_date"
}

if (length(ConcePTION_CDM_EAV_tables)!=0 ){
  for (dom in alldomain) {
    for (i in 1:(length(ConcePTION_CDM_EAV_tables[["Diagnosis"]]))){
      for (ds in append(ConcePTION_CDM_tables[[dom]],ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]])) {
        if (ds==ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]]) {
          if (str_detect(ds,"^SURVEY_OB"))  ConcePTION_CDM_codvar[["Diagnosis"]][[ds]]="so_source_value"
          if (str_detect(ds,"^MEDICAL_OB"))  ConcePTION_CDM_codvar[["Diagnosis"]][[ds]]="mo_source_value"
        }else{
          if (dom=="Medicines") ConcePTION_CDM_codvar[[dom]][[ds]]="medicinal_product_atc_code"
          if (dom=="Diagnosis") ConcePTION_CDM_codvar[[dom]][[ds]]="event_code"
          if (dom=="Diagnosis_free_text") ConcePTION_CDM_codvar[[dom]][[ds]]="event_free_text"
        }
      }
    }
  }
}else{
  for (dom in alldomain) {
    for (ds in ConcePTION_CDM_tables[[dom]]) {
      if (dom=="Medicines") ConcePTION_CDM_codvar[[dom]][[ds]]="medicinal_product_atc_code"
      if (dom=="Diagnosis") ConcePTION_CDM_codvar[[dom]][[ds]]="event_code"
      # TODO check
      if (dom=="Diagnosis_free_text") ConcePTION_CDM_codvar[[dom]][[ds]]="event_free_text"
    }
  }
}

# ConcePTION_CDM_coding_system_cols <-vector(mode="list")
# #coding system
# for (dom in alldomain) {
#   for (ds in ConcePTION_CDM_tables[[dom]]) {
#     if (dom=="Diagnosis") ConcePTION_CDM_coding_system_cols[[dom]][[ds]] = "event_record_vocabulary"
#     #    if (dom=="Medicines") ConcePTION_CDM_coding_system_cols[[dom]][[ds]] = "code_indication_vocabulary"
#   }
#}

#coding system
ConcePTION_CDM_coding_system_cols <-vector(mode="list")
if (length(ConcePTION_CDM_EAV_tables)!=0){
  for (dom in alldomain) {
    for (i in 1:(length(ConcePTION_CDM_EAV_tables[["Diagnosis"]]))){
      for (ds in append(ConcePTION_CDM_tables[[dom]],ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]])) {
        if (ds==ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]]) {
          if (str_detect(ds,"^SURVEY_OB"))  ConcePTION_CDM_coding_system_cols[["Diagnosis"]][[ds]]="so_unit"
          if (str_detect(ds,"^MEDICAL_OB"))  ConcePTION_CDM_coding_system_cols[["Diagnosis"]][[ds]]="mo_record_vocabulary"
        }else{
          # if (dom=="Medicines") ConcePTION_CDM_coding_system_cols[[dom]][[ds]]="product_ATCcode"
          if (dom=="Diagnosis") ConcePTION_CDM_coding_system_cols[[dom]][[ds]]="event_record_vocabulary"
          if (dom=="Diagnosis_free_text") ConcePTION_CDM_coding_system_cols[[dom]][[ds]]="event_record_vocabulary"
          if (dom=="Procedures") ConcePTION_CDM_coding_system_cols[[dom]][[ds]]="procedure_code_vocabulary"
        }
      }
    }
  }
}else{
  for (dom in alldomain) {
    for (ds in ConcePTION_CDM_tables[[dom]]) {
      if (dom=="Diagnosis") ConcePTION_CDM_coding_system_cols[[dom]][[ds]] = "event_record_vocabulary"
      if (dom=="Procedures") ConcePTION_CDM_coding_system_cols[[dom]][[ds]] = "procedure_code_vocabulary"
      #    if (dom=="Medicines") ConcePTION_CDM_coding_system_cols[[dom]][[ds]] = "code_indication_vocabulary"
      # 
      if (dom=="Diagnosis_free_text") ConcePTION_CDM_coding_system_cols[[dom]][[ds]] = "event_record_vocabulary"
    }
  }
}


#coding system
for (dom in alldomain) {
  for (ds in ConcePTION_CDM_tables[[dom]]) {
    if (dom=="Diagnosis") ConcePTION_CDM_coding_system_cols[[dom]][[ds]] = "event_record_vocabulary"
    if (dom=="Diagnosis_free_text") ConcePTION_CDM_coding_system_cols[[dom]][[ds]] = "event_record_vocabulary"
    #    if (dom=="Medicines") ConcePTION_CDM_coding_system_cols[[dom]][[ds]] = "code_indication_vocabulary"
  }
}

# assign 2 more 2-level lists: -id- -date-. They encode from the data model the name of the column(s) of each data table that contain, respectively, the personal identifier and the date. Those 2 lists are to be inputted in the rename_col option of the function. 
#NB: GENERAL  contains the names columns will have in the final datasets

if (length(ConcePTION_CDM_EAV_tables)!=0 ){
  for (dom in alldomain) {
    for (i in 1:(length(ConcePTION_CDM_EAV_tables[[dom]]))){
      for (ds in append(ConcePTION_CDM_tables[[dom]],ConcePTION_CDM_EAV_tables[[dom]][[i]][[1]][[1]])) {
        person_id [[dom]][[ds]] = "person_id"
      }
    }
  }
}else{
  for (dom in alldomain) {
    for (ds in ConcePTION_CDM_tables[[dom]]) {
      person_id [[dom]][[ds]] = "person_id"
    }
  }
}


if (length(ConcePTION_CDM_EAV_tables)!=0 ){
  for (dom in alldomain) {
    for (i in 1:(length(ConcePTION_CDM_EAV_tables[["Diagnosis"]]))){
      for (ds in append(ConcePTION_CDM_tables[[dom]],ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]])) {
        if (ds==ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]]) {
          if (str_detect(ds,"^SURVEY_OB")){
            date[["Diagnosis"]][[ds]]="so_date"
            date[["Diagnosis_free_text"]][[ds]]="so_date"
          }
          if (str_detect(ds,"^MEDICAL_OB")) {
            date[["Diagnosis"]][[ds]]="mo_date"
            date[["Diagnosis_free_text"]][[ds]]="mo_date"
          }
        }else{
          if (dom=="Medicines") { 
            if (thisdatasource %in% datasources_prescriptions ){
              date[[dom]][[ds]]="date_prescription"
            }else{
              date[[dom]][[ds]]="date_dispensing"
            }
          }
          if (dom=="Diagnosis" | dom=="Diagnosis_free_text") date[[dom]][[ds]]="start_date_record"
        }
      }
    }
  }
}else{
  for (dom in alldomain) {
    for (ds in ConcePTION_CDM_tables[[dom]]) { 
      if (dom=="Medicines") { 
        if (thisdatasource_has_prescriptions == TRUE){
          date[[dom]][[ds]]="date_prescription"
        }else{
          date[[dom]][[ds]]="date_dispensing"
        }
      }
      if (dom=="Diagnosis") date[[dom]][[ds]]="start_date_record"
      if (dom=="Diagnosis_free_text") date[[dom]][[ds]]="start_date_record"
    }
  }
}


# #DA CMD_SOURCE
# ConcePTION_CDM_EAV_attributes<-vector(mode="list")
# datasources<-c("ARS","TEST")
# 
# if (length(ConcePTION_CDM_EAV_tables)!=0 ){
#   for (dom in alldomain) {
#     for (i in 1:(length(ConcePTION_CDM_EAV_tables[[dom]]))){
#       for (ds in ConcePTION_CDM_EAV_tables[[dom]][[i]][[1]][[1]]) {
#         for (dat in datasources) {
#           if (dom=="Diagnosis") ConcePTION_CDM_EAV_attributes[[dom]][[ds]][[dat]][["ICD9"]] <-  list(list("RMR","CAUSAMORTE"))
#           ConcePTION_CDM_EAV_attributes[[dom]][[ds]][[dat]][["ICD10"]] <-  list(list("RMR","CAUSAMORTE_ICDX"))
#           ConcePTION_CDM_EAV_attributes[[dom]][[ds]][[dat]][["SNOMED"]] <-  list(list("AP","COD_MORF_1"),list("AP","COD_MORF_2"),list("AP","COD_MORF_3"),list("AP","COD_TOPOG"))
#           #        if (dom=="Medicines") ConcePTION_CDM_EAV_attributes[[dom]][[ds]][[dat]][["ICD9"]] <-  list(list("CAP1","SETTAMEN_ARSNEW"),list("CAP1","GEST_ECO"),list("AP","COD_MORF_1"),list("AP","COD_MORF_2"),list("AP","COD_MORF_3"),list("AP","COD_TOPOG"))
#         }
#       }
#     }
#   }
# }


#NEW ATTRIBUTES DEFINITION
files_par<-sub('\\.RData$', '', list.files(dirpargen))

if(length(files_par)>0){
  for (i in 1:length(files_par)) {
    if (str_detect(files_par[i],"^ConcePTION_CDM_EAV_attributes")) { 
      load(paste0(dirpargen,files_par[i],".RData")) 
      load(paste0(dirpargen,"ConcePTION_CDM_coding_system_list.RData")) 
      print("upload existing EAV_attributes")
    } else {
      print("create EAV_attributes")
      
      ConcePTION_CDM_coding_system_list<-vector(mode="list")
      METADATA<-fread(paste0(dirinput,"METADATA.csv"))
      ConcePTION_CDM_coding_system_list<-unique(unlist(str_split(unique(METADATA[type_of_metadata=="list_of_values" & (columnname=="so_unit" | columnname=="mo_record_vocabulary"),values])," ")))
      
      ConcePTION_CDM_EAV_attributes<-vector(mode="list")
      
      if (length(ConcePTION_CDM_EAV_tables)!=0 ){
        for (i in 1:(length(ConcePTION_CDM_EAV_tables[["Diagnosis"]]))){
          for (ds in ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]]) {
            temp <- fread(paste0(dirinput,ds,".csv"))
            for( cod_syst in ConcePTION_CDM_coding_system_list) {
              if ("mo_source_table" %in% names(temp) ) {
                temp1<-unique(temp[mo_record_vocabulary %in% cod_syst,.(mo_source_table,mo_source_column)])
                if (nrow(temp1)!=0) ConcePTION_CDM_EAV_attributes[["Diagnosis"]][[ds]][[thisdatasource]][[cod_syst]]<-as.list(as.data.table(t(temp1)))
              } else{
                temp1<-unique(temp[so_unit %in% cod_syst,.(so_source_table,so_source_column)])
                if (nrow(temp1)!=0) ConcePTION_CDM_EAV_attributes[["Diagnosis"]][[ds]][[thisdatasource]][[cod_syst]]<-as.list(as.data.table(t(temp1)))
              }
              
            }
          }
        }
      }
      
      ConcePTION_CDM_EAV_attributes_this_datasource<-vector(mode="list")
      
      if (length(ConcePTION_CDM_EAV_attributes)!=0 ){
        for (t in  names(ConcePTION_CDM_EAV_attributes)) {
          for (f in names(ConcePTION_CDM_EAV_attributes[[t]])) {
            for (s in names(ConcePTION_CDM_EAV_attributes[[t]][[f]])) {
              if (s==thisdatasource ){
                ConcePTION_CDM_EAV_attributes_this_datasource[[t]][[f]]<-ConcePTION_CDM_EAV_attributes[[t]][[f]][[s]]
              }
            }
          }
        }
      }
      
      save(ConcePTION_CDM_EAV_attributes_this_datasource, file = paste0(dirpargen,"ConcePTION_CDM_EAV_attributes.RData"))
      save(ConcePTION_CDM_coding_system_list, file = paste0(dirpargen,"ConcePTION_CDM_coding_system_list.RData"))
      
    }
  }
} else {
  
  print("create EAV_attributes")
  
  ConcePTION_CDM_coding_system_list<-vector(mode="list")
  METADATA<-fread(paste0(dirinput,"METADATA.csv"))
  ConcePTION_CDM_coding_system_list<-unique(unlist(str_split(unique(METADATA[type_of_metadata=="list_of_values" & (columnname=="so_unit" | columnname=="mo_record_vocabulary"),values])," ")))
  
  ConcePTION_CDM_EAV_attributes<-vector(mode="list")
  
  if (length(ConcePTION_CDM_EAV_tables)!=0 ){
    for (i in 1:(length(ConcePTION_CDM_EAV_tables[["Diagnosis"]]))){
      for (ds in ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]]) {
        temp <- fread(paste0(dirinput,ds,".csv"))
        for( cod_syst in ConcePTION_CDM_coding_system_list) {
          if ("mo_source_table" %in% names(temp) ) {
            temp1<-unique(temp[mo_record_vocabulary %in% cod_syst,.(mo_source_table,mo_source_column)])
            if (nrow(temp1)!=0) ConcePTION_CDM_EAV_attributes[["Diagnosis"]][[ds]][[thisdatasource]][[cod_syst]]<-as.list(as.data.table(t(temp1)))
          } else{
            temp1<-unique(temp[so_unit %in% cod_syst,.(so_source_table,so_source_column)])
            if (nrow(temp1)!=0) ConcePTION_CDM_EAV_attributes[["Diagnosis"]][[ds]][[thisdatasource]][[cod_syst]]<-as.list(as.data.table(t(temp1)))
          }
          
        }
      }
    }
  }
  
  ConcePTION_CDM_EAV_attributes_this_datasource<-vector(mode="list")
  
  if (length(ConcePTION_CDM_EAV_attributes)!=0 ){
    for (t in  names(ConcePTION_CDM_EAV_attributes)) {
      for (f in names(ConcePTION_CDM_EAV_attributes[[t]])) {
        for (s in names(ConcePTION_CDM_EAV_attributes[[t]][[f]])) {
          if (s==thisdatasource ){
            ConcePTION_CDM_EAV_attributes_this_datasource[[t]][[f]]<-ConcePTION_CDM_EAV_attributes[[t]][[f]][[s]]
          }
        }
      }
    }
  }
  
  save(ConcePTION_CDM_EAV_attributes_this_datasource, file = paste0(dirpargen,"ConcePTION_CDM_EAV_attributes.RData"))
  save(ConcePTION_CDM_coding_system_list, file = paste0(dirpargen,"ConcePTION_CDM_coding_system_list.RData"))
  
}

ConcePTION_CDM_EAV_attributes_this_datasource<-vector(mode="list")

if (length(ConcePTION_CDM_EAV_attributes)!=0 ){
  for (t in  names(ConcePTION_CDM_EAV_attributes)) {
    for (f in names(ConcePTION_CDM_EAV_attributes[[t]])) {
      for (s in names(ConcePTION_CDM_EAV_attributes[[t]][[f]])) {
        if (s==thisdatasource ){
          ConcePTION_CDM_EAV_attributes_this_datasource[[t]][[f]]<-ConcePTION_CDM_EAV_attributes[[t]][[f]][[s]]
        }
      }
    }
  }
}

if (length(ConcePTION_CDM_EAV_tables)!=0 ){
  for (dom in alldomain) {
    for (i in 1:(length(ConcePTION_CDM_EAV_tables[["Diagnosis"]]))){
      for (ds in append(ConcePTION_CDM_tables[[dom]],ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]])) {
        if (ds==ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]]) {
          if (str_detect(ds,"^SURVEY_OB")) ConcePTION_CDM_datevar[["Diagnosis"]][[ds]]="so_date"
          if (str_detect(ds,"^MEDICAL_OB"))  ConcePTION_CDM_datevar[["Diagnosis"]][[ds]]="mo_date"
        }else{
          if (dom=="Medicines") ConcePTION_CDM_datevar[[dom]][[ds]]= list("date_dispensing","date_prescription")
          if (dom=="Diagnosis") ConcePTION_CDM_datevar[[dom]][[ds]]=list("start_date_record","end_date_record")
          if (dom=="Diagnosis_free_text") ConcePTION_CDM_datevar[[dom]][[ds]]=list("start_date_record","end_date_record")
        }
      }
    }
  }
}else{
  for (dom in alldomain) {
    for (ds in ConcePTION_CDM_tables[[dom]]) { 
      if (dom=="Medicines") ConcePTION_CDM_datevar[[dom]][[ds]]= list("date_dispensing","date_prescription")
      if (dom=="Diagnosis") ConcePTION_CDM_datevar[[dom]][[ds]]=list("start_date_record","end_date_record")
      if (dom=="Diagnosis_free_text") ConcePTION_CDM_datevar[[dom]][[ds]]=list("start_date_record","end_date_record")
    }
  }
}


ConcePTION_CDM_datevar_retrieve<-list()

if (length(ConcePTION_CDM_EAV_tables)!=0 ){
  for (i in 1:(length(ConcePTION_CDM_EAV_tables[["Diagnosis"]]))){
    for (ds in ConcePTION_CDM_EAV_tables[["Diagnosis"]][[i]][[1]][[1]]) {
      ConcePTION_CDM_datevar_retrieve = ConcePTION_CDM_datevar [["Diagnosis"]]
    }
  }
}

rm(temp,temp1)