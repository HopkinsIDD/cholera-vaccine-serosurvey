

#dropbox token

token <- read_rds('data/generated_data/dropbox_token/token.rds')



#####################sample selection #################################


#' @param name_start character string, first letters in file name
#' @param path character string, path to folder of interest, end with "/"
#' @param exclude character string, patterns to exclude from the file names of interest
#'
#' @return character string, path to most recent file
#' @export
#'
#' @examples
find_recent_file <- function(name_start, path, exclude=NULL, verbose=T){
    if(substring(path, nchar(path))!="/")
        warning("Path does not end with a '/', problems may ensue.")
    ## view all files of that name at that path
    file_list <- list.files(path=path,
                            pattern=paste0(name_start, "*"))
    ## remove files with unwanted patterns
    if(!is.null(exclude)){
        for(i in 1:length(exclude))
            file_list <- file_list[!grepl(pattern = exclude[i], file_list)]
    }
    if(length(file_list)==0){
        warning('File not found')
        return(NA)
    }
    ## view file info
    file_info <- file.info(paste0(path, file_list))
    ## find most recent file
    most_recent_file <- paste0(path,
                               file_list[which.max(file_info$mtime)])

    if(verbose)
        cat(sprintf("Loaded file: \n %s last updated on \n %s \n",
                    most_recent_file,file_info$mtime[which.max(file_info$mtime)]))

    return(most_recent_file)
}

#' Find ROC curve cutoff
#' Finds the optimal cutoff for probabilities using either the Youden statistic or the point with the maximum sensitivity given the maximum specificity
#'
#' @param preds numeric vector of predictions, each between 0 and 1
#' @param labs numeric or logical vector of outcomes, each equaling 1/0 or T/F
#' @param cutoff_type character, either "youden" or "spec"
#'
#' @return numeric cutoff
#' @export
#'
#' @examples
find_roc_cutoff <- function(preds,
                        labs,
                        cutoff_type=c("youden","spec"),
                        train_spec_cutoff=0.99
                        ){
    cutoff_type <- match.arg(cutoff_type)
    roc_sens_spec <- prediction(preds, labs) %>%
        performance("sens","spec")
    if(cutoff_type=="youden"){
        roc_dat <- data.frame(spec=roc_sens_spec@x.values[[1]] %>% unlist,
                              sens=roc_sens_spec@y.values[[1]] %>% unlist,
                              cutoff=roc_sens_spec@alpha.values[[1]] %>% unlist) %>%
            mutate(youden=sens+spec-1)
        cutoff <- roc_dat$cutoff[which.max(roc_dat$youden)]
    }
    if(cutoff_type=="spec"){
        roc_dat <- data.frame(spec=roc_sens_spec@x.values[[1]] %>% unlist,
                              sens=roc_sens_spec@y.values[[1]] %>% unlist,
                              cutoff=roc_sens_spec@alpha.values[[1]] %>% unlist) %>%
                                    filter(cutoff<=1) %>%
                                # filter(spec==max(spec))
                                filter(spec>=train_spec_cutoff)
        
        cutoff <- roc_dat$cutoff[which.max(roc_dat$sens)]
    }
    return(cutoff)
}

#' Fit ranger model to cohort data and predict survey data, using cross validation to find sensitivity and specificity
#'
#' @param cohort_dat data.frame, cohort study data
#' @param survey_dat data.frame, cross-sectional serosurvey data
#' @param days_back numeric, time frame of interest for seropositivity
#' @param k_folds character or numeric, how should the cross validation be carried out? If "loo" then leave-one-out cross validation is performed. If a character string that matches a column name, then folds are determined by that covariate (e.g. "id" will create a new fold for each id). If a number is given, then that is the number of CV folds, e.g. 10-fold CV.
#' @param cutoff_type character, use the Youden statistic cutoff ("youden") or maximize the specificity ("spec")
#' @param save_results logical, save the individual results within the function?
#' @param ... ranger parameters
#'
#' @return list with ranger objects and diagnostics as well as predictions for the survey data
#' @export
#'
#' @examples
fit_ranger_cv <- function(cohort_dat,
                          days_back,
                          rf_var,
                          k_folds=c("id", "loo"),
                          cutoff_type=c("youden","spec"),
                          save_results=T,
                          ...){
    cutoff_type <- match.arg(cutoff_type)
    ## assign observations to folds
    if(k_folds=="loo"){
        fold <- seq(1:nrow(cohort_dat))
    } else
        if(k_folds %in% colnames(cohort_dat)){
            ids <- unique(cohort_dat[,k_folds,T])
            fold <- rep(NA, nrow(cohort_dat))
            for(i in 1:length(ids)){
                fold[cohort_dat[,k_folds,T]==ids[i]] <- i
            }
        } else
            if(as.numeric(k_folds)){
                fold <- sample(k_folds, nrow(cohort_dat), replace=T)
            } else
                stop("k_fold should be 'loo', the name of a column in cohort_dat, or a number.")

    rf_formula <- paste0("factor(inf_",days_back,") ~",
                         paste(rf_var, collapse="+")) %>% as.formula
    ## we run k-fold cross validation with random forests leaving out
    ## one fold, fitting on the remaining observations, then predicting
    ## the left-out fold
    k_preds <- foreach(i=1:max(fold), .combine=rbind) %dopar%{
        print(paste("Rep", i, "of", length(unique(fold))))
        set.seed(i)
        k_idx <- fold==i
        tmp_forest <- ranger(rf_formula,
                             data=cohort_dat[!k_idx,],
                             probability=T,
                             ...)
        tmp_cutoff <- find_roc_cutoff(tmp_forest$predictions[,2],
                                  cohort_dat[!k_idx,paste0("inf_",
                                                           days_back)])
        # tmp_roc <- prediction(tmp_forest$predictions[,2],
        #                       cohort_dat[!k_idx,paste0("inf_",
        #                                                days_back)]) %>%
        #     performance("sens","spec")
        #
        # if(cutoff_type=="youden"){
        #     tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
        #                              sens=tmp_roc@y.values[[1]] %>% unlist,
        #                              cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
        #         mutate(youden=sens+spec-1)
        #     tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$youden)]
        # }
        # if(cutoff_type=="youden_plus"){
        #     tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
        #                              sens=tmp_roc@y.values[[1]] %>% unlist,
        #                              cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
        #         mutate(youden=sens+1.1*spec-1)
        #     tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$youden)]
        # }
        # if(cutoff_type=="spec"){
        #     tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
        #                              sens=tmp_roc@y.values[[1]] %>% unlist,
        #                              cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
        #         filter(cutoff<=1) %>%
        #         filter(spec==max(spec))
        #     tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$sens)]
        # }
        votes <- predict(tmp_forest,
                         cohort_dat[k_idx,],
                         type="response")$predictions[,2]

        return(cohort_dat[k_idx,] %>%
                   mutate(rf_cutoff=tmp_cutoff,
                          rf_votes=votes,
                          rf_preds=votes>tmp_cutoff))
    }
    if(save_results)
        saveRDS(k_preds, file=paste0("generated_data/cohort-w-rf-preds-",days_back,"-",cutoff_type,"-",Sys.Date(),".rds"))

    ## fit and save full random forest model
    set.seed(1)
    rf_model <- ranger(rf_formula,
                       data=cohort_dat,
                       probability=T,
                       ...)

    if(save_results)
        saveRDS(rf_model, paste0("generated_data/forest_",days_back,".rds"))

    return(list(rf_formula=rf_formula,
                cohort_preds=k_preds,
                rf_model=rf_model,
                days_back=days_back,
                k_folds=fold,
                cutoff_type=cutoff_type))
}

#' Find best sample for predicting the rest of the cohort study
#'
#' @param cohort_dat data frame of cohort study
#' @param rf_var vector of character strings naming columns to be used by random forest
#' @param cutoff_type character, use the Youden statistic cutoff ("youden") or maximize the specificity ("spec")
#' @param search_type character, only "random" is available for now
#' @param n_samples numeric, number of samples to fit the random forest with
#' @param ... other arguments for the ranger package
#'
#' @return
#' @export
#'
#' @examples
find_best_sample <- function(cohort_dat,
                             rf_var,
                             cutoff_type=c("youden", "youden_plus","spec"),
                             search_type=c("random"),
                             n_samples,
                             hh_samples=NULL,
                             balance_test_set=F,
                             ...){
    cutoff_type <- match.arg(cutoff_type)
    search_type <- match.arg(search_type)

    ## create rf_formula using days_back and rf_var
    rf_formula <- paste0("factor(inf_365) ~",
                         paste(rf_var, collapse="+")) %>% as.formula

    ## list adult and child ids
    child_ids <- cohort_dat$id[cohort_dat$age<5] %>% unique
    adult_ids <- cohort_dat$id[cohort_dat$age>=5 & cohort_dat$cx<9] %>% unique
    adult_ids <- adult_ids[-which(adult_ids %in% child_ids)]
    hh_ids <- cohort_dat$id[cohort_dat$cx==9] %>% unique
    if(is.null(hh_samples)){
        adult_ids <- c(adult_ids, hh_ids)
        hh_samples <- 0
    }
    ## choose samples for fitting and predicting
    if(search_type=="random"){
        sample_dat <- foreach(i=1:n_samples, .combine=rbind) %dopar% {
            set.seed(i)
            c(sample(adult_ids, 20-hh_samples),
              sample(hh_ids, hh_samples),
              sample(child_ids, 20))
        }
    }

    ## fit random forest models and make predictions
    oos_preds <- foreach(i=1:n_samples, .combine=rbind) %dopar% {
        if(i %% 1e3==0)
            print(paste("Rep", i, "of", n_samples,";", Sys.time()))
        set.seed(i)
        sample_ids <- sample_dat[i,]
        train_dat <- cohort_dat %>% filter(id %in% sample_ids)
        test_dat <- cohort_dat %>% filter(!(id %in% sample_ids))
        tmp_rf <- ranger(rf_formula,
                         data=train_dat,
                         probability=T,
                         ...)
        tmp_cutoff <- find_roc_cutoff(tmp_rf$predictions[,2],
                                  train_dat$inf_365)
        test_preds <- predict(tmp_rf, test_dat, type="response")$predictions[,2]
        pred_dat <- test_dat %>%
            mutate(rf_votes=test_preds,
                   rf_preds=test_preds>tmp_cutoff,
                   pred_correct=rf_preds==(inf_365==1))
        if(balance_test_set)
            pred_dat <- make_balanced_test_set(pred_dat)
        return(tibble(sample=i,
                      cvAUC=ci.pooled.cvAUC(predictions=pred_dat$rf_votes,
                                            labels=pred_dat$inf_365,
                                            ids=pred_dat$id)$cvAUC,
                      pred_accuracy=mean(pred_dat$pred_correct),
                      pred_sens=mean(pred_dat$pred_correct[pred_dat$inf_365==1]),
                      pred_spec=mean(pred_dat$pred_correct[pred_dat$inf_365==0]),
                      child_accuracy=mean(pred_dat$pred_correct[pred_dat$age<5]),
                      child_sens=mean(pred_dat$pred_correct[pred_dat$age<5 & pred_dat$inf_365==1]),
                      child_spec=mean(pred_dat$pred_correct[pred_dat$age<5 & pred_dat$inf_365==0]),
                      adult_accuracy=mean(pred_dat$pred_correct[pred_dat$age>4]),
                      adult_sens=mean(pred_dat$pred_correct[pred_dat$age>4 & pred_dat$inf_365==1]),
                      adult_spec=mean(pred_dat$pred_correct[pred_dat$age>4 & pred_dat$inf_365==0])))
    }
    return(list(oos_preds=oos_preds,
                sample_dat=sample_dat))
}


make_balanced_test_set <- function(test_dat){
    test_dat %>%
        filter(inf_365==0) %>%
        sample_n(10000, T) %>%
        bind_rows(test_dat %>%
                      filter(inf_365==1, lbday2<45) %>%
                      sample_n(2500, T)) %>%
        bind_rows(test_dat %>%
                      filter(inf_365==1, lbday==90) %>%
                      sample_n(2500, T)) %>%
        bind_rows(test_dat %>%
                      filter(inf_365==1, lbday==180) %>%
                      sample_n(2500, T)) %>%
        bind_rows(test_dat %>%
                      filter(inf_365==1, lbday>=270) %>%
                      sample_n(2500, T))
}


#####################cleaning #################################


#read the excel file in
drop_read_excel <- function(file, dtoken, dest = tempdir(), ...) {
    localfile = paste0(dest, '/', basename(file))
    drop_download(file, localfile, overwrite = TRUE, dtoken = dtoken)
    readxl::read_excel(localfile, ...)
}


#extract data from luminex output (.csv) file and put into a list of dataframes
extractData<- function(file,token){
    
    #read in .CSV
    trial <- drop_read_csv(file,
                           dtoken=token,
                           header=FALSE,
                           col.names=paste("column", 1:26, sep="_"),
                           colClasses = "character"
    )
    
    # #Local
    # trial <- read_csv(file,
    #               # header=FALSE,
    #               col_names=paste("column", 1:26, sep="_"),
    #               col_types  = "c"
    #               )
    
    #find the times the "DataType shows up"
    dataTypes_index <- (1:length(trial$column_1))[trial$column_1=="DataType:"]
    dataTypes_names <- trial$column_2[dataTypes_index]
    
    #find the beginng and end of data types to form data frames
    starts <- dataTypes_index+1
    ends   <- dataTypes_index-1
    ends   <- c(ends[-1],length(trial$column_1))
    
    
    #intitialize a list
    dataList <- list()
    #store meta data
    dataList[["metadata"]] <- trial[1:(dataTypes_index[1]-1),]
    dataList[["batch"]] <- dataList$metadata[5,2]
    #get type, isotype, subclass and layout number
    dataList[["batch_split"]] <-unlist(str_split(dataList$batch,pattern = "_"))
    #get date
    dataList[["date"]] <- dataList$metadata[3,2]
    
    #for each data type make a dataframe to store
    for (i in 1:length(dataTypes_names)){
        
        #create a temporary dataframe
        tmp <- data.frame()
        
        if (starts[i]!=ends[i]){ #only include sections with data in them
            
            #find rows with data
            rows <- (starts[i]+1):ends[i]
            tmp <- trial[rows,]
            
            #add headers
            colnames(tmp) <- trial[starts[i],]
            tmp <- tmp[,colnames(tmp)!=""] #remove columns with no header
            
            #make data long, numeric and analyzeable
            if (dataTypes_names[i] %in% c("Net MFI","Median",
                                          "Std Dev","%CV",
                                          "Count",
                                          "Trimmed % CV of Microspheres"
            )){
                tmp <- tmp %>%
                    # select(-Sample,-`Total Events`) %>% #remove useless columns
                    gather(antigen,value,-c(Location,Sample,`Total Events`)) %>% #make long
                    mutate(value=as.numeric(value)) %>% #make value numeric
                    mutate(location=sapply(Location,FUN = get_loc)) %>% #find location in string
                    mutate(batch=dataList$batch) %>% #record batch name
                    mutate(date=dataList$date) %>% #record date
                    select(-Location) #%>%
                    #antigen name fixes
                    # mutate(antigen=ifelse(antigen=="V. cholerae cytolysin (VCC)","VCC",antigen))
                
                #rename value column
                tmp[dataTypes_names[i]] <- tmp$value
                tmp <- tmp %>% select(-value)
                
                if(length(dataList$batch_split)==4){
                    tmp <- tmp %>%
                        mutate(plateType=dataList$batch_split[1]) %>% #record batch info
                        mutate(isotype=dataList$batch_split[2]) %>%
                        mutate(subclass=dataList$batch_split[3]) %>%
                        mutate(plateNumber=dataList$batch_split[4])
                }
                
                
                
            }
            
        }
        
        # tmp <- tmp[,colnames(tmp)!=""] #remove columns with no header
        
        #add dataframe to the list
        dataList[[dataTypes_names[i]]] <- tmp
        
        
    }
    
    
    return(dataList)
}



#get the plate location from the luminex output string
get_loc <- function(string){
    loc <- str_split(string,",")
    loc <- unlist(loc)[2]
    loc <- str_remove(loc,"\\)")
    
    return(loc)
}
#get the participant ID from sample ID
get_id <- function(string){
    id <- str_split(string," ")
    id <- unlist(id)[1]
}

#get the day from the luminex output string
get_day <- function(string){
    day <- str_split(string," ")
    day <- unlist(day)[2]
    day <- str_remove(day,"d")
}

#bring in plate layout
getLayout <-function(path, ranges, type, token){
    
    df <- data.frame()
    
    for (r in 1:length(ranges)){
        tmp <-   case_layout1 <- drop_read_excel(path,
                                                 dtoken=token, sheet=1,
                                                 # col_names=FALSE,
                                                 col_types="text",
                                                 range=ranges[r])
        
        colnames(tmp)[1] <- "letter"
        
        tmp <- tmp %>% gather(number,sample,-letter) %>%
            mutate(location=paste0(letter,number)) %>%
            filter(!is.na(sample)) %>%
            select(location,sample) %>%
            mutate(control=str_detect(sample,
                                      "Blank|Pos|Pool|PBS"))%>%
            mutate(id=sapply(sample,FUN = get_id)) %>%
            mutate(id=ifelse(control,sample,id)) %>%
            mutate(day=sapply(sample,FUN = get_day)) %>%
            mutate(day=as.numeric(ifelse(control,NA,day))) %>%
            mutate(layout=as.character(r)) %>%
            mutate(plateType=type)
        
        
        df <- bind_rows(df,tmp)
        
    }
    
    #get and add file name
    path_unlisted <- path %>% str_split(.,"/") %>% unlist()
    df <- df %>% mutate(file=path_unlisted[length(path_unlisted)])
    
    
    return(df)
    
}


#figure out the dilution from a positive control
get_dilution <- function(string){
    as.numeric(str_remove(unlist(str_split(string,":"))[2],","))
}

#brings in demographic data for individuals
getWideData <- function(){
    
    #cases and contacts
    weill_data1 <- readstata13::read.dta13("data/raw_data/cohort_data/CIRS PIC SMIC Master DB March 2019.dta") 
    labels <- readstata13::varlabel(weill_data1)
    labels <- data.frame(variable=names(labels),description=as.character(labels))
    weill_data2 <- weill_data1 %>% 
        # filter(studycode %in% c("PIC","SMIC")) %>%
        # mutate(id=ifelse(studycode=="PIC","P","S")) %>%
        mutate(id=recode(studycode,
                         "PIC"="P",
                         "SMIC"="S")) %>%
        mutate(id=paste0(id,sid)) %>%
        select(studycode,id,date_enroll=dten,age=yoa,sex,blood=abo,status=casecon,
               inf_category,culture=cx, 
               dehydration=dehycas,
               `Case Watery Diarrhea Day 7`= wdcas,
               `Contact Watery Diarrhea Day 2`=wdcon2,
               `hospital_hrs`=durstay
               
        ) %>%
        mutate(status=factor(status,labels=c("Case","Contact"))) %>%
        mutate(inf_category=factor(inf_category,
                                   labels= c("Uninfected",
                                             "Symptomatic w. seroconversion",
                                             "Symptomatic, no seroconversion",
                                             "Asymptomatic w. seroconversion",
                                             "Asymptomatic no seroconversion")
                                   
        )) %>%
        mutate(blood=str_remove(blood,"\\+|\\-| group"))
    
    
    
    
    #vaccinees
    
    # get age and sex for RIKS2
    vax_adult_demo1 <-  read_excel("data/raw_data/cohort_data/RIKS2_enrollment_gender_modified.xlsx") %>%
        rename(sex=Sexe,age=Laj,simpleid=id)  %>%
        mutate(sex=factor(sex,levels=c("M","F"),labels=c("male","female"))) %>%
        mutate(id=str_remove(`RIK ID Patisipan`," ")) %>%
        mutate(id=str_replace(id,"RIKS2","R2")) %>%
        mutate(id=str_replace(id,"-","_"))

    # add blood group for RIKS2
    vax_adult_demo2 <-  read_excel("data/raw_data/cohort_data/pntd.0007057.s002 (1).xlsx") %>%
        select(simpleid=`ID#`,Blood) %>%
        mutate(blood=factor(Blood, labels=c("O","A","B","AB"))) %>%
        right_join(vax_adult_demo1, by ="simpleid") %>%
        select(id,age,sex,blood) %>%
        mutate(
            studycode="RIKS2",
            status="Vaccinee")
    
    # get demographics sex for RIKS1
    # # child inventory
    vax_R1_inventory<- read_excel("data/raw_data/cohort_data/RIKS 1 Complete Data Set Haiti Only.xls") %>%
        select(`ID#`,Age,Gender,Blood,HIV) %>%
        mutate(id=paste0("R1_",`ID#`)) %>%
        mutate(age=as.numeric(Age)) %>%
        mutate(sex=factor(Gender,labels=c("male","female"))) %>%
        mutate(blood=factor(Blood,labels=c("O","A","B","AB"))) %>%
        mutate(HIV=factor(HIV,labels=c("no","yes"))) %>%
        select(id,age,sex,blood, HIV)%>%
        mutate(            
            studycode="RIKS1",
            status="Vaccinee")
    
    
    #get demographics for Bangladeshi vaccination
    vaccinee_shipment <- read_excel("data/raw_data/sample-availability/Bangladesh_vaccinee_Shipment_edited.xlsx",
                                    range = "B4:I58") %>%
            filter(!is.na(ID)) %>%
            mutate(Day=str_remove_all(Day,"Day \\(|\\)"))%>%
            rename(Age = `Age (Years)`,
                   Volume = `Volume( µl)`
            ) %>%
            mutate(study=str_extract(ID,"IMS|RB")) %>%
            mutate(ID=str_replace(ID,"RB","RB_")) %>%
            mutate(ID=str_replace(ID,"IMS-","IMS_"))
    
    
    
    #get the days for vaccinee
    day_df <- data.frame()
    for(i in 1:nrow(vaccinee_shipment)){
            day_df <- bind_rows(
                    day_df,
                    data.frame(id=vaccinee_shipment$ID[i],
                               day=str_split(vaccinee_shipment$Day,",")[[i]] %>%
                                       as.numeric()
                               
                    )
                    
            ) 
            
    }
    
    #final roster
    BGD_vax_roster <- vaccinee_shipment %>%
            select(id=ID,
                   sex=Gender,age=Age,study_code=study) %>%
            right_join(day_df) %>%
            #Kian was listing as day 1, confirmed it is day 0
            # mutate(day=ifelse(day==0,1,day)) %>%
            mutate(sample=paste0(id,"_D",day))  %>%
            rename(studycode=study_code) %>%
            distinct(studycode,id,age,sex) %>%
            mutate(status="Vaccinee") %>%
            mutate(sex=tolower(sex))
    
    
    
        
    
    final <- bind_rows(weill_data2,
                       vax_adult_demo2,vax_R1_inventory,
                       BGD_vax_roster
                       ) %>%
                mutate(age_group= cut(age,c(0,5,10,18,1000),right=FALSE) ) %>%
                mutate(age_group=factor(age_group,labels = c("<5 years","5-9 years","10-17 years", "18+ years")))
    
    
    return(final)
}


#brings in timing data for individuals
getTimingData <- function(){
    
        # 2021-09-24 see source/final_code/shared/estimateSampleTiming.Rmd
        out <- readRDS("data/generated_data/imputed_sample_times.rds")
    
    return(out)
    
}

#brings in timing data for individuals
getTimingData2 <- function(){
        
        # 2021-09-24 see source/final_code/shared/estimateSampleTiming.Rmd
        # switch to the code in main_text.Rmd eventually
        case_timing <- readRDS("data/generated_data/imputed_sample_times.rds")
        
        
        #vaccination timing data
        vax_dates_import <- read_excel("data/raw_data/cohort_data/RIK Sero2 Electronic Lab Registry_ Sample inventory.xls")%>%
                select(`RIK ID No.`,starts_with("Date")) %>%
                gather(day,date,-(`RIK ID No.`))%>%
                #get the id
                mutate(id=str_remove(`RIK ID No.`,"RIK S2-|RIKS2-")) %>%
                mutate(id=paste0("R2_",as.numeric(id))) %>%
                #remove missing dates
                filter(!is.na(date)) %>%
                mutate(date=as.Date(date))%>%
                #day
                mutate(day=as.integer(str_remove(day,"Date of Day ")))%>%
                mutate(sample=paste0(id," d",day))%>%
                select(id,sample,day,date)
        
        
        vax_baseline_date <- filter(vax_dates_import,day==0L) %>%
                select(id,baseline_date=date)
        
        #calculate difference
        vax_timing <- vax_dates_import %>% 
                left_join(vax_baseline_date, by="id") %>%
                mutate(day_actual=as.numeric(date-baseline_date)) %>%
                select(id,sample,day,day_actual)
        
        #limit to dates where we actually measured the samples
        studydata <- read_rds("data/generated_data/rau_data/2021-10-18_rau_data.rds") %>%
                        distinct(sample,id,day) %>%
                        select(id,day,sample)
        
        #combine and output
        out <- bind_rows(case_timing,vax_timing) %>%
                right_join(studydata,by=c("id","day","sample")) %>%
                #need to include a time for RIKS1 and RIKS2_27 day 90
                # if we are missing a date, just use the regular day in this place
                #as written now, should only apply to vaccinees
                mutate(day_actual=ifelse(is.na(day_actual),
                                                 day, day_actual))

                
        
        return(out)
        
}


#brings in timing data for individuals
getTimingData3 <- function(){
        
        
        #limit which id's are in our study versus not
        studydata <- read_rds("data/generated_data/rau_data/2021-10-18_rau_data.rds") %>%
                distinct(sample,id,day) %>%
                select(id,day,sample)
        
        # 2021-09-24 see source/final_code/shared/estimateSampleTiming.Rmd
        pic_DATES <- bind_rows(
                
                # file 1
                read.csv(here::here('data','raw_data','cohort_data','dates_data','pic_date_OCT.csv'),as.is=T)%>%
                        select(SID,starts_with("CDAT"))  %>%
                        gather(day,date_visit,-SID) %>%
                        mutate(day=str_remove(day,"CDAT") )%>%
                        mutate(day=case_when(
                                day =="7D" ~ 7,
                                day =="1M" ~ 30,
                                day =="3M" ~ 90,
                                day =="6M" ~ 180,
                                day =="9M" ~ 270,
                                day =="12M" ~ 360)) %>%
                        filter(date_visit!="")%>%
                        mutate(date_visit=as.Date(date_visit,
                                                  format="%d-%b-%Y",origin="1970-1-1")),
                
                
                # #file 2
                read.csv(here::here('data','raw_data','cohort_data','dates_data','pic_final_2012_FromAshraf_6.2015.csv'),as.is=T)%>%
                        select(sid,starts_with("CDAT")) %>%
                        rename(SID=sid) %>%
                        gather(day,date_visit,-SID) %>%
                        mutate(day=str_remove(day,"cdat") %>% toupper) %>%
                        filter(day %in% c("7D","1M","3M","6M",
                                          "9M","12M"
                        ))%>%
                        mutate(day=case_when(
                                day =="7D" ~ 7,
                                day =="1M" ~ 30,
                                day =="3M" ~ 90,
                                day =="6M" ~ 180,
                                day =="9M" ~ 270,
                                day =="12M" ~ 360)) %>%
                        filter(date_visit!="")%>%
                        mutate(date_visit=as.Date(date_visit,
                                                  format="%d/%m/%Y",origin="1970-1-1")) 
                
        )  %>%
                mutate(id=paste0("P",SID)) %>%
                mutate(sample=paste0(id," d",day)) %>%
                distinct(id,day,sample,date_visit)
        
        
        ## to make it easier to get into same format at smic_lab, I will deal with the three visits
        ## seperatley. Not efficient but...
        rename_smic_lab <- function(dat,day=7){
                dat  %>% rename_at(.vars = vars(ends_with(day)),
                                   .funs = funs(sub(paste0(day,"$"), "", .))) %>%
                        rename_all(toupper)
                
        }
        
        ## merge with long term dates for cases
        path_to_excel_file <- here('data','raw_data','cohort_data','dates_data','SMIC follow up Enrollment Dates.xlsx')
        
        smic_post365_dates <- full_join(
                read_excel(path_to_excel_file,sheet = "D540",na=c(NA,"Dropout","--")) %>% rename(SID = `SMIC ID`,
                                                                                                 date_540d =  `Date of Enrollment, D540`),
                read_excel(path_to_excel_file,sheet = "D720",na=c(NA,"Dropout","--")) %>% rename(SID = `SMIC ID`,
                                                                                                 date_720d =  `Date of Enrollment, D720`)) %>% 
                full_join(.,
                          read_excel(path_to_excel_file,sheet = "D900",na=c(NA,"Dropout","--")) %>% rename(SID = `SMIC ID`,
                                                                                                           date_900d =  `Date of Enrollment, D900`)
                ) %>%
                mutate(date_540d = as.Date(date_540d,format = "%d/%m/%y"),
                       date_720d = as.Date(date_720d,format = "%d/%m/%y"),
                       date_900d = as.Date(date_900d,format = "%d/%m/%y"),
                       date_540d = if_else(SID == 29,as.Date("2013-10-01"),date_540d) ## fixing this seemingly clear mistake
                ) %>%
                filter(!(is.na(date_540d) & is.na(date_720d) & is.na(date_900d)))
        
        
        smic_lab2_540 <- read.csv(here('data','raw_data','cohort_data','dates_data','SMICdata.vbx.plasma.Azman.2018.csv')) %>%
                select(sid,matches("7$")) %>%
                rename_smic_lab(day="7") %>%
                mutate(LBDAY = 540) %>%
                select(-starts_with("pm"),-starts_with("sm")) ## getting rid of IgM since we bring in seperateley
        
        
        
        smic_DATES<- bind_rows(
                
                #file 3
                read.csv(here('data','raw_data','cohort_data','dates_data','smic_date_OCT.csv'),as.is=T)%>%
                        select(SID,starts_with("CDAT"))  %>%
                        gather(day,date_visit,-SID) %>%
                        mutate(day=str_remove(day,"CDAT") )%>%
                        mutate(day=case_when(
                                day =="7D" ~ 7,
                                day =="1M" ~ 30,
                                day =="3M" ~ 90,
                                day =="6M" ~ 180,
                                day =="9M" ~ 270,
                                day =="12M" ~ 360))%>%
                        filter(date_visit!="")%>%
                        mutate(date_visit=as.Date(date_visit,
                                                  format="%d-%b-%Y",origin="1970-1-1")),
                
                
                
                #file 4
                read.csv(here('data','raw_data','cohort_data','dates_data','SMICdata.vbx.plasma.Azman.2018.csv')) %>%
                        select(SID=sid,starts_with("CDAT"))%>%
                        mutate(date_visit=as.character(cdat7d)) %>%
                        select(-cdat7d) %>%
                        mutate(day=7) %>%
                        filter(date_visit!="") %>%
                        mutate(date_visit=as.Date(date_visit,
                                                  format="%d-%b-%y",origin="1970-1-1")),
                
                
                #file 5
                smic_post365_dates %>%
                        gather(day,date_visit,-SID) %>%
                        mutate(day=str_remove(day,"date_")) %>%
                        mutate(day=str_remove(day,"d") %>% as.numeric) %>%
                        filter(!is.na(date_visit)),
                
                
                #         #file 6
                smic_lab2_540 %>% select(SID,date_visit=QFHDAT,day=LBDAY) %>%
                        filter(date_visit!="") %>%
                        mutate(date_visit=as.Date(date_visit,
                                                  format="%d-%b-%y",origin="1970-1-1"))
        ) %>% mutate(id=paste0("S",SID)) %>%
                mutate(sample=paste0(id," d",day)) %>%
                distinct(id,day,sample,date_visit)
        
        
        caseDATES <- bind_rows(smic_DATES,pic_DATES)        
        
        master_data <- readstata13::read.dta13("data/raw_data/cohort_data/CIRS PIC SMIC Master DB March 2019.dta") %>%
                filter(studycode %in% c("SMIC","PIC")) %>%
                # #limit to only cases no contacts
                # #contact timing does not really matter
                filter(casecon==0) %>%
                mutate(studycode=substr(studycode,1,1)) %>%
                mutate(id=paste0(studycode,sid)) %>%
                select(id, date_enroll=dten,diadur)
        
        case_timing <- studydata %>%
                inner_join(master_data,by="id") %>%
                left_join(caseDATES) %>%
                #calculate the difference
                mutate(diff=(date_visit-date_enroll)%>% as.numeric) %>%
                mutate(diff=ifelse(is.na(diff),day,diff)) %>%
                # - diarrhea duration
                # - incubation period (1.4 days)
                dplyr::mutate(symptom_lag = ifelse(is.na(diadur),0,diadur/24)) %>%
                #calculate time to enrollment
                mutate(tte=round(symptom_lag+1.4)) %>%
                mutate(day_actual=tte+diff) %>%
                #change all day 2 to 4 days (we dont have dates)
                mutate(day_actual=ifelse(day==2,4,day_actual)) %>%
                select(id,sample,day,day_actual)
        
        
        
        #vaccination timing data
        vax_dates_import <- read_excel("data/raw_data/cohort_data/RIK Sero2 Electronic Lab Registry_ Sample inventory.xls")%>%
                select(`RIK ID No.`,starts_with("Date")) %>%
                gather(day,date,-(`RIK ID No.`))%>%
                #get the id
                mutate(id=str_remove(`RIK ID No.`,"RIK S2-|RIKS2-")) %>%
                mutate(id=paste0("R2_",as.numeric(id))) %>%
                #remove missing dates
                filter(!is.na(date)) %>%
                mutate(date=as.Date(date))%>%
                #day
                mutate(day=as.integer(str_remove(day,"Date of Day ")))%>%
                mutate(sample=paste0(id," d",day))%>%
                select(id,sample,day,date)
        
        
        vax_baseline_date <- filter(vax_dates_import,day==0L) %>%
                select(id,baseline_date=date)
        
        #calculate difference
        vax_timing <- vax_dates_import %>% 
                left_join(vax_baseline_date, by="id") %>%
                mutate(day_actual=as.numeric(date-baseline_date)) %>%
                select(id,sample,day,day_actual) %>%
                #as written now, should only apply to vaccinees
                right_join(studydata %>% filter(str_detect(id,"R2")),
                           by=c("id","sample","day")) %>%
                # if we are missing a date (RIKS2_27 day 90), just use the regular day in this place
                mutate(day_actual=ifelse(is.na(day_actual),
                                         day, day_actual))
        
        
        
        
        #combine and output
        out <- bind_rows(case_timing,vax_timing) %>%
                right_join(studydata,by=c("id","day","sample")) %>%
                #if RIKS1 or contact, assume exact timing of sample collection
                mutate(day_actual=ifelse(str_detect(id,"R1|\\."),day,day_actual))
        
        return(out)
        
}




#loads data from the dropbox and associate plate layouts
#loads data from the dropbox and associate plate layouts
getLongLuminexRaw <- function(reload=TRUE){
        
        if(reload){
                #bring in dropbox token
                token <- read_rds('data/generated_data/dropbox_token/token.rds')
                
                
                # List of all the filenames you want to read in
                file_list <- drop_dir("vc_luminex_longtudinal/raw_plate_data/final/summaryFI",dtoken = token)
                
                #read in the data
                run_list <- list()
                for (i in 1:nrow(file_list)){
                        # for (i in (nrow(file_list)-2):nrow(file_list)){
                        
                        run_list[[file_list$name[i]]] <- extractData(file_list$path_display[i],token)
                }
                
                
                #get data frames ready
                netmfi_df <-  run_list[[1]]$`Net MFI` %>%
                        mutate(path=file_list$path_display[1])
                median_df <-  run_list[[1]]$Median
                count_df <-   run_list[[1]]$Count
                cv_df <-      run_list[[1]]$`%CV`
                trimcv_df <-  run_list[[1]]$`Trimmed % CV of Microspheres`
                
                for (i in 2:length(run_list)){
                        netmfi_df <- bind_rows(netmfi_df,
                                               run_list[[i]]$`Net MFI`%>%
                                                       mutate(path=file_list$path_display[i]))
                        median_df <- bind_rows(median_df,run_list[[i]]$Median)
                        count_df <- bind_rows(count_df,run_list[[i]]$Count)
                        cv_df <- bind_rows(cv_df,run_list[[i]]$`%CV`)
                        trimcv_df <- bind_rows(trimcv_df,run_list[[i]]$`Trimmed % CV of Microspheres`)
                        
                }
                
                
                total_df <- median_df %>% left_join(netmfi_df) %>%
                        left_join(cv_df) %>%
                        left_join(count_df) %>%
                        left_join(trimcv_df)
                
                
                layout_metadata<- drop_read_excel("vc_luminex_longtudinal/layout_lookup/layout_matching.xlsx",
                                                  dtoken=token, sheet=1,
                                                  # col_names=FALSE,
                                                  col_types="text") %>%
                        select(-date)
                
                #get plate layouts
                classes_layout <- getLayout("vc_luminex_longtudinal/layout_lookup/layouts/Final Plate Layout_Classes.xlsx",
                                            ranges=c("A12:Y28",
                                                     "A32:Y48",
                                                     "A52:Y68"),
                                            type="Case",token
                ) %>% mutate(class="Full Isotypes")
                
                subclasses_layout <- getLayout("vc_luminex_longtudinal/layout_lookup/layouts/Final Plate Layout_Subclasses.xlsx",
                                               ranges=c("A11:Y27",
                                                        "A31:Y47",
                                                        "A51:Y67"),
                                               type="Case",token
                ) %>% mutate(class="Non-IgG1/IgG2 subclass")
                
                IGG1N2_layout <- getLayout("vc_luminex_longtudinal/layout_lookup/layouts/Final Plate Layout_Subclasses_IgG1 and IgG2.xlsx",
                                           ranges=c("A11:Y27",
                                                    "A31:Y47",
                                                    "A51:Y67"),
                                           type="Case",token
                ) %>% mutate(class="IgG1/IgG2")
                
                dilutionseries1 <- getLayout("vc_luminex_longtudinal/layout_lookup/layouts/IgG_total_002 IgG_total_003 extended dilution series.xlsx",
                                             ranges=c("A2:Y18",
                                                      "A22:Y38"
                                             ),
                                             type="Case",token
                ) %>% mutate(class="Dilution Series")
                
                dilutionseries2 <- getLayout("vc_luminex_longtudinal/layout_lookup/layouts/IgG_total_002 IgG_total_003 extended dilution series.xlsx",
                                             ranges=c("A42:Y58"),
                                             type="Series",token
                ) %>% mutate(class="Dilution Series") %>%
                        mutate(layout="3")
                
                vaxlayout <- getLayout("vc_luminex_longtudinal/layout_lookup/layouts/Plate Layout_RIKS and SMIC Add-Ons_Clean.xlsx",
                                       ranges=c("A1:Y17",
                                                "A20:Y36",
                                                "A39:Y55"),
                                       type="Vaccine",token
                ) %>% mutate(class="Full Isotypes") 
                
                
                layouts <- bind_rows(classes_layout,subclasses_layout,IGG1N2_layout,dilutionseries1,dilutionseries2,vaxlayout) 
                
                split_on_underscore <- function(string,n){
                        unlist(str_split(string,"_"))[n]
                }
                
                analysis_df <- total_df %>%  #join with metadata to find proper layout
                        left_join(layout_metadata,by="batch") %>%
                        #match with sample ID
                        left_join(layouts,by = c("location", "plateType", "file", "layout"))%>%
                        #make antigen a factor
                        mutate(antigen= factor(antigen,
                                               levels=
                                                       c("Flu","Ogawa OSP:BSA",
                                                         "Inaba OSP:BSA", "O139:BSA",
                                                         "CtxB",  "LT-B" , "CT HT","LTh",
                                                         "VCC" ,  "Sialidase","TcpA"  )))%>%
                        #deal with IGX and subclass X
                        mutate(newisotype=ifelse(isotype=="IGX",
                                                 sapply(sample,split_on_underscore,n=1),
                                                 isotype)) %>% 
                        mutate(newsubclass=ifelse(subclass=="X",
                                                  sapply(sample,split_on_underscore,n=2),
                                                  subclass)) %>% 
                        mutate(sample=ifelse(isotype=="IGX"|subclass=="X",
                                             sapply(sample,split_on_underscore,n=3),
                                             sample)) %>% 
                        mutate(id=ifelse(isotype=="IGX"|subclass=="X",
                                         sapply(id,split_on_underscore,n=3),
                                         id)) %>%
                        #replace the IGX and subclass X with new values and remove extra variables
                        mutate(isotype=newisotype,sublass=newsubclass) %>%
                        select(-newisotype,-newsubclass) %>%
                        #change median to MFI
                        rename(MFI=Median) %>%
                        # find the dilution
                        mutate(dilution=sapply(sample,get_dilution)) %>%
                        #date 
                        mutate(date=as.Date(date,"%m/%d/%Y")) %>%
                        #make full_series a logical
                        mutate(full_series=as.logical(full_series)) %>%
                        #rename isotypes
                        mutate(isotype=str_replace(isotype,"G","g"))
        } else {
                
                analysis_df <- read_rds("data/generated_data/2021-10-25-raw-data.rds")
                
        }
        
        return(analysis_df)
}

# get a long-style database with avgMFI for each sample
getLongLuminexTidy <- function(reload
                              # logical to determine whether to load from
                              # dropbox or from a previous file
                              ){

    if(reload==TRUE){

        df <- getLongLuminexRaw(reload=TRUE) %>%
            #remove those plates that required a rerun
            filter(required_rerun==0) %>%
            #exclude those with a low bead count or missing MFI value
            filter(!is.na(MFI)) %>%
            filter(Count>=30) %>%
            #calculate the average MFI across triplicates
            group_by(batch,date,isotype,subclass,antigen,day,id,sample,control,dilution,full_series) %>%
            summarize(avgMFI=mean(MFI),
                      count=n()
                      ) %>%
            ungroup()

    }

    if(reload==FALSE){
        df <- read_rds("data/generated_data/tidy_data/2021-10-18-tidy_data.rds")
    }

    return(df)

}



#get old ELISA and vibriocidals
getLongSeroTidy <-function(){
    weill_data1 <- readstata13::read.dta13("data/raw_data/cohort_data/CIRS PIC SMIC Master DB March 2019.dta") 
    # labels <- readstata13::varlabel(weill_data1)
    # labels <- data.frame(variable=names(labels),description=as.character(labels))
    weill_data2 <- weill_data1 %>% 
        mutate(id=recode(studycode,
                         "PIC"="P",
                         "SMIC"="S")) %>%
        mutate(id=paste0(id,sid)) %>%
        
        #standardize the elisas
        mutate(
               #anti-lps IgG        
               lgptl2=   lgptl2/lgpp2, 
               lgptl7=   lgptl7/lgpp7, 
               lgptl30=  lgptl30/lgpp30, 
               lgptl90=  lgptl90/lgpp90,
               lgptl180= lgptl180/lgpp180, 
               lgptl270= lgptl270/lgpp270, 
               lgptl360= lgptl360/lgpp360, 
               lgptl540= lgptl540/lgpp540,
               lgptl900= lgptl900/lgpp900, 
               lgptl1080= lgptl1080/lgpp1080,
               #anti-lps IgA
               laptl2=laptl2/lapp2,
               laptl7=laptl7/lapp7,
               laptl30=laptl30/lapp30,
               laptl90=laptl90/lapp90,
               laptl180=laptl180/lapp180,
               laptl270=laptl270/lapp270, 
               laptl360=laptl360/lapp360, 
               laptl540=laptl540/lapp540,
               laptl900=laptl900/lapp900,
               laptl1080=laptl1080/lapp1080,
               
               #anti-lps IgM
               lmptl2=lmptl2/lmpp2,
               lmptl7=lmptl7/lmpp7,
               lmptl30=lmptl30/lmpp30,
               lmptl90=lmptl90/lmpp90,
               lmptl180=lmptl180/lmpp180,
               lmptl270=lmptl270/lmpp270,
               lmptl360=lmptl360/lmpp360,
               lmptl540=lmptl540/lmpp540,
               lmptl900=lmptl900/lmpp900,
               lmptl1080=lmptl1080/lmpp1080,

               #anti-ctxb IgG
               xgptl2=xgptl2/xgpp2, 
               xgptl7=xgptl7/xgpp7,
               xgptl30=xgptl30/xgpp30,
               xgptl90=xgptl90/xgpp90,
               xgptl180=xgptl180/xgpp180,
               xgptl270=xgptl270/xgpp270,
               xgptl360=xgptl360/xgpp360,
               xgptl540=xgptl540/xgpp540,
               xgptl900=xgptl900/xgpp900,
               xgptl1080=xgptl1080/xgpp1080,
               
               #anti-ctxb IgA
               xaptl2=xaptl2/xapp2,
               xaptl7=xaptl7/xapp7,
               xaptl30=xaptl30/xapp30,
               xaptl90=xaptl90/xapp90,
               xaptl180=xaptl180/xapp180,
               xaptl270=xaptl270/xapp270,
               xaptl360=xaptl360/xapp360,
               xaptl540=xaptl540/xapp540,
               xaptl900=xaptl900/xapp900,
               xaptl1080=xaptl1080/xapp1080,
               
               #anti-ctxb IgM
               xmptl2=xmptl2/xmpp2,
               xmptl7=xmptl7/xmpp7,
               xmptl30=xmptl30/xmpp30,
               xmptl90=xmptl90/xmpp90,
               xmptl180=xmptl180/xmpp180,
               xmptl270=xmptl270/xmpp270,
               xmptl360=xmptl360/xmpp360,
               xmptl540=xmptl540/xmpp540,
               xmptl900=xmptl900/xmpp900,
               xmptl1080=xmptl1080/xmpp1080
               
        )  %>%
        
        select(id,
               #vibriocidal ogawa
               vibogaw2=vibogaw,
               vibogaw7, vibogaw30, vibogaw90,
               vibogaw180, vibogaw270, vibogaw360,vibogaw540,
               vibogaw900,
               #vibriocidal inaba
               vibinab2=vibinab,
               vibinab7, vibinab30, vibinab90,
               vibinab180, vibinab270, vibinab360,vibinab540,
               vibinab900,
               
               #anti-lps IgG        
               lgptl2, lgptl7, lgptl30, lgptl90,
               lgptl180, lgptl270, lgptl360, lgptl540,
               lgptl900,lgptl1080, 
               #anti-lps IgA        
               laptl2, laptl7, laptl30, laptl90,
               laptl180, laptl270, laptl360, laptl540,
               laptl900, laptl1080,
               #anti-lps IgM        
               lmptl2, lmptl7, lmptl30, lmptl90,
               lmptl180, lmptl270, lmptl360, lmptl540,
               lmptl900, lmptl1080,
               
               #anti-ctxb IgG        
               xgptl2, xgptl7, xgptl30, xgptl90,
               xgptl180, xgptl270, xgptl360, xgptl540,
               xgptl900, xgptl1080,
               #anti-ctxb IgA        
               xaptl2, xaptl7, xaptl30, xaptl90,
               xaptl180, xaptl270, xaptl360, xaptl540,
               xaptl900,xaptl1080,
               #anti-ctxb IgM        
               xmptl2, xmptl7, xmptl30, xmptl90,
               xmptl180, xmptl270, xmptl360, xmptl540,
               xmptl900,xmptl1080
               
        )  %>%
        # select(-studycode,-fid,-sid)%>%
        gather(code,value,-id) %>%
        mutate(day=as.numeric(str_remove(code,
                                         "vibogaw|vibinab|lgptl|laptl|lmptl|xgptl|xaptl|xmptl"))) %>%
        mutate(full_test=str_remove(code,"1080|900|540|360|270|180|90|30|7|2")) %>%
        mutate(full_test=recode(full_test,
                                "vibogaw"="Vibriocidal_Ogawa",
                                "vibinab"="Vibriocidal_Inaba",
                                
                                "lgptl"="ELISA_LPS_IgG",
                                "laptl"="ELISA_LPS_IgA",
                                "lmptl"="ELISA_LPS_IgM",
                                
                                "xgptl"="ELISA_CtxB_IgG",
                                "xaptl"="ELISA_CtxB_IgA",
                                "xmptl"="ELISA_CtxB_IgM"
        )) %>%
        select(-code) %>%
        mutate(isotype= ifelse(str_detect(full_test,"IgG"),"IgG",NA),
               isotype= ifelse(str_detect(full_test,"IgA"),"IgA",isotype),
               isotype= ifelse(str_detect(full_test,"IgM"),"IgM",isotype),
               
               subclass=ifelse(is.na(isotype),NA,"total"),
               
               antigen = ifelse(str_detect(full_test,"LPS"),"LPS",NA),
               antigen = ifelse(str_detect(full_test,"CtxB"),"CtxB",antigen),
               
               test_type = ifelse(str_detect(full_test,"ELISA"),"ELISA",NA),
               test_type = ifelse(str_detect(full_test,
                                             "Vibriocidal"),"Vibriocidal",test_type),
               
               
               unit = ifelse(test_type=="Vibriocidal","Titer",NA),
               unit = ifelse(test_type=="ELISA","ug per mL?",unit),
               
               
               
        ) %>%
        filter(!is.na(value)) %>%
        ungroup()
        
    return(weill_data2)
}



##### variable adjustment

# function to predict from the standard dose response curve
antiLL4 <- function(y,fit){
    # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0146021
    
    pVec <- as.numeric(fit$parmMat)
    
    b <- pVec[1]
    c <- pVec[2]
    d <- pVec[3]
    e <- pVec[4]
    #inverse of equation 2
    x <- exp(log((d-c)/(y-c)-1)/b+log(e))
    
    return(x)
}

LuminexRecPaperData <- function(vaccinee=FALSE){
    
    
    # bring in demographics and timing data
    demo_id <-getWideData() %>% mutate(blood=factor(blood))
    case_days<- getTimingData3() 
    
    #bring in the Luminex data
    #include both net mfi and RAU
    compare_data<- getLuminexSamplesRAU() %>%
            #pretty up some of the text
                mutate(
                    antigen=str_replace_all(antigen," |:|-",""),
                    test_type="Luminex",
                    full_test = paste(test_type,isotype,antigen,sep="_"),
                ) %>%
            #log transform the values
                mutate(RAU=log10(1/RAU),
                       `Net MFI` =ifelse(`Net MFI`<10,10,`Net MFI`),
                       `Net MFI`=log10(`Net MFI`)
                       ) %>%
            #incorporate demographic and timing data
            left_join(demo_id, by="id") %>%
            filter(!status %in% "Vaccinee"[!vaccinee]) %>%
            left_join(case_days, by=c("id","sample","day")) 
            
    
    #Vibriocidals/ ELISA
    vibELISA_data <-getLongSeroTidy() %>% 
        mutate(sample=paste0(id," d",day)) %>%
        filter(sample %in% unique(compare_data$sample)) %>% # limit to those in Luminex sample
        #large amount of missingness for IGM data
        filter(isotype %in% c("IgA","IgG") | test_type=="Vibriocidal") %>%
        # change vibriocidal to be logged
        mutate(value=ifelse(test_type=="Vibriocidal" & !is.na(value),
                            log2(value/5),value)) %>%
        #censor data to value of 11 to deal with weird decimal values
        mutate(value=ifelse(test_type=="Vibriocidal" & value >11,11,value)) %>%
        mutate(unit=ifelse(test_type=="Vibriocidal","log2(titer/5)",unit)) %>%
        #change ELISA values to make more sense
        mutate(value=ifelse(test_type=="ELISA",log10(value),value)) %>%
        mutate(unit=ifelse(test_type=="ELISA","log10(ug per mL)",unit)) %>%
        #censor variable
        mutate(censor= ifelse(test_type=="Vibriocidal",as.integer(value),0L))
    
    #combine all assay data together
    assay_data_long <- bind_rows(
        vibELISA_data,
        #rau data
        compare_data %>%
            mutate(value=RAU,
                   unit="log10(1/RAU)",
                   #provide censoring values for the RAU
                   censor=recode(RAUcens,
                                 "<100" = 1L, # upper bound
                                 ">1e+05" = -1L, # lower bound
                                 .default = 0L)) %>% # exact value
            select(id,value,day,sample,full_test,isotype,
                   subclass,antigen,test_type,unit,censor),
        #net mfi data
        compare_data %>%
            mutate(value=`Net MFI`,
                   unit="log10(Net MFI)",
                   #censoring for net MFI
                   censor = ifelse(`Net MFI`==1,
                                   -1L,# lower bound (Net MFI <10)
                                   0L) # exact value
                   ) %>%
            select(id,value,day,sample,full_test,isotype,
                   subclass,antigen,test_type,unit,censor)
        ) %>%
        #limit to non-vaccinees
        left_join(demo_id) %>%
        filter(!status %in% "Vaccinee"[!vaccinee]) %>%
        left_join(case_days) %>%
        #define day groups
        mutate(day_contact = ifelse(day==2, "Day 2 / Contact", paste("Day",day)),
               day_contact = ifelse(status=="Contact", "Day 2 / Contact", day_contact),
               day_contact = factor(day_contact, 
                                    levels= c("Day 2 / Contact",
                                              "Day 7","Day 30","Day 90",
                                              "Day 180","Day 270","Day 360",
                                              "Day 540","Day 720","Day 900","Day 1080"))
        ) %>%
        #define infection windows
        mutate(inf_90=factor(ifelse(day_actual<=90 & day_contact!="Day 2 / Contact",
                                    1,0)),
               inf_180=factor(ifelse(day_actual<=180 & day_contact!="Day 2 / Contact",
                                     1,0)),
               inf_270=factor(ifelse(day_actual<=270 & day_contact!="Day 2 / Contact",
                                     1,0)),
               inf_365=factor(ifelse(day_actual<=365 & day_contact!="Day 2 / Contact",
                                     1,0)),
               inf_540=factor(ifelse(day_actual<=540 & day_contact!="Day 2 / Contact",
                                     1,0))) %>%
            #remove samples with individuals who were potentially infected during follow-up
            filter(!sample %in% c(
                          # potential reinfections
                          "P25 d360",
                          "P27 d180",
                          "P27 d270",
                          "P27 d360",
                          #potentially infected vaccinated people
                          # "R2_33 d360", #not selected
                          "R2_4 d180",
                          "R2_4 d360",
                          "R2_46 d180",
                          "R2_46 d360",
                          "R2_50 d360",
                          "R2_58 d360",
                          "R2_63 d360",
                          "R2_7 d180",
                          "R2_7 d360")
            )
            
            
    
    #get all datasets together for RAU output
    RAU_data =list()
    #limit data set to only RAU, vibs and ELISA
    RAU_data$long_data <- filter(assay_data_long,
                                    unit %in% c("log10(1/RAU)","log2(titer/5)","log10(ug per mL)"))

    #get data into form for prediction / decay models
    RAU_data$wide_data <- RAU_data$long_data %>%
        #change to wide
        select(id,sample, day_actual, day,
               inf_90,inf_180,inf_270,inf_365,inf_540,
               full_test,value) %>%
        spread(full_test,value) %>%
        #add demographic factors
        left_join(
            demo_id %>% select(id,status,age,sex,blood,culture)
        ) %>%
        #make dummy variables
        mutate(`Under 10` = ifelse(age<10,1,0)) %>%
        mutate(Female = ifelse(sex=="female",1,0)) %>%
        mutate(`O Blood`=ifelse(blood=="O",1,0)
        ) %>%
        mutate(`Ogawa`=ifelse(culture=="O1-Ogawa",1,0))

    #get censoring into a wide format
    RAU_data$wide_censor <- RAU_data$long_data %>%
        #change to wide
        select(id,sample, day_actual, day,
               inf_90,inf_180,inf_270,inf_365,inf_540,
               full_test,censor) %>%
        spread(full_test,censor) %>%
        left_join(demo_id)

    #get all datasets together for Net MFI output
    netMFI_data <- list()
    netMFI_data$long_data <- filter(assay_data_long,
                                 unit %in% c("log10(Net MFI)","log2(titer/5)","log10(ug per mL)"))

    #get data into form for prediction / decay models
    netMFI_data$wide_data <- netMFI_data$long_data %>%
        #change to wide
        select(id,sample, day_actual, day,
               inf_90,inf_180,inf_270,inf_365,inf_540,
               full_test,value) %>%
        spread(full_test,value) %>%
        #add demographic factors
        left_join(
            demo_id %>% select(id,status,age,sex,blood,culture)
        ) %>%
        #make dummy variables
        mutate(`Under 10` = ifelse(age<10,1,0)) %>%
        mutate(Female = ifelse(sex=="female",1,0)) %>%
        mutate(`O Blood`=ifelse(blood=="O",1,0)
        ) %>%
        mutate(`Ogawa`=ifelse(culture=="O1-Ogawa",1,0))

    #get censoring into a wide format
    netMFI_data$wide_censor <- netMFI_data$long_data %>%
        #change to wide
        select(id,sample, day_actual, day,
               inf_90,inf_180,inf_270,inf_365,inf_540,
               full_test,censor) %>%
        spread(full_test,censor)%>%
        left_join(demo_id)
    
    
    #new way forward?? have RAU and netmfi together
    
    new_long_data <- assay_data_long %>%
            mutate(full_test= case_when(
                    unit == "log10(Net MFI)" ~str_replace(full_test,"Luminex","NetMFI"),
                    unit == "log10(1/RAU)" ~str_replace(full_test,"Luminex","RAU"),
                    TRUE ~ as.character(full_test)
            ))
    
    new_wide_data <- new_long_data%>%
            #change to wide
            select(id,sample, day_actual, day,
                   inf_90,inf_180,inf_270,inf_365,inf_540,
                   full_test,value) %>%
            spread(full_test,value) %>%
            #add demographic factors
            left_join(
                    demo_id %>% select(id,status,age,sex,blood,culture)
            ) %>%
            #make dummy variables
            mutate(`Under 10` = ifelse(age<10,1,0)) %>%
            mutate(Female = ifelse(sex=="female",1,0)) %>%
            mutate(`O Blood`=ifelse(blood=="O",1,0)
            ) %>%
            mutate(`Ogawa`=ifelse(culture=="O1-Ogawa",1,0))
    
    new_wide_censor <- new_long_data  %>%
            #change to wide
            select(id,sample, day_actual, day,
                   inf_90,inf_180,inf_270,inf_365,inf_540,
                   full_test,censor) %>%
            spread(full_test,censor)%>%
            left_join(demo_id) %>%
            #make dummy variables
            mutate(`Under 10` = ifelse(age<10,1,0)) %>%
            mutate(Female = ifelse(sex=="female",1,0)) %>%
            mutate(`O Blood`=ifelse(blood=="O",1,0)
            ) %>%
            mutate(`Ogawa`=ifelse(culture=="O1-Ogawa",1,0))
    
    
    #output files
    output <- list(
                #for comparing NetMFI to RAU
                compare_data = compare_data,
                #RAU data (long and wide)
                RAU_data = RAU_data,
                #NetMFI data (long and wide)
                netMFI_data =netMFI_data,
                
                new_long_data=new_long_data,
                new_wide_data = new_wide_data,
                new_wide_censor = new_wide_censor
    )
           
    return(output)
}


##################### variable adjustment  #################################

find4PL_parameters <- function(){
    
    #only worry about totals not subclass for now
    dat <- getLongLuminexTidy(reload=FALSE) %>%
            filter(subclass=="total")
    
    #frequentist model
    
    # limit to the full dilution series
    full_series_data <- filter(dat,full_series
                               # batch %in% c("Series_IGX_total_001",
                               #              "Case_IGG_total_005",
                               #              "Case_IGG_total_006")
                               ) %>%
        filter(str_detect(sample,"Hi Pos"))
    
    # fit the frequentist model to each and  store parameters
    # for each full series plate, fit a frequentist model
    full_series_freq <- data.frame()
    
    for(i in unique(full_series_data$batch)){
        this_plate1 <- filter(full_series_data,batch==i)
        for(j in unique(this_plate1$isotype)){
            this_plate2 <- filter(this_plate1,isotype==j)
            for(k in unique(this_plate2$antigen)){
                
                cat(i,j,k,"\n")
                
                this_plate3 <- filter(this_plate2,antigen==k)
                
                tmp <- FreqFit(this_plate3) %>%
                    mutate(
                        batch=i,
                        isotype=j,
                        subclass="total",
                        antigen=k
                    )
                
                full_series_freq <- bind_rows(full_series_freq,tmp)
            }
        }
    }
    
    #combine the IGG full dilution series for priors
    full_series_freq_avg <- full_series_freq %>% group_by(antigen,isotype,parameter) %>%
        summarize(value=mean(value)) %>%
        filter(parameter %in% c("B","C","D","E"))
    
    
    
    #bayes model
    
    #open up stan
    library(rstan)
    options(mc.cores = parallel::detectCores())

    #get the reduced data series
    reduced_series_data <- filter(dat, !full_series
                                  # !(batch %in% c("Series_IGX_total_001",
                                  #               "Case_IGG_total_005",
                                  #               "Case_IGG_total_006" ))
    ) %>%
        filter(str_detect(sample,"Hi Pos|Blank")) 

    # fit the bayesian model to each and  store parameters
    # for each reduced series plate, fit a bayesian model
    reduced_series_bayes <- data.frame()
    for(i in unique(reduced_series_data$batch)){
        this_plate1 <- filter(reduced_series_data,batch==i)
        for(j in unique(this_plate1$isotype)){
            this_plate2 <- filter(this_plate1,isotype==j)
            for(k in unique(this_plate2$subclass)){
                this_plate3 <- filter(this_plate2,subclass==k)
                for(m in unique(this_plate3$antigen)){
                    this_plate4 <- filter(this_plate2,antigen==m)

                cat(i,j,k,m,"\n")
                # print(head(this_plate4))
                

                prior_centered<- full_series_freq_avg %>%
                    filter(
                        isotype==j,
                        antigen==m
                    ) %>% ungroup()
                
               
                

                tmp <- BayesFit(this_plate4,
                                B_mu = prior_centered %>%
                                    filter(parameter=="B") %>%
                                    select(value) %>%unlist(),
                                D_mu = prior_centered %>%
                                    filter(parameter=="D") %>%
                                    select(value) %>%unlist(),
                                E_mu = prior_centered %>%
                                    filter(parameter=="E") %>%
                                    select(value) %>%unlist()
                )
                tmp <- tmp[[1]] %>%
                    mutate(
                        batch=i,
                        isotype=j,
                        subclass=k,
                        antigen=m
                    )
                
                print(head(tmp))

                reduced_series_bayes <- bind_rows(reduced_series_bayes,tmp)
            
                }
            }
        }
    }
    
    
    parameters <- bind_rows(
                    full_series_freq,
                    reduced_series_bayes
            ) %>%  filter(parameter %in% c("B","C","D","E"))
    
    
    return(parameters)
    
    
}



getLuminexSamplesRAU <- function(
                                 pred_new=FALSE,
                                 param=NA, #generated from 4PL_parameters()
                                 generate_plots=FALSE
                                 ){
    #this data includes both samples and controls
    #only find the RAU for those that are total
    dat <-getLongLuminexTidy(reload=FALSE) %>%
        filter(subclass=="total")
    
    if (pred_new==TRUE){
        
        #get sample data  ready
        sample_data <- dat %>%
            filter(!control) %>% mutate(RAU=NA) 
        
        sample_data_RAU <- data.frame()

        # df to store dilution curves 
        dilution_curves <- data.frame()
        x_s <- 10^(seq(0,10,0.1))
     
        #limit to only total subclass parameters (find4PL_parameters output)
        parameters <- filter(param, subclass=="total")

        for(i in unique(parameters$batch)){
            this_plate1 <- filter(parameters,batch==i)
            for(j in unique(this_plate1$isotype)){
                this_plate2 <- filter(this_plate1,isotype==j)
                for(k in unique(this_plate2$subclass)){
                    this_plate3 <- filter(this_plate2,subclass==k)
                    for(m in unique(this_plate3$antigen)){
                        this_plate4 <- filter(this_plate2,antigen==m)
                        
                        
                        cat(i,j,k,m,"\n")
                       
                        if(i != "Series_IGX_total_001"){
                            
                        
                        min_dil <- 100
                        max_dil <- 100000
                        
                        ub <- getCurve4(this_plate4,min_dil)
                        lb <- getCurve4(this_plate4,max_dil)
                        
                        sample_data_RAU_plate <- sample_data %>%
                            #first just fit to the curve
                            filter(batch==i & 
                                                  isotype==j &
                                                  subclass==k &
                                                  antigen==m) %>%
                                mutate(RAU= getRAU4(this_plate4,avgMFI))  %>%
                                mutate(RAU=ifelse(avgMFI>ub,min_dil,RAU)) %>%
                                mutate(RAU=ifelse(avgMFI<lb,max_dil,RAU)) %>%
                                #show as a censored value
                                mutate(RAUcens=ifelse(avgMFI>ub,paste0("<",min_dil),as.character(RAU))) %>%
                                mutate(RAUcens=ifelse(avgMFI<lb,paste0(">",max_dil),as.character(RAUcens)))                        
                        
                       
                        
                        
                        sample_data_RAU <- bind_rows(sample_data_RAU,
                                                     sample_data_RAU_plate)
                        
                        
                        
                        
                        }
                        
                        # get dilution curves
                        tmp <- this_plate4 %>%
                            getCurve4(x_s)
                        
                        
                        tmp_df <- data.frame(
                            dilution= x_s,
                            avgMFI=tmp,
                            batch=i,
                            isotype=j,
                            subclass=k,
                            antigen=m
                        )
                        
                        dilution_curves <- bind_rows(dilution_curves,tmp_df)
                        
                       
                        
                        
                        
                        
                    }
                    }
            }
        }  
        
    # plot all of the dilution curves
    if(generate_plots==TRUE){
        
        cat("Now the plots....")
        
        
    for(i in unique(parameters$batch)){
        this_plate1 <- filter(parameters,batch==i)
        # if(i != "Series_IGX_total_001"){
        for(j in unique(this_plate1$isotype)){
            this_plate2 <- filter(this_plate1,isotype==j)
            for(k in unique(this_plate2$subclass)){
                
                cat(i,j,k,m,"\n")
                

                        p <- plotDilutionCurve(MFIdat=dat,
                                               curve=dilution_curves,
                                               plate=i,
                                               iso=j,
                                               sub=k)

                        ggsave(plot=p,
                               filename=paste0("figures/exploratory_figures/dilution_curves/",
                                               Sys.Date(),i,j,k,".pdf"))


                    }}}
      }
        
        # return(sample_data_RAU)
            
        
        
        
        
    }
    
    if (pred_new==FALSE){
        sample_data_RAU <- read_rds("data/generated_data/rau_data/2021-10-18_rau_data.rds") %>%
                select(-blankMFI,-`Net MFI`)

    }
    
    #also calculate the Net MFI
    blanks <- dat %>% ungroup() %>% filter(sample=="Blank") %>%
        distinct(batch,isotype,subclass,antigen,blankMFI=avgMFI)
    final <- sample_data_RAU %>% left_join(blanks,
                                           by = c("batch", "isotype", "subclass", "antigen")
                                           ) %>%
            mutate(`Net MFI`=avgMFI-blankMFI)
    
    return(final)
    
}


FreqFit <- function(data){
    
    fit_ll4 <- drc::drm(data=data,log(avgMFI)~dilution,fct=drc::LL.4())
    
    
    out <- data.frame(
        parameter= names(fit_ll4$coefficients) %>%
            substr(1,1) %>%
            toupper(),
        value = as.numeric(fit_ll4$coefficients)   
    )
    
    return(out)
    
}

BayesFit <- function(data,  warmup =1000, iter = 5000,
                     # full_series = FALSE,
                     B_mu,#=0.7,
                     B_sigma=1,
                     D_mu,#=10,
                     D_sigma=20,
                     E_mu,#=1000,
                     E_sigma=1000,
                     C_sigma=0.5
                     
){
    
    # separate Blank and dilution data
    dilution_data_reduced <- data %>% filter(sample!="Blank")
    lower <- data %>% filter(sample=="Blank") %>% ungroup() %>%
        select(avgMFI) %>% unlist()
    
    
    # shared priors for 4PL
    shared_priors4 <- list(
        avg_B=B_mu,
        sd_B=B_sigma,
        avg_D=D_mu,
        sd_D=D_sigma,
        avg_E=E_mu,
        sd_E=E_sigma
    )
    #stan format for data
    stan_data <- list(
        N = nrow(dilution_data_reduced),
        y= log(dilution_data_reduced$avgMFI),
        x= dilution_data_reduced$dilution,
        avg_C = log(lower),
        sd_C= C_sigma
    )
    
    
    # print(c(c(stan_data, shared_priors4)))
    
    #compile stan code
    # fileName <- "source/exploration/snippets/12-loglogistic4-model.stan"
    # ret <- stanc(fileName) # Check Stan file
    # ret_sm <- stan_model(stanc_ret = ret) # Compile Stan code
    # save(fileName, ret, ret_sm, file="source/exploration/stan/12-loglogistic4-model.stan.Rdata")
    
    stan_model_path <- "source/final_code/luminex_recommendation/stan/loglogistic4-model.stan"
    
    dilution_model <- stan(stan_model_path, warmup=warmup, iter=iter,
                               seed=123,
                               data=c(stan_data, shared_priors4),
                               chains=4#,
                               # control=list(adapt_delta=.95)
                           )
    # ret_sm <- load(file="source/exploration/stan/12-loglogistic4-model.stan.Rdata")
    # 
    # dilution_model <- sampling(ret_sm, warmup=warmup, iter=iter,
    #                            seed=123,
    #                            data=c(stan_data, shared_priors4),
    #                            chains=4,
    #                            control=list(adapt_delta=.95))
    
    output <- summary(dilution_model)$summary %>% data.frame() %>% mutate(parameter=rownames(.)) %>%
        select(parameter,value=X50.)
    
    return(list(output,dilution_model))
    
}

# new getCurve4
getCurve4<- function(df,dilution_x){
    
    B <- df %>% filter(parameter=="B") %>% select(value) %>%unlist()
    C <- df %>% filter(parameter=="C") %>% select(value) %>%unlist()
    D <- df %>% filter(parameter=="D") %>% select(value) %>%unlist()
    E <- df %>% filter(parameter=="E") %>% select(value) %>%unlist()
    
    pred <- C+ (D-C)/(1+exp(B*(log(dilution_x) -log(E))))
    
    return(exp(pred))
    
}

# new  RAU function
getRAU4<- function(df,MFI){
    
    B <- df %>% filter(parameter=="B") %>% select(value) %>%unlist()
    C <- df %>% filter(parameter=="C") %>% select(value) %>%unlist()
    D <- df %>% filter(parameter=="D") %>% select(value) %>%unlist()
    E <- df %>% filter(parameter=="E") %>% select(value) %>%unlist()
    
    RAU = exp(log((D-C)/(log(MFI) - C) -1)/B +log(E))
    
    return(RAU)
    
}



##################### analysis #################################


#complete MW test
getMannWhitney <- function(v1,v2){
    test <- wilcox.test(v1~v2,
                        paired=FALSE,
                        digits.rank=7
    )
    test$p.value
    
}

#complete wilcoxon rank sum test
getWilcoxonRS <- function(v0,v1,v2,v3){
    
    #v0 = id
    #v1 = value
    #v2 = comparison
    #v3 = day
    
    tmp <- data.frame(v0,v1,v2,v3) %>%
        group_by(v0,v2) %>%
        filter(v3==min(v3)) %>%
        group_by(v0)%>%
        mutate(n=n()) %>%
        filter(n>1)
    
    tmp %>% group_by(v2) %>% count() %>% print    
    # colnames(tmp) <- c("id","group1","group2")
    
    test <- wilcox.test(tmp$v1~tmp$v2,
                        paired=TRUE,
                        digits.rank=7
    )
    test$p.value
}



#compare distributions using the MW test
compareDays <- function(data){
    
    
    final <- data.frame()
    
    for (d in unique(data$control_day)){
        print(d)
        if (d != "Control / Day 2") {
            #find out which tests we have data on in both groups
            tmp1 <- data %>% filter(control_day ==d) %>%
                    # filter(!is.na(RAU)) %>%
                    group_by(test) %>%
                    count() 
            
            #f conduct the statistical test
            tmp2 <- data %>% filter(control_day %in% c("Control / Day 2",d)) %>% 
                filter(test %in% tmp1$test)%>%
                group_by(test) %>%
                summarize(
                    count=n(),
                    MW_pvalue=getMannWhitney(RAU,control_day),
                    WR_pvalue=getWilcoxonRS(id,RAU,control_day,day),
                    control_day=d
                )
        
           
            final <- bind_rows(final,tmp2)
        }
        
    }
    
    return(final)
    
}

# getTimeBasedWeights<- function(df,
#                                start_window = 5,
#                                end_window = 180
# ){
#     
#     inside_window <- df %>% 
#         filter(status=="Case") %>%
#         filter(day>=start_window) %>%
#         filter(day<=end_window)
#     
#     dist_day<- inside_window %>%
#         distinct(id,day,status) %>%
#         group_by(status,day) %>%
#         count()%>% 
#         ungroup() %>%
#         arrange(status,day) %>%
#         mutate(weight=NA) %>%
#         mutate(weight=ifelse(day==min(day),
#                              (day-start_window + (lead(day,1)-day)/2)/(end_window-start_window) * sum(n)/n,
#                              weight)) %>%
#         mutate(weight=ifelse(day==max(day),
#                              (end_window-day + (day-lag(day,1))/2)/(end_window-start_window) * sum(n)/n,
#                              weight)) %>%
#         mutate(weight=ifelse(day!=min(day) & day!=max(day),
#                              ((lead(day,1)-lag(day,1))/2)/(end_window-start_window) * sum(n)/n,
#                              weight)) %>%
#         select(-n)
#     
#     out <- df %>% left_join(dist_day,by = c("day", "status")) %>%
#         mutate(weight=ifelse(is.na(weight),1,weight))
#     
#     return(out)
#     
#     
# }


#calculate weights that account for both
# class imbalance
#distribution of infection times within class
getWeight<- function(data, #wide dataset
                     annual_incidence=0.1, #expected annual incidence (assuming constant hazard)
                     start_window = 5, #choose when the window begins
                     end_window, #choose when the window ends
                     decayed=540 #choose the timing when decay is complete
){
        
        if(end_window>1000) stop("Weights cannot be accurately guessed this far out")
        #ensure that decayed is defined later than the end_window
        if(decayed<end_window) decayed=end_window
        
        #categorize by day and window
        new_day_df <- data %>% 
                #ensure that contacts and later cases are in the same bucket
                mutate(new_day=ifelse(day_actual>=decayed | status=="Contact",
                                      decayed,day)) %>%
                mutate(window = ifelse(day_actual>=start_window & 
                                               day_actual <=end_window &
                                               status=="Case", 1,0 )) %>%
                mutate(window_type= ifelse(window==1,"During","After")) %>%
                mutate(window_type= ifelse(window==0 & new_day < start_window,
                                           "Before",window_type)) %>%
                #temporary fix for day 2 people in later buckets
                mutate(window=ifelse(day==2,0,window),
                       window_type=ifelse(day==2,"Before",window_type)
                       )
                
        
        
        
        #deal with the 1 nuissance 720 observation (day_actual >721) for
        # large infection windows (smush together with the 540 window)
        if(end_window>=721){
                new_day_df <- new_day_df %>%      
                                mutate(new_day=ifelse(new_day==720,540,new_day))
        } 
        
        
        #calculate daily hazard
        lambda <- annual_incidence/365.25
        
        #calculate the weights for each category
        weight_df <- new_day_df %>%
                group_by(window,window_type,new_day) %>%
                count() %>% 
                ungroup() %>%
                #weight to adjust for class imbalance
                mutate(class_weight= (window*(sum(n*window)) +
                                              (1-window)*(sum(n*(1-window))))/sum(n),
                       class_weight=0.5/class_weight
                ) %>%
                
                #make bounds for each category
                arrange(window_type,new_day) %>%
                group_by(window_type) %>%
                #find midpoints between categories
                mutate(lb=(new_day-lag(new_day,1))/2+lag(new_day,1)) %>%
                mutate(ub=lead(lb,1)) %>%
                
                #set bounds for data bordering thresholds
                mutate(lb=ifelse(window_type== "Before" & 
                                         new_day==min(new_day),
                                 0,lb)) %>%
                mutate(ub=ifelse(window_type== "Before" & 
                                         new_day==max(new_day),
                                 start_window,ub)) %>%
                mutate(lb=ifelse(window_type== "During" & 
                                         new_day==min(new_day),
                                 start_window,lb)) %>%
                mutate(ub=ifelse(window_type== "During" & 
                                         new_day==max(new_day),
                                 end_window,ub)) %>%
                mutate(lb=ifelse(window_type== "After" & 
                                         new_day==min(new_day),
                                 end_window,lb)) %>%
                mutate(ub=ifelse(window_type== "After" & 
                                         new_day==max(new_day),
                                 Inf,ub))%>%
                
                
                #calculate the probability
                group_by(window) %>%
                mutate(obs_prob=n/sum(n)) %>%
                mutate(exp_prob=pexp(ub,rate=lambda)-pexp(lb,rate=lambda)) %>%
                
                #weight to adjust for timing within class
                mutate(obs_rel=(obs_prob/sum(obs_prob))) %>%
                mutate(exp_rel=(exp_prob/sum(exp_prob))) %>%
                ungroup()
        
        #combine baseline with decayed individuals
        final_weight_df<-  bind_rows(
                                weight_df %>% filter(!((new_day==max(new_day) & window_type=="After") | 
                                                               window_type=="Before")),
                                weight_df %>% filter((new_day==max(new_day) & window_type=="After") | 
                                                              window_type=="Before") %>%
                                        mutate(obs_rel=sum(obs_rel),
                                               exp_rel=sum(exp_rel),
                                        )
                                )%>%
                mutate(time_weight= exp_rel/obs_rel) %>%
                #make final standardized weight
                mutate(final_weight=class_weight*time_weight)
        
        
        if(max(final_weight_df$final_weight)>50) warning(paste0("Some weights are above 10 (end_window=",end_window,")"))
        if(min(final_weight_df$final_weight)<0.05) warning(paste0("Some weights are below 0.05 (end_window=",end_window,")"))
        if(min(final_weight_df$final_weight)<0) stop(paste0("Some weights are below 0 (end_window=",end_window,")"))
        
        #return data with weights
        output <-  new_day_df %>%
                #attach  to original days
                left_join(final_weight_df,by = c("window", "window_type", "new_day")) %>%
                rename(weight=final_weight)
        
        return(output)
        
}



##' Uses a cross validation approach to get the opimimal cutpoints
##' then tests performance on hold-out. Only uses one observation per
##' person
##' @param dat
##' @param my_threshold threshold for postivie/negative (if known), otehrwise will compute based on weighted youden index
##' @param w_spec weight for youden index
##' @param biomarker
##' @param col_of_interest
##' @param nsims
##' @param print_me
##' @param training_prop
##' @return
##' @author
test_thresholds <- function(dat,
                            my_threshold=NA,
                            w_spec=.5,
                            # biomarker='vib',
                            # col_of_interest='inf_100',
                            nsims=1000,
                            print_me=FALSE,
                            training_prop=0.7){
    
    #maybe bring back if you bring in vibriocidal titers
    # if(biomarker=='vib'){
    #     ## get unique values of the titer
    #     unique_titers <- c(dat$vibinab,dat$vibogaw) %>% unique %>% sort
    #     
    #     ## add new column to dat for max_vib
    #     dat <- dat %>% dplyr::mutate(vib=pmax(vibogaw,vibinab))
    #     
    # } else {
    #     if(!biomarker %in% colnames(dat)) stop("biomarker is not a valid name of a column in dat. Only non-column name allowed is 'vib,' which is for max vibriocidal")
    #     
        unique_titers <- round(dat[,"avgMFI"]) %>% unlist %>% unique %>% sort
    # }
    
    unique_ids <- unique(dat$id) ## unique ids in the dataset
    obs_indices <- sapply(1:length(unique_ids),
                          function(x) which(dat$id == unique_ids[x])) ## list of observations per person
    
    ## number of people in training set
    n_train <- (length(unique_ids) * training_prop) %>% ceiling

    truths <- preds <- vector(mode='list', length=nsims)
    testing_perf <- matrix(ncol=5,nrow=nsims) %>% data.frame()
    colnames(testing_perf) <- c('threshold','sens_test','spec_test','sens_train','spec_train')
    
    

    for (i in 1:nsims){

        train_ids <- sample(length(unique_ids),n_train,replace=F) ## get ids of those who we want to sample for training
        test_ids <- setdiff(1:length(unique_ids),train_ids) ## then those for testing

        ## sample one time point for each person
        train_obs_inds <- sapply(train_ids,function(x) sample(obs_indices[[x]],1))
        test_obs_inds <- sapply(test_ids,function(x) sample(obs_indices[[x]],1))

        training_titers <- dat[train_obs_inds,"avgMFI"]  #biomarker]
        training_truths <- dat[train_obs_inds,"inf_wind"]

        testing_titers <- dat[test_obs_inds,"avgMFI"] #biomarker]
        testing_truths <- dat[test_obs_inds,"inf_wind"]
        

        ## now get optimal cut-point (if needed) from the
        ## training set
        if(is.na(my_threshold)){
            ## returns a matrix with columns for each possible cutpoint and rows for observations
            preds_training <- sapply(unique_titers,function(x) training_titers>=x)

            youden_max_id <-  which.max(apply(preds_training,2,function(x)
                youden_func(x,w_spec=w_spec,truths=training_truths)))

            sens_spec_thresh <- c(unique_titers[youden_max_id],
                                  apply(preds_training,2,function(x) sens_func(x,truths=training_truths))[youden_max_id],
                                  apply(preds_training,2,function(x) spec_func(x,truths=training_truths))[youden_max_id])

            my_thresh <- sens_spec_thresh[1]
            sens_train <- sens_spec_thresh[2]
            spec_train <- sens_spec_thresh[3]

        } else {

            my_thresh <- my_threshold
            preds_training <- training_titers>=my_thresh

            sens_train <- sens_func(preds_training,truths=training_truths)
            spec_train <- spec_func(preds_training,truths=training_truths)

        }

        ## now using the testing data, what is the sens and spec from these cut-offs?
        testing_preds <- data.frame(pred=testing_titers >= my_thresh,
                                    truth=testing_truths)

        testing_perf[i,] <- c(threshold=my_thresh,
                              sens_test=sens_func(testing_preds[,1],testing_preds[,2]),
                              spec_test=spec_func(testing_preds[,1],testing_preds[,2]),
                              sens_train=sens_train,
                              spec_train=spec_train)
    }

    ## don't want to keep this loaded
    # library(reshape2)
    thresh_df <- reshape2::melt(testing_perf,id.vars='threshold')
    # detach('package:reshape2')

    if (print_me){
        thresh_df %>% ggplot() +
            geom_histogram(aes(value,fill=factor(threshold)),position='identity') +
            facet_wrap(~variable) +
            labs(title=col_of_interest)  -> gg
        print(gg)
    }

    return(thresh_df)
}




## some helper functions
sens_func <- function(x,truths){
    mean((truths==x)[which(truths==1)])
}

spec_func <- function(x,truths){
    mean((truths==x)[which(truths==0)])
}

youden_func <- function(x,truths,w_spec=.5){
    (1-w_spec)*sens_func(x,truths) + w_spec*spec_func(x,truths) - 1
}






## cross validate single titer measuresment
cv_titer_repeated_measure_auc <- function(dat,k,col_of_interest){
    library(cvAUC)
    
    titers <- dat[,"avgMFI"] %>% unlist
    truths <- dat[,col_of_interest] %>% unlist
    ids <- dat[,'id'] %>% unlist
    
    folds <- cv_folds_rm(dat,k)
    preds <- vector("list",length=k)
    
    for (k in 1:length(folds)){
        perf <- ROCR::performance(ROCR::prediction(titers[-folds[[k]]],truths[-folds[[k]]]), "sens", "spec")
        df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@x.values[[1]], spec = perf@y.values[[1]])
        cutoff <- df[which.max(df$sens + df$spec), "cut"]
        preds[[k]] <- ifelse(titers[folds[[k]]]>=cutoff,1,0) %>% unlist %>% unname
    }
    
    out <- ci.pooled.cvAUC(predictions=preds %>% unlist,
                           labels=truths[unlist(folds)],
                           ids=ids[unlist(folds)],
                           folds=rep(1:k,map(folds,length)),
                           confidence=0.95)
    
    return(out)
}


## creates folds for cross validation stratified by outcome for balance
cv_folds <- function(Y, K){
    
    Y0 <- split(sample(which(Y==0)), rep(1:K, length=length(which(Y==0))))
    Y1 <- split(sample(which(Y==1)), rep(1:K, length=length(which(Y==1))))
    folds <- vector("list", length=K)
    
    for (k in seq(K)) {
        folds[[k]] <- c(Y0[[k]], Y1[[k]])
    }
    
    return(folds)
}

cv_folds_rm <- function(dat, k){
    
    ## get number of observations per person
    unique_ids <- unique(dat$id) ## unique ids in the dataset
    obs_indices <- lapply(1:length(unique_ids),function(x) which(dat$id == unique_ids[x])) ## list of observations per person
    ## get warnings for not having equal values in splits
    person_id_folds <- suppressWarnings(split(unique_ids,rep(1:k,length(unique_ids))))
    obs_id_folds <- lapply(person_id_folds,
                           function(my_fold){which(dat$id %in% my_fold)})
    
    return(obs_id_folds)
}


## cross validated AUC for RF models
cv_rf_model <- function(dat,my_formula,k=20,ntrees=2000,...){
    
    ## creating folds
    folds <- cv_folds_rm(dat,k)
    
    my_preds <- my_truths <- numeric()
    window <- str_extract(my_formula[2] %>% as.character,"inf_[0-9_]+")
    tmp_pred <- var_importance <- vector("list",length=k)
    
    for (f in 1:length(folds)){
        # cat(sprintf("Cross-validating AUC, fold %s \n",f))
    
        my_forest <- randomForest(my_formula,
                                  ntree=ntrees,
                                  data=dat[-c(folds[[f]]),],
                                  importance=TRUE,...)

        my_truths <- c(my_truths,dat[folds[[f]],window] %>% unlist)
        tmp_pred[[f]] <- predict(my_forest,newdata=dat[c(folds[[f]]),],predict.all=TRUE)
        var_importance[[f]] <- importance(my_forest)

        ## get probability from trees
        my_preds <- c(my_preds,
                      apply(tmp_pred[[f]][[2]],1,function(x) mean(x==1)) %>%
                          as.numeric
        )
    }

    out <- vector("list",length=4)

    out[[1]] <- ci.pooled.cvAUC(predictions=my_preds,
                                labels=my_truths,
                                ids=dat[unlist(folds),'id'] %>% unlist,
                                folds=rep(1:k,map(folds,length)),
                                confidence=0.95) %>%
        unlist %>% t %>%
        data.frame

    out[[1]] <- out[[1]] %>% mutate(time_window=window)

    out[[2]] <- bind_cols(truth=my_truths,pred=my_preds,fold=rep(1:k,map(folds,length)))

    out[[3]] <- tmp_pred

    out[[4]] <- var_importance

    names(out) <- c("perf_summary","outs","preds_full","var_imp")
    return(out)
}



make_formula <- function(y,x){
    covariates <- x %>% paste(collapse=" + ")
    paste0(y," ~ ",covariates) %>% as.formula()
}



##' runs the cforest function on a wide dataset
##' @formula formula used for fitting cforest model
##' @data wide dataset including weights with a w
##' @ntrees number of trees
##' @return list containing fit and conditional importance dataframe
cforest_vimp <- function(formula,
                         data,
                         ntrees,
                         weights,
                         importance,
                         # replacement=TRUE,
                         frac
                         ){
    
    #get the number of covariates in model
    m <- formula[3] %>% as.character() %>% str_split("\\+") %>% unlist() %>% length()
    
    # fit conditional random forest model with weights if supplied
    # make sure samples are selected using sample with replacement

        fit <- party::cforest(formula,
                              data,
                              weights=weights,
                              control=cforest_control(teststat = "quad", testtype = "Univ", 
                                                      mincriterion = 0, 
                                                      replace = FALSE,
                                                      fraction = frac,
                                                      ntree=ntrees,
                                                      mtry=ceiling(sqrt(m))) )

        
        
    # fit conditional random forest model without weights
    # samples are selected without replacement
    importance_df <- NULL
    if(importance){
        
        # calculate conditional importance with permutation measure
        vimp <- permimp::permimp(fit, conditional = TRUE, asParty = TRUE)
        
        #put conditional importance into a dataframe
        importance_df <- data.frame(Biomarker=names(vimp$values),
                                    Cond_Importance=as.numeric(vimp$values)) %>%
                                    arrange(desc(Cond_Importance)) %>%
                        mutate(isotype=str_extract(Biomarker,"Ig[A|G|M]")) %>%
                        mutate(antigen=ifelse(is.na(isotype),NA,str_remove(Biomarker,"Luminex_Ig[A|G|M]_"))) 
            
    }
    
    #put everything in a list to return
    final <- list(model = fit,
                  importance= importance_df)
    
    return(final)
}



##' runs k fold cross-validation and produces prediction
##' @formula formula used for fitting cforest model
##' @data wide dataset including weights with a w
##' @ntrees number of trees
##' @return dataframe of predictions
kfold_cv_cforest <- function(formula,
                             data,
                             k=20,
                             ntrees=500,
                             frac=0.632,
                             weighted){
    
    ## creating folds (should equal number of people)
    folds <- cv_folds_rm(data,k)
    
    # start objects to fill
    my_preds <- my_truths <- numeric()
    my_samples <- character()
    window <- str_extract(formula[2] %>% as.character,"inf_[0-9_]+")
    tmp_pred <- vector("list",length=k)
    tmp_truth <- vector("list",length=k)
    var_importance <- data.frame()
    
    # cat("0 of",k,"folds complete\n")
    for (f in 1:length(folds)){
        
        #limit to the fold data
        fold_data <- data[-c(folds[[f]]),] 
        
        #get weights
        w <- NULL
        if(weighted) {
            # find end window
            window_number <- window %>% 
                str_extract("[0-9]+") %>% 
                as.numeric()
            
            # get weights
            w <- getWeight(fold_data,end_window =window_number) %>%
                    pull(weight)
    
            # #fake weights
            # w <- fold_data %>%
            #     mutate(weight=ifelse(day==180,
            #                          2000,0.0000001)) %>%
            #     mutate(inf_180=ifelse(day_actual<=180 & day_actual>5 & status=="Case",
            #                           1,0
            #                           ))%>%
            #     group_by(inf_180) %>%
            #     mutate(weight=weight/sum(weight)) %>%
            #     ungroup() %>%
            #     mutate(weight=n()/2*weight) %>%
            #     pull(weight)

        }    
        
        #run the model
        this_run <- cforest_vimp(
            formula=formula,
            data=fold_data,
            ntrees=ntrees,
            importance=FALSE,
            weights =w,
            frac=frac
            )
        
        #store truth
        tmp_truth[[f]] <-data[c(folds[[f]]),window] %>% as.numeric() -1
        my_truths <- c(my_truths,data[c(folds[[f]]),window] %>% as.numeric() -1 )
        my_samples <- c(my_samples,data[c(folds[[f]]),"sample"])
        
        #get probability
        tmp_pred[[f]] <- predict(this_run$model,newdata=data[c(folds[[f]]),],
                                 type="prob")
        my_preds <- c(my_preds,
                      lapply(tmp_pred[[f]],function(x) x[2] )%>% unlist()
        )
        
        # cat(f,"of",k,"folds complete\n")
        
        
    }
    
    out <- list()

    #calculate cvAUC
    out$cvAUC <- cvAUC::ci.pooled.cvAUC(predictions=my_preds,
                                labels=my_truths,
                                ids=data[unlist(folds),'id'] %>% unlist,
                                folds=rep(1:k,map(folds,length)),
                                confidence=0.95
                                ) %>%
                                unlist %>% t %>%
                                data.frame  %>% mutate(time_window=window)
    
    #store predictions
    out$prediction_df <- bind_cols(sample=my_samples,
                                   truth=my_truths,
                                   pred=my_preds,
                                   fold=rep(1:k,map(folds,length)))
    out$prediction_list <- tmp_pred
    out$truth_list <- tmp_truth

    #store the folds
    out$folds <- folds
    
    return(out)
    
}




##' runs leave one out cross-validation and produces prediction
##' @formula formula used for fitting cforest model
##' @data wide dataset including weights with a w
##' @ntrees number of trees
##' @return dataframe of predictions
# loocv_cforest <- function(formula,data,ntrees){
#     
#     ## creating folds (should equal number of people)
#     k <- length(unique(data$id))
#     folds <- cv_folds_rm(data,k)
#     
#     my_preds <- my_truths <- numeric()
#     window <- str_extract(formula[2] %>% as.character,"inf_[0-9_]+")
#     tmp_pred <- vector("list",length=k)
#     
# 
#     for (f in 1:length(folds)){
#             #run the model
#             this_run <- cforest_vimp(
#                                     formula,
#                                     data[-c(folds[[f]]),],
#                                     ntrees)
#             
#         
#         my_truths <- c(my_truths,data[folds[[f]],window] %>%
#                            as.numeric() -1 )
#         tmp_pred[[f]] <- predict(this_run$model,newdata=data[c(folds[[f]]),],
#                                  type="prob"
#                                  )
# 
#         ## get probability
#         my_preds <- c(my_preds,
#                       lapply(tmp_pred[[f]],function(x) x[2] )%>% unlist()
#         )
#     }
#     
#     out <- vector("list",length=4)
#     
#     out <- list(my_truths,my_preds)
#     
#     
# }





##' Very specific helper function for making ROC curve and variable
##' importance plot from simulations
##' @title
##' @param preds predictions from cross validation
##' @param truths truths (binary)
##' @param imps importance objects
##' @param my_title plot title
##' @return multiplot object
##' @author
make_roc_varimp_plot <- function(preds,truths,imps,my_title,panel_label="",sub=FALSE,ribbon=TRUE,...){
    
    my_roc <- make_rocs(preds,truths,length(preds),title=my_title,ribbon=ribbon,...)
    
    importance_df <- do.call('rbind',imps) %>% as.data.frame() 
    
    importance_df$variable <- rownames(importance_df) %>% pretty_antibody
    
    importance_df <- importance_df %>%
        mutate(variable=str_remove(variable,"[.][0-9]+$"))
    
    levels_vars <- importance_df %>% 
        
        group_by(variable) %>%
        dplyr::summarize(med=median(MeanDecreaseAccuracy)) %>%
        select(med) %>% unlist() %>% order()
    
    
    
    
    
    impplot <- importance_df %>% 
        group_by(variable) %>%
        dplyr::summarize(MeanDecreaseAccuracy=median(MeanDecreaseAccuracy)) %>% 
        dplyr::mutate(variable1=factor(variable,
                                       levels=sort(unique(importance_df$variable))[levels_vars])) %>%
        ggplot() +
        geom_bar(aes(x=variable1,y=MeanDecreaseAccuracy,fill=variable),stat='identity',alpha=0.5) +
        coord_flip() + scale_fill_discrete(guide=F) + xlab('') + ylab('relative importance') +
        theme_classic() + theme(axis.title=element_text(size=8),
                                axis.text.x = element_text(size=6))
    
    if(sub){
        gg <- my_roc[[1]] + annotate("text",x=0,y=0,label=my_title,hjust=0) +
            annotate("text",x=0,y=1,label=panel_label,hjust=0) +
            annotation_custom(ggplotGrob(impplot), xmin=0.165, xmax=1, ymin=.1, ymax=.7)
        
        return(gg)
        
    } else{
        
        return(multiplot(impplot,my_roc[[1]],cols=2))
        
    }
}
## ROCs for cross validated results (multiple)
make_rocs <- function(cvData,k_folds, end_window=120, ribbon=FALSE,title="",annot_auc=TRUE){
        
        library(ROCR)
        
        outcome <- paste0("inf_",end_window)
        
        perf_roc <- sapply(1:k_folds,
                           function(x){
                                   foldData <- filter(cvData,fold==x)
                                   prediction(foldData$pred, 
                                              foldData[outcome] %>% unlist()) %>%
                                           performance(.,"tpr","fpr")        
                                   
                           },
                simplify = F)
        
        
      
        perf_auc <- sapply(1:k_folds,
                           function(x){
                                   foldData <- filter(cvData,fold==x)
                                   prediction(foldData$pred, 
                                              foldData[outcome] %>% unlist()) %>%
                                           performance(.,"auc")      
                                   
                           })

        auc <- sapply(perf_auc,function(x) x@y.values[[1]])
        
        # perf_roc <- sapply(1:k_folds,function(x) {
        #         prediction(votes[[x]],truths[[x]]) %>% performance(.,"tpr","fpr")},
        #         simplify = F)
        
        df <- data.frame(fold = rep(1:k_folds,sapply(perf_roc,function(x) x@x.values %>% unlist %>% length)),
                         fpr=sapply(perf_roc,function(x) x@x.values[[1]],simplify=FALSE) %>% unlist,
                         tpr=sapply(perf_roc,function(x) x@y.values[[1]],simplify=FALSE) %>% unlist)%>% dplyr::mutate(fpr_r = round(fpr,2),tpr_r=round(tpr,2))%>%
                mutate(
                        sens=tpr,
                        spec=1-fpr,
                        youden=sens+spec-1
                )
        
        
        # cat(sprintf("AUC Range: %.2f-%.2f \n",range(auc)[1],range(auc)[2]))
        # cat(sprintf("AUC IQR: %.2f-%.2f \n",quantile(auc,.25),quantile(auc,.75)))
        # cat(sprintf("AUC Median and Mean: %.2f, %.2f \n",median(auc),mean(auc)))
        
        # df <- df %>% dplyr::mutate(fpr_r = round(fpr,2),tpr_r=round(tpr,2))%>%
        #         mutate(
        #                 sens=tpr,
        #                 spec=1-fpr,
        #                 youden=sens+spec-1
        #         )
        

        ci_ribbon <- df %>% group_by(fpr_r) %>% dplyr::summarize(med=median(tpr),low_ci=quantile(tpr,.025),upper_ci=quantile(tpr,.975))

        # get youden index for each curve
        youden <-df  %>%
                group_by(fold) %>%
                arrange(youden)%>%
                filter(youden==max(youden))  %>%
                filter(row_number()==ceiling(n()/2))


        # take the mean of the youden indices
        mean_ss <- youden %>% ungroup()%>%
                summarize(tpr=mean(tpr),
                          fpr=mean(fpr))

        if(ribbon){
                #calculat
                ggplot(data=ci_ribbon) +
                        geom_line(aes(x=fpr_r,y=med),col='royalblue',lwd=1.2)+
                        geom_ribbon(data=ci_ribbon,aes(x=fpr_r,ymin=low_ci,ymax=upper_ci),alpha=.5,fill='royalblue') +
                        xlab('False Positive Rate') + ylab('True Positive Rate')+
                        theme_bw()-> gg

        } else {

                ggplot(data=df,aes(x=fpr,y=tpr)) +
                        geom_line(aes(group=fold),
                                  alpha=0.1,
                                  color=AddAlpha('#2c7fb8',0.9)
                        ) +
                        geom_point(data=youden,alpha=0.5,col="#2c7fb8")+
                        geom_point(data=mean_ss,col="red",shape=8)+
                        scale_x_continuous('False Positive Rate',
                                           breaks=c(0,0.2,0.4,0.6,0.8,1)) +
                        scale_y_continuous('True Positive Rate',
                                           breaks=c(0,0.2,0.4,0.6,0.8,1)) +
                        theme_cowplot()+
                        geom_abline(lty=2)+
                        annotate("text", x = 0.5, y = 0.15, size=2,label = title) -> gg

                if(annot_auc){
                        # gg <- gg + annotate("text", x = 0.5, y = 0.1,size=4, label = sprintf("AUC IQR: %.2f-%.2f \n",quantile(auc,.25),quantile(auc,.75)))
                        # gg <- gg + annotate("text", x = 0.5, y = 0.1,size=4, label = paste0("AUC Range=",round(min(auc),3),"-",round(max(auc),3)))
                        gg <- gg + annotate("text", x = 0.5, y = 0.07,size=4,
                                            label = sprintf("Mean Sensitivity: %.2f\n Mean Specificity: %.2f\n",mean_ss$tpr,1-mean_ss$fpr))
                }

        }

        return(list(gg,df,auc))
        
}


# ## ROCs for cross validated results (multiple)
# make_rocs <- function(preds,truths,k_folds,ribbon=FALSE,title="",annot_auc=TRUE){
#     
#     library(ROCR)
#     
#     votes <- lapply(preds, FUN=function(pred) {
#         lapply(pred,FUN=function(x) x[2]) %>% unlist()        
#     })
#     
#     perf_auc <- sapply(1:k_folds,function(x)
#         prediction(votes[[x]],truths[[x]]) %>% performance(.,"auc"))
#     
#     auc <- sapply(perf_auc,function(x) x@y.values[[1]])
#     
#     perf_roc <- sapply(1:k_folds,function(x) {
#         prediction(votes[[x]],truths[[x]]) %>% performance(.,"tpr","fpr")},
#         simplify = F)
#     
#     df <- data.frame(fold = rep(1:k_folds,sapply(perf_roc,function(x) x@x.values %>% unlist %>% length)),
#                      fpr=sapply(perf_roc,function(x) x@x.values[[1]],simplify=FALSE) %>% unlist,
#                      tpr=sapply(perf_roc,function(x) x@y.values[[1]],simplify=FALSE) %>% unlist)
#     
#     
#     # cat(sprintf("AUC Range: %.2f-%.2f \n",range(auc)[1],range(auc)[2]))
#     # cat(sprintf("AUC IQR: %.2f-%.2f \n",quantile(auc,.25),quantile(auc,.75)))
#     # cat(sprintf("AUC Median and Mean: %.2f, %.2f \n",median(auc),mean(auc)))
#     
#     df <- df %>% dplyr::mutate(fpr_r = round(fpr,2),tpr_r=round(tpr,2))%>%
#         mutate(
#             sens=tpr,
#             spec=1-fpr,
#             youden=sens+spec-1
#         )
#     
#     ci_ribbon <- df %>% group_by(fpr_r) %>% dplyr::summarize(med=median(tpr),low_ci=quantile(tpr,.025),upper_ci=quantile(tpr,.975))
#     
#     # get youden index for each curve
#     youden <-df  %>%
#         group_by(fold) %>%
#         arrange(youden)%>%
#         filter(youden==max(youden))  %>%
#         filter(row_number()==ceiling(n()/2))        
#     
#     
#     # take the mean of the youden indices
#     mean_ss <- youden %>% ungroup()%>%
#         summarize(tpr=mean(tpr),
#                   fpr=mean(fpr))
#     
#     if(ribbon){
#         #calculat
#         ggplot(data=ci_ribbon) +
#             geom_line(aes(x=fpr_r,y=med),col='royalblue',lwd=1.2)+
#             geom_ribbon(data=ci_ribbon,aes(x=fpr_r,ymin=low_ci,ymax=upper_ci),alpha=.5,fill='royalblue') +
#             xlab('False Positive Rate') + ylab('True Positive Rate')+
#             theme_bw()-> gg
#         
#     } else {
#         
#         ggplot(data=df,aes(x=fpr,y=tpr)) +
#             geom_line(aes(group=fold),
#                       alpha=0.1,
#                       color=AddAlpha('#2c7fb8',0.9)
#             ) +
#             geom_point(data=youden,alpha=0.5,col="#2c7fb8")+
#             geom_point(data=mean_ss,col="red",shape=8)+
#             scale_x_continuous('False Positive Rate',
#                                breaks=c(0,0.2,0.4,0.6,0.8,1)) + 
#             scale_y_continuous('True Positive Rate',
#                                breaks=c(0,0.2,0.4,0.6,0.8,1)) +
#             theme_cowplot()+
#             geom_abline(lty=2)+
#             annotate("text", x = 0.5, y = 0.15, size=2,label = title) -> gg
#         
#         if(annot_auc){
#             # gg <- gg + annotate("text", x = 0.5, y = 0.1,size=4, label = sprintf("AUC IQR: %.2f-%.2f \n",quantile(auc,.25),quantile(auc,.75)))
#             # gg <- gg + annotate("text", x = 0.5, y = 0.1,size=4, label = paste0("AUC Range=",round(min(auc),3),"-",round(max(auc),3)))   
#             gg <- gg + annotate("text", x = 0.5, y = 0.07,size=4, 
#                                 label = sprintf("Mean Sensitivity: %.2f\n Mean Specificity: %.2f\n",mean_ss$tpr,1-mean_ss$fpr))
#         }
#         
#     }
#     
#     return(list(gg,df,auc))
#     
# }


#old version of make_rocs
# ## ROCs for cross validated results (multiple)
# make_rocs <- function(preds,truths,k_folds,ribbon=TRUE,title="",annot_auc=TRUE){
#     votes <- sapply(preds,function(pred)
#         apply(pred$individual,1,function(x) mean(x==1)),simplify=FALSE)
#     
#     perf_auc <- sapply(1:k_folds,function(x)
#         prediction(votes[[x]],truths[[x]]) %>% performance(.,"auc"))
#     
#     auc <- sapply(perf_auc,function(x) x@y.values[[1]])
#     
#     perf_roc <- sapply(1:k_folds,function(x) {
#         prediction(votes[[x]],truths[[x]]) %>% performance(.,"tpr","fpr")},
#         simplify = F)
#     
#     df <- data.frame(fold = rep(1:k_folds,sapply(perf_roc,function(x) x@x.values %>% unlist %>% length)),
#                      fpr=sapply(perf_roc,function(x) x@x.values[[1]],simplify=FALSE) %>% unlist,
#                      tpr=sapply(perf_roc,function(x) x@y.values[[1]],simplify=FALSE) %>% unlist)
#     
#     
#     # cat(sprintf("AUC Range: %.2f-%.2f \n",range(auc)[1],range(auc)[2]))
#     # cat(sprintf("AUC IQR: %.2f-%.2f \n",quantile(auc,.25),quantile(auc,.75)))
#     # cat(sprintf("AUC Median and Mean: %.2f, %.2f \n",median(auc),mean(auc)))
#     
#     df <- df %>% dplyr::mutate(fpr_r = round(fpr,2),tpr_r=round(tpr,2))
#     ci_ribbon <- df %>% group_by(fpr_r) %>% dplyr::summarize(med=median(tpr),low_ci=quantile(tpr,.025),upper_ci=quantile(tpr,.975))
#     
#     if(ribbon){
#         
#         ggplot(data=ci_ribbon) +
#             geom_line(aes(x=fpr_r,y=med),col='royalblue',lwd=1.2)+
#             geom_ribbon(data=ci_ribbon,aes(x=fpr_r,ymin=low_ci,ymax=upper_ci),alpha=.5,fill='royalblue') +
#             xlab('False Positive Rate') + ylab('True Positive Rate') -> gg
#         
#     } else {
#         
#         ggplot(data=df,aes(x=fpr,y=tpr)) +
#             geom_line(aes(group=fold),color=AddAlpha('#2c7fb8',0.9)) +
#             xlab('False Positive Rate') + ylab('True Positive Rate') +
#             annotate("text", x = 0.5, y = 0.15, size=2,label = title) -> gg
#         
#         if(annot_auc) gg <- gg + annotate("text", x = 0.5, y = 0.1,size=2, label = paste0("AUC Range=",round(min(auc),3),"-",round(max(auc),3)))
#         
#     }
#     
#     return(list(gg,df,auc))
#     
# }

##' Adds alpha to a set of colors
##' @title
##' @param COLORS
##' @param ALPHA
##' @return
AddAlpha <- function(COLORS, ALPHA){
    if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
    RGB <- col2rgb(COLORS, alpha=TRUE)
    RGB[4,] <- round(RGB[4,]*ALPHA)
    NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
    return(NEW.COLORS)
}


##' helper function to make pretty names for covariates
##' either give a dataframe so it works on column names, otherwise
##' give a vector and it replaces names lnly
##' @param x
##' @param dplyr::rename_cols
##' @return
##' @author asa
pretty_antibody <- function(x){
    
    if(!is.null(ncol(x)) && ncol(x)>1){
        x <- x %>%
            dplyr::rename('plpsa_s',"anti-LPS IgA") %>%
            dplyr::rename('ptctxa_s',"anti-CTB IgA") %>%
            dplyr::rename('vibinab' ,"vibriocidal (Inaba)") %>%
            dplyr::rename('vibogaw' ,"vibriocidal (Ogawa)") %>%
            dplyr::rename('pglsp_s' ,"anti-LPS IgG") %>%
            dplyr::rename('pgctxb_s', "anti-CTB IgG") %>%
            dplyr::rename('plpsm_s',"anti-LPS IgM") %>%
            dplyr::rename('pmctxb_s',"anti-CTB IgM") %>%
            dplyr::rename('lbday2',"Days from Symptom Onset") %>%
            dplyr::rename('wt' ,"weight") %>%
            dplyr::rename('ht' ,"height") %>%
            dplyr::rename('o_group' ,"O blood group")
        
    } else {
        
        x <- x %>%
            str_replace_all('plpsa_s',"anti-LPS IgA") %>%
            str_replace_all('ptctxa_s',"anti-CTB IgA") %>%
            str_replace_all('vibinab' ,"vibriocidal (Inaba)") %>%
            str_replace_all('vibogaw' ,"vibriocidal (Ogawa)") %>%
            str_replace_all('pglsp_s' ,"anti-LPS IgG") %>%
            str_replace_all('pgctxb_s', "anti-CTB IgG") %>%
            str_replace_all('plpsm_s' ,"anti-LPS IgM") %>%
            str_replace_all('pmctxb_s' ,"anti-CTB IgM") %>%
            str_replace_all('lbday2' ,"Days from Symptom Onset") %>%
            str_replace_all('^wt$' ,"weight") %>%
            str_replace_all('^ht$' ,"height") %>%
            str_replace_all('o_group' ,"O blood group")
        
    }
    return(x)
}


## to make a pdf
to.pdf <- function(expr, filename, ...)
    to.dev(expr, pdf, filename, ...)


##' wrapper for sending output from function to device
##' (note: this is used by to pdf to create pdfs)
##' @param expr
##' @param dev
##' @param filename
##' @param ...
##' @param verbose
##' @return
##' @author http://nicercode.github.io/blog/2013-07-09-figure-functions/
to.dev <- function(expr, dev, filename, ..., verbose=TRUE) {
    if ( verbose )
        cat(sprintf("Creating %s\n", filename))
    dev(filename, ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}


##' Very specific helper function for making ROC curve and variable
##' importance plot from simulations
##' @title
##' @param preds predictions from cross validation
##' @param truths truths (binary)
##' @param imps importance objects
##' @param my_title plot title
##' @return multiplot object
##' @author
make_roc_varimp_df <- function(preds,truths,imps,my_title,panel_label="",sub=FALSE,ribbon=TRUE,...){
    
    my_roc <- make_rocs(preds,truths,length(preds),title=my_title,ribbon=ribbon,...)
    
    importance_df <- do.call('rbind',imps) %>% as.data.frame() 
    
    importance_df$variable <- rownames(importance_df) %>% pretty_antibody
    
    importance_df <- importance_df %>%
        mutate(variable=str_remove(variable,"[.][0-9]+$")) 
        
    
    levels_vars <- importance_df %>% 
        
        group_by(variable) %>%
        dplyr::summarize(med=median(MeanDecreaseAccuracy)) %>%
        select(med) %>% unlist() %>% order()
    
    
    imp_df<- importance_df %>% 
        group_by(variable) %>%
        dplyr::summarize(MeanDecreaseAccuracy=median(MeanDecreaseAccuracy)) %>% 
        dplyr::mutate(variable1=factor(variable,
                                       levels=sort(unique(importance_df$variable))[levels_vars])) 
    
    
    return(imp_df)
}



## cross validated AUC for RF models
cv_caret_model <- function(dat, #wide dataset
                           inputs, #predictors
                           outcomes, #infection window
                           model="gbm", #character describing which caret model
                           k=5, # folds
                           ...
                           
                           ){
    
    ## creating folds
    folds <- cv_folds_rm(dat,k)
    
    my_preds <- my_truths <- numeric()
    window <- str_extract(outcomes,"inf_[0-9_]+")
    tmp_pred <- var_importance <- vector("list",length=k)
    
    for (f in 1:length(folds)){
        # cat(sprintf("Cross-validating AUC, fold %s \n",f))
        
        my_ML <- train(dat[-c(folds[[f]]),inputs] %>% as.data.frame(),
                       dat[-c(folds[[f]]),outcomes]  %>% unlist()%>% as.factor(),
              method = model, 
              weights = dat[-c(folds[[f]]),"weigh"]  %>% unlist()
              )
        
        my_truths <- c(my_truths,dat[folds[[f]],outcomes] %>% unlist)
        tmp_pred[[f]] <- predict(my_ML,newdata=dat[c(folds[[f]]),inputs] %>%
                                     as.data.frame(),
                                 type = "prob") #%>%
                        # as.character() %>% as.numeric()
        # var_importance[[f]] <- importance(my_forest)

        ## get probability from trees
        my_preds <- c(my_preds,tmp_pred[[f]][,2])
        
     }
    
    out <- vector("list",length=4)

    out[[1]] <- ci.pooled.cvAUC(predictions=my_preds,
                                labels=my_truths,
                                ids=dat[unlist(folds),'id'] %>% unlist,
                                folds=rep(1:k,map(folds,length)),
                                confidence=0.95) %>%
        unlist %>% t %>%
        data.frame %>% mutate(method=model)

    out[[1]] <- out[[1]] %>% mutate(time_window=window)

    out[[2]] <- bind_cols(truth=my_truths,pred=my_preds,fold=rep(1:k,map(folds,length)))

    out[[3]] <- tmp_pred

    # out[[4]] <- var_importance

    names(out) <- c("perf_summary","outs","preds_full","var_imp")
    return(out)
}

## select stan model to avoid recompiling every time
select_univariate_model <- function(model_type, covariates =NULL){
    
    model_choice <- data.frame(decay=c("exponential","exponential","biphasic","biphasic"),
                               co_variates=c(FALSE,TRUE,FALSE,TRUE),
                               path=c(
                                   "source/final_code/luminex_recommendation/stan/uni_exp_decay_muD.stan",
                                   "source/final_code/luminex_recommendation/stan/uni_exp_decay_muD_covariate.stan",
                                   "source/final_code/luminex_recommendation/stan/uni_biphas_decay_muD.stan",
                                   NA)) %>%
        filter(decay==model_type) %>%
        filter(co_variates == !is.null(covariates))
    
    model <- stan_model(model_choice$path)
    
    return(model)
    
    
}



#takes a wide dataset in and puts it into a format
# acceptable by stan model
stanfit_decay_univariate <- function(data,
                                     cens_df,
                                     marker,
                                     delay=5, #if d>=delay, boost
                                     covariates =NULL, # leave null if you want no covariates
                                     # model_type,
                                     iterations=2000,
                                     # luminex_data_type = "RAU",
                                     model,
                                     chains=4,
                                     warmup=1000
                                     ){
        time <- Sys.time()
        
        #remove rows with missing data for the marker
        missing <- is.na(data[marker])
        data <- data[!missing,] %>%
                #create novel id
                mutate(id_new=factor(id)) %>%
                arrange(id_new,day_actual) %>%
                mutate(id_new=as.numeric(id_new))
        
        cens_df <- cens_df[!missing,] %>%
                #create novel id
                mutate(id_new=factor(id)) %>%
                arrange(id_new,day_actual) %>%
                mutate(id_new=as.numeric(id_new))
        
        
        #check to make sure they are aligned properly
        if(any(data$sample!=cens_df$sample))stop("data and censored data misaligned")
        
        cov_df <- data[c("id_new",covariates)] %>% distinct()
        
        #create list for stan to take in 
        stan_input <- list(
                
                #actual data
                L = length(unique(data$id)), # number of individuals
                N = nrow(data),
                t = data$day_actual, #- avg_t,        #timing of values
                y = data[marker] %>% unlist(),   # select marker for analysis
                x =  cov_df[covariates] %>% unlist(),
                cens = cens_df[marker] %>% unlist(),
                
                id_new =  data$id_new, 
                delay = delay, #- avg_t, # delay for jump in value 
                
                
                #value for priors
                #if exponential
                prior_mean_delta = 0.02,
                prior_sd_delta = .05,
                
                #if biphasic
                prior_mean_theta = c(0.02,0.02),
                prior_cov_theta = matrix(c(.1,-0.01,-0.01,.1),nrow=2),
                
                #store individual level data with object
                id_df = select(data,id,id_new,sample,day_actual) %>%
                        mutate(marker=marker),
                #store stull dataset
                full_data = data
                
        )
        
        #identify measurement type and priors
        if(str_detect(marker,"Luminex|RAU|NetMFI")) {
                stan_input$measure <- 1
                # if(luminex_data_type=="RAU"){
                #         #avg boost and baseline values
                #         stan_input$prior_mean_params_mu <- c(-4,2)
                #         stan_input$prior_cov_params_mu  <- matrix(c(10,-1,-1,10),nrow=2)
                # }
                # if(luminex_data_type=="netMFI"){
                #         #avg boost and baseline values
                #         stan_input$prior_mean_params_mu <- c(2,2)
                #         stan_input$prior_cov_params_mu  <- matrix(c(10,-1,-1,10),nrow=2)
                # }


        }
        if(str_detect(marker,"Vibriocidal")){
                stan_input$measure <- 2
                #avg boost and baseline values
                # stan_input$prior_mean_params_mu <- c(2,4)
                # stan_input$prior_cov_params_mu  <- matrix(c(10,-1,-1,10),nrow=2)
        }
        if(str_detect(marker,"ELISA")) {
                stan_input$measure <- 3
                #avg boost and baseline values
                # stan_input$prior_mean_params_mu <- c(-1,1)
                # stan_input$prior_cov_params_mu  <- matrix(c(10,-1,-1,10),nrow=2)
        }
        
        
        #select stan file for model    
        # model_choice <- data.frame(decay=c("exponential","exponential","biphasic","biphasic"),
        #                             co_variates=c(FALSE,TRUE,FALSE,TRUE),
        #                             path=c(
        #                                 "source/final_code/luminex_recommendation/stan/uni_exp_decay.stan",
        #                                 "source/final_code/luminex_recommendation/stan/uni_exp_decay_covariate.stan",
        #                                 "source/final_code/luminex_recommendation/stan/uni_biphas_decay.stan",
        #                                     NA)) %>%
        #                 filter(decay==model_type) %>%
        #                 filter(co_variates == !is.null(covariates))
        
        uni_stan_fit<- sampling(model,
                                data = stan_input,
                                chains=chains, iter=iterations, warmup = warmup
        )
        
        # uni_stan_fit<- stan(file = model_choice$path,
        #                        data = stan_input,
        #                     chains=4, iter=iterations, warmup = 500
        #   )
        
        #store summary statistics
        summary_df <- summary(uni_stan_fit)$summary %>% data.frame() %>%
                mutate(parameter=rownames(.)) %>%
                mutate(marker=marker) %>%
                
                #match with original id
                mutate(id_new=ifelse(str_detect(parameter,"ind"),
                                     str_extract(parameter,"\\[(.*?)\\]"),
                                     NA))%>%
                mutate(id_new=ifelse(str_detect(id_new,"\\[(.*?),"),
                                     str_extract(id_new,"\\[(.*?),"),
                                     NA))%>%
                mutate(id_new=str_remove_all(id_new,"[\\[\\],]")) %>%
                mutate(id_new=as.numeric(id_new)) %>%
                left_join(distinct(data,id_new,id)) %>%
                
                #match predicted days
                mutate(day_id=ifelse(str_detect(parameter,"trajectory"),
                                     str_extract(parameter,"\\[(.*?)\\]"),
                                     NA)) %>%
                mutate(day_id=ifelse(str_detect(day_id,",(.*?)\\]"),
                                     str_extract(parameter,",(.*?)\\]"),
                                     day_id)) %>%
                mutate(day_id=str_remove_all(day_id,"[\\[\\],]")) 
        
        #get the output together
        out <- list(data =stan_input,
                    fit  =uni_stan_fit,
                    summary= summary_df
        )
        
        print(Sys.time()-time)
        
        return(out)
        
        
}


#Old version doesnt account for net MFI vs. RAU
# #takes a wide dataset in and puts it into a format
# # acceptable by stan model
# stanfit_decay_univariate <- function(data,
#                                      cens_df,
#                                      marker,
#                                      delay=5, #if d>=delay, boost
#                                      covariates =NULL, # leave null if you want no covariates
#                                      # model_type,
#                                      iterations=2000,
#                                      model){
#     time <- Sys.time()
#     
#     #remove rows with missing data for the marker
#     missing <- is.na(data[marker])
#     data <- data[!missing,] %>%
#         #create novel id
#         mutate(id_new=factor(id)) %>%
#         arrange(id_new,day_actual) %>%
#         mutate(id_new=as.numeric(id_new))
#     
#     cens_df <- cens_df[!missing,] %>%
#         #create novel id
#         mutate(id_new=factor(id)) %>%
#         arrange(id_new,day_actual) %>%
#         mutate(id_new=as.numeric(id_new))
#     
#     
#     #check to make sure they are aligned properly
#     if(any(data$sample!=cens_df$sample))stop("data and censored data misaligned")
#     
#     
#     #create list for stan to take in 
#     stan_input <- list(
#         
#         #actual data
#         L = length(unique(data$id)), # number of individuals
#         N = nrow(data),
#         t = data$day_actual, #- avg_t,        #timing of values
#         y = data[marker] %>% unlist(),   # select marker for analysis
#         x = data[covariates] %>% unlist(),
#         cens = cens_df[marker] %>% unlist(),
#         
#         id_new =  data$id_new, 
#         D = delay, #- avg_t, # delay for jump in value 
#         
#         
#         #value for priors
#         #if exponential
#         prior_mean_delta = 0.02,
#         prior_sd_delta = .05,
#         
#         #if biphasic
#         prior_mean_theta = c(0.02,0.02),
#         prior_cov_theta = matrix(c(.1,-0.01,-0.01,.1),nrow=2),
#         
#         #store individual level data with object
#         id_df = select(data,id,id_new,sample,day_actual) %>%
#             mutate(marker=marker)
#         
#     )
#     
#     #identify measurement type and priors
#     if(str_detect(marker,"Luminex")) {
#         stan_input$measure <- 1
#         #avg boost and baseline values
#         stan_input$prior_mean_params_mu <- c(-4,2)
#         stan_input$prior_cov_params_mu  <- matrix(c(10,-1,-1,10),nrow=2)
#     }
#     if(str_detect(marker,"Vibriocidal")){
#         stan_input$measure <- 2
#         #avg boost and baseline values
#         stan_input$prior_mean_params_mu <- c(2,4)
#         stan_input$prior_cov_params_mu  <- matrix(c(10,-1,-1,10),nrow=2)
#     }
#     if(str_detect(marker,"ELISA")) {
#         stan_input$measure <- 3
#         #avg boost and baseline values
#         stan_input$prior_mean_params_mu <- c(-1,1)
#         stan_input$prior_cov_params_mu  <- matrix(c(10,-1,-1,10),nrow=2)
#     }
#             
#     
#     #select stan file for model    
#     # model_choice <- data.frame(decay=c("exponential","exponential","biphasic","biphasic"),
#     #                             co_variates=c(FALSE,TRUE,FALSE,TRUE),
#     #                             path=c(
#     #                                 "source/final_code/luminex_recommendation/stan/uni_exp_decay.stan",
#     #                                 "source/final_code/luminex_recommendation/stan/uni_exp_decay_covariate.stan",
#     #                                 "source/final_code/luminex_recommendation/stan/uni_biphas_decay.stan",
#     #                                     NA)) %>%
#     #                 filter(decay==model_type) %>%
#     #                 filter(co_variates == !is.null(covariates))
#     
#     uni_stan_fit<- sampling(model,
#                             data = stan_input,
#                             chains=4, iter=iterations, warmup = 500
#     )
#     
#     # uni_stan_fit<- stan(file = model_choice$path,
#     #                        data = stan_input,
#     #                     chains=4, iter=iterations, warmup = 500
#     #   )
#     
#     #store summary statistics
#     summary_df <- summary(uni_stan_fit)$summary %>% data.frame() %>%
#         mutate(parameter=rownames(.)) %>%
#         mutate(marker=marker) %>%
#         
#         #match with original id
#         mutate(id_new=ifelse(str_detect(parameter,"ind"),
#                              str_extract(parameter,"\\[(.*?)\\]"),
#                              NA))%>%
#         mutate(id_new=ifelse(str_detect(id_new,"\\[(.*?),"),
#                              str_extract(id_new,"\\[(.*?),"),
#                              NA))%>%
#         mutate(id_new=str_remove_all(id_new,"[\\[\\],]")) %>%
#         mutate(id_new=as.numeric(id_new)) %>%
#         left_join(distinct(data,id_new,id)) %>%
#         
#         #match predicted days
#         mutate(day_id=ifelse(str_detect(parameter,"trajectory"),
#                              str_extract(parameter,"\\[(.*?)\\]"),
#                              NA)) %>%
#         mutate(day_id=ifelse(str_detect(day_id,",(.*?)\\]"),
#                              str_extract(parameter,",(.*?)\\]"),
#                              day_id)) %>%
#         mutate(day_id=str_remove_all(day_id,"[\\[\\],]")) 
#     
#     #get the output together
#     out <- list(data =stan_input,
#                 fit  =uni_stan_fit,
#                 summary= summary_df
#     )
#     
#     print(Sys.time()-time)
#     
#     return(out)
#     
#     
# }




#takes a wide dataset in and puts it into a format
# acceptable by stan model and then fits a multivariate model
stanfit_decay_multivariate <- function(data,markers,
                                       delay, #if d>=delay, boost
                                       iterations=2000
){
    time <- Sys.time()
    
    #remove rows with missing data for the marker
    missing <- rowSums(is.na(wide_analysis[,markers])) !=0
    data <- data[!missing,] %>%
        #create novel id
        mutate(id_new=factor(id)) %>%
        arrange(id_new,day_actual) %>%
        mutate(id_new=as.numeric(id_new))
    
    #calculate average to center the timing of values (necessary?)
    # avg_t <- mean(data$day_actual)
    
    #create list for stan to take in 
    stan_input <- list(
        
        K = length(markers),
        
        #actual data
        L = length(unique(data$id)), # number of individuals
        N= nrow(data),
        t = data$day_actual,       #timing of values
        y = data[markers],   # select marker for analysis
        
        id_new =  data$id_new, 
        D = delay, 
        
        
        #value for priors
        #avg boost and baseline values
        prior_mean_params_mu = c(-4,2),
        prior_cov_params_mu = matrix(c(10,-1,-1,10),nrow=2),
        
        #if exponential
        prior_mean_delta = 0.02,
        prior_sd_delta = .1,
        
        
        #store individual level data with object
        id_df = select(data,id,id_new,sample,day_actual) %>%
            mutate(markers=paste0(markers,collapse=" + "))
        
    )
    
    #fit stan model
    multi_stan_fit<- stan(file = "source/final_code/epidemic_reconstruction/stan/multivariate_decay_model.stan",
                          data = stan_input,
                          chains=4, iter=iterations, warmup = 500
    )
    
    #store summary statistics
    summary_df <- summary(multi_stan_fit)$summary %>% data.frame() %>%
        mutate(parameter=rownames(.)) %>%
        mutate(markers=paste0(markers,collapse=" + ")) #%>%
    
    # #match with original id
    # mutate(id_new=ifelse(str_detect(parameter,"ind"),
    #                      str_extract(parameter,"\\[(.*?)\\]"),
    #                      NA))%>%
    # mutate(id_new=ifelse(str_detect(id_new,"\\[(.*?),"),
    #                      str_extract(id_new,"\\[(.*?),"),
    #                      NA))%>%
    # mutate(id_new=str_remove_all(id_new,"[\\[\\],]")) %>%
    # mutate(id_new=as.numeric(id_new)) %>%
    # left_join(distinct(data,id_new,id)) %>%
    # 
    # #match predicted days
    # mutate(day_id=ifelse(str_detect(parameter,"trajectory"),
    #                      str_extract(parameter,"\\[(.*?)\\]"),
    #                      NA)) %>%
    # mutate(day_id=ifelse(str_detect(day_id,",(.*?)\\]"),
    #                      str_extract(parameter,",(.*?)\\]"),
    #                      day_id)) %>%
    # mutate(day_id=str_remove_all(day_id,"[\\[\\],]")) %>%
    # left_join(data.frame(day_id=as.character(1:length(predict_days)),
    #                      predict_day=predict_days
    # ))
    
    #get the output together
    out <- list(data =stan_input,
                fit  =multi_stan_fit,
                summary= summary_df
    )
    
    print(Sys.time()-time)
    
    return(out)
    
    
}

# get the performance of a ranger model
# determines if test set is seropositive based on
# Youden index
# includes sens, spec, AUC
getPerf_ranger<- function(data, m,  truth_variable){

        # get truths predictions
        truths <- data[truth_variable] %>% unlist() %>%
                as.character() %>% as.numeric()
        pred_obj <- predict(m,data=data,predict.all=TRUE)
        preds <- pred_obj[[1]] %>% rowMeans()-1


        # get performance
        output_df <- NULL
        ss <- NULL
        fr <- NULL
        auc <- NULL

        if(!(all(truths==1)|all(truths==0))){ #if some recent or non-recent cases

                prediction <- ROCR::prediction(predictions=preds,
                                               labels=truths)
                # get sens and spec
                ss <- ROCR::performance(prediction,"sens","spec")
                fr <- ROCR::performance(prediction,"tpr","fpr")
                auc <- ROCR::performance(prediction,"auc")@y.values %>%
                        unlist()

                perf_df <- data.frame(
                        unlist(ss@x.values),
                        unlist(ss@y.values),
                        unlist(ss@alpha.values))
                colnames(perf_df) <- c(ss@x.name,ss@y.name,ss@alpha.name)
                perf_df <- perf_df %>%
                        mutate(Youden=Sensitivity+Specificity-1)
                # find the cutoff for the youden index
                youden_cutoff <- perf_df %>%
                        filter(Youden==max(Youden)) %>%
                        filter(Sensitivity == max(Sensitivity)) %>%
                        filter(Specificity == max(Specificity)) %>%
                        pull(Cutoff)


                # find the cutoff for the specificity
                spec99_cutoff <- perf_df %>%
                        filter(Specificity>=0.99) %>%
                        filter(Sensitivity == max(Sensitivity)) %>%
                        #want the value lowest cutoff to maximize sensitivity
                        filter(Cutoff == min(Cutoff)) %>%
                        pull(Cutoff)
                
                # find the cutoff for the specificity
                spec95_cutoff <- perf_df %>%
                        filter(Specificity>=0.95) %>%
                        filter(Sensitivity == max(Sensitivity)) %>%
                        #want the value lowest cutoff to maximize sensitivity
                        filter(Cutoff == min(Cutoff)) %>%
                        pull(Cutoff)
                
                # find the cutoff for the specificity
                spec90_cutoff <- perf_df %>%
                        filter(Specificity>=0.90) %>%
                        filter(Sensitivity == max(Sensitivity)) %>%
                        #want the value lowest cutoff to maximize sensitivity
                        filter(Cutoff == min(Cutoff)) %>%
                        pull(Cutoff)


                output_df <- data %>%
                        mutate(prediction=preds) %>%
                        #get youden seropositivity
                        mutate(youden_cutoff = youden_cutoff) %>%
                        mutate(youden_seropos = ifelse(prediction>=youden_cutoff,1,0)) %>%
                        #get 99% specificity seropositivity
                        mutate(spec99_cutoff = spec99_cutoff) %>%
                        mutate(spec99_seropos = ifelse(prediction>=spec99_cutoff,1,0)) %>%
                        #get 95% specificity seropositivity
                        mutate(spec95_cutoff = spec95_cutoff) %>%
                        mutate(spec95_seropos = ifelse(prediction>=spec95_cutoff,1,0))%>%
                        #get 90% specificity seropositivity
                        mutate(spec90_cutoff = spec90_cutoff) %>%
                        mutate(spec90_seropos = ifelse(prediction>=spec90_cutoff,1,0))


        }
        
        list(
                #data
                data= output_df,
                #performance metrics
                ss = ss,
                fr = fr,
                auc = auc
        )
}


getPerf_SuperLearner <- function(data, m,  truth_variable,variables){
        
        # get truths predictions
        truths <- data[truth_variable] %>% unlist() %>%
                as.character() %>% as.numeric()
        
        pred_obj <- predict(m, data[variables], onlySL = TRUE) 
        # predict(m,data=data,predict.all=TRUE)
        preds <- pred_obj$pred %>% unlist()
        
        
        # get performance
        output_df <- NULL
        ss <- NULL
        fr <- NULL
        auc <- NULL
        
        if(!(all(truths==1)|all(truths==0))){ #if some recent or non-recent cases
                
                prediction <- ROCR::prediction(predictions=preds,
                                               labels=truths)
                # get sens and spec
                ss <- ROCR::performance(prediction,"sens","spec")
                fr <- ROCR::performance(prediction,"tpr","fpr")
                auc <- ROCR::performance(prediction,"auc")@y.values %>%
                        unlist()
                
                perf_df <- data.frame(
                        unlist(ss@x.values),
                        unlist(ss@y.values),
                        unlist(ss@alpha.values))
                colnames(perf_df) <- c(ss@x.name,ss@y.name,ss@alpha.name)
                perf_df <- perf_df %>%
                        mutate(Youden=Sensitivity+Specificity-1)
                # find the cutoff for the youden index
                youden_cutoff <- perf_df %>%
                        filter(Youden==max(Youden)) %>%
                        filter(Sensitivity == max(Sensitivity)) %>%
                        filter(Specificity == max(Specificity)) %>%
                        pull(Cutoff)
                
                
                # find the cutoff for the specificity
                spec99_cutoff <- perf_df %>%
                        filter(Specificity>=0.99) %>%
                        filter(Sensitivity == max(Sensitivity)) %>%
                        #want the value lowest cutoff to maximize sensitivity
                        filter(Cutoff == min(Cutoff)) %>%
                        pull(Cutoff)
                
                # find the cutoff for the specificity
                spec95_cutoff <- perf_df %>%
                        filter(Specificity>=0.95) %>%
                        filter(Sensitivity == max(Sensitivity)) %>%
                        #want the value lowest cutoff to maximize sensitivity
                        filter(Cutoff == min(Cutoff)) %>%
                        pull(Cutoff)
                
                # find the cutoff for the specificity
                spec90_cutoff <- perf_df %>%
                        filter(Specificity>=0.90) %>%
                        filter(Sensitivity == max(Sensitivity)) %>%
                        #want the value lowest cutoff to maximize sensitivity
                        filter(Cutoff == min(Cutoff)) %>%
                        pull(Cutoff)
                
                
                output_df <- data %>%
                        mutate(prediction=preds) %>%
                        #get youden seropositivity
                        mutate(youden_cutoff = youden_cutoff) %>%
                        mutate(youden_seropos = ifelse(prediction>=youden_cutoff,1,0)) %>%
                        #get 99% specificity seropositivity
                        mutate(spec99_cutoff = spec99_cutoff) %>%
                        mutate(spec99_seropos = ifelse(prediction>=spec99_cutoff,1,0)) %>%
                        #get 95% specificity seropositivity
                        mutate(spec95_cutoff = spec95_cutoff) %>%
                        mutate(spec95_seropos = ifelse(prediction>=spec95_cutoff,1,0))%>%
                        #get 90% specificity seropositivity
                        mutate(spec90_cutoff = spec90_cutoff) %>%
                        mutate(spec90_seropos = ifelse(prediction>=spec90_cutoff,1,0))
                
                
        }
        
        list(
                #data
                data= output_df,
                #performance metrics
                ss = ss,
                fr = fr,
                auc = auc
        )
}


#get the infection window added
addInfectionWindow<- function(data, end_window){
        
        #define infection window outcome
        outcome <- paste0("inf_",end_window)
        data[outcome] <- (data$day_actual <= end_window &
                                  data$day_actual >= 5 &
                                  data$status == "Case" )%>%
                as.numeric() %>%
                as.factor()
        
        return(data)
}


#conduct loocv on wide dataset by id
loocv_ranger <- function(wide_data, weighted=TRUE, 
                         end_window,
                         variables,
                         num.trees=1000
                         ){
        
        #add the outcome variable to data
        wide_data <- wide_data %>%
                addInfectionWindow(end_window=end_window)
        outcome <- paste0("inf_",end_window)
        
        
        #loop through each id to get loocv
        loocv_df <- data.frame()
        auc_list <- list()
        for ( i in unique(wide_data$id)){
                
                # 0. Choose the individual
                # remove the particular id
                fit_data <- wide_data %>%
                        filter(id != i)  
                
                # 1. run 100 models to find a cutoff
                #keep 50% fit and 50% to test 
                cut_data <- data.frame()
                
                for(j in 1:100){
                        
                        n_ind <-unique(fit_data$id) %>% length() 
                                
                        #define who is in and who is out
                        outside <-sample(fit_data$id,round(n_ind*0.5),replace=FALSE)
                        outside_df <- filter(fit_data,id %in% outside)
                        inside_df  <- filter(fit_data,!id %in% outside)%>%
                                addInfectionWindow(end_window=end_window)
                        
                        # cat(nrow(fit_data)," ",nrow(outside_df)," ",nrow(inside_df),"\n")
                        
                        
                        inside_w = NULL
                        
                        if(weighted){
                                inside_w <- inside_df %>%
                                        getWeight(end_window = end_window)  %>%
                                        pull(weight)
                                
                        }
                        
                        #fit rf model
                        inside_fit <- ranger::ranger(formula=make_formula(paste0("inf_",end_window),
                                                                        variables),
                                                   data=inside_df,
                                                   num.trees = 500,
                                                   case.weights = inside_w,
                                                   replace=TRUE
                        )
                        
                        
                        #find the cutoff for each 
                        performance <- 
                                outside_df %>% addInfectionWindow(end_window) %>%
                                getPerf_ranger(inside_fit,paste0("inf_",end_window)#,
                                               # spec=train_spec_cutoff
                                               )
                        cut <- distinct(performance$data,youden_cutoff,
                                        spec99_cutoff,
                                        spec95_cutoff,
                                        spec90_cutoff
                                        )
                        cut_data <- bind_rows(cut_data,cut)
                }
                
                #find cutoffs
                median_cut <- cut_data %>%
                                summarize(n =n(),
                                          youden_cutoff=median(youden_cutoff),
                                          spec99_cutoff=median(spec99_cutoff),
                                          spec95_cutoff=median(spec95_cutoff),
                                          spec90_cutoff=median(spec90_cutoff)
                                          )
                
                cat(i, "\n")
                print(median_cut)
                
                # 2. fit the model, predict, and apply cutoff 
                #fit the model

                #set weights if needed
                w= NULL
                if(weighted){
                        w <- fit_data%>%
                                getWeight(end_window = end_window)  %>%
                                pull(weight)

                }

                this_fit<- ranger::ranger(data=fit_data,
                                          formula=make_formula(outcome,variables),
                                          num.trees= num.trees,
                                          case.weights = w,
                                          replace = TRUE
                )
                
                # keep the particular id
                #calculate the individual's samples seropositivity
                test_data <- wide_data %>%
                        filter(id == i)
                test_data["truth"] <- test_data[outcome]
                out_preds <-predict(this_fit,
                                    test_data,
                                    predict.all=TRUE)

                #store the data
                loocv_df <- bind_rows(loocv_df,
                                      test_data %>%
                                              mutate(prediction=out_preds[[1]] %>% rowMeans()-1) %>%
                                              mutate(youden_cutoff=median_cut$youden_cutoff,
                                                     spec99_cutoff=median_cut$spec99_cutoff,
                                                     spec95_cutoff=median_cut$spec95_cutoff,
                                                     spec90_cutoff=median_cut$spec90_cutoff
                                                     ) %>%
                                              mutate(
                                                      youden_seropos = ifelse(prediction>youden_cutoff,
                                                                              1,0
                                                      ),
                                                      spec99_seropos = ifelse(prediction>spec99_cutoff,
                                                                            1,0),
                                                      spec95_seropos = ifelse(prediction>spec95_cutoff,
                                                                              1,0),
                                                      spec90_seropos = ifelse(prediction>spec90_cutoff,
                                                                              1,0)
                                                      
                                              ))



}

        return(list(cv_df = loocv_df,
                    formula = make_formula(outcome,variables)
        ))
}


fit_tvsens <- function(data,
                       end_window,
                       seropos_type, # youden_seropos, spec99 spec95 spec90
                       last_time=end_window
){
        
        #get data ready
        tvsens_data <- data %>% filter(truth==1) %>%
                mutate(log_day=log(day_actual)-mean(log(day_actual))) %>%
                mutate(day1=log_day) %>%
                mutate(day2=log_day^2)%>%
                mutate(day3=log_day^3)
        #choose seropositivity cutoff type        
        tvsens_data$seropos_model <- tvsens_data[seropos_type] %>% unlist
        
        
        if(!end_window %in% c(45,120,200,300)){
                stop("not a set time window")
        }
        
        
        
        if(end_window==45){
                logitfit<- stan_glmer(seropos_model ~ (1|id), 
                                      data = tvsens_data,
                                      family = binomial(link = "logit"), 
                )
                spread <- logitfit %>% 
                        spread_draws(`(Intercept)`,
                                     b[term,group])
                
                new_df <- data.frame()
                for(t in unique(c(7:9,seq(10,last_time,7),last_time))){
                        
                        time <- log(t) -mean(log(tvsens_data$day_actual))
                        
                        tmp_df <- spread %>%
                                mutate(logit_theta_j=`(Intercept)` + b
                                ) %>%
                                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                                group_by(.draw) %>%
                                summarize(theta=mean(theta_j),n()) %>%
                                mutate(time=t)
                        
                        new_df <- bind_rows(new_df,tmp_df)
                }
                
        }
        if(end_window==120){
                logitfit<- stan_glmer(seropos_model ~ day1 +(1|id), 
                                      data = tvsens_data,
                                      family = binomial(link = "logit"), 
                )
                
                spread <- logitfit %>% 
                        spread_draws(`(Intercept)`,
                                     b[term,group],
                                     day1)
                new_df <- data.frame()
                for(t in unique(c(7:9,seq(10,last_time,7),last_time))){
                        
                        time <- log(t) -mean(log(tvsens_data$day_actual))
                        
                        tmp_df <- spread %>%
                                mutate(logit_theta_j=`(Intercept)` + b+
                                               day1*time
                                ) %>%
                                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                                group_by(.draw) %>%
                                summarize(theta=mean(theta_j),n()) %>%
                                mutate(time=t)
                        
                        new_df <- bind_rows(new_df,tmp_df)
                }
                
        }
        if(end_window==200){
                
                logitfit<- stan_glmer(seropos_model ~ day1+day2 +(1|id), 
                                      data = tvsens_data,
                                      family = binomial(link = "logit"), 
                )
                
                spread <- logitfit %>% 
                        spread_draws(`(Intercept)`,
                                     b[term,group],
                                     day1,day2)
                new_df <- data.frame()
                for(t in unique(c(7:9,seq(10,last_time,7),last_time))){
                        
                        time <- log(t) -mean(log(tvsens_data$day_actual))
                        
                        tmp_df <- spread %>%
                                mutate(logit_theta_j=`(Intercept)` + b+
                                               day1*time +
                                               day2*time^2
                                ) %>%
                                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                                group_by(.draw) %>%
                                summarize(theta=mean(theta_j),n()) %>%
                                mutate(time=t)
                        
                        new_df <- bind_rows(new_df,tmp_df)
                }
                
                
        }
        if(end_window==300){
                
                logitfit<- stan_glmer(
                        seropos_model ~ day1+day2 + day3+ (1|id), 
                        data = tvsens_data,
                        family = binomial(link = "logit"), 
                )
                
                spread <- logitfit %>% 
                        spread_draws(`(Intercept)`,
                                     b[term,group],
                                     day1,day2,day3)
                
                new_df <- data.frame()
                for(t in unique(c(7:9,seq(10,last_time,7),last_time))){
                        
                        time <- log(t) -mean(log(tvsens_data$day_actual))
                        
                        tmp_df <- spread %>%
                                mutate(logit_theta_j=`(Intercept)` + b+
                                               day1*time +
                                               day2*time^2+
                                               day3*time^3
                                       
                                ) %>%
                                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                                group_by(.draw) %>%
                                summarize(theta=mean(theta_j),n()) %>%
                                mutate(time=t)
                        
                        new_df <- bind_rows(new_df,tmp_df)
                }
                
                
        }
        
        summary_df <- new_df %>%
                group_by(time) %>%
                summarize(
                        lb = quantile(theta,0.025),
                        median=median(theta),
                        ub = quantile(theta,0.975),
                ) %>% mutate(end_window=end_window)
        
        list(
                fit=logitfit,
                data=tvsens_data,
                out_df=new_df,
                summary_df =summary_df
                
        )
        
        
}



fit_spec <- function(data,
                     end_window,
                     seropos_type){
        
        
        #get data ready
        spec_data <- data %>% filter(truth==0) 
        
        #choose seropositivity cutoff type        
        spec_data$seropos_model <- spec_data[seropos_type] %>% unlist
        
        #check if correct window is specified
        if(!end_window %in% c(45,120,200,300)){
                stop("not a set time window")
        }
        
        #run stan model
        logitfit<- stan_glmer(seropos_model ~ (1|id), 
                              data = spec_data,
                              family = binomial(link = "logit"), 
        )
        spread <- logitfit %>% 
                spread_draws(`(Intercept)`,
                             b[term,group])
        
        
        new_df <- spread %>%
                mutate(logit_theta_j=`(Intercept)` + b
                ) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                group_by(.draw) %>%
                summarize(theta=mean(theta_j),n()) 
        
        
        summary_df <- new_df %>%
                summarize(
                        lb = 1-quantile(theta,0.975),
                        median=1-median(theta),
                        ub = 1-quantile(theta,0.025), 
                ) %>% mutate(end_window=end_window)
        
        
        list(
                fit=logitfit,
                data=spec_data,
                out_df=new_df,
                summary_df=summary_df
        )
        
        
}


fit_tvfpr <- function(data,
                       end_window,
                       seropos_type, # youden_seropos, spec99 spec95 spec90
                       last_time=360,
                       curve="cubic" #"linear"
                      
){
        
        #get data ready
        fpr_data <- data %>% 
                mutate(new_day=ifelse(day==0,1,day_actual)) %>%
                mutate(log_day=log(new_day)-mean(log(new_day))) %>%
                mutate(day1=log_day) %>%
                mutate(day2=log_day^2)%>%
                mutate(day3=log_day^3)
        #choose seropositivity cutoff type        
        fpr_data$seropos_model <- fpr_data[seropos_type] %>% unlist
        
        if(!end_window %in% c(45,120,200,300)){
                stop("not a set time window")
        }
        
        #return NULL if no misclassification
        if(sum(fpr_data$seropos_model)==0){
                
                return(list(
                        fit=NULL,
                        data=NULL,
                        out_df=NULL,
                        summary_df =NULL
                        
                ))
        }
        
        if(curve=="cubic"){
                #fit the stan model
                logitfit<- stan_glmer(
                        seropos_model ~ day1+day2 + day3+
                                (1|id), 
                        data = fpr_data,
                        family = binomial(link = "logit"), 
                )
                
                #get the draws        
                spread <- logitfit %>% 
                        spread_draws(`(Intercept)`,
                                     b[term,group],
                                     day1,day2,day3
                        )
                
                new_df <- data.frame()
                # for(t in unique(c(1:9,seq(10,last_time,7),last_time))){
                for(t in 1:last_time){
                        
                        time <- log(t) -mean(log(fpr_data$new_day))
                        
                        tmp_df <- spread %>%
                                mutate(logit_theta_j=`(Intercept)` + b+
                                               day1*time +
                                               day2*time^2+
                                               day3*time^3
                                       
                                ) %>%
                                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                                group_by(.draw) %>%
                                summarize(theta=mean(theta_j),n()) %>% #average across individuals here
                                mutate(time=t)
                        
                        new_df <- bind_rows(new_df,tmp_df)
                }
                
        }
        
        if(curve=="linear"){
                #fit the stan model
                logitfit<- stan_glmer(
                        seropos_model ~ day1+ (1|id), 
                        data = fpr_data,
                        family = binomial(link = "logit"), 
                )
                
                #get the draws        
                spread <- logitfit %>% 
                        spread_draws(`(Intercept)`,
                                     b[term,group],
                                     day1
                        )
                
                new_df <- data.frame()
                # for(t in unique(c(1:9,seq(10,last_time,7),last_time))){
                for(t in 1:last_time){
                                
                        
                        time <- log(t) -mean(log(fpr_data$new_day))
                        
                        tmp_df <- spread %>%
                                mutate(logit_theta_j=`(Intercept)` + b+
                                               day1*time 
                                       
                                ) %>%
                                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                                group_by(.draw) %>%
                                summarize(theta=mean(theta_j),n()) %>% #average across individuals here
                                mutate(time=t)
                        
                        new_df <- bind_rows(new_df,tmp_df)
                }
                
        }
                
        
        summary_df <- new_df %>%
                group_by(time) %>%
                summarize(
                        lb = quantile(theta,0.025),
                        median=median(theta),
                        ub = quantile(theta,0.975),
                ) 
        
        list(
                fit=logitfit,
                data=fpr_data,
                out_df=new_df,
                summary_df =summary_df
                
        )
        
        
}

# #conduct loocv on wide dataset by id
# loocv_ranger <- function(wide_data, weighted, 
#                          end_window,
#                          variables,
#                          num.trees=1000,
#                          train_spec_cutoff=0.99
# ){
#         
#         #add the outcome variable to data
#         wide_data <- wide_data %>%
#                 addInfectionWindow(end_window=end_window)
#         outcome <- paste0("inf_",end_window)
#         
#         
#         #loop through each id to get loocv
#         loocv_df <- data.frame()
#         auc_list <- list()
#         for ( i in unique(wide_data$id)){
#                 
#                 # cat(i,"\n")
#                 # remove the particular id
#                 fit_data <- wide_data %>%
#                         filter(id != i)  %>%
#                         getWeight(end_window = end_window) 
#                 #set weights if needed
#                 w= NULL
#                 if(weighted){
#                         w <- fit_data %>% pull(weight)
#                         
#                 }
#                 #fit the model 
#                 this_fit<- ranger::ranger(data=fit_data,
#                                           formula=make_formula(outcome,variables),
#                                           num.trees= num.trees,
#                                           case.weights = w,
#                                           replace = TRUE
#                 )
#                 
#                 #get the predicted probability in TRAINING set to get youden index
#                 tmp_pred_obj <- predict(this_fit,
#                                         fit_data,
#                                         predict.all=TRUE)
#                 tmp_preds <- tmp_pred_obj[[1]] %>% rowMeans()-1
#                 
#                 
#                 #find the regular cutoffs
#                 youden_cutoff<- find_roc_cutoff(tmp_preds,fit_data[outcome] %>% unlist(),
#                                                 cutoff_type="youden")
#                 spec99_cutoff<- find_roc_cutoff(tmp_preds,fit_data[outcome] %>% unlist(),
#                                                 cutoff_type="spec",
#                                                 train_spec_cutoff=train_spec_cutoff
#                 )
#                 
#                 #find the cutoffs using weighting
#                 weighted_cutoffs <- find_roc_cutoff_weighted(
#                         tmp_preds,
#                         fit_data[outcome] %>% unlist(),
#                         weights = fit_data %>% pull(weight),
#                         train_spec_cutoff=train_spec_cutoff
#                 )
#                 
#                 weighted_youden_cutoff <- weighted_cutoffs$youden_cutoff
#                 weighted_spec99_cutoff <- weighted_cutoffs$spec_cutoff
#                 
#                 # keep the particular id
#                 #calculate the individual's samples seropositivity
#                 test_data <- wide_data %>%
#                         filter(id == i)
#                 test_data["truth"] <- test_data[outcome]
#                 out_preds <-predict(this_fit,
#                                     test_data,
#                                     predict.all=TRUE)
#                 
#                 #store the data
#                 loocv_df <- bind_rows(loocv_df,
#                                       test_data %>%
#                                               mutate(prediction=out_preds[[1]] %>% rowMeans()-1) %>%
#                                               mutate(youden_cutoff=youden_cutoff,
#                                                      spec99_cutoff=spec99_cutoff,
#                                                      weighted_youden_cutoff=weighted_youden_cutoff,
#                                                      weighted_spec99_cutoff=weighted_spec99_cutoff) %>%
#                                               
#                                               mutate(
#                                                       youden_seropos = ifelse(prediction>youden_cutoff,
#                                                                               1,0
#                                                       ),
#                                                       spec_seropos = ifelse(prediction>spec99_cutoff,
#                                                                             1,0),
#                                                       weighted_youden_seropos = ifelse(prediction>weighted_youden_cutoff,
#                                                                                        1,0),
#                                                       weighted_spec_seropos = ifelse(prediction>weighted_spec99_cutoff,
#                                                                                      1,0) 
#                                               ))
#                 
#                 
#                 
#         }
#         
#         return(list(cv_df = loocv_df,
#                     formula = make_formula(outcome,variables)
#         ))
# }




find_roc_cutoff_weighted <-function(preds,
                                    labs,
                                    weights,
                                    train_spec_cutoff=0.99
){
        
        #calculate fpr and tpr for all possible cutpoints
        calc_df <- data.frame(
                preds=preds,
                labs=labs,
                weights=weights
        ) %>% 
                arrange(labs,preds) %>%
                group_by(labs,preds) %>%
                summarize(
                        preds=unique(preds),
                        weights=sum(weights)
                        
                ) %>%
                group_by(labs) %>%
                mutate(cusum=cumsum(weights)) %>%
                mutate(pr=1-cusum/max(cusum)) %>%
                mutate(labs=ifelse(labs==1,"tpr","fpr")) %>%
                spread(labs,pr) %>% ungroup() %>%
                bind_rows(
                        data.frame(preds=0, tpr=1,fpr=1),
                        data.frame(preds=Inf, tpr=0,fpr=0)
                ) %>%
                arrange(preds,-fpr,tpr)
        
        #fill in missing spots
        while(any(is.na(c(calc_df$tpr,calc_df$fpr)))){
                
                calc_df <- calc_df %>%
                        mutate(tpr=ifelse(is.na(tpr),lag(tpr,1),tpr)) %>%
                        mutate(fpr=ifelse(is.na(fpr),lag(fpr,1),fpr))
                
        }
        
        #readjust and calculate youden index at each point
        final_df <- calc_df %>%
                mutate(fpr=lag(fpr,1),
                       tpr=lag(tpr,1)
                ) %>%
                filter(preds!=0) %>%
                mutate(youden=tpr-fpr)%>%
                arrange(-preds)
        
        #calculate cutoffs
        youden_cutoff <- final_df$preds[which.max(final_df$youden)]

        spec_final <- final_df %>%
                filter(preds<=1) %>%
                filter(fpr<=1-train_spec_cutoff) 
        
        spec_cutoff <- spec_final$preds[which.max(spec_final$tpr)]
        
        return(list(youden_cutoff=youden_cutoff,
                    spec_cutoff=spec_cutoff))
}


# #conduct loocv on wide dataset by id
# loocv_ranger <- function(wide_data, weighted=FALSE, 
#                          end_window,
#                          variables,
#                          num.trees=1000
# ){
#         
#         #add the outcome variable to data
#         wide_data <- wide_data %>%
#                 addInfectionWindow(end_window=end_window)
#         outcome <- paste0("inf_",end_window)
#         
#         
#         #loop through each id to get loocv
#         loocv_df <- data.frame()
#         auc_list <- list()
#         for ( i in unique(wide_data$id)){
#                 # remove the particular id
#                 fit_data <- wide_data %>%
#                         filter(id != i) 
#                 
#                 #set weights if needed
#                 w= NULL
#                 if(weighted){
#                         fit_data <- fit_data %>%
#                                 getWeight(end_window = end_window) 
#                         w <- fit_data %>% pull(weight)
#                         
#                 }
#                 #fit the model (weighted and unweighted)
#                 this_fit<- ranger::ranger(data=fit_data,
#                                           formula=make_formula(outcome,variables),
#                                           num.trees= num.trees,
#                                           case.weights = w,
#                                           replace = TRUE
#                 )
#                 
#                 # keep the particular id
#                 test_data <- wide_data   %>%
#                         filter(id == i )
#                 
#                 #get the predicted probability (alternatively we could store this to get an cvAUC???)
#                 #calculate the youden index
#                 #calculate the individual's samples seropositivity
#                 pred_obj <- test_data %>%
#                         getPerf_ranger(m=this_fit,truth_variable=outcome)
#                 
#                 #store individual data ONLY if there is a mix of cases and non cases in fold
#                 if(!is.null(pred_obj$data )){
#                         ind_data <- pred_obj$data %>%
#                                 mutate(end_window=end_window,
#                                        weighted=weighted)
#                         colnames(ind_data)[colnames(ind_data)==outcome] <- "truth"
#                         
#                         #store the seropositivity
#                         loocv_df <- bind_rows(loocv_df,ind_data) 
#                         
#                 }
#                 
#                 #store the auc
#                 auc_list[[i]] <-   pred_obj$auc
#                 
#                 
#         }
#         
#         return(list(cv_df = loocv_df,
#                     auc_list = auc_list,
#                     formula = make_formula(outcome,variables)
#         ))
# }



#conduct loocv on wide dataset by id
kfoldcv_ranger <- function(wide_data, k=10,
                           weighted=FALSE, 
                           end_window,
                           variables,
                           num.trees=1000
){
        
        #add the outcome variable to data
        wide_data <- wide_data %>%
                addInfectionWindow(end_window=end_window)
        outcome <- paste0("inf_",end_window)
        
        ## creating folds (should equal number of people)
        folds <- cv_folds_rm(wide_data,k)
        
        #loop through each fold to do cv
        cv_df <- data.frame()
        auc_list <- list()
        # cat("0 of",k,"folds complete\n")
        for (f in 1:length(folds)){
                
                # remove the particular id
                fit_data <- wide_data[-c(folds[[f]]),] 
                
                #set weights if needed
                w= NULL
                if(weighted){
                        fit_data <- fit_data %>%
                                getWeight(end_window = end_window) 
                        w <- fit_data %>% pull(weight)
                        
                }
                #fit the model (weighted and unweighted)
                this_fit<- ranger::ranger(data=fit_data,
                                          formula=make_formula(outcome,variables),
                                          num.trees= num.trees,
                                          case.weights = w,
                                          replace = TRUE
                )
                
                # keep the particular id
                test_data <- wide_data[c(folds[[f]]),]
                
                #get the predicted probability (alternatively we could store this to get an cvAUC???)
                #calculate the youden index
                #calculate the individual's samples seropositivity
                pred_obj <- test_data %>%
                        getPerf_ranger(m=this_fit,
                                       truth_variable=outcome)
                
                #store the seropositivity
                cv_df <- bind_rows(cv_df,
                                   pred_obj$data %>%
                                           mutate(fold=f)
                )
                
                #store the auc
                auc_list[[f]] <-   pred_obj$auc
                
                # cat(f," of",k,"folds complete\n")
                
        }
        
        #calculate cvAUC using influence curve
        cvAUC <- cvAUC::ci.pooled.cvAUC(
                predictions=cv_df$prediction,
                labels=cv_df[outcome] %>% unlist(),
                ids=cv_df$id,
                folds=cv_df$fold,
                confidence=0.95
        )%>% unlist() %>% t() %>%
        as.data.frame() %>%
        mutate(weighted=weighted,
               end_window=end_window
               )
        
        return(list(cv_df = cv_df %>%
                            mutate(weighted=weighted,
                                   end_window=end_window
                            ),
                    auc_list = auc_list,
                    formula = make_formula(outcome,variables),
                    cvAUC= cvAUC
        ))
}




#conduct loocv on wide dataset by id
kfoldcv_SuperLearner <- function(wide_data, k=10,
                           # weighted=FALSE, 
                           end_window,
                           variables,
                           num.trees=1000,
                           libs = c(
                                    "SL.glmnet",
                                    "SL.bartMachine", 
                                    "SL.ranger",
                                    "SL.xgboost"
                           )
                           
){
        
        #add the outcome variable to data
        wide_data <- wide_data %>%
                addInfectionWindow(end_window=end_window)
        outcome <- paste0("inf_",end_window)
        
        ## creating folds (should equal number of people)
        folds <- cv_folds_rm(wide_data,k)
        
        #loop through each fold to do cv
        cv_df <- data.frame()
        auc_list <- list()
        # cat("0 of",k,"folds complete\n")
        for (f in 1:length(folds)){
                
                # remove the particular id
                fit_data <- wide_data[-c(folds[[f]]),] 
                
                #set weights if needed
                # w= NULL
                # if(weighted){
                #         fit_data <- fit_data %>%
                #                 getWeight(end_window = end_window) 
                #         w <- fit_data %>% pull(weight)
                #         
                # }
                #fit the model (weighted and unweighted)
                # this_fit<- ranger::ranger(data=fit_data,
                #                           formula=make_formula(outcome,variables),
                #                           num.trees= num.trees,
                #                           case.weights = w,
                #                           replace = TRUE
                # )
                
                this_fit<- SuperLearner(X=fit_data[variables],
                                        Y=fit_data[outcome] %>% 
                                                unlist()%>%
                                                as.character() %>%
                                                as.numeric() ,
                                        family = binomial(), 
                                        SL.library = libs
                                        )
                
                # keep the particular id
                test_data <- wide_data[c(folds[[f]]),]
                
                #get the predicted probability (alternatively we could store this to get an cvAUC???)
                #calculate the youden index
                #calculate the individual's samples seropositivity
                pred_obj <- test_data %>%
                        getPerf_SuperLearner(m=this_fit,
                                       truth_variable=outcome,
                                       variables=variables
                                       )
                
                #store the seropositivity
                cv_df <- bind_rows(cv_df,
                                   pred_obj$data %>%
                                           mutate(fold=f)
                )
                
                #store the auc
                auc_list[[f]] <-   pred_obj$auc
                
                cat(f," of",k,"folds complete\n")
                
        }
        
        #calculate cvAUC using influence curve
        cvAUC <- cvAUC::ci.pooled.cvAUC(
                predictions=cv_df$prediction,
                labels=cv_df[outcome] %>% unlist(),
                ids=cv_df$id,
                folds=cv_df$fold,
                confidence=0.95
        )%>% unlist() %>% t() %>%
                as.data.frame() %>%
                mutate(#weighted=weighted,
                       end_window=end_window
                )
        
        return(list(cv_df = cv_df %>%
                            mutate(#weighted=weighted,
                                   end_window=end_window
                            ),
                    auc_list = auc_list,
                    formula = make_formula(outcome,variables),
                    cvAUC= cvAUC
        ))
}


exactBin <- function(p,size,prob) qbinom(p=p,size = size,prob=prob)




### caret functions

#puts things into LOOCV folds
groupKFold_fojo <- function (group, k = length(unique(group))) {
        g_unique <- unique(group)
        m <- length(g_unique)
        if (k > m) {
                stop("`k` should be less than ", m)
        }
        g_folds <- sample(k, size = m, replace = FALSE)
        out <- split(seq_along(group), g_folds[match(group, g_unique)])
        names(out) <- paste0("Fold", gsub(" ", "0", format(seq_along(out))))
        lapply(out, function(z) seq_along(group)[-z])
}




##################### figures #################################


plotDilutionCurve<- function(MFIdat,curve,plate,iso,sub){
    
    #get data for the graph
    plate_data <- filter(MFIdat, batch==plate & isotype ==iso & subclass==sub)
    plate_data_dil <- filter(plate_data, !is.na(dilution))
    plate_data_PBS <- filter(plate_data, sample=="PBS")
    
    
    
    plate_data_samples <- filter(plate_data, !control)
    
    curve_data <- curve %>% filter(batch==plate & isotype ==iso & subclass==sub)


    ggplot()+
        geom_hline(data= plate_data_PBS,aes(yintercept=avgMFI),
                   lty=2,col="blue")+
        geom_point(data=plate_data_dil,aes(x=dilution,y=avgMFI))+
        geom_line(data=curve_data,aes(x=dilution,y=avgMFI))+
        geom_rug(data=plate_data_samples,aes(y=avgMFI#,
              #col=outside
              ),alpha=0.4)+
        # scale_color_manual(values = c("black","red"))+
        facet_wrap(.~antigen,#scales="free"
        )+
        scale_x_continuous(trans="log10")+
        scale_y_continuous(trans="log10")+
        theme_bw()+
        ggtitle(plate)+
        theme(legend.position = "none")
    
    # return(plate_data)
}


#Use loo package to compare bayesian model fits
loo_ing <-function(obj){
    extract <- obj$fit %>% loo::extract_log_lik(merge_chains = F)
    r_eff <- loo::relative_eff(exp(extract))
    loo_obj <- loo::loo(extract,r_eff=r_eff,moment_match = TRUE)
    return(loo_obj)
}

#fit exponential and biphasic bayesian models for a given marker and set of delays
compare_decay_delay_form<- function(delays,marker,iterations=2000,data,cens_df){
    
    # list to hold fits  
    delay_form_fits <- list()
    
    #fit exponential models
    exp_model <- select_univariate_model("exponential")
    
    for (d in delays){
        
        exponential_fit <-  stanfit_decay_univariate(
                data= data %>% filter(status == "Case"),
                cens_df= cens_df %>% filter(status == "Case"),
                marker=marker,
                delay=d,
                covariates=NULL,
                iterations=iterations,
                model= exp_model#,
                # luminex_data_type = "netMFI"
            )
        
        delay_form_fits[[paste0("exponential_",d)]] <- exponential_fit
        
        
    }
    
    
    #fit biphasic models
    bi_model <- select_univariate_model("biphasic")
    
    for (d in delays){
        biphasic_fit <- stanfit_decay_univariate(
                data= data %>% filter(status == "Case"),
                cens_df= cens_df %>% filter(status == "Case"),
                marker=marker,
                delay=d,
                covariates=NULL,
                iterations=iterations,
                model= bi_model#,
                # luminex_data_type = "netMFI"
            )
        #store fits      
        delay_form_fits[[paste0("biphasic_",d)]] <- biphasic_fit
        
    }
    
    #use loo package to compare
    loo_list <- lapply(delay_form_fits,FUN=loo_ing)
    comparison <- loo::loo_compare(loo_list)
    
    
    out <- list(marker=marker,
                delays=delays,
                fits=delay_form_fits,
                loo=comparison)
    return(out)
    
}




# # calculate specificity
# loocv_sens_spec <- function(data,end_window, variables,
#                            weighted=TRUE, 
#                            seropos_type= "youden_seropos" #youden vs spec
#                                    ){
#         # 1. Conduct LOOCV 
#         cv_obj <- loocv_ranger(data,
#                                weighted=weighted,
#                                end_window=end_window,
#                                variables= variables,
#                                num.trees=1000)
# 
#         # 2. Run stan model to get specificity estimates 
#         loocv_df_spec <- cv_obj$cv_df %>%
#                 filter(truth==0)
#         #input function to stan function
#         stan_input_reg <- list(
#                 N = nrow(loocv_df_spec),
#                 J = length(unique(loocv_df_spec$id)),
#                 seropos = loocv_df_spec[seropos_type] %>% unlist(),
#                 stan_id = loocv_df_spec$id %>% as.factor()%>%
#                         as.numeric(),
#                 match = data.frame(id = loocv_df_spec$id,
#                                    stan_id =loocv_df_spec$id %>%
#                                            as.factor()%>%
#                                            as.numeric()
#                 )
#         )
# 
#         filename = "source/final_code/luminex_recommendation/stan/positivity-rate.stan"
#         specfit <- stan(file = filename,data=stan_input_reg,
#                         iter=10000
#         )
#         
#         # 3. Run stan model to get time-carying estimates 
#         loocv_df_sens <- cv_obj$cv_df %>%
#                 filter(truth==1)
#         #input function to stan function
#         stan_input_tv<- list(
#                 N = nrow(loocv_df_sens),
#                 J = length(unique(loocv_df_sens$id)),
#                 seropos = loocv_df_sens[seropos_type] %>% unlist(),
#                 stan_id = loocv_df_sens$id %>% as.factor()%>%
#                         as.numeric(),
#                 log_t = log(loocv_df_sens$day_actual)-mean(log(loocv_df_sens$day_actual)),
#                 match = data.frame(id = loocv_df_sens$id,
#                                    stan_id =loocv_df_sens$id %>%
#                                            as.factor()%>%
#                                            as.numeric()
#                 )
#         )
#         filename = "source/final_code/luminex_recommendation/stan/positivity-rate-time-varying.stan"
#         tvsens_fit <- stan(file = filename,data=stan_input_tv,
#                            iter=10000#, control = list(adapt_delta=0.9)
#         )
#         
#         
#         # 4. Get estimates from both fits
#         #first specificity
#         spec_estimate <-  specfit %>% spread_draws(rate) %>%
#                 summarize(
#                         specificity=1-median(rate),
#                         lb=1-quantile(rate,0.975),
#                         ub=1-quantile(rate,0.025)
#                 ) %>%
#                 mutate( weighted=weighted,
#                         end_window=end_window)
#         
#         #then time-varying sensitivity
#         tvsens_fit_draws <-tvsens_fit %>% spread_draws(alpha_j[stan_id],`beta[1]`,`beta[2]`,`beta[3]`)
#         tvsens_fit_thin <- tvsens_fit_draws %>%
#                 filter(.draw %in% round(seq(1,nrow(.),length.out=2000)))
# 
#         tvsens_df <- data.frame()
#         for(t in c(7:9,seq(10,end_window,7))){
# 
#                 time <- log(t) - mean(log(loocv_df_sens$day_actual))
# 
#                 tmp_df <- tvsens_fit_thin %>%
#                         mutate(logit_theta_j=alpha_j +
#                                        `beta[1]`*time +
#                                        `beta[2]`*time^2+
#                                        `beta[3]`*time^3) %>%
#                         mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
#                         group_by(.draw) %>%
#                         summarize(theta=mean(theta_j),n()) %>%
#                         mutate(time=t)
#                 tvsens_df <- bind_rows(tvsens_df,tmp_df)
#         }
#         
#         # 5. return output
#         return(
#                 list( data=cv_obj$cv_df,
#                       estimates = list(
#                               spec= spec_estimate,
#                               tv_sens=tvsens_df %>%
#                                   mutate(end_window=end_window) %>%
#                                   mutate(weighted=weighted)
#                       ),
#                       fits = list(
#                               spec=specfit,
#                               tv_sens=tvsens_fit
#                               
#                       ))
#         )
# }




# # calculate specificity
# calc_spec <- function(data,end_window, variables,weighted=TRUE){
#         
#         cv_obj <- loocv_ranger(data, 
#                                weighted=weighted, 
#                                end_window=end_window,
#                                variables= variables,
#                                num.trees=1000)
#         #first try spec
#         loocv_df <- cv_obj$cv_df %>%
#                 filter(truth==0)
#         
#         #input function to stan function
#         stan_input_reg <- list(
#                 N = nrow(loocv_df),
#                 J = length(unique(loocv_df$id)),
#                 seropos = loocv_df$seropos,
#                 stan_id = loocv_df$id %>% as.factor()%>%
#                         as.numeric(),
#                 match = data.frame(id = loocv_df$id,
#                                    stan_id =loocv_df$id %>% 
#                                            as.factor()%>%
#                                            as.numeric()
#                 )
#         )
#         
#         filename = "source/final_code/luminex_recommendation/stan/positivity-rate.stan"
#         
#         specfit <- stan(file = filename,data=stan_input_reg,
#                         iter=10000
#         )
#         
#         
#         # specfit %>% shinystan::launch_shinystan()
#         
#         
#         estimate= specfit %>% spread_draws(rate) %>%
#                 summarize(
#                         specificity=1-median(rate),
#                         lb=1-quantile(rate,0.975),
#                         ub=1-quantile(rate,0.025)
#                 ) %>%
#                 mutate( weighted=weighted, 
#                         end_window=end_window)
#         
#         return(list(data=loocv_df,
#                     estimate=estimate
#         ))
#         
#         
# }
# 
# calc_TVsens <- function(data,end_window, variables,weighted=TRUE){
#         
#         #fit model
#         cv_obj <- loocv_ranger(data, 
#                                weighted=weighted, 
#                                end_window=end_window,
#                                variables= variables,
#                                num.trees=1000)
#         
#         loocv_df <- cv_obj$cv_df %>%
#                 filter(truth==1)
#         
#         #input function to stan function
#         stan_input_tv<- list(
#                 N = nrow(loocv_df),
#                 J = length(unique(loocv_df$id)),
#                 seropos = loocv_df$seropos,
#                 stan_id = loocv_df$id %>% as.factor()%>%
#                         as.numeric(),
#                 log_t = log(loocv_df$day_actual)-mean(log(loocv_df$day_actual)),
#                 match = data.frame(id = loocv_df$id,
#                                    stan_id =loocv_df$id %>% 
#                                            as.factor()%>%
#                                            as.numeric()
#                 )
#         )
#         
#         
#         filename = "source/final_code/luminex_recommendation/stan/positivity-rate-time-varying.stan"
#         
#         tvsens_fit <- stan(file = filename,data=stan_input_tv,
#                            iter=10000#, control = list(adapt_delta=0.9)
#         )
#         
#         
#         
#         tvsens_fit_draws <-tvsens_fit %>% spread_draws(alpha_j[stan_id],`beta[1]`,`beta[2]`,`beta[3]`) 
#         
#         tvsens_fit_thin <- tvsens_fit_draws %>%
#                 filter(.iteration %in% round(seq(1,5000,length.out=1000)))
#         
#         avg_df <- data.frame()
#         for(t in c(7:9,seq(10,end_window,7))){
#                 
#                 time <- log(t) - mean(log(loocv_df$day_actual))
#                 
#                 tmp_df <- tvsens_fit_thin %>% 
#                         mutate(logit_theta_j=alpha_j + 
#                                        `beta[1]`*time + 
#                                        `beta[2]`*time^2+  
#                                        `beta[3]`*time^3) %>%
#                         mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
#                         group_by(.draw) %>%
#                         summarize(theta=mean(theta_j),n()) %>%
#                         mutate(time=t)
#                 avg_df <- bind_rows(avg_df,tmp_df)
#                 
#                 
#                 
#                 
#         }
#         
#         return(list(estimate=avg_df %>%
#                             mutate(end_window=end_window) %>%
#                             mutate(weighted=weighted)
#                     ,
#                     data=loocv_df))
#         
# }

# calc_TVfpr <- function(fit_data,end_window, variables,weighted=TRUE,
#                        test_data
# ){
#         
#         #declare infection window
#         fit_data <- addInfectionWindow(fit_data,end_window=end_window)
#         
#         #get weights
#         w <- NULL
#         if(weighted){
#                 fit_data <- fit_data %>% getWeight(end_window = end_window)
#                 w <- fit_data %>% pull(weight) 
#         } 
#         
#         #fit rf model
#         full_fit <- ranger::ranger(formula=make_formula(paste0("inf_",end_window),
#                                                         variables),
#                                    data=fit_data,
#                                    num.trees = 10000,
#                                    case.weights = w
#         )
#         performance <- getPerf_ranger(fit_data,full_fit,"inf_180")
#         cut <- unique(performance$data$cutoff)
#         
#         
#         #get predictions for vaccines
#         pred_obj <- predict(full_fit,data=test_data,predict.all=TRUE)
#         
#         loocv_df <- bind_cols(
#                 test_data,
#                 preds= pred_obj[[1]] %>% rowMeans()-1) %>%
#                 mutate(cutoff=cut) %>%
#                 mutate(seropos= ifelse(preds>cutoff,
#                                        1,0)) %>%
#                 #add 1 to include day 0
#                 mutate(day_actual=ifelse(day==0,1,day))%>%
#                 #for the vaccinees, limit to those above 18
#                 filter(age>=18)
#         
#         
#         #input function to stan function
#         stan_input_tv<- list(
#                 N = nrow(loocv_df),
#                 J = length(unique(loocv_df$id)),
#                 seropos = loocv_df$seropos,
#                 stan_id = loocv_df$id %>% as.factor()%>%
#                         as.numeric(),
#                 log_t = log(loocv_df$day_actual)-mean(log(loocv_df$day_actual)),
#                 match = data.frame(id = loocv_df$id,
#                                    stan_id =loocv_df$id %>% 
#                                            as.factor()%>%
#                                            as.numeric()
#                 )
#         )
#         #fit stan model
#         filename = "source/final_code/luminex_recommendation/stan/positivity-rate-time-varying.stan"
#         tvspec_fit <- stan(file = filename,data=stan_input_tv,
#                            iter=10000#, control = list(adapt_delta=0.9)
#         )
#         
#         #thin it out
#         tvspec_fit_draws <-tvspec_fit %>% spread_draws(alpha_j[stan_id],`beta[1]`,`beta[2]`,`beta[3]`) 
#         tvspec_fit_thin <- tvspec_fit_draws %>%
#                 filter(.iteration %in% round(seq(1,5000,length.out=1000)))
#         
#         #calculate the fpr
#         avg_df <- data.frame()
#         for(t in c(1:9,seq(10,360,7))){
#                 
#                 time <- log(t) - mean(log(loocv_df$day_actual))
#                 
#                 tmp_df <- tvspec_fit_thin %>% 
#                         mutate(logit_theta_j=alpha_j + 
#                                        `beta[1]`*time + 
#                                        `beta[2]`*time^2+  
#                                        `beta[3]`*time^3) %>%
#                         mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
#                         group_by(.draw) %>%
#                         summarize(theta=mean(theta_j),n()) %>%
#                         mutate(time=t)
#                 avg_df <- bind_rows(avg_df,tmp_df)
#         }
#         
#         return(list(
#                 estimate=avg_df %>%
#                         mutate(end_window=end_window) %>%
#                         mutate(weighted=weighted),
#                 data=loocv_df))
# }     



# ggplot geom to add a 2D HPD plot
## https://github.com/mjskay/ggdist/issues/89

stat_hpd_2d <- function(mapping = NULL, data = NULL, geom = "polygon",
                        position = "identity", na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE, n = 100, prob = 0.95, ...) {
        ggplot2::layer(
                stat = StatHPDContour, data = data, mapping = mapping, geom = geom,
                position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, n = n, prob = prob, ...)
        )
}

StatHPDContour <- ggplot2::ggproto(
        "hpd_2d"
        , Stat
        , compute_group = function (data, scales, na.rm = FALSE, h = NULL,
                                    n = 100, prob = 0.95)
        {
                if (is.null(h)) {
                        h <- c(MASS::bandwidth.nrd(data$x), MASS::bandwidth.nrd(data$y))
                }
                dens <- MASS::kde2d(data$x, data$y, h = h, n = n,
                                    lims = c(scales$x$dimension(), scales$y$dimension()))
                df <- data.frame(expand.grid(x = dens$x, y = dens$y), z = as.vector(dens$z))
                df$group <- data$group[1]
                
                dx <- diff(dens$x[1:2])
                dy <- diff(dens$y[1:2])
                sz <- sort(dens$z)
                c1 <- cumsum(sz) * dx * dy
                
                breaks <- sapply(prob, function(x) {
                        withCallingHandlers(
                                stats::approx(c1, sz, xout = 1 - x)$y
                                , warning = function(w) {
                                        if (grepl("collapsing to unique 'x' values", w$message))
                                                invokeRestart("muffleWarning")
                                }
                        )
                })
                
                ggplot2::StatContour$compute_panel(df, scales, breaks = breaks)
        }
        , required_aes = c("x", "y")
)







#need  a different square root function
# for ggplot
mysqrt_trans <- function() {
        trans_new("mysqrt", 
                  transform = base::sqrt,
                  inverse = function(x) ifelse(x<0, 0, x^2),
                  domain = c(0, Inf))
}

