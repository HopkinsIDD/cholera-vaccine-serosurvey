
#extract data from luminex output (.csv) file and put into a list of dataframes
extractExponentcsv <- function(file,token){
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
                                        gather(antigen,value,-c(Location,Sample,`Total Events`)) %>% #make long
                                        mutate(value=as.numeric(value)) %>% #make value numeric
                                        mutate(location=sapply(Location,FUN = get_loc)) %>% #find location in string
                                        mutate(batch=dataList$batch) %>% #record batch name
                                        mutate(date=dataList$date) %>% #record date
                                        mutate(date= as.Date(date,"%m/%d/%Y")) %>%
                                        select(-Location,-`Total Events`) 
                                
                                #rename value column
                                tmp[dataTypes_names[i]] <- tmp$value
                                tmp <- tmp %>% select(-value)
                                
                                
                                
                                
                        }
                        
                }
                
                # tmp <- tmp[,colnames(tmp)!=""] #remove columns with no header
                
                #add dataframe to the list
                dataList[[dataTypes_names[i]]] <- tmp
                
                
        }
        
        
        dataList[["out_df"]] <- left_join(dataList$Median,dataList$Count,
                                          by = c("Sample", "antigen", "location", "batch", "date")) %>% 
                rename(MFI=Median,
                       sample_exponent=Sample,
                       antigen_exponent=antigen
                ) %>%
                left_join(dataList$Units[1,] %>%
                                  select(-`Analyte:`) %>%
                                  gather(antigen_exponent,BeadID),
                          by = "antigen_exponent"
                          )
        
        
        
        return(dataList)
}


drop_excel_sheets <- function(file, dtoken, dest = tempdir(), ...) {
        localfile = paste0(dest, '/', basename(file))
        drop_download(file, localfile, overwrite = TRUE, dtoken = dtoken)
        readxl::excel_sheets(localfile, ...)
}


#extract data from layout file

extractLayout <- function(file,token){
        
        # find the name
        sheets <- drop_excel_sheets(file,dtoken=token)
        
        if("version" %in% sheets){
                #first read in the file to get the version
                version <- drop_read_excel(file,sheet="version",
                                           range=c("A1:A2"),
                                           dtoken = token)
                
        } else version = "1.0"
        
        
        if(version=="1.0") out <- extractLayout_v1.0(file=file,token=token)
        if(version=="1.1") out <- extractLayout_v1.1(file=file,token=token)
        
        return(out)
        
}

extractLayout_v1.0  <- function(file,token){
        layout_1 <- drop_read_excel(file,sheet=1,
                                    range=c("A1:Y17"),
                                    dtoken = token) 
        
        colnames(layout_1)[1] <-"Letter"
        
        description_df <- drop_read_excel(file,sheet=2,
                                          col_types="text",range=c("A2:B4"),
                                          col_names=FALSE,
                                          dtoken = token) 
        
        # need to account for semicolons eventually
        antibody_df <- drop_read_excel(file,sheet=2,
                                       range=c("A7:C16"),
                                       col_names=TRUE,
                                       dtoken = token) %>%
                filter(!is.na(wells)) %>%
                mutate(min=sapply(wells,
                                  FUN= function(x) str_split_fixed(x,"\\:",n=2)[1]))%>%
                mutate(max=sapply(wells,
                                  FUN= function(x) str_split_fixed(x,"\\:",n=2)[2]))%>%
                mutate(min_letter=str_extract(min,"[A-Z]{1,2}"),
                       max_letter=str_extract(max,"[A-Z]{1,2}"),
                       min_number=str_extract(min,"[0-9]{1,2}"),
                       max_number=str_extract(max,"[0-9]{1,2}")
                )%>%
                mutate(letter_range=paste0("[",min_letter,"-",max_letter,"]"))
        
        
        #do a for loop for each row and bind to get this
        ind_antibody_df <- data.frame()
        for(i in 1:nrow(antibody_df)){
                tmp_antibody <-    expand.grid(
                        letter=LETTERS[str_detect(LETTERS,antibody_df$letter_range[i])],
                        number=antibody_df$min_number[i]:antibody_df$max_number[i]) %>%
                        mutate(location=paste0(letter,number),
                               isotype=antibody_df$isotype[i],
                               subclass=antibody_df$subclass[i]
                        )%>%
                        select(-letter,-number)
                
                ind_antibody_df <- bind_rows(ind_antibody_df,tmp_antibody)
                
        }
        
        
        antigen_df <- drop_read_excel(file,sheet=2,
                                      range=c("A19:D110"),
                                      col_names=TRUE,
                                      col_types=c("text","date","text","text"),
                                      dtoken = token) %>%
                filter(!is.na(antigen))%>%
                mutate(coupling_date=as.Date(coupling_date)) %>%
                filter(!is.na(wells)) %>%
                mutate(min=sapply(wells,
                                  FUN= function(x) str_split_fixed(x,"\\:",n=2)[1]))%>%
                mutate(max=sapply(wells,
                                  FUN= function(x) str_split_fixed(x,"\\:",n=2)[2]))%>%
                mutate(min_letter=str_extract(min,"[A-Z]{1,2}"),
                       max_letter=str_extract(max,"[A-Z]{1,2}"),
                       min_number=str_extract(min,"[0-9]{1,2}"),
                       max_number=str_extract(max,"[0-9]{1,2}")
                )%>%
                mutate(letter_range=paste0("[",min_letter,"-",max_letter,"]"))%>%
                rename(antigen_name=antigen)
        
        if(any(str_detect(c(antibody_df$wells, antigen_df$wells),";"))) stop(paste0("semicolon in wells column:",description_df[1,2]))
        
        ind_antigen_df <- data.frame()
        for(i in 1:nrow(antigen_df)){
                tmp_antigen <- expand.grid(letter=LETTERS[str_detect(LETTERS,antigen_df$letter_range[i])],
                                           number=antigen_df$min_number[i]:antigen_df$max_number[i]) %>%
                        mutate(location=paste0(letter,number),
                               lot_number=antigen_df$lot_number[i],
                               antigen_name=antigen_df$antigen_name[i]
                        ) %>%
                        select(-letter,-number)
                
                ind_antigen_df <- bind_rows(ind_antigen_df,tmp_antigen)
                
        }
        
        
        final <- layout_1 %>% gather(Number,sample,-Letter) %>%
                mutate(location=paste0(Letter,Number)) %>%
                select(location,sample) %>%
                filter(!is.na(sample)) %>%
                mutate(batch= description_df[1,2] %>% unlist()) %>%
                left_join(ind_antibody_df, by = "location") %>%
                #label layout version
                mutate(layout_version="1.0")
                #no bead region id in these data
                # # left_join(ind_antigen_df) %>%
                # group_by(location) %>%
                # arrange(antigen_name)%>%
                # mutate(bead_group=paste0(antigen_name,lot_number,collapse="+"))  %>%
                # ungroup()%>%
                # mutate(bead_group=as.numeric(factor(bead_group))) 
        
        
        return(list(summary=list(
                version="1.0",
                layout=layout_1,
                description=description_df,
                antibody=antibody_df,
                antigen= antigen_df),
                out_df = final  ))     
}


extractLayout_v1.1 <- function(file,token){
        layout_1 <- drop_read_excel(file,sheet=1,
                                    range=c("A1:Y17"),
                                    dtoken = token) 
        
        colnames(layout_1)[1] <-"Letter"
        
        
        description_df <- drop_read_excel(file,sheet=2,
                                          col_types="text",range=c("A2:D4"),
                                          col_names=FALSE,
                                          dtoken = token) 
        
        
        # need to account for semicolons eventually
        antibody_df <- drop_read_excel(file,sheet=2,
                                       range=c("A7:D16"),
                                       col_names=TRUE,
                                       dtoken = token) %>%
                filter(!is.na(wells)) %>%
                mutate(min=sapply(wells,
                                  FUN= function(x) str_split_fixed(x,"\\:",n=2)[1]))%>%
                mutate(max=sapply(wells,
                                  FUN= function(x) str_split_fixed(x,"\\:",n=2)[2]))%>%
                mutate(min_letter=str_extract(min,"[A-Z]{1,2}"),
                       max_letter=str_extract(max,"[A-Z]{1,2}"),
                       min_number=str_extract(min,"[0-9]{1,2}"),
                       max_number=str_extract(max,"[0-9]{1,2}")
                ) %>%
                mutate(letter_range=paste0("[",min_letter,"-",max_letter,"]"))
        
        #do a for loop for each row and bind to get this
        ind_antibody_df <- data.frame()
        for(i in 1:nrow(antibody_df)){
                tmp_antibody <-    expand.grid(
                        letter=LETTERS[str_detect(LETTERS,antibody_df$letter_range[i])],
                        number=antibody_df$min_number[i]:antibody_df$max_number[i]) %>%
                        mutate(location=paste0(letter,number),
                               isotype=antibody_df$isotype[i],
                               subclass=antibody_df$subclass[i]
                        )%>%
                        select(-letter,-number)
                
                ind_antibody_df <- bind_rows(ind_antibody_df,tmp_antibody)
                
        }
        
        
        
        antigen_df <- drop_read_excel(file,sheet=2,
                                      range=c("A19:E110"),
                                      col_names=TRUE,
                                      col_types=c("text","text","date","text","text"),
                                      dtoken = token) %>%
                filter(!is.na(antigen_name))%>%
                mutate(coupling_date=as.Date(coupling_date)) %>%
                filter(!is.na(wells)) %>%
                mutate(min=sapply(wells,
                                  FUN= function(x) str_split_fixed(x,"\\:",n=2)[1]))%>%
                mutate(max=sapply(wells,
                                  FUN= function(x) str_split_fixed(x,"\\:",n=2)[2]))%>%
                mutate(min_letter=str_extract(min,"[A-Z]{1,2}"),
                       max_letter=str_extract(max,"[A-Z]{1,2}"),
                       min_number=str_extract(min,"[0-9]{1,2}"),
                       max_number=str_extract(max,"[0-9]{1,2}")
                )%>%
                mutate(letter_range=paste0("[",min_letter,"-",max_letter,"]"))
        
        #do a for loop for each row and bind to get this
        ind_antigen_df <- data.frame()
        for(i in 1:nrow(antigen_df)){
                tmp_antigen <- expand.grid(letter=LETTERS[str_detect(LETTERS,antigen_df$letter_range[i])],
                                           number=antigen_df$min_number[i]:antigen_df$max_number[i]) %>%
                        mutate(location=paste0(letter,number),
                               lot_number=antigen_df$lot_number[i],
                               antigen_name=antigen_df$antigen_name[i],
                               BeadID=antigen_df$BeadID[i]
                        ) %>%
                        select(-letter,-number)
                
                ind_antigen_df <- bind_rows(ind_antigen_df,tmp_antigen)
                
        }
        
        if(any(str_detect(c(antibody_df$wells, antigen_df$wells),";"))) stop(paste0("semicolon in wells column:",description_df[1,2]))
        
        
        final <- layout_1 %>% gather(Number,sample,-Letter) %>%
                mutate(location=paste0(Letter,Number)) %>%
                select(location,sample) %>%
                filter(!is.na(sample)) %>%
                mutate(batch= description_df[1,2] %>% unlist()) %>%
                left_join(ind_antibody_df, by = "location") %>%
                left_join(ind_antigen_df, by = "location") %>%
                group_by(location) %>%
                arrange(antigen_name)%>%
                mutate(bead_group=paste0(antigen_name,lot_number,collapse="+"))  %>%
                ungroup()%>%
                mutate(bead_group=as.numeric(factor(bead_group))) %>%
                #label layout version
                mutate(layout_version="1.1")
        
        
        return(list(summary=list(
                version="1.1",
                layout=layout_1,
                description=description_df,
                antibody=antibody_df,
                antigen= antigen_df),
                out_df = final  ))   
        
}



combineExtraction<-function(batch_id,register,token){
        
        #find the row with the batch id
        this_batch <- filter(register,batch==batch_id)
        
        #set paths of where csv and layouts are
        csv_path <- "vc_luminex_longtudinal/raw_plate_data/final/summaryFI/"
        layout_path <-"vc_luminex_longtudinal/layout_lookup/layouts/new_layouts/"
        csv_file  <-paste0(csv_path,this_batch$data_file %>% unlist())
        xlsx_file <-paste0(layout_path,this_batch$file%>% unlist())
        
        
        #extract the data AND merge
        csv <- extractExponentcsv(csv_file,token)
        layout <- extractLayout(xlsx_file,token)
        
        if (layout$summary$version=="1.0"){
                final <- left_join(csv$out_df,layout$out_df,
                                   by = c("location", "batch")
                                   ) 
        }
        if (layout$summary$version!="1.0"){
                final <- left_join(csv$out_df,layout$out_df,
                                   by = c("location", "batch", "BeadID")
                                   ) 
        }
        
        #run the checks
        #do the number of rows match up between the sets
        if(nrow(csv$out_df)!=nrow(layout$out_df) ) {
                if(nrow(csv$out_df)!=nrow(layout$out_df)*length(unique(layout$summary$antigen$antigen_name))){
                        stop(paste0(csv$batch,":Layout rows (",nrow(layout$out_df),") doesnt match .csv rows (",nrow(csv$out_df),")"))
                }}
        
        #does every sample have 
        # sample id
        # batch id
        # bead count
        # isotype
        
        if(select(final,batch,Count,sample,isotype) %>% is.na() %>% any()){
                final %>% filter(is.na(batch)|is.na(Count)|is.na(sample)|is.na(isotype)) %>%
                        select(batch,location,sample_exponent,Count,sample,isotype)
                stop("Something went wrong with the merge for the above data")
        }
        
        return(final)
        
        
}


fitLL_model <- function(data, title="") {
        
        if(length(unique(data$dilution))<=5){
                
                # data <- data %>% mutate(netmfi_value=log10(netmfi_value))
                #consider adding weights or log transforming
                fit <-  drc::drm(
                        formula=netmfi_value~dilution,
                        data=data,
                        fct=drc::LL.4(
                                fixed = c(NA, 10, NA, NA), 
                                names = c("b", "c", "d", "e")))
                
                
                model_type="Net MFI, 3P-log-logistic"
                
                
                #get the R^2 value
                # avg <- mean(data$netmfi_value)
                # preds <- fit$predres[,1]
                # R.squared <- sum((preds-avg)^2)/sum((data$netmfi_value-avg)^2)
                
                avg <- mean(data$netmfi_value)
                residuals <- fit$predres[,2]
                R.squared <- 1-sum(residuals^2)/sum((data$netmfi_value-avg)^2)
                
                
                if(R.squared>1) warning("R-squared over 1")
                
        } else {
                # data <- data %>% mutate(MFI=log10(MFI))
                
                
                #consider adding weights or log transforming
                fit <-  drc::drm(MFI~dilution,
                                 data=data,
                                 fct=drc::LL.5())
                model_type="MFI, 5P-log-logistic"
                
                #get the R^2 value
                # avg <- mean(data$MFI)
                # preds <- fit$predres[,1]
                # R.squared <- sum((preds-avg)^2)/sum((data$MFI-avg)^2)
                
                avg <- mean(data$MFI)
                residuals <- fit$predres[,2]
                R.squared <- 1-sum(residuals^2)/sum((data$MFI-avg)^2)
                
                
                if(R.squared>1) warning("R-squared over 1")
                
                
        }
        
        #plot the curve
        # plot(fit,type="all")
        # title(combo_df[r,] %>% unlist %>% paste(collapse=" "))
        # text(x=500, y=150,col="red",
        #      labels=paste("R^2: ",round(R.squared,2))
        # )
        # 
        
        
        
        
        out <- list( data=data,
                     fit=fit,
                     model_type=model_type,
                     R.squared = R.squared#,
                     # plot=recordPlot()
                     
        )
        
        return(out)
        
}


# new  RAU function
get3PLL_RAU<- function(value,value_range,dilution_range,fit,lb=10){
        
        if(value>=max(value_range)) RAU <- min(dilution_range)
        if(value<=min(value_range)) RAU <- max(dilution_range)
        
        if(value <max(value_range) & value > min(value_range)){
                B <- fit$coefficients[1]
                C <- lb
                D <- fit$coefficients[2]
                E <- fit$coefficients[3]
                
                RAU <- exp(log((D-C)/(value - C) -1)/B +log(E))
        }
        
        return(RAU)
        
}

get5PLL_RAU<- function(value,fit){
        
        value_range <- fit$data$MFI
        dilution_range <- fit$data$dilution
        
        if(value>=max(value_range)) RAU <- min(dilution_range)
        if(value<=min(value_range)) RAU <- max(dilution_range)
        
        if(value <max(value_range) & value > min(value_range)){
                B <- fit$coefficients[1]
                C <- fit$coefficients[2]
                D <- fit$coefficients[3]
                E <- fit$coefficients[4]
                F_param <- fit$coefficients[5]
                
                Q <- (D-C)/(value - C)
                
                RAU <- exp(log(Q^(1/F_param)-1)/B +log(E))
        }
        
        return(RAU)
        
}

