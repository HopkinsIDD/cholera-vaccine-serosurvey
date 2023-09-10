
#extract data from luminex output (.csv) file and put into a list of dataframes
extractExponentcsv <- function(file){
        
        #get the plate location from the luminex output string
        get_loc <- function(string){
                loc <- str_split(string,",")
                loc <- unlist(loc)[2]
                loc <- str_remove(loc,"\\)")
                
                return(loc)
        }
        
        
        #read in .CSV
        # trial <- drop_read_csv(file,
        #                        dtoken=token,
        #                        header=FALSE,
        #                        col.names=paste("column", 1:26, sep="_"),
        #                        colClasses = "character"
        # )
        
        trial <- read.csv(file,
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

extractLayout <- function(file){
        
        # find the name
        sheets <- readxl::excel_sheets(file)
        
        if("version" %in% sheets){
                #first read in the file to get the version
                version <- readxl::read_excel(file,sheet="version",
                                           range=c("A1:A2"))
                
        } else version = "1.0"
        
        
        if(version=="1.0") out <- extractLayout_v1.0(file=file)
        if(version=="1.1") out <- extractLayout_v1.1(file=file)
        
        return(out)
        
}

extractLayout_v1.0  <- function(file){
        layout_1 <- readxl::read_excel(file,sheet=1,
                                    range=c("A1:Y17")) 
        
        colnames(layout_1)[1] <-"Letter"
        
        description_df <- readxl::read_excel(file,sheet=2,
                                          col_types="text",range=c("A2:B4"),
                                          col_names=FALSE) 
        
        # need to account for semicolons eventually
        antibody_df <- readxl::read_excel(file,sheet=2,
                                       range=c("A7:C16"),
                                       col_names=TRUE) %>%
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
        
        
        antigen_df <- readxl::read_excel(file,sheet=2,
                                      range=c("A19:D110"),
                                      col_names=TRUE,
                                      col_types=c("text","date","text","text")) %>%
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
        layout_1 <- readxl::read_excel(file,sheet=1,
                                    range=c("A1:Y17")) 
        
        colnames(layout_1)[1] <-"Letter"
        
        
        description_df <- readxl::read_excel(file,sheet=2,
                                          col_types="text",range=c("A2:D4"),
                                          col_names=FALSE) 
        
        
        # need to account for semicolons eventually
        antibody_df <- readxl::read_excel(file,sheet=2,
                                       range=c("A7:D16"),
                                       col_names=TRUE) %>%
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
        
        
        
        antigen_df <- readxl::read_excel(file,sheet=2,
                                      range=c("A19:E110"),
                                      col_names=TRUE,
                                      col_types=c("text","text","date","text","text")
                                      ) %>%
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

fitLL_model <- function(data) {
        
        #check on 
        missing_blank<- any(is.na(data$blank))
        available_dilutions<- length(unique(data$dilution))
        
        if(!missing_blank & available_dilutions<=5 &available_dilutions>3){
                
                blank_value <- unique(data$blank)
                
                #consider adding weights or log transforming
                fit <-  drc::drm(
                        formula=avgMFI~dilution,
                        data=data,
                        fct=drc::LL.4(
                                fixed = c(NA, blank_value, NA, NA), 
                                names = c("b", "c", "d", "e")))
                
                
                model_type="Net MFI, 3P-log-logistic"
                
                avg <- mean(data$avgMFI)
                residuals <- fit$predres[,2]
                R.squared <- 1-(sum(residuals^2)/sum((data$avgMFI-avg)^2))
                
                
                if(R.squared>1) warning("R-squared over 1")
                
                out <- list( data=data,
                             fit=fit,
                             model_type=model_type,
                             R.squared = R.squared
                )
        } 
        if(available_dilutions>5){
                
                #consider adding weights or log transforming
                fit <-  drc::drm(avgMFI~dilution,
                                 data=data,
                                 fct=drc::LL.5())
                model_type="MFI, 5P-log-logistic"
                
                avg <- mean(data$avgMFI)
                residuals <- fit$predres[,2]
                R.squared <- 1-sum(residuals^2)/sum((data$avgMFI-avg)^2)
                
                
                if(R.squared>1) warning("R-squared over 1")
                
                out <- list( data=data,
                             fit=fit,
                             model_type=model_type,
                             R.squared = R.squared
                )
                
        }
        
        
        if(!(!missing_blank & available_dilutions<=5 &available_dilutions>3) & !(available_dilutions>5)){
                
                out <- list(model_type="Error",
                            missing_blank=missing_blank,
                            available_dilutions=available_dilutions
                            )
        }
                
        
        
        return(out)
        
}


# new  RAU function
get3PLL_RAU<- function(value,fit, max_dilution=10^5,min_dilution=10^2){
        
        
        max_value<- predict(fit$fit,newdata = data.frame(dilution=min_dilution))
        min_value<- predict(fit$fit,newdata = data.frame(dilution=max_dilution))
        

        if(value>=max_value) RAU <- min_dilution
        if(value<=min_value) RAU <- max_dilution

        if(value <max_value& value > min_value){
                B <- fit$fit$coefficients[1]
                C <- unique(fit$data$blank)
                D <- fit$fit$coefficients[2]
                E <- fit$fit$coefficients[3]

                RAU <- exp(log((D-C)/(value - C) -1)/B +log(E)) %>%
                        as.numeric()

        }
        
        return(RAU)
        
}

get5PLL_RAU<- function(value,fit, max_dilution=10^5,min_dilution=10^2){
        
        max_value<- predict(fit$fit,newdata = data.frame(dilution=min_dilution))
        min_value<- predict(fit$fit,newdata = data.frame(dilution=max_dilution))
        
        if(value>=max_value) RAU <- min_dilution
        if(value<=min_value) RAU <- max_dilution
        
        if(value <max_value & value > min_value){
                B <- fit$fit$coefficients[1]
                C <- fit$fit$coefficients[2]
                D <- fit$fit$coefficients[3]
                E <- fit$fit$coefficients[4]
                F_param <- fit$fit$coefficients[5]
                
                Q <- (D-C)/(value - C)
                
                RAU <- exp(log(Q^(1/F_param)-1)/B +log(E))
        }
        
        return(RAU)
        
}

#combine lapply and brind_rows
bind_lapply <-function(your_list,
                       your_function,
                       your_id){
        lapply(your_list,FUN = your_function) %>%
                bind_rows(.id=your_id)
} 


#takes in a antigen object from standard curves and gives a data.frame
pred_curve <- function(standard_antigen,
                       pred_diluctions= 10^(seq(0.1,7,0.1))
) {
        
        
        predicted_value <- predict(
                standard_antigen$fit,
                newdata = data.frame(dilution=pred_diluctions)
        )
        
        data.frame(prediction=predicted_value,
                   dilution=pred_diluctions
        ) %>%
                mutate(model_type=standard_antigen$model_type) %>%
                mutate(r_squared = standard_antigen$R.squared)
}

# # new  RAU function
# get3PLL_RAU<- function(value,fit, max_dilution=10^5){
# 
#         value_range <- fit$data$netmfi_value
#         dilution_range <- fit$data$dilution
#         
#         min_value<- predict(this_fit$fit,newdata = data.frame(dilution=max_dilution))
#         
#         if(value>=max(value_range)) RAU <- min(dilution_range)
#         if(value<=min_value) RAU <- max_dilution
#         
#         if(value <max(value_range) & value > min_value){
#                 B <- fit$coefficients[1]
#                 C <- 0
#                 D <- fit$coefficients[2]
#                 E <- fit$coefficients[3]
#                 
#                 RAU <- exp(log((D-C)/(value - C) -1)/B +log(E))
#         }
#         
#         return(RAU)
#         
# }
# 
# get5PLL_RAU<- function(value,fit){
#         
#         value_range <- fit$data$MFI
#         dilution_range <- fit$data$dilution
#         
#         if(value>=max(value_range)) RAU <- min(dilution_range)
#         if(value<=min(value_range)) RAU <- max(dilution_range)
#         
#         if(value <max(value_range) & value > min(value_range)){
#                 B <- fit$coefficients[1]
#                 C <- fit$coefficients[2]
#                 D <- fit$coefficients[3]
#                 E <- fit$coefficients[4]
#                 F_param <- fit$coefficients[5]
#                 
#                 Q <- (D-C)/(value - C)
#                 
#                 RAU <- exp(log(Q^(1/F_param)-1)/B +log(E))
#         }
#         
#         return(RAU)
#         
# }

