file <- paste0(
        "data/raw_data/layout_lookup/layouts/new_layouts/",
        matched_layout$file)

layout_1 <- read_excel(file,sheet=1,
                            range=c("A1:Y17")) 

colnames(layout_1)[1] <-"Letter"

description_df <- read_excel(file,sheet=2,
                                  col_types="text",range=c("A2:B4"),
                                  col_names=FALSE) 

# need to account for semicolons eventually
antibody_df <- read_excel(file,sheet=2,
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


antigen_df <- read_excel(file,sheet=2,
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

beadID_df<- extract$Units[1,-1] %>% gather(antigen,BeadID)

full_data <- extract$Median %>% left_join(extract$Count) %>%
        left_join(beadID_df)%>%
        left_join(ind_antigen_df)%>%
        left_join(ind_antibody_df) %>%
        mutate(location_letter= factor(str_extract(location,"[A-Z]")),
               location_number= as.numeric(str_extract(location,"[0-9]{1,}")),
        ) %>%
        mutate(location_letter=fct_rev(location_letter)) %>%
        mutate(antigen =ifelse(str_detect(antigen,"VCC"),"VCC",antigen))# %>%
        # left_join(final)

# if(length(unique(full_data$antigen))>length(unique(full_data$lot_number))){
#         full_data <- full_data %>% mutate(antigen=paste0(antigen,
#                                                          " (",lot_number,")"))
#         
# }


# if (nrow(antibody_df==1)) full_data <- full_data %>% mutate(isotype=antibody_df$isotype)
# if (nrow(antibody>1))  full_data
#         mutate(hi=str_detect(location_letter,"[B-G]") & location_number<13)


MFI_df <- full_data  %>%
        filter(Count>=30)


control_MFI <- MFI_df %>% filter(str_detect(Sample,"Control")) %>% 
        mutate(Sample=ifelse(str_detect(Sample,"Control"),
                             Sample,"Non-Control"
        )) %>% 
        mutate(Sample= str_remove(Sample,"Control: ")) %>%
        mutate(dilution=sapply(Sample,get_dilution))


sample_MFI <- filter(MFI_df,!str_detect(Sample,"Control|Hi Pos"))