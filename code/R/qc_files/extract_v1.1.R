## code to extract data from template file

layout_1 <- read_excel(paste0(
  "data/raw_data/layout_lookup/layouts/new_layouts/",
  matched_layout$file),sheet=1,
  range=c("A1:Y17")) 
colnames(layout_1)[1] <-"Letter"

description_df <- read_excel(paste0(
  "data/raw_data/layout_lookup/layouts/new_layouts/",
  matched_layout$file),sheet=2,
  col_types="text",range=c("A2:D4"),
  col_names=FALSE) 

# need to account for semicolons eventually
antibody_df <- read_excel(paste0(
  "data/raw_data/layout_lookup/layouts/new_layouts/",
  matched_layout$file),sheet=2,
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
  tmp_antibody <- expand.grid(letter=LETTERS[str_detect(LETTERS,antibody_df$letter_range[i])],
                              number=antibody_df$min_number[i]:antibody_df$max_number[i]) %>%
    mutate(location=paste0(letter,number),
           isotype=antibody_df$isotype[i],
           subclass=antibody_df$subclass[i]
    )%>%
    select(-letter,-number)
  
  ind_antibody_df <- bind_rows(ind_antibody_df,tmp_antibody)
  
}

antigen_df <- read_excel(
  paste0(
  "data/raw_data/layout_lookup/layouts/new_layouts/",
    matched_layout$file),sheet=2,
  range=c("A19:E110"),
  col_names=TRUE,
  col_types=c("text","text","date","text","text")) %>%
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
  ) %>%
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


beadID_df<- extract$Units[1,-1] %>% gather(antigen,BeadID)

full_data <- extract$Median %>% left_join(extract$Count) %>%
  left_join(beadID_df)%>%
  left_join(ind_antigen_df)%>%
  left_join(ind_antibody_df) %>% 
  mutate(location_letter= factor(str_extract(location,"[A-Z]")),
         location_number= as.numeric(str_extract(location,"[0-9]{1,}")),
  ) %>%
  mutate(location_letter=fct_rev(location_letter)) %>%
  mutate(antigen =ifelse(str_detect(antigen,"VCC"),"VCC",antigen))

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


############ Calculate RAU


standards <- control_MFI %>% filter(str_detect(Sample,"O1 Standard C")) %>%
                group_by(Sample, dilution,antigen) %>%
                summarize(MFI=(mean(Median)),
                       n())

RAU_data <- data.frame()

for(anti in unique(standards$antigen)){
        
        fitLL <- fitLL_model(filter(standards,antigen==anti))
        
        tmp_sample_data <- sample_MFI %>% filter(antigen==anti) %>%
                mutate(RAU=sapply(Median,FUN=get5PLL_RAU,fit=fitLL$fit))
        
        
        
        RAU_data <- bind_rows(RAU_data,tmp_sample_data)

}









