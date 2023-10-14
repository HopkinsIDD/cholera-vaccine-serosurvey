source("code/R/packages.R")
source("code/R/utils.R")

#load cleaned data
clean_sample_data <- read_rds('data/generated_data/quality_control/clean_sample_data.rds')


#check the different study codes for trajectories to see if there are any bad ones
study_codes <- unique(clean_sample_data$studycode)
s = 4


tmp_trajectory <- clean_sample_data %>%
        filter(subclass=="total")%>%
        filter(studycode ==study_codes[s]) %>%
        filter(str_detect(antigen_exponent,"CtxB|Ogawa"))

tmp_ids <- unique(tmp_trajectory$id)
select_ids <-1:10

tmp_trajectory %>%
        filter(id %in% tmp_ids[select_ids]) %>%
        ggplot(aes(x=day,y=1/RAU_value,col=antigen_exponent))+
        geom_line()+
        geom_point()+
        scale_y_continuous(trans = "log10")+
        facet_grid(isotype~id)+
        cowplot::theme_cowplot()+
        theme(axis.text.x = element_text(angle=45,hjust=1))

select_ids <- select_ids +10

#Maybe remove
#IMS_C004 Day 17
#IMS_B004 Day 17
#maybe remove those with a meteoric rise?


#get samples ready for analysis
#output for analysis
ok_samples <- clean_sample_data %>%
        filter(subclass=="total")%>%
        mutate(marker=paste(isotype,antigen_exponent)) %>%
        #identify those with at least 31 markers (not including IgA and IgM O139)
        filter(!marker %in% c(
                "IgA O139:BSA",
                "IgM O139:BSA"
        )) %>%
        group_by(sample) %>%
        summarize(n=n()) %>%
        filter(n==31) %>%
        pull(sample)


#bring in the demographic data
demographic_data <- getWideData() %>%
        select(id, status,age,age_group,culture,sex,blood)

#get the long dataset together
final_df <- clean_sample_data %>%
                #limit to those only with enough measurements or challenge study (did not of 0139 IgG)
                filter(subclass=="total")%>%
                filter(studycode %in% c("IMS","P","S","R2","E01JH")) %>%
                filter(sample %in% ok_samples | studycode =="E01JH") %>%
                #get RAU value in a log10 inverse form
                mutate(RAU_value=log10(1/RAU_value)) %>%
                select(studycode,id,sample,day,isotype,antigen_exponent,RAU_value) %>%
                #add demographic data
                left_join(demographic_data,by="id") %>%
                mutate(marker=glue::glue("RAU_{isotype}_{str_remove_all(antigen_exponent,' |:')}")) %>%
                mutate(antigen_pretty = case_when(
                        antigen_exponent=="Inaba OSP:BSA" ~ "Inaba OSP",
                        antigen_exponent=="Ogawa OSP:BSA" ~ "Ogawa OSP",
                        antigen_exponent=="O139:BSA" ~ "O139 OSP",
                        antigen_exponent=="CT HT" ~ "CT-H",
                        antigen_exponent=="CtxB" ~ "CT-B",
                        antigen_exponent=="LTh" ~ "LT-H",
                        TRUE ~ antigen_exponent
                )) %>%
                mutate(cohort=case_when(
                        studycode == "R2" ~ "HTI Vaccinee",
                        studycode %in% c("S","P") ~ "SMIC/PIC",
                        studycode=="IMS" ~ "BGD Vaccinee",
                        studycode=="E01JH" ~ "ETEC Challenge"
                )) %>%
                #maybe a temporary or permanent fix depending on needs
                mutate(day_actual=day)


final_wide <- final_df %>% select(studycode,cohort,status,id, day,
                                    sample,age,sex,marker,RAU_value,culture) %>%
          spread(marker,RAU_value) %>%
        #maybe a temporary or permanent fix depending on needs
          mutate(day_actual=day)


write_rds(final_df,"data/generated_data/analysis_data/final_df.rds")
write_rds(final_wide,"data/generated_data/analysis_data/final_wide.rds")
