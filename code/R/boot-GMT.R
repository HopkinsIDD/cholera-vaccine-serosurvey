#boot strap GM MFI ratio
samples_vax_child <- final_df %>%
        filter(status=="Vaccinee",age<10,day==0) %>%
        distinct(sample,status) %>%
        mutate(`Age Group`="Children")

samples_vax_adult <- final_df%>%
        filter(status=="Vaccinee",age>=10,day==0) %>%
        distinct(sample,status) %>%
        mutate(`Age Group`="Adults")

samples_case_child <- final_df %>%
        filter(status=="Case",age<10,day==2) %>%
        distinct(sample,status) %>%
        mutate(`Age Group`="Children")

samples_case_adult <- final_df %>%
        filter(status=="Case",age>=10,day==2) %>%
        distinct(sample,status) %>%
        mutate(`Age Group`="Adults")


baseline_sero_data <- final_df %>% 
        filter((status=="Vaccinee" & day==0)|status=="Case" & day==2) %>%
        filter(antigen %in% c("CtxB","OgawaOSPBSA","InabaOSPBSA","O139BSA")) %>%
        mutate(`Age Group`=ifelse(age>=10,"Adults","Children")) 


bootDF <- data.frame()
for(i in 1:1000){
        
        cat(i,"\n")
        #sample bootstrap
        tmp_sample_vax_child <- samples_vax_child[sample(1:nrow(samples_vax_child),
                                                         nrow(samples_vax_child),replace=TRUE),]
        tmp_sample_vax_adult <- samples_vax_adult[sample(1:nrow(samples_vax_adult),
                                                         nrow(samples_vax_adult),replace=TRUE),]
        tmp_sample_case_child<- samples_case_child[sample(1:nrow(samples_case_child),
                                                          nrow(samples_case_child),replace=TRUE),]
        tmp_sample_case_adult <- samples_case_adult[sample(1:nrow(samples_case_adult),
                                                           nrow(samples_case_adult),replace=TRUE),]
        
        #calculate the GMT ratio
        tmp_ratio <- bind_rows(
                tmp_sample_vax_child,
                tmp_sample_vax_adult,
                tmp_sample_case_child,
                tmp_sample_case_adult
        ) %>% 
                left_join(baseline_sero_data,by=c("sample", "status", "Age Group")) %>%
                group_by(isotype,antigen,`Age Group`,status) %>%
                summarize(avg_value=mean(RAU_value)) %>%
                group_by(isotype,antigen,`Age Group`) %>%
                summarize(difference=avg_value[status=="Vaccinee"]-avg_value[status=="Case"]) %>%
                mutate(mfi_ratio=10^difference)
        
        #store value
        bootDF <- bind_rows(bootDF,
                            tmp_ratio %>% mutate(bootid=i)
        )
        
        
}


write_rds(bootDF,"source/final_code/vax_strategy/generated_rds/bootDF.rds")

