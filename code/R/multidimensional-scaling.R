
source("code/R/packages.R")
source("code/R/utils.R")
source("code/R/load-data.R")

#generate the mds 
generate_MDS <- function(data,regex_string="Ogawa|Inaba|TcpA|CtxB|IgG_O139"){
        
        mds_data <- data %>% 
                select(status,matches(regex_string))
        # select(vaxinf_class,matches("Ogawa|Inaba|TcpA|IgG_O139"))
        
        mds_dist = dist(x = scale(mds_data[,-1]))
        mds_analysis <- smacof::mds(delta = mds_dist , ndim = 2 , type = "ratio")
        
        out <- bind_cols(select(mds_data,status),
                         as.data.frame(mds_analysis$conf))
        
        return(out)
        
}



# compare_data_wide <- final_df %>%
#         # only individuals from bandladesh
#         filter(cohort %in% c("BGD Vaccinee","SMIC/PIC")) %>%
#         # recent
#         filter(day>3,day<50) %>%
#         #only cases and vaccinees
#         filter(status %in% c("Case","Vaccinee"))%>%
#         select(sample,id,isotype,antigen_pretty,marker,RAU_value,day,age,status) %>% 
#         select(-isotype,-antigen_pretty) %>%
#         # mutate(RAU_value=log10(1/RAU_value)) %>%
#         spread(marker,RAU_value)
# mds_obj <- bind_rows(
#         compare_data_wide %>% 
#                 generate_MDS() %>% mutate(variables="All"),
#         compare_data_wide %>% 
#                 generate_MDS(regex_string="IgG_Ogawa|IgG_Inaba|IgG_CtxB") %>% mutate(variables="Original"),
#         compare_data_wide %>% 
#                 generate_MDS(regex_string="CtxB") %>%  mutate(variables="CtxB")
# )


compare_data_wide <- final_df %>%
        # only individuals from bandladesh
        filter(cohort %in% c("BGD Vaccinee","SMIC/PIC", "HTI Vaccinee")) %>%
        #only cases and vaccinees
        filter(status %in% c("Case","Vaccinee"))%>%
        select(sample,id,isotype,antigen_pretty,marker,RAU_value,day,age,status) %>% 
        select(-isotype,-antigen_pretty) %>%
        # mutate(RAU_value=log10(1/RAU_value)) %>%
        spread(marker,RAU_value)


mds_obj <- bind_rows(
        compare_data_wide %>% filter(day<=45) %>%
                generate_MDS() %>% mutate(time_window="<=45 days"),
        compare_data_wide %>% filter(day>45 & day<=120)%>%
                generate_MDS() %>% mutate(time_window="46-120 days"),
        compare_data_wide %>% filter(day>120 & day<=200)%>%
                generate_MDS() %>% mutate(time_window="121-200 days")
)%>%
        mutate(time_window=factor(time_window,
                                  levels=c("<=45 days",
                                           "46-120 days",
                                           "121-200 days"
                                  )
        ))


write_rds(mds_obj,"data/generated_data/mds_obj.rds")




        