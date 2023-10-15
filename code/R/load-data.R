final_df <- read_rds("data/generated_data/analysis_data/final_df.rds")
final_wide <- read_rds("data/generated_data/analysis_data/final_wide.rds")


# final_df <- read_rds("data/generated_data/sample_measurements.rds") %>%
#         #remove haiti children
#         filter(study_code!="R1") %>%
#         #remove weird cohort bangladeshi
#         filter(!str_detect(id,"RB"))%>%
#         mutate(day_actual=day) %>%
#         mutate(RAU_value=log10(1/RAU_value)) %>%
#         #create comparison groups
#         mutate(Vaccinee=case_when(
#                 cohort == "BGD Vaccinee" & age<10 ~ "Bangladeshi <10 years",
#                 cohort == "BGD Vaccinee" & age>=10 ~  "Bangladeshi 10+ years",
#                 cohort == "HTI Vaccinee" & age>=10 ~  "Haitian 18+ years"
#         )) %>%
#         mutate(Under18=ifelse(age<18,"<18 years","18+ years"))
# 
# 
# 
# final_wide <- final_df %>% select(study_code,cohort,status,id, day,
#                                   sample,age,Under18,sex,marker,RAU_value,culture) %>%
#         spread(marker,RAU_value) %>%
#         #remove those missing any of the three key antigens for IgG
#         filter(!if_any(.cols= c("RAU_IgG_OgawaOSPBSA",
#                                 "RAU_IgG_InabaOSPBSA",
#                                 "RAU_IgG_CtxB"),
#                        .fns = is.na )) %>%
#         mutate(day_actual=day)
# 
# #limit only to those individuals who have these three antigens
# final_df <- final_df %>% filter(sample %in% final_wide$sample)
# 
# 
# write_rds(final_df,"data/generated_data/analysis_data/final_df.rds")
# write_rds(final_wide,"data/generated_data/analysis_data/final_wide.rds")

