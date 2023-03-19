analysisData<- LuminexRecPaperData(vaccinee = TRUE)
casecon_analysis <- analysisData$RAU_data$long_data
final_df <- read_rds("data/generated_data/sample_measurements.rds")




#list of antigens associated with O1 cholera infection
O1_antigens <- c("CTHT","CtxB","InabaOSPBSA","OgawaOSPBSA","Sialidase","TcpA","VCC")

#data frame to join to make antigen names look nice for plots
antigen_df <- distinct(casecon_analysis,antigen) %>%
        filter(!is.na(antigen)) %>%
        mutate(O1_antigen= antigen %in% O1_antigens) %>%
        mutate(antigen_pretty = recode(antigen,
                                       "InabaOSPBSA" = "Inaba OSP",
                                       "OgawaOSPBSA" = "Ogawa OSP",
                                       "O139BSA" = "O139 OSP",
                                       "CTHT" = "CT-H",
                                       "CtxB" = "CT-B",
                                       "LTh" = "LT-H",
                                       "LTB"  = "LT-B"
        )) %>%
        mutate(antigen_pretty = factor(antigen_pretty,
                                       levels=c("CT-B","Ogawa OSP","Inaba OSP","O139 OSP","TcpA","VCC","Sialidase","CT-H","LT-B","LT-H","LPS","Flu")
        ))



vax_long_RAU <- analysisData$RAU_data$long_data %>%
        filter(status=="Vaccinee") %>%
        left_join(antigen_df)%>%
        mutate(isotype=factor(isotype,levels=c("IgG","IgA","IgM")))

vax_wide_RAU <- analysisData$RAU_data$wide_data %>%
        filter(status=="Vaccinee")


vax_long_NetMFI<- analysisData$netMFI_data$long_data%>%
        filter(status=="Vaccinee")%>%
        left_join(antigen_df)%>%
        mutate(isotype=factor(isotype,levels=c("IgG","IgA","IgM")))


casecon_analysis <- casecon_analysis %>% filter(status!="Vaccinee")%>%
        left_join(antigen_df) %>%
        mutate(isotype=factor(isotype,levels=c("IgG","IgA","IgM")))

wide_analysis <- analysisData$RAU_data$wide_data %>% filter(status!="Vaccinee")


casecon_analysisMFI <- analysisData$netMFI_data$long_data %>% filter(status!="Vaccinee")%>%
        left_join(antigen_df)%>%
        mutate(isotype=factor(isotype,levels=c("IgG","IgA","IgM")))

wide_analysisMFI <- analysisData$netMFI_data$wide_data %>% filter(status!="Vaccinee")


#create dataset for distinguishing classes of cases and vaccienes
multi_class_data <- analysisData$netMFI_data$long_data %>%
        filter(test_type=="Luminex")%>% 
        filter(age>=18) %>%
        #limit vaccinees to only baseline and "recent samples"
        filter(!(status=="Vaccinee" & (!day %in% c(0,7,21,44)))) %>%
        mutate(vaxinf_class =case_when(
                status=="Case" & day !=2 & day<=120 ~ "Recently Infected",
                status=="Vaccinee" & day %in% c(7,21,44) ~ "Recently Vaccinated",
                TRUE ~ "Neither"
        )) %>%
        left_join(antigen_df) 
