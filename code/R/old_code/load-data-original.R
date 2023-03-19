analysisData<- LuminexRecPaperData()

# casecon_analysis <- analysisData$RAU_data$long_data
# wide_analysis <- analysisData$RAU_data$wide_data
# wide_censor <- analysisData$RAU_data$wide_censor

casecon_analysis <- analysisData$new_long_data
wide_analysis <- analysisData$new_wide_data
wide_censor <- analysisData$new_wide_censor

wide_analysis_nomissing <- filter(wide_analysis,
                                  !is.na(Vibriocidal_Ogawa),
                                  !is.na(ELISA_LPS_IgG)
)

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
        ))