
final_df <- read_rds("data/generated_data/sample_measurements.rds")  %>%
        mutate(day_actual=day)



df_list <-list(
        
        cases_child_BGD = final_df %>%
                filter(
                        status=="Case",
                        age<10
                ) ,
        
        
        cases_adult_BGD = final_df %>%
                filter(
                        status=="Case",
                        age>=10
                ) ,
        
        vax_child_BGD = final_df %>%
                filter(
                        cohort=="BGD Vaccinee",
                        age<10
                ),
        
        vax_adult_BGD = final_df %>%
                filter(
                        cohort=="BGD Vaccinee",
                        age<10
                ),
        
        vax_adult_HTI = final_df %>%
                filter(
                        cohort=="HTI Vaccinee",
                        age>=18
                )
)


exp_model <- select_univariate_model("exponential")
markers <- final_df %>% distinct(marker) %>%
                filter(str_detect(marker,"Ogawa|Inaba|O139|TcpA|CtxB")) %>%
                pull(marker)

path <- "source/final_code/vax_strategy/generated_rds/exponential_fits/group_fits/"


for(group in names(df_list)){
        for(m in markers[2:15]) {
                
                tmp_value <-df_list[[group]] %>% select(cohort,status,id, sample, day_actual,age,sex,marker,RAU_value) %>%
                        mutate(RAU_value=log10(1/RAU_value)) %>%
                        spread(marker,RAU_value)
                tmp_censor <-df_list[[group]] %>% select(cohort,status,id, sample, day_actual,age,sex,marker,RAU_censor) %>%
                        spread(marker,RAU_censor)
                
                tmp_fit <- stanfit_decay_univariate(
                        data = tmp_value,
                        cens_df = tmp_censor,
                        marker =m,
                        delay=5, #if d>=delay, boost
                        # model_type,
                        iterations=2000,
                        # luminex_data_type = "RAU",
                        model=exp_model,
                        chains=4,
                        warmup=1000
                )
                
                
                
                write_rds(tmp_fit,paste0(path,paste(group,m,"expmodel.rds",sep="_")))
                
                
        }
}



summary_df <- data.frame()

for(group in names(df_list)){
        for(m in markers) {
                
                cat(group,m,"\n")
                
                tmp_fit <- read_rds(paste0(path,paste(group,m,"expmodel.rds",sep="_")))

                
                summary_df <- bind_rows(summary_df,
                                        tmp_fit$summary %>% 
                                                mutate(group=group)
                )
                
                
                
        }
}





# #########---------------
# 
# #vaccinees
# adultvax_wide_netMFI <- analysisData$netMFI_data$wide_data %>%
#         filter(status=="Vaccinee") %>%
#         filter(age>=18) 
# 
# adultvax_wide_netMFIcens <- analysisData$netMFI_data$wide_censor %>%
#         filter(status=="Vaccinee") %>%
#         filter(age>=18)
# 
# childvax_wide_netMFI <- analysisData$netMFI_data$wide_data %>%
#         filter(status=="Vaccinee") %>%
#         filter(age<18) 
# 
# childvax_wide_netMFIcens <- analysisData$netMFI_data$wide_censor %>%
#         filter(status=="Vaccinee") %>%
#         filter(age<18)
# 
# #cases
# adultcase_wide_netMFI <- analysisData$netMFI_data$wide_data %>%
#         filter(status=="Case") %>%
#         filter(age>=18) 
# 
# adultcase_wide_netMFIcens <- analysisData$netMFI_data$wide_censor %>%
#         filter(status=="Case") %>%
#         filter(age>=18)
# 
# childcase_wide_netMFI <- analysisData$netMFI_data$wide_data %>%
#         filter(status=="Case") %>%
#         filter(age<18)
# 
# childcase_wide_netMFIcens <- analysisData$netMFI_data$wide_censor %>%
#         filter(status=="Case") %>%
#         filter(age<18)
# 
# 
# 
# 
# ##########--------------
# 
# 
# 
# markers_of_interest <-analysisData$netMFI_data$long_data %>%
#         # filter(antigen %in% c("CtxB","OgawaOSPBSA","InabaOSPBSA","O139BSA")) %>%
#         # filter(antigen %in% c("TcpA")) %>%
#         # filter(isotype!="IgG")%>%
#         filter(test_type=="Luminex") %>%
#         distinct(full_test)%>%
#         pull(full_test)
# 
# 
# populations <- c("Child Vaccinee","Adult Vaccinee","Child Case","Adult Case")
# 
# options(mc.cores = parallel::detectCores())
# 
# exp_model <- select_univariate_model("exponential")
# 
# fit_list <- list()
# fit_df_1000 <- data.frame()
# for(i in markers_of_interest){
#         
#         cat(i,"\n")
#         
#         #child vaccinees        
#         tmp_fit1 <- stanfit_decay_univariate(
#                 data= childvax_wide_netMFI,
#                 cens_df= childvax_wide_netMFIcens,
#                 marker=i,
#                 delay=5,
#                 covariates=NULL,
#                 iterations=3000,
#                 model= exp_model,
#                 luminex_data_type = "netMFI"
#         )
#         
#         tmp_df1 <- tmp_fit1$fit %>%
#                 spread_draws(foldchange,halflife) %>%
#                 mutate(status=populations[1]) %>%
#                 mutate(full_test=i)
#         
#         fit_list[[i]][[populations[1]]] <- tmp_fit1
#         
#         
#         #adult vaccinees        
#         tmp_fit2 <- stanfit_decay_univariate(
#                 data= adultvax_wide_netMFI,
#                 cens_df= adultvax_wide_netMFIcens,
#                 marker=i,
#                 delay=5,
#                 covariates=NULL,
#                 iterations=3000,
#                 model= exp_model,
#                 luminex_data_type = "netMFI"
#         )
#         
#         tmp_df2 <- tmp_fit2$fit %>%
#                 spread_draws(foldchange,halflife) %>%
#                 mutate(status=populations[2]) %>%
#                 mutate(full_test=i)
#         
#         fit_list[[i]][[populations[2]]] <- tmp_fit2
#         
#         #child cases        
#         tmp_fit3 <- stanfit_decay_univariate(
#                 data= childcase_wide_netMFI,
#                 cens_df= childcase_wide_netMFIcens,
#                 marker=i,
#                 delay=5,
#                 covariates=NULL,
#                 iterations=3000,
#                 model= exp_model,
#                 luminex_data_type = "netMFI"
#         )
#         
#         tmp_df3 <- tmp_fit3$fit %>%
#                 spread_draws(foldchange,halflife) %>%
#                 mutate(status=populations[3]) %>%
#                 mutate(full_test=i)
#         
#         fit_list[[i]][[populations[3]]] <- tmp_fit3
#         
#         #adult cases        
#         tmp_fit4 <-  stanfit_decay_univariate(
#                 data= adultcase_wide_netMFI,
#                 cens_df= adultcase_wide_netMFIcens,
#                 marker=i,
#                 delay=5,
#                 covariates=NULL,
#                 iterations=3000,
#                 model= exp_model,
#                 luminex_data_type = "netMFI"
#         )
#         
#         tmp_df4 <- tmp_fit4$fit %>%
#                 spread_draws(foldchange,halflife) %>%
#                 mutate(status=populations[4]) %>%
#                 mutate(full_test=i)
#         
#         fit_list[[i]][[populations[4]]] <- tmp_fit4
#         
#         
#         #store 1000 iterations of each fit in a dataframe
#         fit_df_1000 <- bind_rows(fit_df_1000,
#                                  tmp_df1[round(seq(1,nrow(tmp_df1),length.out=1000)),],
#                                  tmp_df2[round(seq(1,nrow(tmp_df2),length.out=1000)),],
#                                  tmp_df3[round(seq(1,nrow(tmp_df3),length.out=1000)),],
#                                  tmp_df4[round(seq(1,nrow(tmp_df4),length.out=1000)),],
#         )
#         
#         
#         
#         
# }
# 
# write_rds(fit_list, "source/final_code/vax_strategy/generated_rds/exponential_fits/fit_list_netMFI.rds")
# write_rds(fit_df_1000, "source/final_code/vax_strategy/generated_rds/exponential_fits/fit_df_1000_netMFI.rds")
# 
# 



