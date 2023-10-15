rm(list=ls())

source("code/R/packages.R")
source("code/R/utils.R")
source("code/R/load-data.R")


#variable list
var_list <- list(
        `Reduced IgG Panel` = c("RAU_IgG_CtxB","RAU_IgG_InabaOSPBSA","RAU_IgG_OgawaOSPBSA", "age"),
        `All Variables`= data.frame(marker=colnames(final_wide)) %>%
                filter(str_detect(marker,"CtxB|Ogawa|Inaba|TcpA|IgG_O139")) %>%
                unlist() %>% c("age") %>% unname()
)



#get data sets for the analysis
data_list <- list(
        smicpic_fitdata =  final_wide %>%
                filter(cohort %in% c("SMIC/PIC")),
        all_fitdata = final_wide %>%
                filter(cohort %in% c("SMIC/PIC","BGD Vaccinee","HTI Vaccinee")),
        bgd_fitdata = final_wide %>%
                filter(cohort %in% c("SMIC/PIC","BGD Vaccinee"))
)



#important variables
tree_param <- 1000
t <-200


loocv_list <- list()
raw_df <- data.frame()

for(v in 1:length(var_list)){
        for (q in names(data_list)){
        
        base_data <- data_list[[q]]
        ids <- unique(base_data$id) 
                
        for(i in ids){
                
                
                marks <- var_list[[v]]
                
                cat("Data: ",q, "Variables: ",names(var_list)[v],"Individiual: ",i,"\n")        
                
                fit_data <- base_data %>%
                        filter(id!=i) %>%
                        addInfectionWindow(end_window=t) 
                
                outcome <- paste0("inf_",t)
                
                # 1. run 100 models to find a cutoff
                #keep 50% fit and 50% to test 
                cut_data <- data.frame()
                
                for(j in 1:100){
                        

                        #define who is in and who is out
                        outside <-fit_data %>%
                                distinct(cohort,id) %>%
                                group_by(cohort) %>% 
                                slice_sample(prop=0.5,replace = FALSE) %>%
                                pull(id)
                        outside_df <- filter(fit_data,id %in% outside)
                        inside_df  <- filter(fit_data,!id %in% outside)%>%
                                addInfectionWindow(end_window=t)
                        
                        
                        inside_w = NULL
                        
                        # if(weighted){
                        # inside_w <- inside_df %>%
                        #         getWeight(end_window = t)  %>%
                        #         pull(weight)
                        
                        # }
                        
                        #fit rf model
                        inside_fit <- ranger::ranger(formula=make_formula(paste0("inf_",t),
                                                                          var_list[[v]]),
                                                     data=inside_df,
                                                     num.trees = 500,
                                                     case.weights = inside_w,
                                                     replace=TRUE
                        )
                        
                        
                        #find the cutoff for each 
                        performance <-  outside_df %>% addInfectionWindow(end_window=t) %>%
                                getPerf_ranger(inside_fit,paste0("inf_",t)#,
                                               # spec=train_spec_cutoff
                                )
                        cut <- distinct(performance$data,youden_cutoff,
                                        spec99_cutoff,
                                        spec95_cutoff,
                                        spec90_cutoff
                        )
                        cut_data <- bind_rows(cut_data,cut)
                }
                
                #find cutoffs
                median_cut <- cut_data %>%
                        summarize(n =n(),
                                  youden_cutoff=median(youden_cutoff),
                                  spec99_cutoff=median(spec99_cutoff),
                                  spec95_cutoff=median(spec95_cutoff),
                                  spec90_cutoff=median(spec90_cutoff)
                        )
                
                print(median_cut)
                
                # 2. fit the model, predict, and apply cutoff 
                #fit the model
                
                #set weights if needed
                
                w <- NULL
                
                # w <- fit_data%>%
                #         getWeight(end_window = t)  %>%
                #         pull(weight)
                
                
                
                this_fit<- ranger::ranger(data=fit_data,
                                          formula=make_formula(outcome,var_list[[v]]),
                                          num.trees= tree_param,
                                          case.weights = w,
                                          replace = TRUE
                )
                
                #
                #calculate the individual's samples seropositivity
                # test_data <- vax_wide_RAU %>%
                #         addInfectionWindow(end_window=t)
                test_data <- base_data %>%
                        filter(id==i) %>%
                        addInfectionWindow(end_window=t)
                
                test_missing <- rowSums(is.na(test_data[var_list[[v]]])) >0
                
                test_data <- test_data[!test_missing,]  
                
                
                test_data["truth"] <- test_data[outcome]
                out_preds <-predict(this_fit,
                                    test_data,
                                    predict.all=TRUE)
                
                #store the data
                fpr_df <- test_data %>%
                        mutate(prediction=out_preds[[1]] %>% rowMeans()-1) %>%
                        mutate(youden_cutoff=median_cut$youden_cutoff,
                               spec99_cutoff=median_cut$spec99_cutoff,
                               spec95_cutoff=median_cut$spec95_cutoff,
                               spec90_cutoff=median_cut$spec90_cutoff
                        ) %>%
                        mutate(
                                youden_seropos = ifelse(prediction>youden_cutoff,
                                                        1,0
                                ),
                                spec99_seropos = ifelse(prediction>spec99_cutoff,
                                                        1,0),
                                spec95_seropos = ifelse(prediction>spec95_cutoff,
                                                        1,0),
                                spec90_seropos = ifelse(prediction>spec90_cutoff,
                                                        1,0),
                                
                        ) %>%
                        mutate(variables=names(var_list)[v]) %>%
                        mutate(end_window=as.character(t)) %>%
                        mutate(base_data_name=q)
                
                out_tmp <- list(
                        cutoffs=cut_data,
                        train=fit_data,
                        fit=this_fit,
                        preds=fpr_df
                )
                
                loocv_list[[names(var_list)[v]]][[q]][[i]] <- out_tmp
                raw_df <- bind_rows(raw_df,fpr_df)
                
                
                write_rds(raw_df, "data/generated_data/analysis_objects/loocv_df.rds")
                
                
                
}}}



write_rds(raw_df, "data/generated_data/analysis_objects/loocv_df.rds")


# #get fit data
# bangladeshi_fitdata <- final_wide %>%
#         filter(cohort %in% c("BGD Vaccinee" , "SMIC/PIC"))  
# missing <-bangladeshi_fitdata[,var_list$`All Variables`] %>% is.na() %>% rowSums()
# bangladeshi_fitdata <- bangladeshi_fitdata[missing==0,]
# 
# 
# tree_param <- 1000
# t <-200
# ids <- unique(bangladeshi_fitdata$id) 
# 
# 
# loocv_list <- list()
# raw_df <- data.frame()
# 
# for(v in 1:length(var_list)){
#         for(i in ids){
# 
#                 
#                 marks <- var_list[[v]]
#                 
#                 cat("Window: ",t, "Variables: ",names(var_list)[v],"Individiual: ",i,"\n")        
#                 
#                 fit_data <- bangladeshi_fitdata %>%
#                         filter(id!=i) %>%
#                         addInfectionWindow(end_window=t) 
#                 
#                 outcome <- paste0("inf_",t)
#                 
#                 # 1. run 100 models to find a cutoff
#                 #keep 50% fit and 50% to test 
#                 cut_data <- data.frame()
#                 
#                 for(j in 1:100){
#                         
#                         n_ind <-unique(fit_data$id) %>% length() 
#                         
#                         #define who is in and who is out
#                         outside <-sample(fit_data$id,round(n_ind*0.5),replace=FALSE)
#                         outside_df <- filter(fit_data,id %in% outside)
#                         inside_df  <- filter(fit_data,!id %in% outside)%>%
#                                 addInfectionWindow(end_window=t)
#                         
#                         
#                         inside_w = NULL
#                         
#                         # if(weighted){
#                         # inside_w <- inside_df %>%
#                         #         getWeight(end_window = t)  %>%
#                         #         pull(weight)
#                         
#                         # }
#                         
#                         #fit rf model
#                         inside_fit <- ranger::ranger(formula=make_formula(paste0("inf_",t),
#                                                                           var_list[[v]]),
#                                                      data=inside_df,
#                                                      num.trees = 500,
#                                                      case.weights = inside_w,
#                                                      replace=TRUE
#                         )
#                         
#                         
#                         #find the cutoff for each 
#                         performance <- 
#                                 outside_df %>% addInfectionWindow(end_window=t) %>%
#                                 getPerf_ranger(inside_fit,paste0("inf_",t)#,
#                                                # spec=train_spec_cutoff
#                                 )
#                         cut <- distinct(performance$data,youden_cutoff,
#                                         spec99_cutoff,
#                                         spec95_cutoff,
#                                         spec90_cutoff
#                         )
#                         cut_data <- bind_rows(cut_data,cut)
#                 }
#                 
#                 #find cutoffs
#                 median_cut <- cut_data %>%
#                         summarize(n =n(),
#                                   youden_cutoff=median(youden_cutoff),
#                                   spec99_cutoff=median(spec99_cutoff),
#                                   spec95_cutoff=median(spec95_cutoff),
#                                   spec90_cutoff=median(spec90_cutoff)
#                         )
#                 
#                 print(median_cut)
#                 
#                 # 2. fit the model, predict, and apply cutoff 
#                 #fit the model
#                 
#                 #set weights if needed
#                 
#                 w <- NULL
#                 
#                 # w <- fit_data%>%
#                 #         getWeight(end_window = t)  %>%
#                 #         pull(weight)
#                 
#                 
#                 
#                 this_fit<- ranger::ranger(data=fit_data,
#                                           formula=make_formula(outcome,var_list[[v]]),
#                                           num.trees= tree_param,
#                                           case.weights = w,
#                                           replace = TRUE
#                 )
#                 
#                 #
#                 #calculate the individual's samples seropositivity
#                 # test_data <- vax_wide_RAU %>%
#                 #         addInfectionWindow(end_window=t)
#                 test_data <- bangladeshi_fitdata %>%
#                         filter(id==i) %>%
#                         addInfectionWindow(end_window=t)
#                 
#                 test_missing <- rowSums(is.na(test_data[var_list[[v]]])) >0
#                 
#                 test_data <- test_data[!test_missing,]  
#                 
#                 
#                 test_data["truth"] <- test_data[outcome]
#                 out_preds <-predict(this_fit,
#                                     test_data,
#                                     predict.all=TRUE)
#                 
#                 #store the data
#                 fpr_df <- test_data %>%
#                         mutate(prediction=out_preds[[1]] %>% rowMeans()-1) %>%
#                         mutate(youden_cutoff=median_cut$youden_cutoff,
#                                spec99_cutoff=median_cut$spec99_cutoff,
#                                spec95_cutoff=median_cut$spec95_cutoff,
#                                spec90_cutoff=median_cut$spec90_cutoff
#                         ) %>%
#                         mutate(
#                                 youden_seropos = ifelse(prediction>youden_cutoff,
#                                                         1,0
#                                 ),
#                                 spec99_seropos = ifelse(prediction>spec99_cutoff,
#                                                         1,0),
#                                 spec95_seropos = ifelse(prediction>spec95_cutoff,
#                                                         1,0),
#                                 spec90_seropos = ifelse(prediction>spec90_cutoff,
#                                                         1,0),
#                                 
#                         ) %>%
#                         mutate(variables=names(var_list)[v]) %>%
#                         mutate(end_window=as.character(t))
#                 
#                 out_tmp <- list(
#                         cutoffs=cut_data,
#                         train=fit_data,
#                         fit=this_fit,
#                         preds=fpr_df
#                 )
# 
#                 loocv_list[[names(var_list)[v]]][[i]] <- out_tmp
#                 raw_df <- bind_rows(raw_df,fpr_df)
#                 
#                 
#                 
#                 
#                 
# }}
# 
# 
# 
# write_rds(raw_df, "data/generated_data/loocv_df.rds")






# for(v in 1:length(var_list)){
#         
#         marks <- var_list[[v]]
#         
#         cat("Window: ",t, "Variables: ",names(var_list)[v],"\n")        
#         
#         fit_data <- final_wide %>%
#                 filter(cohort %in% c("BGD Vaccinee" , "SMIC/PIC"))%>% 
#                 mutate(day_actual=day)%>%
#                 addInfectionWindow(end_window=200) 
#         
#         outcome <- paste0("inf_",t)
#         
#         # 1. run 100 models to find a cutoff
#         #keep 50% fit and 50% to test 
#         cut_data <- data.frame()
#         
#         for(j in 1:100){
#                 
#                 n_ind <-unique(fit_data$id) %>% length() 
#                 
#                 #define who is in and who is out
#                 outside <-sample(fit_data$id,round(n_ind*0.5),replace=FALSE)
#                 outside_df <- filter(fit_data,id %in% outside)
#                 inside_df  <- filter(fit_data,!id %in% outside)%>%
#                         addInfectionWindow(end_window=t)
#                 
#                 
#                 # inside_w = NULL
#                 
#                 # if(weighted){
#                 inside_w <- inside_df %>%
#                         getWeight(end_window = t)  %>%
#                         pull(weight)
#                 
#                 # }
#                 
#                 #fit rf model
#                 inside_fit <- ranger::ranger(formula=make_formula(paste0("inf_",t),
#                                                                   var_list[[v]]),
#                                              data=inside_df,
#                                              num.trees = 500,
#                                              case.weights = inside_w,
#                                              replace=TRUE
#                 )
#                 
#                 
#                 #find the cutoff for each 
#                 performance <- 
#                         outside_df %>% addInfectionWindow(end_window=t) %>%
#                         getPerf_ranger(inside_fit,paste0("inf_",t)#,
#                                        # spec=train_spec_cutoff
#                         )
#                 cut <- distinct(performance$data,youden_cutoff,
#                                 spec99_cutoff,
#                                 spec95_cutoff,
#                                 spec90_cutoff
#                 )
#                 cut_data <- bind_rows(cut_data,cut)
#         }
#         
#         #find cutoffs
#         median_cut <- cut_data %>%
#                 summarize(n =n(),
#                           youden_cutoff=median(youden_cutoff),
#                           spec99_cutoff=median(spec99_cutoff),
#                           spec95_cutoff=median(spec95_cutoff),
#                           spec90_cutoff=median(spec90_cutoff)
#                 )
#         
#         print(median_cut)
#         
#         # 2. fit the model, predict, and apply cutoff 
#         #fit the model
#         
#         #set weights if needed
#         w <- fit_data%>%
#                 getWeight(end_window = t)  %>%
#                 pull(weight)
#         
#         
#         
#         this_fit<- ranger::ranger(data=fit_data,
#                                   formula=make_formula(outcome,var_list[[v]]),
#                                   num.trees= tree_param,
#                                   case.weights = w,
#                                   replace = TRUE
#         )
#         
#         #
#         #calculate the individual's samples seropositivity
#         # test_data <- vax_wide_RAU %>%
#         #         addInfectionWindow(end_window=t)
#         test_data <- final_wide %>%
#                 filter(cohort %in% c("BGD Vaccinee" , "SMIC/PIC"))%>% 
#                 addInfectionWindow(end_window=t)
#         
#         test_missing <- rowSums(is.na(test_data[var_list[[v]]])) >0
#         
#         test_data <- test_data[!test_missing,]  
#         
#         
#         test_data["truth"] <- test_data[outcome]
#         out_preds <-predict(this_fit,
#                             test_data,
#                             predict.all=TRUE)
#         
#         #store the data
#         fpr_df <- test_data %>%
#                 mutate(prediction=out_preds[[1]] %>% rowMeans()-1) %>%
#                 mutate(youden_cutoff=median_cut$youden_cutoff,
#                        spec99_cutoff=median_cut$spec99_cutoff,
#                        spec95_cutoff=median_cut$spec95_cutoff,
#                        spec90_cutoff=median_cut$spec90_cutoff
#                 ) %>%
#                 mutate(
#                         youden_seropos = ifelse(prediction>youden_cutoff,
#                                                 1,0
#                         ),
#                         spec99_seropos = ifelse(prediction>spec99_cutoff,
#                                                 1,0),
#                         spec95_seropos = ifelse(prediction>spec95_cutoff,
#                                                 1,0),
#                         spec90_seropos = ifelse(prediction>spec90_cutoff,
#                                                 1,0),
#                         
#                 ) %>%
#                 mutate(variables=names(var_list)[v]) %>%
#                 mutate(end_window=as.character(t))
#         
#         out_tmp <- list(
#                 cutoffs=cut_data,
#                 train=fit_data,
#                 fit=this_fit,
#                 preds=fpr_df
#         )
#         
#         fpr_fit_list[[as.character(t)]][[names(var_list)[v]]] <- out_tmp
#         raw_df <- bind_rows(raw_df,fpr_df)
#         
# }












# LOOCV

# 
# 
# ### all markers
# 
# formula_all <- make_formula("inf_200",all_var)
# fit_all <- runCaretLOOCV(formula_all, dat=bangladeshi_fitdata)
# 
# 
# 
# 
#         ggplot()+ 
#         plotROC::geom_roc(aes(m = value, d = truth),labels=FALSE,n.cuts=0)+ 
#         geom_abline(lty=3,col="grey50")+
#         xlab("False Positivity Rate")+
#         ylab("True Positivity Rate")+
#         # scale_color_manual("Markers",values=c("red","darkblue"))+
#         theme_cowplot()+
#         theme(legend.position = "bottom")
#         
# #three markers
# 
# three_var <- data.frame(marker=colnames(final_wide)) %>%
#         filter(str_detect(marker,"CtxB|Ogawa|Inaba")) %>%
#         filter(str_detect(marker,"IgG")) %>% 
#         unlist() %>% unname()
# 
# 
# formula_three <- make_formula("inf_200",three_var)
# 
# fit_three <- runCaretLOOCV(formula_three, dat=bangladeshi_fitdata)
# 
# 
# 
# all_roc_df <- fit_all$trainingData %>% 
#         mutate(rowIndex=1:n()) %>%
#         left_join(fit_all$pred, by = "rowIndex")%>%
#         select(obs,RecentlyInfected,NotRecentlyInfected)  %>%
#         gather(Class,value,-c(obs)) %>%
#         mutate(truth=ifelse(obs==Class,1,0)) %>%
#         mutate(Class = str_replace(Class,"Recently","Recently ")) %>% 
#         mutate(markers="all")
# 
# 
# 
# three_roc_df <- fit_three$trainingData %>% 
#         mutate(rowIndex=1:n()) %>%
#         left_join(fit_three$pred, by = "rowIndex")%>%
#         select(obs,RecentlyInfected,NotRecentlyInfected)  %>%
#         gather(Class,value,-c(obs)) %>%
#         mutate(truth=ifelse(obs==Class,1,0)) %>%
#         mutate(Class = str_replace(Class,"Recently","Recently ")) %>% 
#         mutate(markers="three")
# 
# 
# roc_df <- bind_rows(
#         all_roc_df,three_roc_df) 





