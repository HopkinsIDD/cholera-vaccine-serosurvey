rm(list=ls())

source("code/R/packages.R")
source("code/R/utils.R")
source("code/R/load-data.R")


#variable list
var_list <- list(
        `Reduced IgG Panel` = c("RAU_IgG_CtxB","RAU_IgG_InabaOSPBSA","RAU_IgG_OgawaOSPBSA"),
        `All Variables`= data.frame(marker=colnames(final_wide)) %>%
                filter(str_detect(marker,"CtxB|Ogawa|Inaba|TcpA|IgG_O139")) %>%
                unlist()  %>% unname()
)

formula_list <- list(
        reduced_2class = make_formula("inf_200",var_list$`Reduced IgG Panel`),
        all_2class = make_formula("inf_200",var_list$`All Variables`),
        
        reduced_3class= make_formula("three_class",var_list$`Reduced IgG Panel`),
        all_3class= make_formula("three_class",var_list$`All Variables`)
        
)



#get data sets for the analysis
data_list <- list(
        all_fitdata = final_wide %>%
                filter(cohort %in% c("SMIC/PIC","BGD Vaccinee","HTI Vaccinee"))%>%
                mutate(three_class=case_when(
                        status== "Case" & day_actual <= 200 & day_actual >= 5 ~ "RecentlyInfected",
                        status== "Vaccinee" & day_actual <= 200 & day_actual >= 5 ~ "RecentlyVaccinated",
                        TRUE ~ "Neither"
                )) %>%
                mutate(three_class=factor(three_class))%>%
                addInfectionWindow(end_window=200),
        case_fitdata = final_wide %>%
                filter(cohort %in% c("SMIC/PIC"))%>%
                addInfectionWindow(end_window=200)
)



#important variables
tree_param <- 1000
t <-200


loocv_list <- list()
loocv_df <- data.frame()

for(v in 1:length(formula_list)){
        for (q in names(data_list)){

        #dont do 3 class models when no vaccinees are in the dataset
         if(q=="case_fitdata" & names(formula_list)[v] %in% c("reduced_3class","all_3class")) next
                
        base_data <- data_list[[q]]
        ids <- unique(base_data$id)

        for(i in ids){


                marks <- formula_list[[v]]

                cat("Data: ",q, "Variables: ",names(formula_list)[v],"Individiual: ",i,"\n")

                fit_data <- base_data %>%
                        filter(id!=i) %>%
                        addInfectionWindow(end_window=t)

                outcome <- paste0("inf_",t)

                # 1. run 100 models to find a cutoff
                #keep 50% fit and 50% to test
                cut_data <- data.frame()

                for(j in 1:5){ #normally 1:100

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


                        #fit rf model
                        inside_fit <- ranger::ranger(formula=formula_list[[v]],
                                                     data=inside_df,
                                                     num.trees = 500,
                                                     case.weights = inside_w,
                                                     replace=TRUE
                        )

                        class_levels <- inside_fit$predictions %>% levels()

                        if(length(class_levels)==2){

                                #find the cutoff for each
                                performance <-  outside_df %>%
                                        filter(cohort %in% c("SMIC/PIC")) %>%
                                        addInfectionWindow(end_window=t) %>%
                                        getPerf_ranger(inside_fit,paste0("inf_",t))
                        }

                        if(length(class_levels)==3){

                                tmp_performance_df <- outside_df %>%
                                        filter(cohort %in% c("SMIC/PIC")) %>%
                                        addInfectionWindow(end_window=t)

                                tmp_predict <- predict(inside_fit,data=tmp_performance_df,predict.all=TRUE)

                                performance <- getPerf_ranger_new2(tmp_performance_df,
                                                                   preds= rowMeans(tmp_predict[[1]]==2),#decimal (numeric)
                                                                   truths= tmp_performance_df$inf_200#1 0 (numeric)
                                                                   )

                        }



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
                                          formula=formula_list[[v]],
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

                # test_missing <- rowSums(is.na(test_data[formula_list[[v]]])) >0
                #
                # test_data <- test_data[!test_missing,]

                if(length(class_levels)==2){


                        test_data["truth"] <- test_data[outcome]
                        out_preds <-predict(this_fit,
                                            test_data,
                                            predict.all=TRUE)

                        test_pred <- out_preds[[1]] %>% rowMeans()-1

                }
                if(length(class_levels)==3){

                        test_data["truth"] <- test_data[outcome]
                        out_preds <-predict(this_fit,
                                            test_data,
                                            predict.all=TRUE)

                        test_pred <- rowMeans(out_preds[[1]]==2)

                        specific_preds <- out_preds[[1]] %>% t() %>%
                                as.data.frame() %>%
                                gather(V_id,value) %>%
                                count(V_id,value) %>%
                                mutate(value=factor(value,levels=c(1,2,3),labels=class_levels)) %>%
                                spread(value,n)

                }


                #store the data
                fpr_df <- test_data %>%
                        mutate(prediction=test_pred) %>%
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
                        mutate(variables=names(formula_list)[v]) %>%
                        mutate(end_window=as.character(t)) %>%
                        mutate(base_data_name=q)


                if(length(class_levels)==3){

                        fpr_df <- bind_cols(fpr_df,specific_preds)

                }

                out_tmp <- list(
                        cutoffs=cut_data,
                        train=fit_data,
                        fit=this_fit,
                        preds=fpr_df
                )

                loocv_list[[names(formula_list)[v]]][[q]][[i]] <- out_tmp
                loocv_df <- bind_rows(loocv_df,fpr_df)


                write_rds(loocv_df, "data/generated_data/analysis_objects/loocv/loocv_df_new.rds")

                
                
}}}



loocv_df %>%
        filter(status=="Vaccinee")%>%
        count(truth,spec95_seropos,variables,day) %>%
        spread(spec95_seropos,n) %>%
        rename(positive=`1`,
               negative=`0`) %>%
        mutate(positive=ifelse(is.na(positive),0,positive),
               negative=ifelse(is.na(negative),0,negative),
        ) %>%
        mutate(perpos=positive/(positive+negative)) %>%
        ggplot(aes(x=factor(day),y=perpos,fill=variables)
        )+
        geom_col(position=position_dodge2(width=0.1))


loocv_df %>%
        filter(str_detect(variables,"3class")) %>%
        mutate(across(Neither:RecentlyVaccinated,.fns = ~ifelse(is.na(.x),
                                                                0,.x
                                                                ))) %>%
        mutate(prediction=case_when(
                Neither >= RecentlyInfected & Neither >=RecentlyVaccinated ~"Neither",
                RecentlyInfected > Neither & RecentlyInfected >= RecentlyVaccinated ~"RecentlyInfected",
                RecentlyVaccinated > Neither & RecentlyVaccinated > RecentlyInfected ~"RecentlyVaccinated"
        )) %>%
        mutate(three_class=factor(three_class),
                  prediction=factor(prediction)
        ) %>%
        count(variables,three_class,prediction,.drop=FALSE) %>%
        group_by(variables,three_class) %>%
        mutate(proportion=n/sum(n)) %>%
        mutate(three_class=str_replace(three_class,"Recently","Recently\n"))%>%
        mutate(three_class=glue::glue("{three_class}\n(n={sum(n)})"))%>%
        mutate(prediction=str_replace(prediction,"Recently","Recently "))%>%
        ggplot(aes(x=three_class,y=prediction))+
        geom_tile(aes(fill=100*proportion))+
        facet_wrap(.~variables)+
        geom_text(aes(label=paste0(100*round(proportion,2),"%")))+
        scale_fill_distiller("Percent\nClassified",palette = "Oranges",
                             direction = 1
        )+
        xlab("True Class")+
        ylab("Predicted Class")+
        theme_cowplot()


        
loocv_df %>%
        filter(str_detect(variables,"3class")) %>%
        mutate(across(Neither:RecentlyVaccinated,.fns = ~ifelse(is.na(.x),
                                                                0,.x
        ))) %>%
        mutate(prediction=case_when(
                Neither >= RecentlyInfected & Neither >=RecentlyVaccinated ~"Neither",
                RecentlyInfected > Neither & RecentlyInfected >= RecentlyVaccinated ~"RecentlyInfected",
                RecentlyVaccinated > Neither & RecentlyVaccinated > RecentlyInfected ~"RecentlyVaccinated"
        )) %>%
        filter(three_class %in% c("RecentlyVaccinated","RecentlyInfected")) %>%
        group_by(day,variables,three_class)%>%
        summarize(n=n(),
                  vax_pred=mean(prediction=="RecentlyVaccinated"),
                  inf_pred=mean(prediction=="RecentlyInfected"),
                  neither_pred=mean(prediction=="Neither")
                  ) %>%
        gather(pred_class,proportion,-c(day,variables,three_class,n))%>%
        ggplot(aes(x=factor(day),y=proportion))+
        geom_col(aes(fill=pred_class))+
        facet_wrap(variables~three_class,scales = "free")


# loocv_df %>%
#         filter(str_detect(variables,"3class")) %>%
#         mutate(across(Neither:RecentlyVaccinated,.fns = ~ifelse(is.na(.x),
#                                                                 0,.x
#         ))) %>%
#         mutate(prediction=case_when(
#                 Neither >= RecentlyInfected & Neither >=RecentlyVaccinated ~"Neither",
#                 RecentlyInfected > Neither & RecentlyInfected >= RecentlyVaccinated ~"RecentlyInfected",
#                 RecentlyVaccinated > Neither & RecentlyVaccinated > RecentlyInfected ~"RecentlyVaccinated"
#         )) %>%
#         filter(three_class=="RecentlyInfected") %>%
#         group_by(day,variables)%>%
#         summarize(n=n(),
#                   inf_pred=mean(prediction=="RecentlyVaccinated")
#         ) %>%
#         ggplot(aes(x=factor(day),y=vax_pred))+
#         geom_col(aes(fill=variables),position=position_dodge2(0.1))


loocv_df <- read_rds("data/generated_data/analysis_objects/loocv/loocv_df_new.rds")
options(mc.cores = 1)

# ### model time-varying fpr
# 
# tvfpr_fit <- list()
# tvfpr_df <- data.frame()
# 
# vars <- loocv_df$variables %>% unique()
# 
# 
# for(i in vars){
#         
#         cat("Window: ",200, " Population: ",i,"\n")
#         
#         # estimate time-varying sensitivity
#         tvfpr_tmp <- fit_tvfpr(data = loocv_df %>%
#                                        filter(status=="Vaccinee") %>%
#                                        filter(variables==i),
#                                end_window=200,
#                                seropos_type= "spec95_seropos",
#                                last_time = 365,
#                                curve="cubic"
#         )
#         
#         tvfpr_fit[[as.character(200)]][[i]] <- tvfpr_tmp
#         
#         tmp_df<- tvfpr_tmp$summary_df 
#         
#         if(!is.null(tmp_df)){
#                 tmp_df<- tmp_df %>%
#                         # mutate(variables="Reduced IgG Panel")%>%
#                         mutate(seropos_type="spec95_seropos") %>%
#                         mutate(end_window=200) %>%
#                         mutate(population=i)
#                 
#         }
#         
#         
#         tvfpr_df <-bind_rows(tvfpr_df,tmp_df)
#         
#         
#         
# }
# 
# 
# 
# ggplot()+
#         geom_line(data=tvfpr_df,
#                   aes(x=time,y=median,col=population))


#model the the time varying estimates using loocv data
loocv_tvfit_list <- list()
loocv_tvfit_df <- data.frame()

models_df <- distinct(loocv_df, variables,base_data_name,status) %>%
        filter(status %in% c("Case","Vaccinee"))


for(m in 1:nrow(models_df)){

                this_model <- models_df[m,]
                this_model_name <- models_df[m,c("variables","base_data_name")] %>% unlist() %>% paste0(collapse=" with ")
                this_data <- left_join(this_model,loocv_df, by = c("status", "variables", "base_data_name"))                 
                
                cat("Model: ",this_model_name, " Population: ",this_model$status,"\n")  
                
                # estimate time-varying sensitivity
                tvfpr_tmp <- fit_tvfpr(data = this_data,
                                       end_window=200,
                                       seropos_type= "spec95_seropos",
                                       last_time = if(this_model$status=="Case") 900 else 365,
                                       curve="cubic"
                )
                
                loocv_tvfit_list[[this_model$variables]][[this_model$base_data_name]][["200-day"]][[this_model$status]] <- tvfpr_tmp
                
                tmp_df<- tvfpr_tmp$summary_df 
                
                if(!is.null(tmp_df)){
                        tmp_df<- tmp_df %>%
                                mutate(variables=this_model$variables)%>%
                                mutate(base_data_name=this_model$base_data_name)%>%
                                mutate(seropos_type="spec95_seropos") %>%
                                mutate(end_window=200) %>%
                                mutate(population=this_model$status)
                        
                }
                
                
                loocv_tvfit_df <-bind_rows(loocv_tvfit_df,tmp_df)
                
                
}


write_rds(loocv_tvfit_list,
          paste0("data/generated_data/analysis_objects/loocv/","loocv_tvfit_list_new.rds")
)
write_rds(loocv_tvfit_df,
          paste0("data/generated_data/analysis_objects/loocv/","loocv_tvfit_df_new.rds")
)

#model the specificity 
loocv_df <- read_rds("data/generated_data/analysis_objects/loocv/loocv_df_new.rds")
options(mc.cores = 1)


loocv_spec_list <- list()
loocv_spec_df <- data.frame()

models_spec_df <- distinct(loocv_df, variables,base_data_name)
data_spec_list <- list(
                `Outside Window` = filter(loocv_df, cohort=="SMIC/PIC",inf_200==0),
                Vaccinee = filter(loocv_df, status=="Vaccinee")
                )

for(m in 1:nrow(models_spec_df)){
        for(set in names(data_spec_list)){
                this_model <- models_spec_df[m,]
                this_model_name <- models_spec_df[m,] %>% unlist() %>% paste0(collapse=" with ")
                this_data <- left_join(this_model,data_spec_list[[set]],by = c("variables", "base_data_name"))  
                
                if(nrow(this_data)>1){
                        
                        cat("Model: ",this_model_name, " Population: ",set,"\n")  
                        
                        # estimate specificity
                        spec_tmp <- fit_spec(this_data,
                                             end_window=200,
                                             seropos_type= "spec95_seropos")
                        
                        
                        loocv_spec_list[[this_model$variables]][[this_model$base_data_name]][["200-day"]][[set]] <- spec_tmp
                        
                        tmp_df<- spec_tmp$summary_df %>%
                                mutate(variables=this_model$variables)%>%
                                mutate(base_data_name=this_model$base_data_name)%>%
                                mutate(seropos_type="spec95_seropos") %>%
                                mutate(end_window=200) %>%
                                mutate(population=set)
                        
                        loocv_spec_df <-bind_rows(loocv_spec_df,tmp_df)
                        
                        
                }                
        }
}

write_rds(loocv_spec_list,
          paste0("data/generated_data/analysis_objects/loocv/","loocv_spec_list_new.rds")
)

write_rds(loocv_spec_df,
          paste0("data/generated_data/analysis_objects/loocv/","loocv_spec_df_new.rds")
)






