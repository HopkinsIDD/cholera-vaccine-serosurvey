path <- "source/final_code/vax_strategy/generated_rds/fpr-estimates/"
tree_param <- 1000
library(rstanarm)
library(tidybayes)

# #get variables into a list
# form_Luminex <- casecon_analysis %>%
#         filter(test_type=="Luminex") %>%
#         filter(!(antigen %in% c("CTHT","LTh","LTB","Flu","O139BSA"))) %>%
#         distinct(full_test) %>% pull(full_test) %>%
#         c("age","sex","blood")


var_list <- list(
        `Reduced IgG Panel` = c("RAU_IgG_CtxB","RAU_IgG_InabaOSPBSA","RAU_IgG_OgawaOSPBSA")
)


# var_list <- list(
#                  "All MBA (All isotypes)" = form_Luminex,
#                  "All MBA IgG"= form_Luminex[str_detect(form_Luminex,"IgG|age|sex|blood")],
#                  "All MBA IgA"= form_Luminex[str_detect(form_Luminex,"IgA|age|sex|blood")],
#                  "All MBA IgM"= form_Luminex[str_detect(form_Luminex,"IgM|age|sex|blood")]
# )
# 
# var_list[["Reduced Panel (All isotypes)"]]<- form_Luminex[str_detect(form_Luminex,"CtxB|OSP|age|sex|blood")]
# var_list[["Reduced Panel IgG"]]<- var_list$`All MBA IgG`[str_detect(var_list$`All MBA IgG`,"CtxB|OSP|age|sex|blood")]
# var_list[["Reduced Panel IgA"]]<- var_list$`All MBA IgA`[str_detect(var_list$`All MBA IgA`,"CtxB|OSP|age|sex|blood")]
# var_list[["Reduced Panel IgM"]]<- var_list$`All MBA IgM`[str_detect(var_list$`All MBA IgM`,"CtxB|OSP|age|sex|blood")]
# 
# var_list[["CtxB (All isotypes)"]]<- form_Luminex[str_detect(form_Luminex,"CtxB|age|sex|blood")]
# var_list[["CtxB IgG"]]<- var_list$`All MBA IgG`[str_detect(var_list$`All MBA IgG`,"CtxB|age|sex|blood")]
# var_list[["CtxB IgA"]]<- var_list$`All MBA IgA`[str_detect(var_list$`All MBA IgA`,"CtxB|age|sex|blood")]
# var_list[["CtxB IgM"]]<- var_list$`All MBA IgM`[str_detect(var_list$`All MBA IgM`,"CtxB|age|sex|blood")]
# 
# var_list[["OSPs (All isotypes)"]]<- form_Luminex[str_detect(form_Luminex,"OSP|age|sex|blood")]
# var_list[["OSPs IgG"]]<- var_list$`All MBA IgG`[str_detect(var_list$`All MBA IgG`,"OSP|age|sex|blood")]
# var_list[["OSPs IgA"]]<- var_list$`All MBA IgA`[str_detect(var_list$`All MBA IgA`,"OSP|age|sex|blood")]
# var_list[["OSPs IgM"]]<- var_list$`All MBA IgM`[str_detect(var_list$`All MBA IgM`,"OSP|age|sex|blood")]

#####  Fit RF models and make predictions -------------

fpr_fit_list <- list()
raw_df <- data.frame()

for(t in c(45,120,200,300)){
        # for(v in 1:length(var_list)){
                v <- 1        
                marks <- var_list[[v]]
                
                cat("Window: ",t, "Variables: ",names(var_list)[v],"\n")        
                
                
                # fit_data <- wide_analysis %>%
                #         addInfectionWindow(end_window=t)
                fit_data <- final_wide %>%
                        filter(cohort =="SMIC/PIC")%>%
                        addInfectionWindow(end_window=t)
                
                missing <- rowSums(is.na(fit_data[var_list[[v]]])) >0
                  
                fit_data <- fit_data[!missing,]      
                
                outcome <- paste0("inf_",t)
                
                # 1. run 100 models to find a cutoff
                #keep 50% fit and 50% to test 
                cut_data <- data.frame()
                
                for(j in 1:100){
                        
                        n_ind <-unique(fit_data$id) %>% length() 
                        
                        #define who is in and who is out
                        outside <-sample(fit_data$id,round(n_ind*0.5),replace=FALSE)
                        outside_df <- filter(fit_data,id %in% outside)
                        inside_df  <- filter(fit_data,!id %in% outside)%>%
                                addInfectionWindow(end_window=t)
                        
                        
                        # inside_w = NULL
                        
                        # if(weighted){
                        inside_w <- inside_df %>%
                                getWeight(end_window = t)  %>%
                                pull(weight)
                                
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
                        performance <- 
                                outside_df %>% addInfectionWindow(end_window=t) %>%
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
                # w= NULL
                # if(weighted){
                w <- fit_data%>%
                        getWeight(end_window = t)  %>%
                        pull(weight)
                        
                # }
                
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
                test_data <- final_wide %>%
                        filter(str_detect(cohort,"Vaccinee"))%>%
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
                        mutate(end_window=as.character(t))
                
                out_tmp <- list(
                        cutoffs=cut_data,
                        train=fit_data,
                        fit=this_fit,
                        preds=fpr_df
                )
                
                fpr_fit_list[[as.character(t)]][[names(var_list)[v]]] <- out_tmp
                raw_df <- bind_rows(raw_df,fpr_df)
                
        # }
}

raw_df <- raw_df %>% mutate(variables=factor(variables,levels=names(var_list)))

write_rds(fpr_fit_list,
          paste0(path,"fpr_fit_list.rds")
)

write_rds(raw_df,
          paste0(path,"raw_df.rds")
)

##### Estimate time-varying sensitivity  and specificity-------------

raw_df <- read_rds(
          paste0(path,"raw_df.rds")
)


options(mc.cores = 1)


tvfpr_fit <- list()
tvfpr_df <- data.frame()

for(t in c(45,120,200,300)){
        for(v in unique(raw_df$variables)){
                for(s in c(#"youden_seropos",
                           "spec90_seropos",
                           "spec95_seropos","spec99_seropos")){
                        
                        cat("Window: ",t, " Variables: ",v," Seropositivity: ",s,"\n")  
                        
                        #fit FPR to the adults
                        
                        # limit to particular window and variables
                        adult_reg_data <- raw_df %>%
                                filter(end_window==t) %>%
                                filter(variables==v) %>%
                                filter(age>=18) 
                        
                        # estimate time-varying sensitivity
                        adult_tvfpr_tmp <- fit_tvfpr(data = adult_reg_data,
                                                  end_window=t,
                                                  seropos_type= s,
                                                  last_time = 365,
                                                  curve="cubic"
                        )
                        
                        tvfpr_fit[[as.character(t)]][[v]][[s]][["adult"]] <- adult_tvfpr_tmp
                        
                        tmp_df<- adult_tvfpr_tmp$summary_df 
                        
                        if(!is.null(tmp_df)){
                                tmp_df<- tmp_df %>%
                                        mutate(variables=v)%>%
                                        mutate(seropos_type=s) %>%
                                        mutate(end_window=t) %>%
                                        mutate(age_group="Adult")
                                
                        }
                        
                        
                        tvfpr_df <-bind_rows(tvfpr_df,tmp_df)
                        
                        
                        # #fit FPR to the children
                        # 
                        # # limit to particular window and variables
                        # child_reg_data <- raw_df %>%
                        #         filter(end_window==t) %>%
                        #         filter(variables==v) %>%
                        #         filter(age<18) 
                        # 
                        # # estimate time-varying sensitivity
                        # child_tvfpr_tmp <- fit_tvfpr(data = child_reg_data,
                        #                              end_window=t,
                        #                              seropos_type= s,
                        #                              last_time = 21,
                        #                              curve="linear"
                        # )
                        # 
                        # tvfpr_fit[[as.character(t)]][[v]][[s]][["child"]] <- child_tvfpr_tmp
                        # 
                        # tmp_df<- child_tvfpr_tmp$summary_df 
                        # 
                        # if(!is.null(tmp_df)){
                        #         tmp_df<- tmp_df %>%
                        #                 mutate(variables=v)%>%
                        #                 mutate(seropos_type=s) %>%
                        #                 mutate(end_window=t) %>%
                        #                 mutate(age_group="Child")
                        #         
                        # }
                        # 
                        
                        # tvfpr_df <-bind_rows(tvfpr_df,tmp_df)
                        
                        
                        
                        
                }}}





write_rds(tvfpr_fit,
          paste0(path,"tvfpr_fit.rds")
)
write_rds(tvfpr_df,
          paste0(path,"tvfpr_df.rds")
)


# 
# vax_wide_RAU %>% distinct(id,day) %>%
#         ggplot(aes(x=factor(day),y=id))+
#                 geom_point()+
#                 geom_line(aes(group=id))
# 
# 
# 
# 
# 

# tvfpr_df <- read_rds(
#           paste0(path,"tvfpr_df.rds")
# )
# plot_variables <- c("Reduced Panel (All isotypes)", "Reduced Panel IgG","OSPs IgG","CtxB IgG")
# 
# tvfpr_df %>% filter(age_group=="Adult") %>%
#                 filter(variables %in% plot_variables) %>%
#                 mutate(variables=factor(variables,
#                                         levels=plot_variables))%>%
#                 mutate(seropos_type = case_when(
#                         # seropos_type == "youden_seropos" ~ "Youden Index",
#                         seropos_type == "spec99_seropos" ~ "99% Specificity",
#                         seropos_type == "spec95_seropos" ~ "95% Specificity",
#                         seropos_type == "spec90_seropos" ~ "90% Specificity"
#                 )) %>%
#                 mutate(seropos_type=factor(seropos_type,
#                                            levels=c(#"Youden Index",
#                                                     "90% Specificity",
#                                                     "95% Specificity",
#                                                     "99% Specificity"))) %>%
#                 ggplot(aes(x=time,y=median,col=variables,fill=variables))+
#                         geom_line()+
#                         geom_ribbon(aes(ymin=lb,ymax=ub),
#                                     alpha=0.1,col=NA
#                         )+
#                         geom_rug(data=raw_df %>% filter(age>=18) %>%
#                                 distinct(sample,day_actual),
#                         sides="b", alpha=0.1,
#                         aes(x=day_actual),inherit.aes = FALSE
#                         ) +
#                         geom_hline(data=data.frame(
#                                 line= c(0.1,0.05,0.01
#                                 ),
#                                 seropos_type=c(#"Youden Index",
#                                         "90% Specificity",
#                                         "95% Specificity",
#                                         "99% Specificity"
#                                 )
#                         )%>%
#                                 mutate(seropos_type=factor(seropos_type,
#                                                            levels=c(#"Youden Index",
#                                                                    "90% Specificity",
#                                                                    "95% Specificity",
#                                                                    "99% Specificity"))) ,
#                         aes(yintercept=line),lty=2
#                         )+
#         
#                 facet_grid(seropos_type~end_window)+
#        
#                 scale_y_continuous("False positivity rate")+
#                 scale_x_continuous("Days since first dose of vaccination",
#                                    # breaks=c(0,90,180,270,365)
#                                    breaks=c(0,45,120,200,300,365)
#                 )+
#                 scale_color_brewer(palette="Set2","Serological Markers")+
#                 scale_fill_brewer(palette="Set2","Serological Markers")+
#                 theme_bw()+
#                 theme(  panel.grid.minor=element_blank(),
#                         axis.text.x = element_text(angle=45,hjust=1),
#                         legend.position = "bottom"
#                         )
# 
# #for adults remove the OSP (could add ctxB IgA if you want to boost sensitivity)
# 
# 
# tvfpr_df %>% filter(age_group=="Child") %>%
#         filter(variables %in% plot_variables) %>%
#         mutate(variables=factor(variables,
#                                 levels=plot_variables))%>%
#         mutate(seropos_type = case_when(
#                 seropos_type == "youden_seropos" ~ "Youden Index",
#                 seropos_type == "spec99_seropos" ~ "99% Specificity",
#                 seropos_type == "spec95_seropos" ~ "95% Specificity",
#                 seropos_type == "spec90_seropos" ~ "90% Specificity"
#         )) %>%
#         mutate(seropos_type=factor(seropos_type,
#                                    levels=c("Youden Index",
#                                             "90% Specificity",
#                                             "95% Specificity",
#                                             "99% Specificity"))) %>%
#         ggplot(aes(x=time,y=median,col=variables,fill=variables))+
#         geom_line()+
#         geom_ribbon(aes(ymin=lb,ymax=ub),
#                     alpha=0.1,col=NA
#         )+
#         geom_rug(data=raw_df %>% filter(age<18) %>%
#                          distinct(sample,day),
#                  sides="b", alpha=0.1,
#                  aes(x=day),inherit.aes = FALSE
#         ) +
#         scale_y_continuous("False positivity rate")+
#         
#         scale_x_continuous("Days since first dose of vaccination",
#                            breaks=c(0,7,21)
#         )+
#         scale_color_brewer(palette="Set2","Markers")+
#         scale_fill_brewer(palette="Set2","Markers")+
#         theme_bw()+
#         theme(  panel.grid.minor=element_blank())+
#         
#         facet_grid(seropos_type~end_window)+
#         geom_hline(data=data.frame(
#                 line= c(NA,0.1,0.05,0.01),
#                 seropos_type=c("Youden Index",
#                                "90% Specificity",
#                                "95% Specificity",
#                                "99% Specificity")
#         )%>%
#                 mutate(seropos_type=factor(seropos_type,
#                                            levels=c("Youden Index",
#                                                     "90% Specificity",
#                                                     "95% Specificity",
#                                                     "99% Specificity"))) ,
#         aes(yintercept=line),lty=2
#         )







        