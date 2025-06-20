rm(list=ls())
source("code/R/packages.R")
source("code/R/utils.R")

#load the wide dataset
final_wide <- read_rds("data/generated_data/analysis_data/final_wide.rds")


# path <- "source/final_code/vax_strategy/generated_rds/fpr-estimates/"
tree_param <- 1000
library(rstanarm)
library(tidybayes)


var_list <- list(
        # `Reduced IgG Panel and Age` = c("RAU_IgG_CtxB","RAU_IgG_InabaOSPBSA","RAU_IgG_OgawaOSPBSA", "age"),
        `Reduced IgG Panel` = c("RAU_IgG_CtxB","RAU_IgG_InabaOSPBSA","RAU_IgG_OgawaOSPBSA"),
        `All Variables`= data.frame(marker=colnames(final_wide)) %>%
                filter(str_detect(marker,"CtxB|Ogawa|Inaba|TcpA|IgG_O139")) %>%
                unlist()  %>% unname()
)



fpr_fit_list <- list()
raw_df <- data.frame()

for(t in c(45,120,200,300)){
        for(v in 1:length(var_list)){
        marks <- var_list[[v]]
        
        cat("Window: ",t, "Variables: ",names(var_list)[v],"\n")        
        
        fit_data <- final_wide %>%
                filter(cohort =="SMIC/PIC")%>% 
                mutate(day_actual=day)%>%
                addInfectionWindow(end_window=t) 
        
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
        w=NULL
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
        test_data <- final_wide %>%
                filter(cohort %in% c("BGD Vaccinee","HTI Vaccinee","ETEC Challenge"))

        test_missing <- rowSums(is.na(test_data[var_list[[v]]])) >0
        
        test_data <- test_data[!test_missing,]  
        
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
        
        }
}

write_rds(raw_df, "data/generated_data/analysis_objects/misclassification/misclassify_df.rds")


#model the the estimates using only "reduced IgG Panel"
options(mc.cores = 1)

tvfpr_fit <- list()
tvfpr_df <- data.frame()

raw_df_list <- list(
        
        # `Bangladeshi <10 years` = filter(raw_df, cohort=="BGD Vaccinee") %>%
        #         filter(age<10)%>% filter(variables=="Reduced IgG Panel"),
        # `Bangladeshi 10+ years`= filter(raw_df, cohort=="BGD Vaccinee") %>%
        #         filter(age>=10)%>% filter(variables=="Reduced IgG Panel"),
        # `Haitian 18+ years` = filter(raw_df, cohort=="HTI Vaccinee")%>%
        #         filter(variables=="Reduced IgG Panel") ,
        `All Vaccinees (Reduced IgG Panel)` = filter(raw_df, cohort %in% c("HTI Vaccinee","BGD Vaccinee"))%>%
                filter(variables=="Reduced IgG Panel"),
        `All Vaccinees (All Variables)` = filter(raw_df, cohort %in% c("HTI Vaccinee","BGD Vaccinee"))%>%
                filter(variables=="All Variables")
        # `All Vaccinees (Reduced IgG Panel [NO AGE])` = filter(raw_df, cohort %in% c("HTI Vaccinee","BGD Vaccinee"))%>%
        #         filter(variables=="Reduced IgG Panel (NO AGE)")
        
)

for(t in c(45,120,200,300)){
        for(dat in names(raw_df_list)){
                
                
                
                
                cat("Window: ",t, " Population: ",dat,"\n")  
                
                # estimate time-varying sensitivity
                tvfpr_tmp <- fit_tvfpr(data = raw_df_list[[dat]] %>%
                                               filter(end_window==t),
                                       end_window=t,
                                       seropos_type= "spec95_seropos",
                                       last_time = 365,
                                       curve="cubic"
                )
                
                tvfpr_fit[[as.character(t)]][[dat]] <- tvfpr_tmp
                
                tmp_df<- tvfpr_tmp$summary_df 
                
                if(!is.null(tmp_df)){
                        tmp_df<- tmp_df %>%
                                # mutate(variables="Reduced IgG Panel")%>%
                                mutate(seropos_type="spec95_seropos") %>%
                                mutate(end_window=t) %>%
                                mutate(population=dat)
                        
                }
                
                
                tvfpr_df <-bind_rows(tvfpr_df,tmp_df)
                
                
        }}





write_rds(tvfpr_fit,
          paste0("data/generated_data/analysis_objects/misclassification/","tvfpr_fit.rds")
)
write_rds(tvfpr_df,
          paste0("data/generated_data/analysis_objects/misclassification/","tvfpr_df.rds")
)












