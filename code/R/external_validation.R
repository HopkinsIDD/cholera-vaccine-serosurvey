
source("code/R/packages.R")
source("code/R/utils.R")
source("code/R/load-data.R")

path <- "source/final_code/vax_strategy/generated_rds/fpr-estimates/"
tree_param <- 1000
library(rstanarm)
library(tidybayes)


var_list <- list(
        `Reduced IgG Panel` = c("RAU_IgG_CtxB","RAU_IgG_InabaOSPBSA","RAU_IgG_OgawaOSPBSA", "age"),
        `All Variables`= data.frame(marker=colnames(final_wide)) %>%
                filter(str_detect(marker,"CtxB|Ogawa|Inaba|TcpA|IgG_O139")) %>%
                unlist() %>% c("age") %>% unname()
)

#get fit data
bangladeshi_fitdata <- final_wide %>%
        filter(cohort %in% c("BGD Vaccinee" , "SMIC/PIC"))  
missing <-bangladeshi_fitdata[,var_list$`All Variables`] %>% is.na() %>% rowSums()
bangladeshi_fitdata <- bangladeshi_fitdata[missing==0,]




fpr_fit_list <- list()
raw_df <- data.frame()
t <- 200

for(d in c("No vax included","Vax included")){
        for(v in 1:length(var_list)){
        marks <- var_list[[v]]
        
        
        cat("Window: ",t, "Variables: ",names(var_list)[v],"\n")        
        
        fit_data <- bangladeshi_fitdata %>%
                addInfectionWindow(end_window=t) 
        
        if(d == "No vax included"){
                fit_data <- fit_data %>%
                        filter(cohort== "SMIC/PIC") 
                
        } 
        
        
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
                
        
                
                #fit rf model
                inside_fit <- ranger::ranger(formula=make_formula(paste0("inf_",t),
                                                                  var_list[[v]]),
                                             data=inside_df,
                                             num.trees = 500,
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
        
              this_fit<- ranger::ranger(data=fit_data,
                                  formula=make_formula(outcome,var_list[[v]]),
                                  num.trees= tree_param,
                                  replace = TRUE
        )
        
        #
        #calculate the individual's samples seropositivity
        # test_data <- vax_wide_RAU %>%
        #         addInfectionWindow(end_window=t)
        test_data <- final_wide %>%
                filter(cohort %in% c("HTI Vaccinee","ETEC Challenge")) %>%
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
                mutate(end_window=as.character(t))%>%
                mutate(training_set=d)
        
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

write_rds(raw_df, "data/generated_data/external_df.rds")
