
#generate the questionairre status given sens spec
getQuestionairre<- function(x, sens, spec){
        
        if(x==1) out <- rbinom(1,1,sens)
        if(x==0) out <- rbinom(1,1,1-spec)
        
        return(out)
}

#generate seropositivity given sens spec
getSeroPos<- function(v,spec,sens){
        
        if(v==1) out <- rbinom(1,1,sens)
        if(v==0) out <- rbinom(1,1,1-spec)
        
        return(out)
}



#estimate an adjust seroincidence
adjSeropos <- function(sens,spec,n_seropos,n_total){
        max((n_seropos+n_total*(spec-1))/(sens+spec-1),0)
}



#simulate serosurveys in a partially vaccinated population 
# with a questionairre

simVaxSurvey <- function(n_sim=1,
                         n_pop = 2000,
                         unvaxxed_inc,
                         cov,
                         effectiveness =0,
                         sens_quest = 1,
                         spec_quest = 1, # specificity of questionnaire
                         seropos_param
){
        
        #set seropositivity params
        sens_novax = seropos_param$sens_novax
        sens_vax = seropos_param$sens_vax
        spec_vax = seropos_param$spec_vax
        spec_novax = seropos_param$spec_novax
        
        simulation_df <- data.frame()
        for(i in 1:n_sim){
                #simulate vax and recent infection status
                # this simulates the selection process
                vax <- rbinom(n_pop,1,cov) %>% sort()
                recent <- c(rbinom(length(vax[vax==0]),1,unvaxxed_inc), 
                            rbinom(length(vax[vax==1]),1,unvaxxed_inc*(1-effectiveness)))
                
                
                tmp_df  <- data.frame(
                        vax_status= vax,
                        infected_recently =recent
                ) %>%   
                        #get the simulation number
                        # mutate(simulation=i) %>%
                        mutate(simulation=paste0("sim_",str_pad(i,width=5,pad="0",side="left"))) %>%
                        #individual id
                        mutate(individual_id=paste0("id_",str_pad(1:n_pop,width=5,pad="0",side="left")))  %>%
                        #get questionnaire status
                        mutate(questionnaire=sapply(vax_status, 
                                                    FUN=getQuestionairre,
                                                    sens = sens_quest,
                                                    spec = spec_quest
                        )) %>%
                        #get seropositivity status
                        mutate(seropos= ifelse(vax_status ==0,
                                               sapply(infected_recently,
                                                      FUN=getSeroPos,
                                                      sens=sens_novax,
                                                      spec= spec_novax),
                                               sapply(infected_recently,
                                                      FUN=getSeroPos,
                                                      sens=sens_novax,
                                                      spec= spec_vax))
                        )
                
                
                simulation_df <- bind_rows(simulation_df,tmp_df)        
        }
        
        out <- list(
                ind_data = simulation_df %>% 
                        mutate(variables=seropos_param$title),
                n_pop = n_pop,
                unvaxxed_inc = unvaxxed_inc,
                cov = cov,
                effectiveness =effectiveness,
                sens_quest = sens_quest,
                spec_quest = spec_quest, # specificity of questionnaire
                seropos_param = seropos_param
                
        )
        
        return(out)
}



estimateSeroInc_vaxpop <- function(data,model,params,iter=1000){
        #introduce uncertainty in values with a multiplier
        multi=200
        
        stan_input <- list(
                
                #data
                N=nrow(data),
                S=data$seropos,
                Q=data$questionnaire,
                
                #priors
                
                # no prior on seroincidence
                
                #coverage (uniform prior can be changed)
                gamma_p1 = 1,
                gamma_p2 = 1,
                
                # sensitivity 
                # most people will have been infected before vaccination
                # vaccination might boost recently infected people to be more likely
                # to show positive
                
                alpha_p1=(params$sens_novax)*multi,
                alpha_p2=(1-params$sens_novax)*multi,
                
                #specificity
                beta_p1=(params$spec_novax)*multi,
                beta_p2=(1-params$spec_novax)*multi,
                
                epsilon_p1 = (params$spec_vax)*multi,
                epsilon_p2 = (1-params$spec_vax)*multi#,
                
                
                # phi_p1=0.9*multi,
                # phi_p2=0.1*multi,
                # 
                # psi_p1=0.9*multi,
                # psi_p2=0.1*multi,
                # 
                # #vaccine efficacy
                # nu_p1 = 50,
                # nu_p2 = 50
                
                
        )
        
        
        fit <- sampling(model, iter=iter,
                    data=stan_input)
        
        return(fit)
        
        
        
        
}


estimateSeroInc_vaxpop_simple <- function(data,model,params,
                                          multi=200,
                                          iter=1000){
        #introduce uncertainty in values with a multiplier
        
        stan_input <- list(
                
                #data
                S0Q0=data$S0Q0,
                S1Q0=data$S1Q0,
                S0Q1=data$S0Q1,
                S1Q1=data$S1Q1,
                
                #priors
                
                # no prior on seroincidence
                
                #coverage (uniform prior can be changed)
                gamma_p1 = 1,
                gamma_p2 = 1,
                
                # sensitivity 
                # most people will have been infected before vaccination
                # vaccination might boost recently infected people to be more likely
                # to show positive
                
                alpha_p1=(params$sens_novax)*multi,
                alpha_p2=(1-params$sens_novax)*multi,
                
                #specificity
                beta_p1=(params$spec_novax)*multi,
                beta_p2=(1-params$spec_novax)*multi,
                
                epsilon_p1 = (params$spec_vax)*multi,
                epsilon_p2 = (1-params$spec_vax)*multi#,
                
        
                
        )
        
        
        fit <- sampling(model, iter=iter,
                        data=stan_input)
        
        out <-list(
                data=data,
                fit=fit
                )
        
        return(out)
        
        
}


vaxSimSeroInc<- function(incidence, #different incidences
                         coverage, # different levels coverage
                         spec_vaxdf, #df describing specificity among vaccinated people at times,
                         sims
){
        
        
        
        reduced_param <- list(
                title= "IgG Reduced Panel",
                sens_novax = 0.38,
                sens_vax = 0.38,
                spec_novax = 0.95
                
        )
        
        vax_model <- stan_model("source/final_code/vax_strategy/stan/vax-model-simplified2.stan")
        
        #for each time
        for(i in 1:nrow(spec_vaxdf)){
                
                cat("time:",i,"\n")
                
                #get sensitivitiy and specificity parameters
                tmp_param <- reduced_param
                tmp_param$spec_vax <- spec_vaxdf$spec[i]
                
                #generate serosurvey data
                tmp_survey <-  simVaxSurvey(
                        n_sim=sims,
                        unvaxxed_inc = incidence,
                        cov = coverage,
                        sens_quest=1,
                        spec_quest = 1,
                        effectiveness =0,
                        seropos_param=tmp_param
                )
                
                #sum up serosurvey data
                sum_df <- tmp_survey$ind_data %>%
                        mutate(since2=spec_vaxdf$newtime[i],
                               cov=coverage,
                               unvaxxed_inc=incidence,
                        ) %>% 
                        group_by(cov,unvaxxed_inc,variables,simulation,since2) %>%
                        summarize(S1Q1=sum(seropos*questionnaire),
                                  S0Q1=sum((1-seropos)*questionnaire),
                                  S1Q0=sum(seropos*(1-questionnaire)),
                                  S0Q0=sum((1-seropos)*(1-questionnaire)),
                                  study_seropos=mean(seropos),
                                  study_inc=mean(infected_recently)
                        ) %>%
                        mutate(n=S1Q1+S0Q1+S1Q0+S0Q0)
                
                for(j in sum_df$simulation){
                        
                        tmp_df <- sum_df %>%
                                filter(simulation==j)
                        
                        #estimate seroincidence without questionairre
                        tmp_fit_NoQfit<- estimateSeroInc_vaxpop_simple(tmp_df %>%
                                                                               mutate(S1Q0=S1Q1+S1Q0,
                                                                                      S0Q0=S0Q1+S0Q0,
                                                                                      S1Q1 =0,
                                                                                      S0Q1 =0
                                                                               ),
                                                                       model=vax_model,
                                                                       params=tmp_param,
                                                                       multi=200
                        )
                        
                        #estimate seroincidence with questionairre
                        tmp_fit_Qfit<- estimateSeroInc_vaxpop_simple(tmp_df ,
                                                                     model=vax_model,
                                                                     params=tmp_param,
                                                                     multi=200
                        )
                        
                        
                        out_df <- 
                                bind_rows(out_df,
                                          bind_cols(tmp_df,
                                                    tmp_fit_NoQfit$fit %>% spread_draws(lambda) %>% median_qi() %>%
                                                            select(lambda,.lower,.upper) %>%
                                                            mutate(Adjustment="Without Questionairre")),
                                          bind_cols(tmp_df,
                                                    tmp_fit_Qfit$fit %>% spread_draws(lambda) %>% median_qi() %>%
                                                            select(lambda,.lower,.upper) %>%
                                                            mutate(Adjustment="With Questionairre") )
                                          
                                )
                        
                        
                }
                
                
                
        }
        
        
        return(out_df)
        
}



vaxSimSeroInc2<- function(incidence, #different incidences
                         coverage, # different levels coverage
                         # spec_vaxdf, #df describing specificity among vaccinated people at times,
                         sero_params,
                         quest_params,
                         sims
){

        
        vax_model <- stan_model("source/final_code/vax_strategy/stan/vax-model-simplified2.stan")
        
        #for each time
        # for(i in 1:nrow(spec_vaxdf)){
                
                cat("time:",i,"\n")
                
                #get sensitivitiy and specificity parameters
                tmp_param <- reduced_param
                tmp_param$spec_vax <- spec_vaxdf$spec[i]
                
                #generate serosurvey data
                tmp_survey <-  simVaxSurvey(
                        n_sim=sims,
                        unvaxxed_inc = incidence,
                        cov = coverage,
                        sens_quest=1,
                        spec_quest = 1,
                        effectiveness =0,
                        seropos_param=tmp_param
                )
                
                #sum up serosurvey data
                sum_df <- tmp_survey$ind_data %>%
                        mutate(since2=spec_vaxdf$newtime[i],
                               cov=coverage,
                               unvaxxed_inc=incidence,
                        ) %>% 
                        group_by(cov,unvaxxed_inc,variables,simulation,since2) %>%
                        summarize(S1Q1=sum(seropos*questionnaire),
                                  S0Q1=sum((1-seropos)*questionnaire),
                                  S1Q0=sum(seropos*(1-questionnaire)),
                                  S0Q0=sum((1-seropos)*(1-questionnaire)),
                                  study_seropos=mean(seropos),
                                  study_inc=mean(infected_recently)
                        ) %>%
                        mutate(n=S1Q1+S0Q1+S1Q0+S0Q0)
                
                for(j in sum_df$simulation){
                        
                        tmp_df <- sum_df %>%
                                filter(simulation==j)
                        
                        #estimate seroincidence without questionairre
                        tmp_fit_NoQfit<- estimateSeroInc_vaxpop_simple(tmp_df %>%
                                                                               mutate(S1Q0=S1Q1+S1Q0,
                                                                                      S0Q0=S0Q1+S0Q0,
                                                                                      S1Q1 =0,
                                                                                      S0Q1 =0
                                                                               ),
                                                                       model=vax_model,
                                                                       params=tmp_param,
                                                                       multi=200
                        )
                        
                        #estimate seroincidence with questionairre
                        tmp_fit_Qfit<- estimateSeroInc_vaxpop_simple(tmp_df ,
                                                                     model=vax_model,
                                                                     params=tmp_param,
                                                                     multi=200
                        )
                        
                        
                        out_df <- 
                                bind_rows(out_df,
                                          bind_cols(tmp_df,
                                                    tmp_fit_NoQfit$fit %>% spread_draws(lambda) %>% median_qi() %>%
                                                            select(lambda,.lower,.upper) %>%
                                                            mutate(Adjustment="Without Questionairre")),
                                          bind_cols(tmp_df,
                                                    tmp_fit_Qfit$fit %>% spread_draws(lambda) %>% median_qi() %>%
                                                            select(lambda,.lower,.upper) %>%
                                                            mutate(Adjustment="With Questionairre") )
                                          
                                )
                        
                        
                # }
                
                
                
        }
        
        
        return(out_df)
        
}





simVaxSurvey_new <- function(n_sim=1,
                         n_pop = 2000,
                         unvaxxed_inc,
                         cov,
                         effectiveness =0,
                         sens_quest = 1,
                         spec_quest = 1, # specificity of questionnaire
                         sero_params,
                         quest_params
){
        
        #set seropositivity params for questionairre stuff
        sens_novax = quest_params$sens_novax_mu
        sens_vax = quest_params$sens_vax_mu
        spec_vax = quest_params$spec_vax_mu
        spec_novax = quest_params$spec_novax_mu
        
        simulation_df <- data.frame()
        for(i in 1:n_sim){
                #simulate vax and recent infection status
                # this simulates the selection process
                vax <- rbinom(n_pop,1,cov) %>% sort()
                recent <- c(rbinom(length(vax[vax==0]),1,unvaxxed_inc), 
                            rbinom(length(vax[vax==1]),1,unvaxxed_inc*(1-effectiveness)))
                
                
                tmp_df  <- data.frame(
                        vax_status= vax,
                        infected_recently =recent
                ) %>%   
                        #get the simulation number
                        # mutate(simulation=i) %>%
                        mutate(simulation=paste0("sim_",str_pad(i,width=5,pad="0",side="left"))) %>%
                        #individual id
                        mutate(individual_id=paste0("id_",str_pad(1:n_pop,width=5,pad="0",side="left")))  %>%
                        #get questionnaire status
                        mutate(questionnaire=sapply(vax_status, 
                                                    FUN=getQuestionairre,
                                                    sens = sens_quest,
                                                    spec = spec_quest
                        )) %>%
                        #get seropositivity status
                        mutate(seropos= ifelse(vax_status ==0,
                                               sapply(infected_recently,
                                                      FUN=getSeroPos,
                                                      sens=sens_novax,
                                                      spec= spec_novax),
                                               sapply(infected_recently,
                                                      FUN=getSeroPos,
                                                      sens=sens_novax,
                                                      spec= spec_vax))
                        ) %>%
                        #get 3 class seroposiviity
                        mutate(obs_class3 = case_when(
                                
                                infected_recently ==1 ~ "RecentlyInfected",
                                infected_recently ==0 & vax_status ==1 ~ "RecentlyVaccinated", 
                                infected_recently ==0 & vax_status ==0 ~ "Neither" )) %>%
                        mutate(sero_class3=sapply(obs_class3,FUN=getSeroPos3,params=sero_params))
                        
                
                
                simulation_df <- bind_rows(simulation_df,
                                           tmp_df %>%
                                                mutate(coverage=cov,
                                                       incidence=unvaxxed_inc
                                                       )
                                           )        
        }
        
        out <- list(
                ind_data = simulation_df,
                n_pop = n_pop,
                unvaxxed_inc = unvaxxed_inc,
                cov = cov,
                effectiveness =effectiveness,
                sens_quest = sens_quest,
                spec_quest = spec_quest, # specificity of questionnaire
                sero_params,
                quest_params
        )
        
        return(out)
}


#generate seropositivity given sens spec
getSeroPos3<- function(class,params){
        
        if(class=="Neither") {
                out <- sample(c("Neither","RecentlyInfected","RecentlyVaccinated"),
                              1,
                              prob = params$Neither[[1]])
        }
        
        if(class=="RecentlyInfected") {
                out <- sample(c("Neither","RecentlyInfected","RecentlyVaccinated"),
                              1,
                              prob = params$RecentlyInfected[[1]])
        }
        
        if(class=="RecentlyVaccinated") {
                out <- sample(c("Neither","RecentlyInfected","RecentlyVaccinated"),
                              1,
                              prob = params$RecentlyVaccinated[[1]])
        }
                
                
        
        return(out)
}


runCaretLOOCV <- function(form, dat){
        
        train(form,
              data = dat,
              method = "ranger",
              num.trees = 1000,
              metric="logLoss",
              trControl = trainControl(method="LOOCV",
                                       index = groupKFold_fojo(dat$id),
                                       sampling="up",
                                       classProbs = TRUE,
                                       summaryFunction=mnLogLoss,
                                       savePredictions = "final"
              ),
              importance = 'impurity'
        )
}

estimateAlphaConstant <- function(model){
        
        pred_data <- distinct(model$data,id)
        
        linpred_m1 <- posterior_linpred(model,newdata =pred_data,
                                        ndraws=1000)
        
        
        if(nrow(fixef(model))==2){
                
                m1_infected_mat <- t(linpred_m1[,,1]) %>% as.data.frame()
                m1_vaxxed_mat <- t(linpred_m1[,,2])%>% as.data.frame()
                
                m1_combo <- bind_rows(
                        bind_cols(pred_data,
                                  m1_infected_mat) %>%
                                mutate(pred="RecentlyInfected"), 
                        bind_cols(pred_data,
                                  m1_vaxxed_mat)%>%
                                mutate(pred="RecentlyVaccinated")
                ) %>% gather(.draw,lin_pred,-c(id,pred)) 
                
                
                m1_draws <- m1_combo %>% spread(pred,lin_pred) %>%
                        mutate(prob_RecentlyInfected=exp(RecentlyInfected)/(1+exp(RecentlyInfected) + exp(RecentlyVaccinated)),
                               prob_RecentlyVaccinated=exp(RecentlyVaccinated)/(1+exp(RecentlyInfected) + exp(RecentlyVaccinated)),
                               prob_Neither=1/(1+exp(RecentlyInfected) + exp(RecentlyVaccinated)),
                               check= prob_RecentlyInfected + prob_RecentlyVaccinated  + prob_Neither) %>%
                        group_by(.draw) %>%
                        summarize(prob_RecentlyInfected=mean(prob_RecentlyInfected),
                                  prob_RecentlyVaccinated=mean(prob_RecentlyVaccinated),
                                  prob_Neither=mean(prob_Neither)
                        )  %>% 
                        select(-.draw)
                
                
        } 
        
        if(nrow(fixef(model))==1){
                
                
                m1_infected_mat <- t(linpred_m1)%>% as.data.frame()
                
                m1_combo <- bind_cols(pred_data,
                                      m1_infected_mat)%>%
                        gather(.draw,lin_pred,-c(id)) 
                
                
                m1_draws <- m1_combo %>% 
                        mutate(prob_RecentlyInfected=exp(lin_pred)/(1+exp(lin_pred)),
                               prob_Neither=1/(1+exp(lin_pred))
                        ) %>%
                        group_by(.draw) %>%
                        summarize(prob_RecentlyInfected=mean(prob_RecentlyInfected),
                                  prob_Neither=mean(prob_Neither)
                        )  %>%
                        select(-.draw) %>%
                        mutate(prob_RecentlyVaccinated=0)
                
        }
        
        
        return( list(m1_draws %>% 
                             select(prob_Neither,prob_RecentlyVaccinated,prob_RecentlyInfected) %>%
                             colMeans(),
                     m1_draws %>%
                             select(prob_Neither,prob_RecentlyVaccinated,prob_RecentlyInfected) %>%
                             cov()))
}

estimateAlphaTV <- function(model){
        
        
        tv_data <- model$data %>%
                distinct(id) %>%
                expand_grid(day_actual=1:200)
        
        linpred_m2 <- posterior_linpred(model,newdata = tv_data,ndraws=1000)
        
        
        if(nrow(fixef(model))==4){
                
                infected_mat <- t(linpred_m2[,,1]) %>% as.data.frame()
                vaxxed_mat <- t(linpred_m2[,,2])%>% as.data.frame()
                
                colnames(infected_mat) <- 1:ncol(infected_mat)
                colnames(vaxxed_mat) <- 1:ncol(vaxxed_mat)
                
                m2_combo <- bind_rows(
                        bind_cols(tv_data,infected_mat) %>%
                                mutate(pred="RecentlyInfected"), 
                        bind_cols(tv_data,vaxxed_mat)%>%
                                mutate(pred="RecentlyVaccinated")
                ) %>% gather(.draw,lin_pred,-c(id,day_actual,pred)) %>%
                        mutate(.draw=as.numeric(.draw))
                
                m2_draws <- m2_combo %>% spread(pred,lin_pred) %>%
                        mutate(prob_RecentlyInfected=exp(RecentlyInfected)/(1+exp(RecentlyInfected) + exp(RecentlyVaccinated)),
                               prob_RecentlyVaccinated=exp(RecentlyVaccinated)/(1+exp(RecentlyInfected) + exp(RecentlyVaccinated)),
                               prob_Neither=1/(1+exp(RecentlyInfected) + exp(RecentlyVaccinated)),
                               check= prob_RecentlyInfected + prob_RecentlyVaccinated  + prob_Neither) %>%
                        group_by(day_actual,.draw) %>%
                        summarize(prob_RecentlyInfected=mean(prob_RecentlyInfected),
                                  prob_RecentlyVaccinated=mean(prob_RecentlyVaccinated),
                                  prob_Neither=mean(prob_Neither)
                        ) %>%
                        gather(prob,value,-c(day_actual,.draw)) 
                
                
        }
        
        
        if(nrow(fixef(model))==2){
                
                infected_mat <- t(linpred_m2) %>% as.data.frame()
                
                colnames(infected_mat) <- 1:ncol(infected_mat)
                
                m2_combo <- 
                        bind_cols(tv_data,infected_mat) %>%
                        gather(.draw,lin_pred,-c(id,day_actual)) %>%
                        mutate(.draw=as.numeric(.draw))
                
                m2_draws <- m2_combo %>% 
                        mutate(prob_RecentlyInfected=exp(lin_pred)/(1+exp(lin_pred)),
                               prob_Neither=1/(1+exp(lin_pred)),
                               check= prob_RecentlyInfected  + prob_Neither) %>%
                        group_by(day_actual,.draw) %>%
                        summarize(prob_RecentlyInfected=mean(prob_RecentlyInfected),
                                  prob_Neither=mean(prob_Neither)
                        ) %>%
                        mutate(prob_RecentlyVaccinated=0) %>%
                        gather(prob,value,-c(day_actual,.draw)) 
                
                
        }
        
        list( m2_draws  %>%
                      group_by(prob,.draw) %>%
                      summarize(value=mean(value)) %>%
                      spread(prob,value) %>%
                      select(prob_Neither,prob_RecentlyVaccinated,prob_RecentlyInfected,-.draw) %>%
                      colMeans(),
              
              m2_draws  %>%
                      group_by(prob,.draw) %>%
                      summarize(value=mean(value)) %>%
                      spread(prob,value) %>%
                      select(prob_Neither,prob_RecentlyVaccinated,prob_RecentlyInfected,-.draw) %>%
                      cov()
        )
        
}





################# in development....






#' 
#' split_dataset_id <- function(data,proportion) {
#'         #' Data splitter by ID
#'         #'
#'         #' @description splits the data set by a particular proportion by ID
#'         #'
#'         #' @param data: wide dataset
#'         #' @param proportion: the proportion in the training set
#'         
#'         return(list(train,test))
#' 
#' }
#' 
#' 
#' find_spec_cutoff <- function(data, formula, specificity, true_neg, n_splits,proportion) {
#'         #' get cutoff estimate for a particular training set
#'         #'
#'         #' @description Identifies the cut-off needed to get desired specificity
#'         #'
#'         #' @param data: wide dataset
#'         #' @param specificity: desired specificity in proportion
#'         #' @param true_neg: string for the column with the true negatives
#'         #' @param n_splits: number of times the dataset is split to estimate cutoff 
#'         #' @param proportion: the proportion in the training set
#'         
#' 
#'         for(splt in 1:n_splits){
#'                 
#'                 tmp_datalist<- split_dataset_id(data=data,
#'                                  proportion = proportion
#'                                  )
#'                 # tmp_model <- caret::train()
#'                 # tmp_preds <- predict()
#'         
#'                 
#'                 
#'         }
#' 
#'         
#'         
#'         
#'         return(cut_off)
#' }
#' 
#' 
#' determine_serostatus <-function(data, cutoff_col){
#'         #' get cutoff estimate for a particular training set
#'         #'
#'         #' @description Identifies the cut-off needed to get desired specificity
#'         
#'         #' @param data: wide dataset
#'         #' @param pred_neg: the column where negative is considered
#'         #' @param pred_1: the first non-negative column 
#'         #' @param pred_2: the second non-negative column 
#'         #' @param cutoff: string for the column containing the cutoff for this sample
#'         
#' }
#'         
#' 
#' getLOOCV_cutoffs<- function(data,formula, specificity,true_neg,n_splits){
#'         #' find all the cutoffs
#'         #'
#'         #' @description Identifies the cut-offs for every training set
#'         
#'         id_vec <- data$id %>%  unique() 
#'         
#'         cutoff_df <- data.frame()
#'         for(i in id_vec){
#'                 train_dat <- filter(data, id != i)
#'                 train_cut <- find_spec_cutoff(data=train_dat,
#'                                               formula = formula,
#'                                               specificity=specificity,
#'                                               true_neg=true_neg,
#'                                               n_splits=n_splits) 
#'                 
#'                 cutoff_df <- bind_rows(
#'                         cutoff_df,
#'                         data.frame(id=i,
#'                                    cutoff=train_cut
#'                         )
#'                 )
#'         }
#'         
#'         
#'         return(cutoff_df)
#' }        
#' 
#' 
#' 
#' loocv_caret_cutoff <- function(data, formula, 
#'                                specificity=0.95,
#'                                n_splits=100,
#'                                proportion=0.5,
#'                                true_neg){
#'         
#'         #'
#'         #' @description Conducts LOOCV and chooses the cutoff for each individual to predict
#'         #' serostatus that would come with a desired specificity 
#'         #'
#'         #' @param data: wide dataset,
#'         #' @param formula: formula used for the RF model
#'         #' @param specificity: desired specificity in proportion
#'         #' @param true_neg: string for the class with true negative
#'         
#'         
#'         #### LOOCV is conducted using the caret package
#'         loocv_simple <-runCaretLOOCV(form=formula,dat=data)
#'         
#'         #### Cutoffs are found for each individual to ensure desired specificity
#'         loocv_cutoffs <- getLOOCV_cutoffs(data=data,
#'                                           formula=formula,
#'                                           specificity=specificity,
#'                                           n_splits=n_splits,
#'                                           proportion=proportion,
#'                                           true_neg=true_neg)
#'         
#'         #### Join LOOCV results with cutoffs and find serostatus
#'         loocv_joined <- loocv_simple %>%
#'                 left_join(loocv_cutoffs, by="id") %>%
#'                 determine_serostatus()
#'                 
#'         
#'         #### create final output object
#'         return(loocv_joined)
#' }







