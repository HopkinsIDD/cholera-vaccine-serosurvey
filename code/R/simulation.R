rm(list=ls())

library(tidyverse)
library(tidybayes)

source("code/R/packages.R")
source("code/R/utils.R")

options(mc.cores = 1)


#load loocv objects
loocv_spec_list <- read_rds(
        paste0("data/generated_data/analysis_objects/loocv/","loocv_spec_list.rds")
)
loocv_sens_list <- read_rds("data/generated_data/analysis_objects/loocv/loocv_tvfit_list.rds")

loocv_df <- read_rds("data/generated_data/analysis_objects/loocv/loocv_df.rds")

# 0. Stan model list

# Original model:

# 1. Time-varying sensitivity
# 2. Specificity among non-vaccine
# 3. Time-varying vaccine misclasification

# Alternative model:
# 1. Time-varying sensitivity
# 2. Specificity among vax & non-vax (covariate)


# 1. Get original model characteristics-----------

#### Sensitivity

#load stan fit
original_sens_fitdata <-loocv_sens_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$Case$data
original_sens_spread <- loocv_sens_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$Case$fit %>% 
        spread_draws(`(Intercept)`,
                     b[term,group],
                     day1,day2,day3
        )

#for loop through potential days since infection
original_sens_obj <- data.frame()
original_sens_beta_draws <-data.frame()

for(t in 1:200){
        
        cat(t,"\n")
        time <- log(t) -mean(log(original_sens_fitdata$new_day))
        
        #get individual sensitivities for simulation
        tmp_df <- original_sens_spread %>%
                mutate(logit_theta_j=`(Intercept)` + b+
                               day1*time +
                               day2*time^2+
                               day3*time^3
                       
                ) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                group_by(group) %>%
                summarize(theta=mean(theta_j),n()) %>% #average across individuals here
                mutate(days_ago=t)%>%
                rename(pos_ID=group,
                       ind_sens=theta
                )
        
        original_sens_obj <- bind_rows(original_sens_obj,tmp_df)
        
        #get sensitivity distribution for seroincidence estimation
        tmp_dist <- original_sens_spread %>%
                mutate(logit_theta_j=`(Intercept)` + b+
                               day1*time +
                               day2*time^2+
                               day3*time^3
                       
                ) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                group_by(.draw) %>%
                summarize(theta=mean(theta_j),n()) %>%
                mutate(days_ago=t)
        
        original_sens_beta_draws <-bind_rows(
                original_sens_beta_draws,
                tmp_dist
        )
}


#fit beta distribution
original_sens_beta <- original_sens_beta_draws %>% group_by(.draw)%>%
        summarize(sens=mean(theta)) %>%
        pull(sens) %>%
        MASS::fitdistr(
                densfun = "beta", list(shape1=1,shape2=1))

#### Specificity (non-vaccinees)

#get individual specificities for simulation
original_spec_obj<- loocv_spec_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$`Outside Window`$fit%>% 
        spread_draws(`(Intercept)`,
                     b[term,group]) %>%
        mutate(logit_theta_j=`(Intercept)` + b
        ) %>%
        mutate(theta_j = 1/(1+exp(-logit_theta_j))) %>%
        group_by(group)%>%
        summarize(theta=mean(theta_j),n()) %>%
        rename(neg_ID=group
        ) %>%
        mutate(ind_spec=1-theta)

#get specificity distribution for seroincidence estimation
original_spec_dist <- loocv_spec_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$`Outside Window`$fit%>% 
        spread_draws(`(Intercept)`,
                     b[term,group]) %>%
        mutate(logit_theta_j=`(Intercept)` + b
        ) %>%
        mutate(theta_j = 1/(1+exp(-logit_theta_j)))  %>%
        group_by(.draw) %>%
        summarize(spec=1-mean(theta_j))

#fit beta distribution
original_spec_beta <- MASS::fitdistr(original_spec_dist$spec,
                            densfun = "beta", list(shape1=1,shape2=1))

#### Specificity (vaccinees)

#load stan fit
original_tvfpr_fit <- read_rds("data/generated_data/analysis_objects/misclassification/tvfpr_fit.rds")

original_vax_fitdata <-original_tvfpr_fit$`200`$`All Vaccinees (Reduced IgG Panel)`$data
original_vax_spread <- original_tvfpr_fit$`200`$`All Vaccinees (Reduced IgG Panel)`$fit%>% 
        spread_draws(`(Intercept)`,
                     b[term,group],
                     day1,day2,day3
        )

#for loop through potential campaign timings
campaign_timepoints <- c(21,120)#c(21,45,90,120,180)

original_vax_tv_obj <- data.frame()
original_vax_beta_draws <- data.frame()
original_vax_beta <- data.frame()


for(t in campaign_timepoints){
        
        cat(t,"\n")
        time <- log(t) -mean(log(original_vax_fitdata$new_day))
        
        #get individual specificities for simulation
        tmp_df <- original_vax_spread %>%
                mutate(logit_theta_j=`(Intercept)` + b+
                               day1*time +
                               day2*time^2+
                               day3*time^3
                       
                ) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                group_by(group) %>%
                summarize(theta=mean(theta_j),n()) %>% #average across individuals here
                mutate(days_ago=t)%>%
                rename(vax_ID=group) %>%
                mutate(ind_spec=1-theta)
        
        original_vax_tv_obj <- bind_rows(original_vax_tv_obj,tmp_df)

        #get specificity distribution for seroincidence estimation
        tmp_dist <- original_vax_spread %>%
                mutate(logit_theta_j=`(Intercept)` + b+
                               day1*time +
                               day2*time^2+
                               day3*time^3
                       
                ) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                group_by(.draw) %>%
                summarize(theta=1-mean(theta_j),n()) %>% #average across individuals here
                mutate(days_ago=t)
        
        original_vax_beta_draws <- bind_rows(original_vax_beta_draws,tmp_dist)
        
        #fit beta distribution
        tmp_beta <- tmp_dist$theta %>%
                MASS::fitdistr(
                        densfun = "beta", list(shape1=1,shape2=1))
        
        original_vax_beta <- bind_rows(original_vax_beta,
                              tmp_beta$estimate %>% as.matrix()%>%
                                      t() %>% data.frame() %>%
                                      mutate(days_ago=t)
        )

}



# 2. Get alternative model characteristics -----------

#### Sensitivity

#load stan fit
alternative_sens_fitdata <-loocv_sens_list$`Reduced IgG Panel`$all_fitdata$`200-day`$Case$data
alternative_sens_spread <- loocv_sens_list$`Reduced IgG Panel`$all_fitdata$`200-day`$Case$fit %>% 
        spread_draws(`(Intercept)`,
                     b[term,group],
                     day1,day2,day3
        )

#for loop through potential days since infection
alternative_sens_obj <- data.frame()
alternative_sens_beta_draws <-data.frame()

for(t in 1:200){
        
        cat(t,"\n")
        time <- log(t) -mean(log(alternative_sens_fitdata$new_day))
        
        #get individual sensitivities for simulation
        tmp_df <- alternative_sens_spread %>%
                mutate(logit_theta_j=`(Intercept)` + b+
                               day1*time +
                               day2*time^2+
                               day3*time^3
                       
                ) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                group_by(group) %>%
                summarize(theta=mean(theta_j),n()) %>% #average across individuals here
                mutate(days_ago=t)%>%
                rename(pos_ID=group,
                       ind_sens=theta
                )
        
        alternative_sens_obj <- bind_rows(alternative_sens_obj,tmp_df)
        
        #get sensitivity distribution for seroincidence estimation
        tmp_dist <- alternative_sens_spread %>%
                mutate(logit_theta_j=`(Intercept)` + b+
                               day1*time +
                               day2*time^2+
                               day3*time^3
                       
                ) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                group_by(.draw) %>%
                summarize(theta=mean(theta_j),n()) %>%
                mutate(days_ago=t)
        
        alternative_sens_beta_draws <-bind_rows(
                alternative_sens_beta_draws,
                tmp_dist
        )
}


#fit beta distribution
alternative_sens_beta <- alternative_sens_beta_draws %>% group_by(.draw)%>%
        summarize(sens=mean(theta)) %>%
        pull(sens) %>%
        MASS::fitdistr(
                densfun = "beta", list(shape1=1,shape2=1))





# fit new stan model        
new_fit <- bind_rows(
        filter(loocv_df,cohort=="SMIC/PIC",inf_200==0),
        filter(loocv_df,status=="Vaccinee")        
) %>%
        filter(variables=="Reduced IgG Panel",end_window==200,base_data_name=="all_fitdata") %>%
        fit_spec_byvax(
                end_window=200,
                seropos_type= "spec95_seropos")

#### Specificity (non-vaccinees)

# #get individual specificities for simulation
alternative_spec_obj<- new_fit$fit%>% 
        spread_draws(`(Intercept)`,
                     b[term,group]
        ) %>%
        mutate(id=str_remove(group,"id\\:")) %>%
        left_join(distinct(new_fit$data,id,status)) %>%
        filter(status!="Vaccinee")%>%
        mutate(logit_theta_j=`(Intercept)` + b ) %>%
        mutate(theta_j = 1/(1+exp(-logit_theta_j))) %>%
        group_by(group)%>%
        summarize(theta=mean(theta_j),n()) %>%
        rename(neg_ID=group
        ) %>%
        mutate(ind_spec=1-theta)
        

# #### Specificity (vaccinees)
# 
# #load stan fit
alternative_vax_obj <- new_fit$fit%>%
        spread_draws(`(Intercept)`,
                     b[term,group],
                     vax
        ) %>%
                mutate(id=str_remove(group,"id\\:")) %>%
                left_join(distinct(new_fit$data,id,status)) %>%
                filter(status=="Vaccinee")%>%
                mutate(logit_theta_j=`(Intercept)` + b + vax) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j))) %>%
                group_by(group)%>%
                summarize(theta=mean(theta_j),n()) %>%
                rename(vax_ID=group
                ) %>%
                mutate(ind_spec=1-theta)

#probably want something where we can estimate the specificity with a covariate given the two differences
#can leave this be for now
# see fit_spec_covariate for start on this


new_draws <- new_fit$fit%>% 
        spread_draws(`(Intercept)`,
                     b[term,group],
                     vax
                     ) %>%
        mutate(id=str_remove(group,"id\\:")) %>%
        left_join(distinct(new_fit$data,id,status)) %>%
        mutate(new_status=ifelse(status=="Vaccinee","Vaccinee","Outside Window"))%>%
        mutate(logit_theta_j=`(Intercept)` + b + vax*(new_status=="Vaccinee")) %>%
        mutate(theta_j = 1/(1+exp(-logit_theta_j)))  %>%
        group_by(.draw,new_status) %>%
        summarize(spec=1-mean(theta_j)) %>%
        spread(new_status,spec)




# 3. Run simulations and stan models to get seroincidence estimates -----------

simulation_df <-data.frame()
sims_original <- list()
lambda_df <- data.frame()

simulations <- 1000
n_survey <- 1000
coverage_values <- c(0, 0.25, 0.5, 0.75)
campaign_timepoints <- c(21,120)

stan_model_compiled <- stan_model("code/stan/vax-model-simplified.stan")

true_seroinc <- 0.1
# for(true_seroinc in c(0.05,0.1,0.2)){ #choose several values of true seroincidence to look at
for(cov in coverage_values){
        
        
        #simulate the true values
        tmp_sim <- sim_truth(simulations,n_survey,coverage = cov, incidence = true_seroinc)
        
        #simulate seropositivity
        for(t in campaign_timepoints){
                
                cat("Coverage: ",cov, " Timing: ", t,"\n")
                
                sims_original[[as.character(t)]] <- tmp_sim %>% lapply(make_seropos,model_name = "Original Model",
                                                    sens_df= original_sens_obj,
                                                    spec_df= original_spec_obj,
                                                    vax_df = original_vax_tv_obj,
                                                    vax_day=t
                ) %>%  bind_rows(.id="simulation") 
                
                
                simulation_df <- bind_rows(
                        simulation_df,
                        sims_original[[as.character(t)]] %>%
                                group_by(simulation, coverage, campaign_day,incidence) %>%
                                summarize(truth=mean(R),seropos=mean(`Original Model`)) %>% 
                                mutate(model="Original")
                )
                
        }
                
        sims_alternative <- tmp_sim %>% lapply(make_seropos,model_name = "Alternative Model",
                                               sens_df= alternative_sens_obj,
                                               spec_df= alternative_spec_obj,
                                               vax_df = alternative_vax_obj,
                                               vax_day=NA
        ) %>%  bind_rows(.id="simulation")
        
        simulation_df <- bind_rows(
                simulation_df,
                sims_alternative %>%
                        group_by(simulation, coverage, campaign_day,incidence) %>%
                        summarize(truth=mean(R),seropos=mean(`Alternative Model`)) %>% 
                        mutate(model="Alternative")
        )
        
             for(i in names(tmp_sim)){
                     for(t in campaign_timepoints){
                             
                     cat("Simulation: ",i , "Timepoint: ", t, "\n")
                             
                        ### Strategy 1: Ignore
                        this_count_vec <- c(0,0,0,0)
                        
                        this_count <- sims_original[[as.character(t)]] %>%
                                filter(simulation==i) %>%
                                # filter(campaign_day==t)%>%
                                count(simulation,`Original Model`,V,coverage,incidence,campaign_day) %>%
                                mutate(category=glue::glue("S{`Original Model`}Q{V}")) %>%
                                select(-`Original Model`,-V) %>%
                                spread(category,n) 
                        
                        if("S0Q0" %in% colnames(this_count)) this_count_vec[1] <- this_count$S0Q0
                        if("S1Q0" %in% colnames(this_count)) this_count_vec[2] <- this_count$S1Q0
                        if("S0Q1" %in% colnames(this_count)) this_count_vec[3] <- this_count$S0Q1
                        if("S1Q1" %in% colnames(this_count)) this_count_vec[4] <- this_count$S1Q1
                        
                        #no adjustment
                        biased_data <- list(
                                S0Q0 = this_count_vec[1] + this_count_vec[3],
                                S1Q0 = this_count_vec[2] + this_count_vec[4],
                                S0Q1 = 0,
                                S1Q1 = 0,
                                
                                alpha_p1 = original_sens_beta$estimate[1],
                                alpha_p2 = original_sens_beta$estimate[2],
                                beta_p1 = original_spec_beta$estimate[1],
                                beta_p2 = original_spec_beta$estimate[2],
                                epsilon_p1 = original_vax_beta %>%
                                        filter(days_ago==t) %>%
                                        pull(shape1),
                                epsilon_p2 =original_vax_beta %>%
                                        filter(days_ago==t) %>%
                                        pull(shape2)
                        )
                        
                        biased_fit<- sampling(stan_model_compiled,
                                              data = biased_data
                                              )
                
                        ### Strategy 2: Questionnaire adjustment
                        questionnaire_data <- list(
                                S0Q0 = this_count_vec[1],
                                S1Q0 = this_count_vec[2],
                                S0Q1 = this_count_vec[3],
                                S1Q1 = this_count_vec[4],
                                
                                alpha_p1 = original_sens_beta$estimate[1],
                                alpha_p2 = original_sens_beta$estimate[2],
                                beta_p1 = original_spec_beta$estimate[1],
                                beta_p2 = original_spec_beta$estimate[2],
                                epsilon_p1 = original_vax_beta %>%
                                        filter(days_ago==t) %>%
                                        pull(shape1),
                                epsilon_p2 =original_vax_beta %>%
                                        filter(days_ago==t) %>%
                                        pull(shape2)
                        )
                        
                        questionnaire_fit <- sampling(stan_model_compiled,
                                                      data = questionnaire_data)
                        
                
                        lambda_df <- bind_rows(
                                lambda_df,
                                bind_rows(
                                        spread_draws(biased_fit,lambda) %>% mean_qi() %>%
                                                select( lambda,.lower,.upper) %>%
                                                mutate(adjustment="Ignore"),
                                        spread_draws(questionnaire_fit,lambda) %>% mean_qi() %>%
                                                select( lambda,.lower,.upper) %>%
                                                mutate(adjustment="Questionnaire")
                                ) %>%
                                        mutate(coverage = cov,
                                               campaign_day=t,
                                               simulation=i
                                        )
                        )
                        
                }
                     #use the alternative model
                     this_count_vec <- c(0,0,0,0)
                     
                     this_count <- sims_alternative %>%
                             filter(simulation==i) %>%
                             count(simulation,`Alternative Model`,V,coverage,incidence,campaign_day) %>%
                             mutate(category=glue::glue("S{`Alternative Model`}Q{V}")) %>%
                             select(-`Alternative Model`,-V) %>%
                             spread(category,n)
                     
                     if("S0Q0" %in% colnames(this_count)) this_count_vec[1] <- this_count$S0Q0
                     if("S1Q0" %in% colnames(this_count)) this_count_vec[2] <- this_count$S1Q0
                     if("S0Q1" %in% colnames(this_count)) this_count_vec[3] <- this_count$S0Q1
                     if("S1Q1" %in% colnames(this_count)) this_count_vec[4] <- this_count$S1Q1
                     
                     
                     ### Strategy 3: Alternative model
                        # run the rapid vax coverage survey
                        rapid_cov <- rbinom(1,500,cov)/500
                     
                        #calculate the beta distribution for specificity distribution
                        alternative_spec_beta <- new_draws %>%
                                mutate(combined=rapid_cov*Vaccinee+(1-rapid_cov)*`Outside Window`) %>%
                                pull(combined) %>%
                                MASS::fitdistr(
                                        densfun = "beta", list(shape1=1,shape2=1))
                        
                     
                     alternative_data <- list(
                             S0Q0 = this_count_vec[1] + this_count_vec[3],
                             S1Q0 = this_count_vec[2] + this_count_vec[4],
                             S0Q1 = 0,
                             S1Q1 = 0,
                             
                             alpha_p1 = alternative_sens_beta$estimate[1],
                             alpha_p2 = alternative_sens_beta$estimate[2],
                             beta_p1 = alternative_spec_beta$estimate[1],
                             beta_p2 = alternative_spec_beta$estimate[2],
                             epsilon_p1 = alternative_spec_beta$estimate[1],
                             epsilon_p2 =alternative_spec_beta$estimate[2]
                     )
                     
                     alternative_fit <- sampling(stan_model_compiled,
                                                   data = alternative_data)
                     
                     
                     lambda_df <- bind_rows(
                             lambda_df,
                             spread_draws(alternative_fit,lambda) %>% mean_qi() %>%
                                     select( lambda,.lower,.upper) %>%
                                     mutate(adjustment="Alternative") %>%
                             mutate(coverage = cov,
                                    campaign_day=NA,
                                    simulation=i
                                    )
                                        )
                     
                     write_rds(simulation_df, paste0("data/generated_data/analysis_objects/simulation/simulation_df_",
                                                     true_seroinc*100,
                                                     ".rds"))
                     write_rds(lambda_df, paste0("data/generated_data/analysis_objects/simulation/lambda_df_",
                               true_seroinc*100,
                               ".rds"))
                     write_rds(sims_original,paste0("data/generated_data/analysis_objects/simulation/sims_original_",
                               true_seroinc*100, 
                               ".rds"))
                     
             }
}
# }
                


# simulation_df <- read_rds("data/generated_data/analysis_objects/simulation/simulation_df_5.rds")
# lambda_df <- read_rds("data/generated_data/analysis_objects/simulation/lambda_df_5.rds")
# sims_original <- read_rds("data/generated_data/analysis_objects/simulation/sims_original_5.rds")
# 
# 
# 
# lambda_df  %>%
#                 ggplot(aes(x=adjustment,y=lambda))+
#                 ggbeeswarm::geom_beeswarm(alpha=0.5)+
#                 geom_hline(yintercept=0.05, lty=2, col="red")+
#                 facet_grid(coverage~campaign_day)+
#                 cowplot::theme_cowplot()+
#                 theme(axis.text.x = element_text(angle=45,hjust = 1))+
#                 ylab("200-day seroincidence")+
#                 xlab("Method")
# 
# 
# table(simulation_df$incidence)
# 
# sims_viz <- filter(simulation_df,campaign_day==21, coverage==0) %>%
#         filter(incidence==0.10)%>%
#                 distinct(simulation,truth) %>%
#                 mutate(`Ignore 21`=truth,
#                        `Ignore 120`=truth,
#                        `Questionnaire 21`=truth,
#                        `Alternative NA`=truth,
#                        ) %>%
#                 select(-truth) %>%
#                 gather(new_strategy, truth,-simulation) %>%
#         mutate(new_strategy = factor(new_strategy,
#                                      levels=c("Ignore 21", "Ignore 120",
#                                               "Questionnaire 21", "Alternative NA"
#                                      ),
#                                      labels=c("Original Model\n(21 days)", "Original Model\n(120 days)",
#                                               "Original Model +\nQuestionnaire\n(21 days)", "Alternative Model 1\n(21 days)")
#         )
#         ) %>%
#         expand_grid(coverage=c("0% Coverage",
#                         "25% Coverage",
#                         "50% Coverage",
#                         "75% Coverage"
#         ))%>%
#         mutate(new_coverage=factor(coverage,
#                                    labels=c("0% Coverage",
#                                             "25% Coverage",
#                                             "50% Coverage",
#                                             "75% Coverage"
#                                    )
#         )) 
#                 
# lambda_df  %>%
#         filter(campaign_day == 21 |
#                (adjustment == "Ignore" & campaign_day == 120)|
#                adjustment == "Alternative"
#                        )%>%
#         mutate(new_strategy=paste(adjustment,campaign_day)) %>%
#         mutate(new_strategy = factor(new_strategy,
#                 levels=c("Ignore 21", "Ignore 120",
#                          "Questionnaire 21", "Alternative NA"
#                          ),
#                 labels=c("Original Model\n(21 days)", "Original Model\n(120 days)",
#                          "Original Model +\nQuestionnaire\n(21 days)", "Alternative Model 1\n(21 days)")
#                 )
#                        ) %>%
#         mutate(new_coverage=factor(coverage,
#                                    labels=c("0% Coverage",
#                                             "25% Coverage",
#                                             "50% Coverage",
#                                             "75% Coverage"
#                                             )
#                                    )) %>%
#         group_by(new_strategy,coverage)%>%
#         mutate(avg_lambda=mean(lambda)) %>%
#         ungroup()%>%
#         select(simulation,new_strategy,new_coverage,lambda,avg_lambda) %>%
#         left_join(sims_viz,by = c("simulation", "new_strategy", "new_coverage")) %>%
#         gather(type,value,-c(simulation,new_strategy, new_coverage, avg_lambda,coverage))%>%
#         ggplot(aes(x=new_strategy))+
#         geom_violin(aes(y=value,col=type,alpha=type),position="identity",fill="grey",
#                     trim=FALSE,scale="area")+
#         scale_color_manual(values=c("blue",NA),na.value = NA)+
#         scale_alpha_manual(values=c(0,0.3))+
#         geom_point(aes(y=avg_lambda),col="blue")+
#         facet_wrap(new_coverage~.)+
#         cowplot::theme_cowplot()+
#         theme(axis.text.x = element_text(angle=45,hjust = 1),
#               legend.position = "none"
#               )+
#         ylab("200-day seroincidence")+
#         xlab("Strategy")
# 
# 
# lambda_df  %>%
#         filter(adjustment== "Ignore", campaign_day == 21) %>%
#         mutate(new_coverage=factor(coverage,
#                                    labels=c("0% Coverage",
#                                             "25% Coverage",
#                                             "50% Coverage",
#                                             "75% Coverage"
#                                    )
#         )) %>%
#         group_by(new_coverage) %>%
#         summarize(lambda=glue::glue("{round(mean(lambda),2)} {round(min(lambda),2)}-{round(max(lambda),2)}")
#         )
# 
# 
# lambda_df  %>%
#         filter(campaign_day == 21 |
#                        (adjustment == "Ignore" & campaign_day == 120)|
#                        adjustment == "Alternative"
#         )%>%
#         mutate(new_strategy=paste(adjustment,campaign_day)) %>%
#         mutate(new_strategy = factor(new_strategy,
#                                      levels=c("Ignore 21", "Ignore 120",
#                                               "Questionnaire 21", "Alternative NA"
#                                      ),
#                                      labels=c("Original Model (21 days)", "Original Model (120 days)",
#                                               "Original Model + Questionnaire (21 days)", 
#                                               "Alternative Model 1 (21 days)")
#         )
#         ) %>%
#         filter(coverage==0.5) %>%
#         group_by(new_strategy) %>%
#         summarize(lambda=glue::glue("{round(mean(lambda),2)} {round(min(lambda),2)}-{round(max(lambda),2)}")
#                   )
# 

# 3. Simulate truth and seropositivity -----------
# sims <- sim_truth(10,1000,0.5,0.05)
# 
# sims_original <- sims %>% lapply(make_seropos,model_name = "Original",
#                 sens_df= original_sens_obj,
#                 spec_df= original_spec_obj,
#                 vax_df = original_vax_tv_obj,
#                 vax_day=21
#                 )
# 
# 
# sims_alternative <- sims %>% lapply(make_seropos,model_name = "Alternative Model",
#                          sens_df= alternative_sens_obj,
#                          spec_df= alternative_spec_obj,
#                          vax_df = alternative_vax_obj,
#                          vax_day=NA
#                          )

# %>%
#         bind_rows(lambda_df %>% select(-.lower,-.upper))%>%
#         # mutate(adjustment = ifelse(adjustment=="Unadjusted",
#         #                            "Adjusted",adjustment
#         # ))%>%
#         # mutate(adjustment=factor(adjustment,levels=c("Incidence","Seropositivity","Adjusted","Questionnaire")))%>%
#         filter(campaign_day<180)%>%
#         mutate(campaign_day=paste(campaign_day,"days")) %>%
#         mutate(campaign_day=factor(campaign_day, 
#                                    levels=c("21 days",
#                                             "45 days",
#                                             "90 days",
#                                             "120 days"
#                                    )))%>%
#         ggplot(aes(x=adjustment,y=lambda))+
#         ggbeeswarm::geom_beeswarm(alpha=0.5)+
#         # geom_boxplot()+
#         geom_hline(yintercept=0.05, lty=2, col="red")+
#         facet_grid(coverage~campaign_day)+
#         cowplot::theme_cowplot()+
#         theme(axis.text.x = element_text(angle=45,hjust = 1))+
#         ylab("200-day seroincidence")+
#         xlab("Method")





### OLD




# #run simulations
# 
# simulation_df <-data.frame()
# lambda_df <- data.frame()
# 
# simulations <- 10
# coverage_values <- c(0, 0.25, 0.5, 0.75)
# 
# stan_model_compiled <- stan_model("code/stan/vax-model-simplified.stan")
# 
# for(cov in coverage_values[3]){
#         for(days in campaign_timepoints[2]){
#                 
#                 cat(cov," "  ,days, "\n")
#                 
#                 new_vax_obj <- new_vax_tv_obj %>%
#                         filter(days_ago==days)%>%
#                         select(-days_ago)
#                 
#                 tmp_vax_run <- simulate_serosurvey_vax(n_sim = simulations,
#                                                         n_survey = 1000,
#                                                         sens_df= new_sens_obj,
#                                                         spec_df= new_spec_obj,
#                                                         vax_df = new_vax_obj,
#                                                         coverage=cov,
#                                                         incidence=0.05
#                 ) %>% bind_rows(.id="simulation") %>%
#                         mutate(campaign_day=days)
#                         
#                 
#                 simulation_df <- bind_rows(
#                         simulation_df,tmp_vax_run
#                 )
#                 
#                 tmp_counts <- tmp_vax_run %>%
#                         mutate(campaign_day=days) %>% 
#                         ungroup()%>%
#                         count(simulation,seropos,V,coverage,incidence,campaign_day) %>%
#                         mutate(category=glue::glue("S{seropos}Q{V}")) %>%
#                         select(-seropos,-V) %>%
#                         spread(category,n)
#                 
#                 for(i in tmp_counts$simulation){
#                         
#                         this_count_vec <- c(0,0,0,0)
#                         
#                         this_count <- tmp_counts %>%
#                                 filter(simulation==i)
#                         
#                         if("S0Q0" %in% colnames(this_count)) this_count_vec[1] <- this_count$S0Q0
#                         if("S1Q0" %in% colnames(this_count)) this_count_vec[2] <- this_count$S1Q0
#                         if("S0Q1" %in% colnames(this_count)) this_count_vec[3] <- this_count$S0Q1
#                         if("S1Q1" %in% colnames(this_count)) this_count_vec[4] <- this_count$S1Q1
#                         
#                         #no adjustment
#                         biased_data <- list(
#                                 S0Q0 = this_count_vec[1] + this_count_vec[3],
#                                 S1Q0 = this_count_vec[2] + this_count_vec[4],
#                                 S0Q1 = 0,
#                                 S1Q1 = 0,
#                                 
#                                 alpha_p1 = sens_beta$estimate[1],
#                                 alpha_p2 = sens_beta$estimate[2],
#                                 beta_p1 = spec_beta$estimate[1],
#                                 beta_p2 = spec_beta$estimate[2],
#                                 epsilon_p1 = vax_beta %>%
#                                         filter(days_ago==days) %>%
#                                         pull(shape1),
#                                 epsilon_p2 =vax_beta %>%
#                                         filter(days_ago==days) %>%
#                                         pull(shape2)
#                         )
#                         
#                         biased_fit<- sampling(stan_model_compiled,
#                                              data = biased_data)
#                         
#                         #questionnaire adjustment
#                         questionnaire_data <- list(
#                                 S0Q0 = this_count_vec[1],
#                                 S1Q0 = this_count_vec[2],
#                                 S0Q1 = this_count_vec[3],
#                                 S1Q1 = this_count_vec[4],
#                                 
#                                 alpha_p1 = sens_beta$estimate[1],
#                                 alpha_p2 = sens_beta$estimate[2],
#                                 beta_p1 = spec_beta$estimate[1],
#                                 beta_p2 = spec_beta$estimate[2],
#                                 epsilon_p1 = vax_beta %>%
#                                         filter(days_ago==days) %>%
#                                         pull(shape1),
#                                 epsilon_p2 =vax_beta %>%
#                                         filter(days_ago==days) %>%
#                                         pull(shape2)
#                         )
#                         
#                         questionnaire_fit <- sampling(stan_model_compiled,
#                                                       data = questionnaire_data)
#                         
#                         
#                         lambda_df <- bind_rows(
#                                 
#                                 lambda_df,
#                                 bind_rows(
#                                 spread_draws(biased_fit,lambda) %>% mean_qi() %>%
#                                         select( lambda,.lower,.upper) %>%
#                                         mutate(adjustment="Adjusted"),
#                                 spread_draws(questionnaire_fit,lambda) %>% mean_qi() %>%
#                                         select( lambda,.lower,.upper) %>%
#                                         mutate(adjustment="Questionnaire")
#                                 ) %>%
#                                         mutate(coverage = cov,
#                                                campaign_day=days
#                                                )
#                         )
#                         
#                 }
#                 
# }}
# 
# simulation_df %>% 
#                 group_by(simulation,coverage,campaign_day) %>%
#                 summarize(Seropositivity=mean(seropos),
#                           Incidence = mean(R)
#                           ) %>%
#         gather(adjustment,lambda,-simulation,-coverage,-campaign_day) %>%
#         bind_rows(lambda_df %>% select(-.lower,-.upper))%>%
#         mutate(adjustment = ifelse(adjustment=="Unadjusted",
#                                    "Adjusted",adjustment
#                                    ))%>%
#         mutate(adjustment=factor(adjustment,levels=c("Incidence","Seropositivity","Adjusted","Questionnaire")))%>%
#         filter(campaign_day<180)%>%
#         mutate(campaign_day=paste(campaign_day,"days")) %>%
#         mutate(campaign_day=factor(campaign_day, 
#                                    levels=c("21 days",
#                                             "45 days",
#                                             "90 days",
#                                             "120 days"
#                                             )))%>%
#         ggplot(aes(x=adjustment,y=lambda))+
#         ggbeeswarm::geom_beeswarm(alpha=0.5)+
#         # geom_boxplot()+
#         geom_hline(yintercept=0.05, lty=2, col="red")+
#         facet_grid(coverage~campaign_day)+
#         cowplot::theme_cowplot()+
#         theme(axis.text.x = element_text(angle=45,hjust = 1))+
#         ylab("200-day seroincidence")+
#         xlab("Method")
# 
# 
# 
# #run it
# 
# cat_counts<- vax21_run_25$sim1 %>% 
#         ungroup()%>%
#         count(seropos,V) %>%
#         mutate(category=glue::glue("S{seropos}Q{V}")) %>%
#         select(n,category) %>%
#         spread(category,n)
# 
# stan_obj <- list(
#           S0Q0 = cat_counts$S0Q0+cat_counts$S0Q1,
#           S1Q0 = cat_counts$S1Q0+cat_counts$S1Q1,
#           S0Q1 = 0,
#           S1Q1 = 0,
#         
#          alpha_p1 = sens_beta$estimate[1],
#          alpha_p2 = sens_beta$estimate[2],
#          beta_p1 = spec_beta$estimate[1],
#          beta_p2 = spec_beta$estimate[2],
#          epsilon_p1 = vax_beta %>%
#                  filter(days_ago==21) %>%
#                  pull(shape1),
#          epsilon_p2 =vax_beta %>%
#                  filter(days_ago==21) %>%
#                  pull(shape2)
# )
# 
# library(rstan)
# 
# second_model <- stan_model("code/stan/vax-model-simplified.stan")
# fit <- fit_second_model <- sampling(second_model, data = stan_obj)
# 
# shinystan::launch_shinystan(fit)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### OLD
# 
# # #get time varying sensitivity
# # sens_obj <- data.frame(
# #         pos_ID = paste0("Z",rep(1:50,200)),
# #         diff=rep(rnorm(50,sd=1),200)
# # ) %>% arrange(pos_ID) %>%
# #         mutate(days_ago=rep(1:200,50),
# #                avg_sens = rep(seq(log(0.75/(1-0.75)),log(0.25/(1-0.25)),length.out=200),50)
# #         ) %>%
# #         # mutate(ind_sens = log(avg_sens/(1-avg_sens))) %>%
# #         mutate(ind_sens = avg_sens+diff) %>%
# #         mutate(ind_sens= 1/(1+exp(-ind_sens))) %>%
# #         mutate(avg_sens= 1/(1+exp(-avg_sens)))
# # 
# # sens_estim <- sens_obj %>%
# #         group_by(days_ago)%>%
# #         summarize(sens=mean(ind_sens))%>%
# #         pull(sens) %>%
# #         mean()
# # 
# # #get specificity (not vaxxed)
# # spec_obj <- data.frame(
# #         neg_ID= paste0("Y",1:50),
# #         avg_spec = 0.95,
# #         diff=rnorm(50,sd=1)
# # )  %>%
# #         mutate(ind_spec = log(avg_spec/(1-avg_spec))) %>%
# #         mutate(ind_spec = ind_spec+diff) %>%
# #         mutate(ind_spec= 1/(1+exp(-ind_spec)))
# # spec_estim <- mean(spec_obj$ind_spec)
# # 
# # 
# # #get specificity (vaxxed)
# # vax_obj <- data.frame(
# #         vax_ID= paste0("X",1:50),
# #         avg_spec = 0.75,
# #         diff=rnorm(50,sd=1)
# # )  %>%
# #         mutate(ind_spec = log(avg_spec/(1-avg_spec))) %>%
# #         mutate(ind_spec = ind_spec+diff) %>%
# #         mutate(ind_spec= 1/(1+exp(-ind_spec)))
# # 
# # 
# # 
# # run_1 <- simulate_serosurvey_novax(n_sim = 10,
# #                                    n_survey = 1000,
# #                                    sens_df=sens_obj,spec_df=spec_obj,
# #                                    incidence=0.01)
# # 
# # run_5 <- simulate_serosurvey_novax(n_sim = 10,
# #                                    n_survey = 1000,
# #                                    sens_df=sens_obj,spec_df=spec_obj,
# #                                    incidence=0.05)
# # 
# # run_10 <- simulate_serosurvey_novax(n_sim = 10,
# #                                     n_survey = 1000,
# #                                     sens_df=sens_obj,spec_df=spec_obj,
# #                                     incidence=0.10)
# # run_20 <- simulate_serosurvey_novax(n_sim = 10,
# #                                     n_survey = 1000,
# #                                     sens_df=sens_obj,spec_df=spec_obj,
# #                                     incidence=0.20)
# # 
# # 
# # bind_rows(
# #         bind_rows(run_1,.id="simulation"),
# #         bind_rows(run_5,.id="simulation"),
# #         bind_rows(run_10,.id="simulation"),
# #         bind_rows(run_20,.id="simulation")
# # ) %>%
# #         group_by(simulation,incidence) %>%
# #         summarize(
# #                 crude_seropos=mean(seropos),
# #                 true_inc=mean(R)
# #         ) %>%
# #         mutate(adjusted_seropos=(crude_seropos + spec_estim-1)/(sens_estim + spec_estim - 1)) %>%
# #         mutate(adjusted_seropos=ifelse(adjusted_seropos<0,0,adjusted_seropos)) %>%
# #         gather(type,percent,-c(simulation,incidence)) %>%
# #         ggplot(aes(y=percent,x=type,col=type))+
# #         ggbeeswarm::geom_beeswarm()+
# #         geom_hline(aes(yintercept=incidence),lty=2)+
# #         facet_wrap(.~incidence,nrow=1)
# # 
# # 
# # run_1 <- simulate_serosurvey_novax(n_sim = 10,
# #                                    n_survey = 1000,
# #                                    sens_df=sens_obj,spec_df=spec_obj,
# #                                    incidence=0.01)
# # 
# # run_5 <- simulate_serosurvey_novax(n_sim = 10,
# #                                    n_survey = 1000,
# #                                    sens_df=sens_obj,spec_df=spec_obj,
# #                                    incidence=0.05)
# # 
# # run_10 <- simulate_serosurvey_novax(n_sim = 10,
# #                                     n_survey = 1000,
# #                                     sens_df=sens_obj,spec_df=spec_obj,
# #                                     incidence=0.10)
# # run_20 <- simulate_serosurvey_novax(n_sim = 10,
# #                                     n_survey = 1000,
# #                                     sens_df=sens_obj,spec_df=spec_obj,
# #                                     incidence=0.20)
# # 
# # 
# # bind_rows(
# #         bind_rows(run_1,.id="simulation"),
# #         bind_rows(run_5,.id="simulation"),
# #         bind_rows(run_10,.id="simulation"),
# #         bind_rows(run_20,.id="simulation")
# # ) %>%
# #         group_by(simulation,incidence) %>%
# #         summarize(
# #                 crude_seropos=mean(seropos),
# #                 true_inc=mean(R)
# #         ) %>%
# #         mutate(adjusted_seropos=(crude_seropos + spec_estim-1)/(sens_estim + spec_estim - 1)) %>%
# #         mutate(adjusted_seropos=ifelse(adjusted_seropos<0,0,adjusted_seropos)) %>%
# #         gather(type,percent,-c(simulation,incidence)) %>%
# #         ggplot(aes(y=percent,x=type,col=type))+
# #         ggbeeswarm::geom_beeswarm()+
# #         geom_hline(aes(yintercept=incidence),lty=2)+
# #         facet_wrap(.~incidence,nrow=1)
# # 
# # 
# # #vax simulations
# # 
# # vax_run_25 <- simulate_serosurvey_vax(n_sim = 10,
# #                                       n_survey = 1000,
# #                                       sens_df=sens_obj,
# #                                       spec_df=spec_obj,
# #                                       vax_df = vax_obj,
# #                                       coverage=0.25,
# #                                       incidence=0.05
# # )
# # 
# # vax_run_50 <- simulate_serosurvey_vax(n_sim = 10,
# #                                       n_survey = 1000,
# #                                       sens_df=sens_obj,
# #                                       spec_df=spec_obj,
# #                                       vax_df = vax_obj,
# #                                       coverage=0.5,
# #                                       incidence=0.05
# # )
# # 
# # vax_run_75 <- simulate_serosurvey_vax(n_sim = 10,
# #                                       n_survey = 1000,
# #                                       sens_df=sens_obj,
# #                                       spec_df=spec_obj,
# #                                       vax_df = vax_obj,
# #                                       coverage=0.75,
# #                                       incidence=0.05
# # )
# # 
# # 
# # 
# # bind_rows(
# #         bind_rows(vax_run_25,.id="simulation"),
# #         bind_rows(vax_run_50,.id="simulation"),
# #         bind_rows(vax_run_75,.id="simulation")
# # )%>%
# #         group_by(simulation,incidence,coverage) %>%
# #         summarize(
# #                 crude_seropos=mean(seropos),
# #                 true_inc=mean(R)
# #         ) %>%
# #         mutate(adjusted_seropos=(crude_seropos + spec_estim-1)/(sens_estim + spec_estim - 1)) %>%
# #         mutate(adjusted_seropos=ifelse(adjusted_seropos<0,0,adjusted_seropos)) %>%
# #         gather(type,percent,-c(simulation,incidence,coverage)) %>%
# #         ggplot(aes(y=percent,x=type,col=type))+
# #         ggbeeswarm::geom_beeswarm()+
# #         geom_hline(aes(yintercept=incidence),lty=2)+
# #         facet_wrap(.~coverage,nrow=1)
# 
# 
# # new_vax_tv_obj %>%
# #         ggplot(aes(x=days_ago,y=ind_spec,group=vax_ID))+
# #         geom_line()+
# #         ylab("Specificity")+
# #         xlab("Days since vaccination")
# 
# # new_vax_21 <- new_vax_tv_obj %>%
# #         filter(days_ago==21)%>%
# #         select(-days_ago)
# # 
# # new_vax_45 <- new_vax_tv_obj %>%
# #         filter(days_ago==45)%>%
# #         select(-days_ago)
# # 
# # new_vax_90 <- new_vax_tv_obj %>%
# #         filter(days_ago==90)%>%
# #         select(-days_ago)
# # 
# # new_vax_120 <- new_vax_tv_obj %>%
# #         filter(days_ago==120)%>%
# #         select(-days_ago)
# # 
# # new_vax_180 <- new_vax_tv_obj %>%
# #         filter(days_ago==180)%>%
# #         select(-days_ago)
# 
# # vax21_run_25 <- simulate_serosurvey_vax(n_sim = 100,
# #                                         n_survey = 1000,
# #                                         sens_df=new_sens_obj,
# #                                         spec_df=new_spec_obj,
# #                                         vax_df = new_vax_21,
# #                                         coverage=0.25,
# #                                         incidence=0.05
# # ) 
# # 
# # vax21_run_50 <- simulate_serosurvey_vax(n_sim = 10,
# #                                         n_survey = 1000,
# #                                         sens_df=new_sens_obj,
# #                                         spec_df=new_spec_obj,
# #                                         vax_df = new_vax_21,
# #                                         coverage=0.5,
# #                                         incidence=0.05
# # )
# # 
# # vax21_run_75 <- simulate_serosurvey_vax(n_sim = 10,
# #                                         n_survey = 1000,
# #                                         sens_df=new_sens_obj,
# #                                         spec_df=new_spec_obj,
# #                                         vax_df = new_vax_21,
# #                                         coverage=0.75,
# #                                         incidence=0.05
# # )
# # 
# # vax120_run_25 <- simulate_serosurvey_vax(n_sim = 100,
# #                                          n_survey = 1000,
# #                                          sens_df=new_sens_obj,
# #                                          spec_df=new_spec_obj,
# #                                          vax_df = new_vax_120,
# #                                          coverage=0.25,
# #                                          incidence=0.05
# # ) 
# # 
# # vax120_run_50 <- simulate_serosurvey_vax(n_sim = 100,
# #                                          n_survey = 1000,
# #                                          sens_df=new_sens_obj,
# #                                          spec_df=new_spec_obj,
# #                                          vax_df = new_vax_120,
# #                                          coverage=0.5,
# #                                          incidence=0.05
# # )
# # 
# # vax120_run_75 <- simulate_serosurvey_vax(n_sim = 100,
# #                                          n_survey = 1000,
# #                                          sens_df=new_sens_obj,
# #                                          spec_df=new_spec_obj,
# #                                          vax_df = new_vax_120,
# #                                          coverage=0.75,
# #                                          incidence=0.05
# # )
# # 
# # 
# # bind_rows(
# #         bind_rows(vax21_run_25,.id="simulation") %>% mutate(campaign_day="21 days prior"),
# #         bind_rows(vax21_run_50,.id="simulation")%>% mutate(campaign_day="21 days prior"),
# #         bind_rows(vax21_run_75,.id="simulation")%>% mutate(campaign_day="21 days prior"),
# #         bind_rows(vax120_run_25,.id="simulation")%>% mutate(campaign_day="120 days prior"),
# #         bind_rows(vax120_run_50,.id="simulation")%>% mutate(campaign_day="120 days prior"),
# #         bind_rows(vax120_run_75,.id="simulation")%>% mutate(campaign_day="120 days prior")
# # )%>%
# #         group_by(simulation,incidence,coverage,campaign_day) %>%
# #         summarize(
# #                 crude_seropos=mean(seropos),
# #                 true_inc=mean(R)
# #         ) %>%
# #         mutate(adjusted_seropos=(crude_seropos + new_spec_estim-1)/(new_sens_estim + new_spec_estim - 1)) %>%
# #         mutate(adjusted_seropos=ifelse(adjusted_seropos<0,0,adjusted_seropos)) %>%
# #         mutate(adjusted_seropos=ifelse(adjusted_seropos>1,1,adjusted_seropos)) %>%
# #         gather(type,percent,-c(simulation,incidence,coverage,campaign_day)) %>%
# #         ggplot(aes(y=percent,x=type,col=type))+
# #         ggbeeswarm::geom_beeswarm()+
# #         geom_hline(aes(yintercept=incidence),lty=2)+
# #         facet_grid(campaign_day~coverage)+
# #         theme(axis.text.x = element_text(angle=45,hjust=1))+
# #         ylab("Seroprevalence")
# 
# 
# # new_spec_obj %>%
# #         ggplot(aes(x=ind_spec))+
# #         geom_histogram()+
# #         xlab("specificity")
# # 
# # new_spec_estim <- mean(new_spec_obj$ind_spec)
# # 
# # 
# # new_spec_obj$ind_spec %>% hist()
# # new_sens_obj %>%
# #         ggplot(aes(x=days_ago,y=ind_sens,group=pos_ID))+
# #                 geom_line() +
# #         xlab("Days since infection")+
# #         ylab("sensitivity")
# # 
# # new_sens_estim <- new_sens_obj %>%
# #         group_by(days_ago)%>%
# #         summarize(sens=mean(ind_sens))%>%
# #         pull(sens) %>%
# #         mean()
