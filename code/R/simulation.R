rm(list=ls())

library(tidyverse)
library(tidybayes)

source("code/R/packages.R")
source("code/R/utils.R")

options(mc.cores = 1)


#load loocv objects
loocv_spec_list <- read_rds(
        paste0("data/generated_data/analysis_objects/loocv/","loocv_spec_list_new.rds")
)
loocv_sens_list <- read_rds("data/generated_data/analysis_objects/loocv/loocv_tvfit_list_new.rds")

loocv_df <- read_rds("data/generated_data/analysis_objects/loocv/loocv_df_new.rds")

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
original_sens_fitdata <-loocv_sens_list$reduced_2class$case_fitdata$`200-day`$Case$data
original_sens_spread <- loocv_sens_list$reduced_2class$case_fitdata$`200-day`$Case$fit %>% 
        spread_draws(`(Intercept)`,
                     b[term,group],
                     day1,day2,day3
        )


# original_sens_fitdata <-loocv_sens_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$Case$data
# original_sens_spread <- loocv_sens_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$Case$fit %>% 
#         spread_draws(`(Intercept)`,
#                      b[term,group],
#                      day1,day2,day3
#         )

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

original_spec_obj<- loocv_spec_list$reduced_2class$case_fitdata$`200-day`$`Outside Window`$fit%>% 
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

# original_spec_obj<- loocv_spec_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$`Outside Window`$fit%>% 
#         spread_draws(`(Intercept)`,
#                      b[term,group]) %>%
#         mutate(logit_theta_j=`(Intercept)` + b
#         ) %>%
#         mutate(theta_j = 1/(1+exp(-logit_theta_j))) %>%
#         group_by(group)%>%
#         summarize(theta=mean(theta_j),n()) %>%
#         rename(neg_ID=group
#         ) %>%
#         mutate(ind_spec=1-theta)

#get specificity distribution for seroincidence estimation

original_spec_dist <- loocv_spec_list$reduced_2class$case_fitdata$`200-day`$`Outside Window`$fit%>% 
        spread_draws(`(Intercept)`,
                     b[term,group]) %>%
        mutate(logit_theta_j=`(Intercept)` + b
        ) %>%
        mutate(theta_j = 1/(1+exp(-logit_theta_j)))  %>%
        group_by(.draw) %>%
        summarize(spec=1-mean(theta_j))

# original_spec_dist <- loocv_spec_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$`Outside Window`$fit%>% 
#         spread_draws(`(Intercept)`,
#                      b[term,group]) %>%
#         mutate(logit_theta_j=`(Intercept)` + b
#         ) %>%
#         mutate(theta_j = 1/(1+exp(-logit_theta_j)))  %>%
#         group_by(.draw) %>%
#         summarize(spec=1-mean(theta_j))

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
alternative_sens_fitdata <-loocv_sens_list$reduced_2class$all_fitdata$`200-day`$Case$data
alternative_sens_spread <- loocv_sens_list$reduced_2class$all_fitdata$`200-day`$Case$fit %>% 
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


#### Specificity (non-vaccinees)

#get individual specificities for simulation
# alternative_spec_obj<- new_fit$fit%>% 
#         spread_draws(`(Intercept)`,
#                      b[term,group]
#         ) %>%
#         mutate(id=str_remove(group,"id\\:")) %>%
#         left_join(distinct(new_fit$data,id,status)) %>%
#         filter(status!="Vaccinee")%>%
#         mutate(logit_theta_j=`(Intercept)` + b ) %>%
#         mutate(theta_j = 1/(1+exp(-logit_theta_j))) %>%
#         group_by(group)%>%
#         summarize(theta=mean(theta_j),n()) %>%
#         rename(neg_ID=group
#         ) %>%
#         mutate(ind_spec=1-theta)


alternative_spec_obj <- loocv_spec_list$reduced_2class$all_fitdata$`200-day`$`Outside Window`$fit%>% 
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
alternative_spec_dist <- loocv_spec_list$reduced_2class$all_fitdata$`200-day`$`Outside Window`$fit%>% 
        spread_draws(`(Intercept)`,
                     b[term,group]) %>%
        mutate(logit_theta_j=`(Intercept)` + b
        ) %>%
        mutate(theta_j = 1/(1+exp(-logit_theta_j)))  %>%
        group_by(.draw) %>%
        summarize(spec=1-mean(theta_j))


# #### Specificity (vaccinees)
# 
# #load stan fit
# alternative_vax_obj <- new_fit$fit%>%
#         spread_draws(`(Intercept)`,
#                      b[term,group],
#                      vax
#         ) %>%
#         mutate(id=str_remove(group,"id\\:")) %>%
#         left_join(distinct(new_fit$data,id,status)) %>%
#         filter(status=="Vaccinee")%>%
#         mutate(logit_theta_j=`(Intercept)` + b + vax) %>%
#         mutate(theta_j = 1/(1+exp(-logit_theta_j))) %>%
#         group_by(group)%>%
#         summarize(theta=mean(theta_j),n()) %>%
#         rename(vax_ID=group
#         ) %>%
#         mutate(ind_spec=1-theta)


#load stan fit
alternative_vax_fitdata <-loocv_sens_list$reduced_2class$all_fitdata$`200-day`$Vaccinee$data
alternative_vax_spread <- loocv_sens_list$reduced_2class$all_fitdata$`200-day`$Vaccinee$fit %>% 
        spread_draws(`(Intercept)`,
                     b[term,group],
                     day1,day2,day3
        )

# #for loop through potential campaign timings
# campaign_timepoints <- c(21,120)#c(21,45,90,120,180)

alternative_vax_tv_obj <- data.frame()
alternative_vax_beta_draws <- data.frame()
# alternative_vax_beta <- data.frame()


for(t in campaign_timepoints){
        
        cat(t,"\n")
        time <- log(t) -mean(log(alternative_vax_fitdata$new_day))
        
        #get individual specificities for simulation
        tmp_df <- alternative_vax_spread %>%
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
        
        alternative_vax_tv_obj <- bind_rows(alternative_vax_tv_obj,tmp_df)
        
        #get specificity distribution for seroincidence estimation
        tmp_dist <- alternative_vax_spread %>%
                mutate(logit_theta_j=`(Intercept)` + b+
                               day1*time +
                               day2*time^2+
                               day3*time^3
                       
                ) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                group_by(.draw) %>%
                summarize(theta=1-mean(theta_j),n()) %>% #average across individuals here
                mutate(days_ago=t)
        
        alternative_vax_beta_draws <- bind_rows(alternative_vax_beta_draws,tmp_dist)
        
}




# 3. Run simulations and stan models to get seroincidence estimates -----------

simulation_df <-data.frame()
sims_original <- list()
sims_alternative <- list()

lambda_df <- data.frame()

simulations <- 100
n_survey <- 1000
coverage_values <- c(0, 0.25, 0.5, 0.75)


stan_model_compiled <- stan_model("code/stan/vax-model-simplified.stan")

true_seroinc <- 0.1
# for(true_seroinc in c(0.05,0.1,0.2)){ #choose several values of true seroincidence to look at
for(cov in coverage_values){
        
        
        #simulate the true values
        tmp_sim <- sim_truth(simulations,n_survey,coverage = cov, incidence = true_seroinc)
        
        #simulate seropositivity
        for(t in campaign_timepoints){
                
                cat("Coverage: ",cov, " Timing: ", t,"\n")
                
                ### simulate the original model seropositivity
                sims_original[[as.character(t)]] <- tmp_sim %>% lapply(make_seropos,model_name = "Original Model",
                                                                       sens_df= original_sens_obj,
                                                                       spec_df= original_spec_obj,
                                                                       vax_df = original_vax_tv_obj,
                                                                       vax_day=t
                ) %>%  bind_rows(.id="simulation") 
                
                
                ### simulate the alternative model seropositivity
                sims_alternative[[as.character(t)]]  <- tmp_sim %>% lapply(make_seropos,model_name = "Alternative Model",
                                                                           sens_df= alternative_sens_obj,
                                                                           spec_df= alternative_spec_obj,
                                                                           vax_df = alternative_vax_tv_obj,
                                                                           vax_day=t
                ) %>%  bind_rows(.id="simulation")
                
                
                #store the seropositivity
                simulation_df <- bind_rows(
                        simulation_df,
                        sims_original[[as.character(t)]] %>%
                                group_by(simulation, coverage, campaign_day,incidence) %>%
                                summarize(truth=mean(R),seropos=mean(`Original Model`)) %>% 
                                mutate(model="Original"),
                        
                        sims_alternative[[as.character(t)]] %>%
                                group_by(simulation, coverage, campaign_day,incidence) %>%
                                summarize(truth=mean(R),seropos=mean(`Alternative Model`)) %>% 
                                mutate(model="Alternative")
                        
                )
                
                
                
        }
        
        
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
                        
                        
                        ### Strategy 3: Alternative model adjustment
                        
                        this_count_vec <- c(0,0,0,0)
                        
                        this_count <- sims_alternative[[as.character(t)]] %>%
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
                        alternative_spec_beta <- bind_cols(
                                alternative_vax_beta_draws %>%
                                        filter(days_ago==t) %>%
                                        select(Vaccinee=theta),
                                alternative_spec_dist %>%
                                        select(`Outside Window`=spec)
                        ) %>%
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
                                bind_rows(
                                        spread_draws(biased_fit,lambda) %>% mean_qi() %>%
                                                select( lambda,.lower,.upper) %>%
                                                mutate(adjustment="Ignore"),
                                        spread_draws(questionnaire_fit,lambda) %>% mean_qi() %>%
                                                select( lambda,.lower,.upper) %>%
                                                mutate(adjustment="Questionnaire"),
                                        spread_draws(alternative_fit,lambda) %>% mean_qi() %>%
                                                select( lambda,.lower,.upper) %>%
                                                mutate(adjustment="Alternative") 
                                        
                                ) %>%
                                        mutate(coverage = cov,
                                               campaign_day=t,
                                               simulation=i
                                        )
                        )
                        
                        
                }
                
                
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






