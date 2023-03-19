
# ANALYSIS FUNCTIONS ------------------------



# Wrapper function that:
# stores arguments
# creates an incidence and case epidemic curve
# generates serological data from a decay model fit
# fits a random forest model to training set
# applies RF model to get predictions from serosurvey
doSerosurvey<- function(
        n_pop, # size of the serosurvey
        serosurvey_date= as.Date("2020-01-01"), #date of serosurvey
        
        N_pop = 600000, #population size,
        years_prior = 10,# number of years back you want to simulate
        true_inc_rate = 0.2, #annual incidence rate
        detect_rate =0.1, #detection rate
        epidemic_type = "flat", #specify flat or sin epidemic
        offset=-90,
        freq = 1,
        window_prop = 0.5,#proportion of samples within the window to train RF on
        variance_factor = 1
        
        
){
        
        cat("Store arguments \n")
        argument_list <- list(
                
                n_pop=n_pop, # size of the serosurvey
                serosurvey_date= serosurvey_date,
                
                N_pop = N_pop, #population size,
                years_prior = years_prior,#dates
                true_inc_rate = true_inc_rate, #annual incidence rate
                detect_rate =detect_rate, #detection rate
                epidemic_type = epidemic_type, #specify flat or sin epidemic
                offset=offset,
                freq = freq,
                window_prop = window_prop,#proportion of samples within the window to train RF on
                variance_factor = variance_factor
        )
        
        #generate epi curve
        cat("Generate epi curve \n")
        inc_df <- getIncidence( N_pop=N_pop,
                                serosurvey_date=serosurvey_date,
                                years_prior=years_prior,
                                true_inc_rate=true_inc_rate,
                                detect_rate=detect_rate,
                                epidemic_type=epidemic_type,
                                offset=offset,
                                freq=freq
        )
        
        #generate serological data for serosurvey and training set
        cat("Generate serological data \n")
        inf_time_training_source<- getInfectionTimes(n_pop=2000,
                                              inc_df=inc_df,
                                              true_inc_rate=true_inc_rate)
        
        
        over52weeks <- filter(inf_time_training_source,day>(52*7)) %>% pull(day) %>%
                                sample(round(1000*(1-window_prop)),replace = FALSE)

        #over samplewithin year 1 for training set
        inf_time_training <- data.frame(
                        day = c(runif(round(1000*window_prop),7*(0:51),7*(1:52)) %>%
                                  ceiling(),over52weeks)
                        ) %>%
                      mutate(ss_id= 1:n())  
        
        inf_time_test  <-getInfectionTimes(n_pop=n_pop,
                                           inc_df=inc_df,
                                           true_inc_rate=true_inc_rate)
        
        sero_data <- list()
        sero_data$training <- getSeroData(inf_time_df=inf_time_training,
                                          variance_factor=variance_factor)%>%
                                                #calculate week as well
                                                mutate(week=ceiling(day/7)) %>%
                                                mutate(week=ifelse(week<53,week,53))
        sero_data$test <-     getSeroData(inf_time_df=inf_time_test,
                                          variance_factor=variance_factor)
        
        
        #fit random forest model
        cat("Fit RF model \n")
        rf_fit <- randomForrest_fit_predict(train=sero_data$train,
                                            test=sero_data$test)
        
        #object outputting everything
        out <- list(
                arguments = argument_list,
                incidence=inc_df,
                serology =sero_data,
                fit =  rf_fit$fit,
                predictions= rf_fit$predictions
                
        )
        
        return(out)
        
}


# create shape of epidemic
generateDailyInfection <- function(x,
                                   true_inc_rate,
                                   N_pop,
                                   epidemic_type="flat",
                                   offset=-90,
                                   freq = 1
){
        if(epidemic_type=="flat") {
                haz <- true_inc_rate/365.25 *N_pop * rep(1,length(x))
                rpois(length(haz),haz)
        } else if(epidemic_type=="sin"){
                # haz <- true_inc_rate/365.25 *N_pop * (sin(2*pi*(x+offset)/(365.25*freq))+1)
                haz <- true_inc_rate/365.25 *N_pop * (sin(2*pi*(x/(365.25*freq)+offset/(365.25)))+1)
                rpois(length(haz),haz)
        } else{
                stop("flat or sin not selected") 
        }
        
        
}



#create the infection epidemic curve
getIncidence<- function(
        N_pop = 600000, #population size,
        serosurvey_date= as.Date("2020-01-01"),
        years_prior = 10,#dates
        true_inc_rate  = 0.05, #annual incidence rate
        detect_rate =0.1, #detection rate
        epidemic_type = "flat", #specify flat or sin epidemic
        offset=-90, #in days
        freq = 1 #number of years per peak
){
        
        day1 <- serosurvey_date-1
        last_day <- as.Date(
                paste0(lubridate::year(serosurvey_date)-years_prior,"-",
                       format(serosurvey_date,"%m-%d")
                ))
        
        data.frame(
                date=seq.Date(day1,last_day,by=-1) 
        ) %>% 
                #year/days/weeks prior to serosurvey
                mutate(year=lubridate::year(date)) %>%
                mutate(years_before_survey=max(year)-year)%>%
                mutate(days_before_survey=1:nrow(.)) %>%
                mutate(weeks_before_survey=ceiling(days_before_survey/7)) %>%
                #generate infections
                mutate(infections=generateDailyInfection(days_before_survey,
                                                         true_inc_rate = true_inc_rate,
                                                         N_pop=N_pop,
                                                         epidemic_type = epidemic_type,
                                                         offset=offset,
                                                         freq=freq
                )) %>%
                # #generate cases (may want to include delays and misclassification)
                mutate(cases=rpois(nrow(.),infections*detect_rate))
        
        
}

#get infection times from incidence curve
getInfectionTimes <- function(n_pop, inc_df, true_inc_rate){
        # figure out last year of infection
        year_lastinfect <- rgeom(n_pop,true_inc_rate)
        max_year <- max(inc_df$years_before_survey)
        year_lastinfect[year_lastinfect>max_year] <- max_year
        
        #get specific day of infection within the year
        times_df <- data.frame()
        for (i in (1:n_pop)){
                
                #limit to these year
                tmp_df <- inc_df %>% filter(years_before_survey==year_lastinfect[i])
                #assign a time of last infection
                poss_times <- tmp_df %>% pull(days_before_survey)
                time_weights <-  tmp_df %>% pull(infections)
                time_weights <- time_weights/sum(time_weights) #normalize
                time_tmp <- sample(poss_times,1,prob=time_weights,replace = FALSE)
                
                # get predicted serosurvey data
                times_df <- bind_rows(times_df,
                                      data.frame(day=time_tmp) %>%
                                              mutate(ss_id=i)
                )
        }
        
        return(times_df)
        
}



# #generate serological data for dataframe with infection times
# getSeroData <- function(inf_time_df,variance_factor){
#         
#         #bring in multivariate fit
#         multi_fit <- read_rds(paste0("source/final_code/epidemic_reconstruction/generated_rds/2021-11-15-multi_fit.rds"))
#         
#         
#         #bring in multivariate model
#         extract <- extract(multi_fit$fit) %>% data.frame()
#         
#         #select parameter values for each individual
#         sample_params <-  bind_cols(inf_time_df, #zero uncertainty in parameter value
#                                     extract[2431,]) 
#         # sample_params <- bind_cols(inf_time_df, #uncertainty in parameter value from sampling
#         # extract[sample(1:nrow(extract),size=nrow(inf_time_df),replace = TRUE),])
#         
#         
#         #make each row be a id marker pair
#         marker_params <- sample_params[,str_detect(colnames(sample_params),
#                                                    "mu|sigma|delta|ss_id|day")] %>%
#                 gather(parameter,value,-c(ss_id,day)) %>%
#                 mutate(marker_number=str_extract(parameter,".[0-9]")) %>%
#                 mutate(marker_number=str_remove(marker_number,".")) %>%
#                 mutate(new_parameter=str_remove(parameter,".[0-9]")) %>%
#                 select(-parameter) %>%
#                 spread(new_parameter,value)
#         
#         #estimate individual baseline/boost variation for each marker
#         ind_matrix <- data.frame()
#         for(i in 1:nrow(marker_params)){
#                 
#                 #get covariance matrix for a particular marker's boost and baseline together
#                 var_mat <- matrix( #0,# set at zero temporarily
#                         c(marker_params$params_sigma.1.1[i],
#                           marker_params$params_sigma.1.2[i],
#                           marker_params$params_sigma.2.1[i],
#                           marker_params$params_sigma.2.2[i]),
#                         nrow=2,
#                         ncol=2)
#                 
#                 #reduce variance for test case
#                 var_mat <- var_mat * variance_factor
#                 
#                 #generate boost and baseline values for the marker
#                 value_vec <- mvtnorm::rmvnorm(1, mean=c(marker_params$mu_omega[i],
#                                                         marker_params$mu_lambda[i]),
#                                               sigma=var_mat) #%>% as.data.frame()
#                 #store in a dataframe
#                 tmp_df <- data.frame(ss_id=marker_params$ss_id[i],
#                                      marker_number=marker_params$marker_number[i],
#                                      value_vec
#                 )
#                 #store in a dataframe
#                 ind_matrix <- bind_rows(ind_matrix,tmp_df)
#                 
#         }
#         colnames(ind_matrix) <- c("ss_id","marker_number","omega_ind","lambda_ind")
#         
#         
#         #estimate the true value for an individual
#         marker_params <- left_join(marker_params,ind_matrix) %>%
#                 mutate(MU=omega_ind+lambda_ind*exp(-delta*(day-5)))%>%
#                 mutate(expected=paste0("MU.",marker_number)) %>%
#                 select(expected,MU,ss_id,day) %>%
#                 spread(expected,MU) %>%
#                 #bring in measurement error
#                 left_join(sample_params[,str_detect(colnames(sample_params),"SIGMA|ss_id")])
#         
#         
#         
#         #connect MU and SIGMA values to generate values
#         MU_string <- str_detect(colnames(marker_params),"MU")
#         SIGMA_string <- str_detect(colnames(marker_params),"SIGMA")
#         n_markers <- multi_fit$data$K
#         
#         #generate values incorporating measurement error
#         new_values <- data.frame()
#         for (i in 1:nrow(marker_params)){
#                 
#                 #get covariance matrix for an individual's marker together
#                 var_mat <-matrix(nrow=n_markers,ncol=n_markers,
#                                  marker_params[i,SIGMA_string] %>% unlist()
#                 )
#                 #reduce variance for test case
#                 var_mat <- var_mat * variance_factor
#                 
#                 #generate marker values
#                 value_vec <- mvtnorm::rmvnorm(1,marker_params[i,MU_string] %>% unlist(),
#                                               sigma=var_mat)
#                 
#                 #store in a dataframe
#                 tmp_df <- data.frame(ss_id=marker_params$ss_id[i],
#                                      day=marker_params$day[i],
#                                      value_vec
#                 )
#                 new_values <- bind_rows(new_values,tmp_df)
#                 
#         }
#         
#         # colnames(new_values) <- c("ss_id","day","Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA")
#         # colnames(new_values) <- c("ss_id","day","Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA",
#         #                           "Luminex_IgA_CtxB","Luminex_IgA_InabaOSPBSA","Luminex_IgA_OgawaOSPBSA",
#         #                           "Luminex_IgM_OgawaOSPBSA","Luminex_IgM_InabaOSPBSA")
#         colnames(new_values) <- c("ss_id","day",
#                                   "Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA",
#                                   "Luminex_IgA_CtxB","Luminex_IgA_OgawaOSPBSA",
#                                   "Luminex_IgM_OgawaOSPBSA")
#         
#         
#         return(new_values)
# }        
        


#get parameters for a set number of children and adults
getSeroParams <- function(n_over10, n_under10){
        
        #markers
        markers <- c(
                "Luminex_IgG_OgawaOSPBSA",
                "Luminex_IgM_OgawaOSPBSA",
                "Luminex_IgA_OgawaOSPBSA",
                "Luminex_IgG_CtxB",
                "Luminex_IgA_CtxB"
        )
        
        #bring in multivariate fit
        original <- read_rds("source/final_code/epidemic_reconstruction/generated_rds/5markermodel.rds")
        
        #CVC matrix
        cvcMAT<- original %>% spread_draws(SIGMA[r,c])%>%
                median_qi() %>%
                select(r,c,SIGMA) %>%
                spread(c,SIGMA) %>%
                select(-r) %>%
                as.matrix()
        
        #over 10 parameters
        Over10_params<- original %>% gather_draws(mu_omega[marker],
                                                  mu_loglambda[marker]
        ) %>%
                median_qi() %>%
                mutate(`Under 10`=0L)%>%
                mutate(.variable=factor(.variable,levels=c("mu_omega","mu_loglambda"))) %>%
                arrange(.variable,marker)
        
        #under 10 parameters
        Under10_params<- original %>% spread_draws(mu_omega[marker],
                                                   mu_loglambda[marker],
                                                   beta_omega[marker],
                                                   beta_loglambda[marker]
        ) %>%
                mutate(mu_omega=mu_omega+beta_omega,
                       mu_loglambda=mu_loglambda+beta_loglambda
                )  %>%
                select(.iteration,mu_omega,mu_loglambda,marker) %>%
                gather(.variable,.value,-c(.iteration,marker))%>%
                group_by(.variable,marker) %>%
                summarize(.value=median(.value))%>%
                mutate(`Under 10`=1L)%>%
                mutate(.variable=factor(.variable,levels=c("mu_omega","mu_loglambda")))%>%
                arrange(.variable,marker)
        
        #sigma parameter
        sigma <-  original %>% spread_draws(sigma[marker]) %>%
                median_qi() %>%
                select(marker,sigma)
        
        #delta parameters
        deltas <- original %>%
                gather_draws(delta_x0[marker],delta_x1[marker]
                ) %>%
                median_qi() %>%
                mutate(`Under 10`=as.integer(str_remove(.variable,"delta_x"))) %>%
                select(marker,delta=.value,`Under 10`)
        
        # V data.frame
        v_df <- data.frame(variable=paste0("V",1:(2*length(markers))),
                           marker= rep(1:length(markers),2),
                           parameter= c(rep("omega_ind",length(markers)),rep("loglambda_ind",length(markers)))
        )
        
        
        #get parameters together
        sim_param <- bind_rows(
                #generate individual parameters for adults
                mvtnorm::rmvnorm(n_over10,Over10_params$.value,cvcMAT) %>%
                        as.data.frame() %>%
                        mutate(`Under 10`=0L),
                #generate individual parameters for children
                mvtnorm::rmvnorm(n_under10,Under10_params$.value,cvcMAT) %>%
                        as.data.frame() %>%
                        mutate(`Under 10`=1L)
        ) %>%
                #provide a simulation id
                mutate(ss_id=1:n())%>%
                #provide names for difference variables
                gather(variable,value,-c(ss_id,`Under 10`)) %>%
                left_join(v_df) %>%
                select(-variable)%>%
                spread(parameter,value)%>%
                #get other parameters ready
                mutate(lambda_ind=exp(loglambda_ind)) %>%
                left_join(deltas) %>%
                mutate(D=ifelse(`Under 10`==1 & marker==1,85L,5L)) %>%
                left_join(sigma) %>%
                #add the marker names
                mutate(full_test=factor(marker,labels=markers))
        
        return(sim_param)

}

#generate serological data from a dataframe that includes infection times
getSeroData <- function(inf_time_df,
                        variance_factor=1,
                        propUnder10
){
        
        little_n <- length(unique(inf_time_df$ss_id))
        
        params <- getSeroParams(
                floor(little_n* (1-propUnder10)),
                ceiling(little_n*propUnder10)
        )
        
        inf_time_df %>% 
                left_join(params) %>% 
                mutate(pred=case_when(
                        day<D & D==5 ~ omega_ind,
                        day<D & D==85 ~ omega_ind + day*lambda_ind/D,
                        day>=D ~omega_ind+lambda_ind*exp(-delta*(day-D))
                )) %>%
                mutate(generated=rnorm(nrow(.),pred,sigma*variance_factor)) %>%
                mutate(generated=ifelse(generated>5.5,5.5,generated))%>%
                mutate(value=ifelse(generated<1,1,generated)) %>%
                select(ss_id,`Under 10`,day,full_test,value) %>%
                spread(full_test,value)
        
}



# #generate serological data for dataframe with infection times
# getSeroData <- function(inf_time_df,variance_factor){
#         
#         #bring in multivariate fit
#         # multi_fit <- read_rds(paste0("source/final_code/epidemic_reconstruction/generated_rds/2021-12-11-multi_fit_censored_9markers.rds"))
#         multi_fit <- read_rds(paste0("source/final_code/epidemic_reconstruction/generated_rds/2021-11-15-multi_fit.rds"))
#         
#         
#         #bring in multivariate model
#         extract <- extract(multi_fit$fit) %>% data.frame()
#         
#         #select parameter values for each individual
#         sample_params <-  bind_cols(inf_time_df, #zero uncertainty in parameter value
#                                     extract[2431,]) 
#         # sample_params <- bind_cols(inf_time_df, #uncertainty in parameter value from sampling
#         # extract[sample(1:nrow(extract),size=nrow(inf_time_df),replace = TRUE),])
#         
#         
#         #make each row be a id marker pair
#         marker_params <- sample_params[,str_detect(colnames(sample_params),
#                                                    "mu|sigma|delta|ss_id|day")] %>%
#                 gather(parameter,value,-c(ss_id,day)) %>%
#                 mutate(marker_number=str_extract(parameter,".[0-9]")) %>%
#                 mutate(marker_number=str_remove(marker_number,".")) %>%
#                 mutate(new_parameter=str_remove(parameter,".[0-9]")) %>%
#                 select(-parameter) %>%
#                 spread(new_parameter,value)
#         
#         #estimate individual baseline/boost variation for each marker
#         ind_matrix <- data.frame()
#         for(i in 1:nrow(marker_params)){
#                 
#                 #get covariance matrix for a particular marker's boost and baseline together
#                 var_mat <- matrix( #0,# set at zero temporarily
#                         c(marker_params$params_sigma.1.1[i],
#                           marker_params$params_sigma.1.2[i],
#                           marker_params$params_sigma.2.1[i],
#                           marker_params$params_sigma.2.2[i]),
#                         nrow=2,
#                         ncol=2)
#                 
#                 #reduce variance for test case
#                 var_mat <- var_mat * variance_factor
#                 
#                 #generate boost and baseline values for the marker
#                 value_vec <- mvtnorm::rmvnorm(1, mean=c(marker_params$mu_omega[i],
#                                                         marker_params$mu_lambda[i]),
#                                               sigma=var_mat) #%>% as.data.frame()
#                 #store in a dataframe
#                 tmp_df <- data.frame(ss_id=marker_params$ss_id[i],
#                                      marker_number=marker_params$marker_number[i],
#                                      value_vec
#                 )
#                 #store in a dataframe
#                 ind_matrix <- bind_rows(ind_matrix,tmp_df)
#                 
#         }
#         colnames(ind_matrix) <- c("ss_id","marker_number","omega_ind","lambda_ind")
#         
#         
#         #estimate the true value for an individual
#         marker_params <- left_join(marker_params,ind_matrix) %>%
#                 mutate(MU=omega_ind+lambda_ind*exp(-delta*(day-5)))%>%
#                 mutate(expected=paste0("MU.",marker_number)) %>%
#                 select(expected,MU,ss_id,day) %>%
#                 spread(expected,MU) %>%
#                 #bring in measurement error
#                 left_join(sample_params[,str_detect(colnames(sample_params),"SIGMA|ss_id")])
#         
#         
#         
#         #connect MU and SIGMA values to generate values
#         MU_string <- str_detect(colnames(marker_params),"MU")
#         SIGMA_string <- str_detect(colnames(marker_params),"SIGMA")
#         n_markers <- multi_fit$data$K
#         
#         #generate values incorporating measurement error
#         new_values <- data.frame()
#         for (i in 1:nrow(marker_params)){
#                 
#                 #get covariance matrix for an individual's marker together
#                 var_mat <-matrix(nrow=n_markers,ncol=n_markers,
#                                  marker_params[i,SIGMA_string] %>% unlist()
#                 )
#                 #reduce variance for test case
#                 var_mat <- var_mat * variance_factor
#                 
#                 #generate marker values
#                 value_vec <- mvtnorm::rmvnorm(1,marker_params[i,MU_string] %>% unlist(),
#                                               sigma=var_mat)
#                 
#                 #store in a dataframe
#                 tmp_df <- data.frame(ss_id=marker_params$ss_id[i],
#                                      day=marker_params$day[i],
#                                      value_vec
#                 )
#                 new_values <- bind_rows(new_values,tmp_df)
#                 
#         }
#         
#         # colnames(new_values) <- c("ss_id","day","Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA")
#         # colnames(new_values) <- c("ss_id","day","Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA",
#         #                           "Luminex_IgA_CtxB","Luminex_IgA_InabaOSPBSA","Luminex_IgA_OgawaOSPBSA",
#         #                           "Luminex_IgM_OgawaOSPBSA","Luminex_IgM_InabaOSPBSA")
#         colnames(new_values) <- c("ss_id","day",
#                                   "Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA",
#                                   "Luminex_IgA_CtxB","Luminex_IgA_OgawaOSPBSA",
#                                   "Luminex_IgM_OgawaOSPBSA")
#         
#         
#         return(new_values)
#         
#         
#         
# }


#fit random forest model and apply to test set
randomForrest_fit_predict <- function(train,test){
        
        #set up model formula
        # markers<- c("Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA")
        # markers <-c("Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA",
        #           "Luminex_IgA_CtxB","Luminex_IgA_InabaOSPBSA","Luminex_IgA_OgawaOSPBSA",
        #           "Luminex_IgM_OgawaOSPBSA","Luminex_IgM_InabaOSPBSA")
        markers <-c(
        "Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA",
        "Luminex_IgA_CtxB","Luminex_IgA_OgawaOSPBSA",
        "Luminex_IgM_OgawaOSPBSA")
        
        timing_formula <- paste0("log_day~",paste(markers,collapse=" + ")) %>%
                as.formula
        
        #fit a randomforest model
        fit_rf <- train %>% mutate(log_day=log(day)) %>%
                randomForest(             formula=timing_formula,
                                          data=.,
                                          ntree=1000
                )
        
        
        #get predictions applied to serosurvey data
        full_preds_rf <- predict(fit_rf,newdata=test,predict.all=TRUE)
        
        # day_cuts <- c(7*(0:52)+1,Inf) #break days into weeks
        
        full_preds_rf_df<- data.frame(exp(full_preds_rf$individual),
                                      day_actual=test$day,
                                      ss_id = test$ss_id
        ) %>%
                gather(tree,prediction,-c(day_actual,ss_id)) %>%
                
                #new method without cuts
                #predicted weeks
                mutate(week=ceiling(prediction/7)) %>%
                mutate(week=ifelse(week<53,week,53)) %>%
                #true week
                mutate(week_actual=ceiling(day_actual/7)) %>%
                mutate(week_actual=ifelse(week_actual<53,week_actual,53)) 
                
                # #old method with cuts
                # #predicted weeks
                # mutate(week=cut(prediction,day_cuts,right=FALSE,
                # )) %>%
                # mutate(week=as.numeric(week)) %>%
                # #true week
                # mutate(week_actual=cut(day_actual,day_cuts,right=FALSE,
                # )) %>%
                # mutate(week_actual=as.numeric(week_actual)) %>%
                
        
        
        out <- list (fit =fit_rf,
                     predictions=full_preds_rf_df
        )
        
        return(out)
}




fitReconstruction <- function(generated_serosurvey, 
                              prior_vector=NULL,
                              model=NULL,
                              autocorrelation=0.2
){
        
        
        daily_incidence <- generated_serosurvey$incidence
        predictions <- generated_serosurvey$predictions
        n_pop <- generated_serosurvey$arguments$n_pop
        N_pop <- generated_serosurvey$arguments$N_pop
        
        #set up stan_input file
        stan_input<- list(
                n=n_pop, #serosurvey size
                N=N_pop, #size of the overall population
                t_max = 53, #number of periods
                #priors on the detection probability
                gamma= 1, 
                delta= 1
        ) 
        
        #take in clinical case data
        stan_input$X <-  daily_incidence %>% group_by(weeks_before_survey) %>%
                summarize(cases=sum(cases)) %>%
                filter(weeks_before_survey<=52) %>%
                pull()
        
        # take in prior or just use training data distribution
        if(is.null(prior_vector)){
                #will only work if there is at least one sample in each catgegory
                stan_input$p_prior <- generated_serosurvey$serology$training %>%
                        group_by(week) %>%
                        count() %>% ungroup() %>%
                        mutate(prop=n/sum(n)) %>%
                        arrange(week) %>%
                        pull(prop)
                        
        } else {
                stan_input$p_prior <- prior_vector
        }
        #variance multiplier for prior infections outside recon period
        stan_input$auto_corr_sigma <-autocorrelation
        
        #use random forest predictions to create probability matrix
        # day_cuts <- c(7*(0:52)+1,Inf) #break days into weeks
        
        #turn the number of trees into probability matrix
        prob_matrix <-  predictions %>% 
                group_by(ss_id,week) %>%
                count() %>%
                group_by(ss_id) %>%
                mutate(prop=n/sum(n))%>%
                select(-n) %>% ungroup()%>%
                spread(ss_id,prop) %>% 
                #get every week included
                right_join(data.frame(week=1:53)) %>% 
                arrange(week) %>%
                select(-week) %>%
                as.matrix()
        prob_matrix[is.na(prob_matrix)] <- 0 # include zero probabilty for trees 
        stan_input$p_mat <- prob_matrix
        
        
        #fit stan model
        if(is.null(model)){
                fileName <- "source/final_code/epidemic_reconstruction/stan/reconstruction_model.stan"
        } else if(model=="serology"){
                        fileName <- "source/final_code/epidemic_reconstruction/stan/reconstruction_model_serologyONLY.stan"
        }
        
        # new stan model that normalizes prior time of infection
        # fileName <- "source/final_code/epidemic_reconstruction/stan/reconstruction_model_new.stan"
        
        recon_fit<- stan(file = fileName, data = stan_input,
                         chains=4, iter=2000
        )
        
        #outputs stan fit and summary statistics
        summary_fit <- summary(recon_fit)$summary %>% data.frame() %>%
                mutate(parameter=rownames(.)) #overall summary statistics
        
        #limits to incidence of infection parameters
        lambdas <- summary_fit %>% filter(str_detect(parameter,"lambda")) %>%
                mutate(week = str_extract(parameter,"\\[(.*?)\\]")) %>%
                mutate(week=str_remove_all(week,"\\[|\\]")) %>%
                mutate(weeks_before_survey=as.numeric(week))
        
        out <- list(
                input=stan_input,
                fit=recon_fit,
                summary = summary_fit,
                lambdas = lambdas
        )
        return(out)
}




# PLOT FUNCTIONS ------------------------


plotEpiCurve <- function(generated,weeks_back=52){
        
        #also get serosurvey samples
        # day_cuts <- c(7*(0:52)+1,max(generated$serology$test$day)+1) #break days into weeks
        # serosamples <- generated$serology$test %>%
        #         #true week
        #         mutate(weeks_before_survey=cut(day,day_cuts,right=FALSE,
        #                              )) %>%
        #         mutate(weeks_before_survey=as.numeric(weeks_before_survey)) 
        
        
        #plot weekly infections and cases for last 10 years
        generated$incidence %>%
                filter(weeks_before_survey<=weeks_back) %>%
                group_by(weeks_before_survey) %>%
                summarize(infections=sum(infections),
                          cases=sum(cases)
                ) %>%
                ggplot(aes(x=weeks_before_survey))+
                geom_point(aes(y=infections),col="red")+
                geom_col(aes(y=cases),fill="blue")+
                ylab("Infections (red) and cases (blue)")+
                scale_x_reverse("Weeks before serosurvey")+
                theme_bw()+
                geom_rug(data=generated$serology$test%>%
                                 filter(day<=weeks_back*7),
                         aes(x=day/7),alpha=0.5
                )
}

plotSerology <- function(truth_wide,generated) {
        
        #what does the data look like
        truth_wide %>% filter(status=="Case") %>%
                select(id,day_actual,
                       Luminex_IgG_CtxB,Luminex_IgG_InabaOSPBSA,Luminex_IgG_OgawaOSPBSA,
                       Luminex_IgA_CtxB,Luminex_IgA_OgawaOSPBSA,
                       Luminex_IgM_OgawaOSPBSA
                ) %>%
                rename(day=day_actual) %>%
                mutate(data_source="original") %>%
                bind_rows(
                        generated$serology$training[sample(1:nrow(generated$serology$training),
                                                           nrow(.),replace = FALSE),] %>% 
                                mutate(data_source="training"),
                        generated$serology$test[sample(1:nrow(generated$serology$test),
                                                       nrow(.),replace = FALSE),] %>%
                                mutate(data_source="test")) %>% 
                select(id,ss_id,day,data_source,
                       Luminex_IgG_CtxB,Luminex_IgG_InabaOSPBSA,Luminex_IgG_OgawaOSPBSA,
                       Luminex_IgA_CtxB,Luminex_IgA_OgawaOSPBSA,
                       Luminex_IgM_OgawaOSPBSA,
                )%>%
                gather(marker,value,-c(id,ss_id,day,data_source)) %>%
                ggplot(aes(x=day,y=value,col=data_source,shape=data_source))+
                geom_hline(yintercept = c(-2,-5),lty=2)+
                geom_point(alpha=0.3)+
                geom_smooth(se=FALSE,method="loess")+
                facet_wrap(.~marker)+
                scale_x_continuous("Day (sqrt transformed)",
                                   trans="sqrt",
                                   breaks=c(7,365,1825,3650)
                )+
                scale_y_continuous("log10(1/RAU)")+
                cowplot::theme_cowplot()
        
}

plotPredictionInd <- function(generated, n=15, firstyear=TRUE){
        
        pred <- generated$predictions
        if(firstyear)  filt <- filter(pred,day_actual<=365)
        if(!firstyear)  filt <- filter(pred,day_actual>365)
        
        samp <- sample(unique(filt$ss_id),n,replace = TRUE)
        
        
        #turn the number of trees into probability mass function
        pred %>% filter(ss_id %in% samp)  %>%
                #get counts
                group_by(ss_id,week,week_actual,day_actual) %>%
                count() %>%
                #convert to proportions
                group_by(ss_id,week_actual) %>%
                mutate(prop=n/sum(n))%>%
                select(-n) %>% ungroup()%>%
                #make bar graphs
                ggplot(aes(x=week,y=prop))+
                geom_col()+
                facet_wrap(.~ss_id)+
                geom_vline(aes(xintercept=week_actual),col="red",lty=2)+
                geom_vline(xintercept=52,col="blue",lty=2)+
                cowplot::theme_cowplot()
        
}

plotPredictionAll <- function(generated,wks=c(1,15,30,50,53)){
        
        #new prediction plot
        generated$predictions %>% 
                filter(week_actual %in% wks) %>%
                group_by(week_actual,week) %>%
                count() %>%
                group_by(week_actual) %>%
                mutate(prop=n/sum(n)) %>%
                ggplot(aes(x=week,y=prop,fill=as.character(week_actual)))+
                geom_col()+
                scale_fill_viridis_d("Week of Infection")+
                scale_x_continuous("Predicted Week",breaks=c(0,16,32,52))+
                ylab("Proportion of Trees")+
                facet_wrap(.~week_actual,nrow=1)+
                # geom_hline(yintercept=c(1-generated$arguments$true_inc_rate,
                #                         generated$arguments$true_inc_rate/(max(generated$predictions$week)-1)),
                #            lty=2)+
                theme_cowplot()+
                theme(legend.position = "bottom")
        
}




plotReconstruction <- function(generated,fit){
        
        #get the plotting dataframe together
        plot_df <- generated$incidence %>% group_by(weeks_before_survey) %>%
                summarize(infections=sum(infections)) %>%
                filter(weeks_before_survey<=52) %>%
                mutate(true_incidence=infections/fit$input$N) %>% 
                left_join(fit$lambdas)
        
        #plot the image
        plot_df %>%     ggplot(aes(x=weeks_before_survey))+
                geom_point(aes(y=true_incidence))+
                geom_line(aes(y=`X50.`))+
                geom_ribbon(aes(ymin=`X2.5.`,ymax=`X97.5.`),col=NA,alpha=0.2)+
                theme_bw()+
                ylab("Incidence of Infection")+
                # scale_y_continuous("Incidence of Infection",
                #                    limits = c(0,max(c(fit$lambdas$`X97.5.`[1:52],
                #                plot_df$true_incidence,
                #                fit$input$p_prior[1]))))+
                scale_x_reverse("Weeks before serosurvey")+
                geom_line(inherit.aes=FALSE,
                          aes(x=1:52,y = fit$input$p_prior[1:52]),col="red",lty=2)
        
}





### extras


doSerosurveyDF <- function(df){
        
        tmp <- doSerosurvey(
                n_pop=df$n_pop, # size of the serosurvey
                true_inc_rate = df$true_inc_rate, #annual incidence rate
                detect_rate =df$detect_rate, #detection rate
                epidemic_type = df$epidemic_type, #specify flat or sin epidemic
                window_prop = df$window_prop#proportion of samples within the window to train RF on
        ) 
        
        return(tmp)
}

calcPriorVec <- function(prior_annual_inc, prior_epidemic_type){
        
        if(prior_epidemic_type=="sin"){
                prior_vec <- sin(2*pi*(1:52-90/7)/(52))+1
                prior_vec <- prior_annual_inc*prior_vec/sum(prior_vec)
                prior_vec <- c(prior_vec,1-sum(prior_vec))
                
        }
        
        if(prior_epidemic_type=="flat"){
                prior_vec <- rep(1,52)
                prior_vec <- prior_annual_inc*prior_vec/sum(prior_vec)
                prior_vec <- c(prior_vec,1-sum(prior_vec))
                
        }
        
        return(prior_vec)
}



#fit random forest model and apply to test set
randomForrest_fit_predict2 <- function(train,test){
        
        #set up model formula
        # markers<- c("Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA")
        # markers <-c("Luminex_IgG_CtxB","Luminex_IgG_InabaOSPBSA","Luminex_IgG_OgawaOSPBSA",
        #           "Luminex_IgA_CtxB","Luminex_IgA_InabaOSPBSA","Luminex_IgA_OgawaOSPBSA",
        #           "Luminex_IgM_OgawaOSPBSA","Luminex_IgM_InabaOSPBSA")
        markers <-c("Luminex_IgM_OgawaOSPBSA",
                    "Luminex_IgG_CtxB",
                    "Luminex_IgG_OgawaOSPBSA",
                    "Under10"
                    )
        
        timing_formula <- paste0("log_day~",paste(markers,collapse=" + ")) %>%
                as.formula
        
        #fit a randomforest model
        fit_rf <- train %>% mutate(log_day=log(day)) %>%
                randomForest(             formula=timing_formula,
                                          data=.,
                                          ntree=2000
                )
        
        
        #get predictions applied to serosurvey data
        full_preds_rf <- predict(fit_rf,newdata=test,predict.all=TRUE)
        
        # day_cuts <- c(7*(0:52)+1,Inf) #break days into weeks
        
        full_preds_rf_df<- data.frame(exp(full_preds_rf$individual),
                                      day_actual=test$day,
                                      ss_id = test$ss_id
        ) %>%
                gather(tree,prediction,-c(day_actual,ss_id)) %>%
                
                #new method without cuts
                #predicted weeks
                mutate(week=ceiling(prediction/7)) %>%
                mutate(week=ifelse(week<53,week,53)) %>%
                #true week
                mutate(week_actual=ceiling(day_actual/7)) %>%
                mutate(week_actual=ifelse(week_actual<53,week_actual,53)) 
        
        # #old method with cuts
        # #predicted weeks
        # mutate(week=cut(prediction,day_cuts,right=FALSE,
        # )) %>%
        # mutate(week=as.numeric(week)) %>%
        # #true week
        # mutate(week_actual=cut(day_actual,day_cuts,right=FALSE,
        # )) %>%
        # mutate(week_actual=as.numeric(week_actual)) %>%
        
        
        
        out <- list (fit =fit_rf,
                     predictions=full_preds_rf_df
        )
        
        return(out)
}



