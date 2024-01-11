
library(tidyverse)
library(tidybayes)

simulate_serosurvey_novax <- function(n_sim, n_survey=1000, sens_df,spec_df,incidence){
      
      out_list <- list()
      
      for(i in 1:n_sim){
            
            #lets do a serosurvey in a place with no vaccination but 5% incidence
            true_data <- data.frame(
                  sero_id=paste0("A",1:n_survey),
                  R=rbinom(n_survey, 1,incidence),
                  days_ago = ceiling(runif(n_survey,0,200))
            ) %>%
                  mutate(days_ago=ifelse(R==0,NA,days_ago))
            
            TP_df <- filter(true_data,R==1) %>%
                  group_by(sero_id)%>%
                  mutate(pos_ID=sample(unique(sens_df$pos_ID),1)) %>%
                  left_join(sens_df, by =c("days_ago", "pos_ID")) %>%
                  mutate(seropos=rbinom(1,1,ind_sens)) #%>%
            # select(sero_id,R,days_ago,pos_ID,seropos)
            
            TN_df <- filter(true_data,R==0) %>%
                  group_by(sero_id)%>%
                  mutate(neg_ID=sample(unique(spec_df$neg_ID),1)) %>%
                  left_join(spec_df,by ="neg_ID") %>%
                  mutate(seropos=1-rbinom(1,1,ind_spec))#%>%
            # select(sero_id,R,days_ago,neg_ID,seropos)
            
            
            sero_data <- bind_rows(TP_df,TN_df) %>%
                  arrange(sero_id) %>%
                  mutate(
                         incidence=incidence)
            
            
            out_list[[paste0("sim",i)]] <-sero_data
      }
      
      
      return(out_list)
}


simulate_serosurvey_vax <- function(n_sim, n_survey=1000,
                                    sens_df,spec_df,vax_df,
                                    coverage, incidence
){
      
      out_list <- list()
      
      for(i in 1:n_sim){
            
            #lets do a serosurvey in a place with no vaccination but 5% incidence
            true_data <- data.frame(
                  sero_id=paste0("A",1:n_survey),
                  R=rbinom(n_survey, 1,incidence),
                  V=rbinom(n_survey, 1,coverage),
                  days_ago = ceiling(runif(n_survey,0,200))
            ) %>%
                  mutate(days_ago=ifelse(R==0,NA,days_ago))
            
            TP_df <- filter(true_data,R==1) %>%
                  group_by(sero_id)%>%
                  mutate(pos_ID=sample(unique(sens_df$pos_ID),1)) %>%
                  left_join(sens_df, by =c("days_ago", "pos_ID")) %>%
                  mutate(seropos=rbinom(1,1,ind_sens)) %>%
                  select(sero_id,R,V,days_ago,pos_ID,seropos)
            
            TN_df_novax <- filter(true_data,R==0,V==0) %>%
                  group_by(sero_id)%>%
                  mutate(neg_ID=sample(unique(spec_df$neg_ID),1)) %>%
                  left_join(spec_df,by ="neg_ID") %>%
                  mutate(seropos=1-rbinom(1,1,ind_spec))%>%
                  select(sero_id,R,V,days_ago,neg_ID,seropos)
            
            TN_df_vax <- filter(true_data,R==0,V==1) %>%
                  group_by(sero_id)%>%
                  mutate(vax_ID=sample(unique(vax_df$vax_ID),1)) %>%
                  left_join(vax_df,by="vax_ID") %>%
                  mutate(seropos=1-rbinom(1,1,ind_spec))%>%
                  select(sero_id,R,V,days_ago,vax_ID,seropos)
            
            sero_data <- bind_rows(TP_df,TN_df_novax,TN_df_vax) %>%
                  arrange(sero_id) %>%
                  mutate(coverage=coverage,
                         incidence=incidence)
            
            out_list[[paste0("sim",i)]] <-sero_data
      }
      
      
      return(out_list)
}


#get time varying sensitivity
sens_obj <- data.frame(
      pos_ID = paste0("Z",rep(1:50,200)),
      diff=rep(rnorm(50,sd=1),200)
) %>% arrange(pos_ID) %>%
      mutate(days_ago=rep(1:200,50),
             avg_sens = rep(seq(log(0.75/(1-0.75)),log(0.25/(1-0.25)),length.out=200),50)
      ) %>%
      # mutate(ind_sens = log(avg_sens/(1-avg_sens))) %>%
      mutate(ind_sens = avg_sens+diff) %>%
      mutate(ind_sens= 1/(1+exp(-ind_sens))) %>%
      mutate(avg_sens= 1/(1+exp(-avg_sens)))

sens_estim <- sens_obj %>%
      group_by(days_ago)%>%
      summarize(sens=mean(ind_sens))%>%
      pull(sens) %>%
      mean()

#get specificity (not vaxxed)
spec_obj <- data.frame(
      neg_ID= paste0("Y",1:50),
      avg_spec = 0.95,
      diff=rnorm(50,sd=1)
)  %>%
      mutate(ind_spec = log(avg_spec/(1-avg_spec))) %>%
      mutate(ind_spec = ind_spec+diff) %>%
      mutate(ind_spec= 1/(1+exp(-ind_spec)))
spec_estim <- mean(spec_obj$ind_spec)


#get specificity (vaxxed)
vax_obj <- data.frame(
      vax_ID= paste0("X",1:50),
      avg_spec = 0.75,
      diff=rnorm(50,sd=1)
)  %>%
      mutate(ind_spec = log(avg_spec/(1-avg_spec))) %>%
      mutate(ind_spec = ind_spec+diff) %>%
      mutate(ind_spec= 1/(1+exp(-ind_spec)))



run_1 <- simulate_serosurvey_novax(n_sim = 10,
                                   n_survey = 1000,
                                   sens_df=sens_obj,spec_df=spec_obj,
                                   incidence=0.01)

run_5 <- simulate_serosurvey_novax(n_sim = 10,
                      n_survey = 1000,
                      sens_df=sens_obj,spec_df=spec_obj,
                      incidence=0.05)

run_10 <- simulate_serosurvey_novax(n_sim = 10,
                                   n_survey = 1000,
                                   sens_df=sens_obj,spec_df=spec_obj,
                                   incidence=0.10)
run_20 <- simulate_serosurvey_novax(n_sim = 10,
                                   n_survey = 1000,
                                   sens_df=sens_obj,spec_df=spec_obj,
                                   incidence=0.20)


bind_rows(
      bind_rows(run_1,.id="simulation"),
      bind_rows(run_5,.id="simulation"),
      bind_rows(run_10,.id="simulation"),
      bind_rows(run_20,.id="simulation")
      ) %>%
      group_by(simulation,incidence) %>%
      summarize(
            crude_seropos=mean(seropos),
            true_inc=mean(R)
      ) %>%
      mutate(adjusted_seropos=(crude_seropos + spec_estim-1)/(sens_estim + spec_estim - 1)) %>%
      mutate(adjusted_seropos=ifelse(adjusted_seropos<0,0,adjusted_seropos)) %>%
      gather(type,percent,-c(simulation,incidence)) %>%
      ggplot(aes(y=percent,x=type,col=type))+
            ggbeeswarm::geom_beeswarm()+
            geom_hline(aes(yintercept=incidence),lty=2)+
      facet_wrap(.~incidence,nrow=1)


run_1 <- simulate_serosurvey_novax(n_sim = 10,
                                   n_survey = 1000,
                                   sens_df=sens_obj,spec_df=spec_obj,
                                   incidence=0.01)

run_5 <- simulate_serosurvey_novax(n_sim = 10,
                      n_survey = 1000,
                      sens_df=sens_obj,spec_df=spec_obj,
                      incidence=0.05)

run_10 <- simulate_serosurvey_novax(n_sim = 10,
                                   n_survey = 1000,
                                   sens_df=sens_obj,spec_df=spec_obj,
                                   incidence=0.10)
run_20 <- simulate_serosurvey_novax(n_sim = 10,
                                   n_survey = 1000,
                                   sens_df=sens_obj,spec_df=spec_obj,
                                   incidence=0.20)


bind_rows(
      bind_rows(run_1,.id="simulation"),
      bind_rows(run_5,.id="simulation"),
      bind_rows(run_10,.id="simulation"),
      bind_rows(run_20,.id="simulation")
      ) %>%
      group_by(simulation,incidence) %>%
      summarize(
            crude_seropos=mean(seropos),
            true_inc=mean(R)
      ) %>%
      mutate(adjusted_seropos=(crude_seropos + spec_estim-1)/(sens_estim + spec_estim - 1)) %>%
      mutate(adjusted_seropos=ifelse(adjusted_seropos<0,0,adjusted_seropos)) %>%
      gather(type,percent,-c(simulation,incidence)) %>%
      ggplot(aes(y=percent,x=type,col=type))+
            ggbeeswarm::geom_beeswarm()+
            geom_hline(aes(yintercept=incidence),lty=2)+
      facet_wrap(.~incidence,nrow=1)


#vax simulations

vax_run_25 <- simulate_serosurvey_vax(n_sim = 10,
                                 n_survey = 1000,
                                 sens_df=sens_obj,
                                 spec_df=spec_obj,
                                 vax_df = vax_obj,
                                 coverage=0.25,
                                 incidence=0.05
                               )

vax_run_50 <- simulate_serosurvey_vax(n_sim = 10,
                                      n_survey = 1000,
                                      sens_df=sens_obj,
                                      spec_df=spec_obj,
                                      vax_df = vax_obj,
                                      coverage=0.5,
                                      incidence=0.05
)

vax_run_75 <- simulate_serosurvey_vax(n_sim = 10,
                                      n_survey = 1000,
                                      sens_df=sens_obj,
                                      spec_df=spec_obj,
                                      vax_df = vax_obj,
                                      coverage=0.75,
                                      incidence=0.05
)



bind_rows(
      bind_rows(vax_run_25,.id="simulation"),
      bind_rows(vax_run_50,.id="simulation"),
      bind_rows(vax_run_75,.id="simulation")
)%>%
      group_by(simulation,incidence,coverage) %>%
      summarize(
            crude_seropos=mean(seropos),
            true_inc=mean(R)
      ) %>%
      mutate(adjusted_seropos=(crude_seropos + spec_estim-1)/(sens_estim + spec_estim - 1)) %>%
      mutate(adjusted_seropos=ifelse(adjusted_seropos<0,0,adjusted_seropos)) %>%
      gather(type,percent,-c(simulation,incidence,coverage)) %>%
      ggplot(aes(y=percent,x=type,col=type))+
      ggbeeswarm::geom_beeswarm()+
      geom_hline(aes(yintercept=incidence),lty=2)+
      facet_wrap(.~coverage,nrow=1)


##### get new sens_df objects and spec_obj


loocv_spec_list <- read_rds(
          paste0("data/generated_data/analysis_objects/loocv/","loocv_spec_list.rds")
)

new_spec_obj<- loocv_spec_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$`Outside Window`$fit%>% 
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

new_spec_obj %>%
        ggplot(aes(x=ind_spec))+
        geom_histogram()+
        xlab("specificity")

new_spec_estim <- mean(new_spec_obj$ind_spec)


new_spec_obj$ind_spec %>% hist()


loocv_sens_list <- read_rds("data/generated_data/analysis_objects/loocv/loocv_tvfit_list_SMICPIC.rds")

sens_fitdata <-loocv_sens_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$Case$data
sens_spread <- loocv_sens_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$Case$fit %>% 
        spread_draws(`(Intercept)`,
                     b[term,group],
                     day1,day2,day3
        )

new_sens_obj <- data.frame()

for(t in 1:200){
        
        time <- log(t) -mean(log(sens_fitdata$new_day))
        
        tmp_df <- sens_spread %>%
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
        
        new_sens_obj <- bind_rows(new_sens_obj,tmp_df)
}

new_sens_obj %>%
        ggplot(aes(x=days_ago,y=ind_sens,group=pos_ID))+
                geom_line() +
        xlab("Days since infection")+
        ylab("sensitivity")

new_sens_estim <- new_sens_obj %>%
        group_by(days_ago)%>%
        summarize(sens=mean(ind_sens))%>%
        pull(sens) %>%
        mean()




run_1 <- simulate_serosurvey_novax(n_sim = 100,
                                   n_survey = 1000,
                                   sens_df=new_sens_obj,spec_df=new_spec_obj,
                                   incidence=0.01)

run_5 <- simulate_serosurvey_novax(n_sim = 100,
                                   n_survey = 1000,
                                   sens_df=new_sens_obj,spec_df=new_spec_obj,
                                   incidence=0.05)

run_10 <- simulate_serosurvey_novax(n_sim = 100,
                                    n_survey = 1000,
                                    sens_df=new_sens_obj,spec_df=new_spec_obj,
                                    incidence=0.10)
run_20 <- simulate_serosurvey_novax(n_sim = 100,
                                    n_survey = 1000,
                                    sens_df=new_sens_obj,spec_df=new_spec_obj,
                                    incidence=0.20)


bind_rows(
        bind_rows(run_1,.id="simulation"),
        bind_rows(run_5,.id="simulation"),
        bind_rows(run_10,.id="simulation"),
        bind_rows(run_20,.id="simulation")
) %>%
        group_by(simulation,incidence) %>%
        summarize(
                crude_seropos=mean(seropos),
                true_inc=mean(R)
        ) %>%
        mutate(adjusted_seropos=(crude_seropos + new_spec_estim-1)/(new_sens_estim + new_spec_estim - 1)) %>%
        mutate(adjusted_seropos=ifelse(adjusted_seropos<0,0,adjusted_seropos)) %>%
        gather(type,percent,-c(simulation,incidence)) %>%
        ggplot(aes(y=percent,x=type,col=type))+
        ggbeeswarm::geom_beeswarm()+
        geom_hline(aes(yintercept=incidence),lty=2)+
        facet_wrap(.~incidence,nrow=1)+
        theme(axis.text.x = element_text(angle=45,hjust=1))


#### new vax_obj
tvfpr_fit <- read_rds("data/generated_data/analysis_objects/misclassification/tvfpr_fit.rds")

vax_fitdata <-tvfpr_fit$`200`$`All Vaccinees`$data
vax_spread <- tvfpr_fit$`200`$`All Vaccinees`$fit%>% 
        spread_draws(`(Intercept)`,
                     b[term,group],
                     day1,day2,day3
        )


new_vax_tv_obj <- data.frame()

for(t in 1:365){
        
        time <- log(t) -mean(log(vax_fitdata$new_day))
        
        tmp_df <- vax_spread %>%
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
        
        new_vax_tv_obj <- bind_rows(new_vax_tv_obj,tmp_df)
}


new_vax_tv_obj %>%
        ggplot(aes(x=days_ago,y=ind_spec,group=vax_ID))+
        geom_line()+
        ylab("Specificity")+
        xlab("Days since vaccination")


new_vax_21 <- new_vax_tv_obj %>%
        filter(days_ago==21)%>%
        select(-days_ago)

new_vax_120 <- new_vax_tv_obj %>%
        filter(days_ago==120)%>%
        select(-days_ago)

vax21_run_25 <- simulate_serosurvey_vax(n_sim = 100,
                                      n_survey = 1000,
                                      sens_df=new_sens_obj,
                                      spec_df=new_spec_obj,
                                      vax_df = new_vax_21,
                                      coverage=0.25,
                                      incidence=0.05
) 

vax21_run_50 <- simulate_serosurvey_vax(n_sim = 10,
                                      n_survey = 1000,
                                      sens_df=new_sens_obj,
                                      spec_df=new_spec_obj,
                                      vax_df = new_vax_21,
                                      coverage=0.5,
                                      incidence=0.05
)

vax21_run_75 <- simulate_serosurvey_vax(n_sim = 10,
                                      n_survey = 1000,
                                      sens_df=new_sens_obj,
                                      spec_df=new_spec_obj,
                                      vax_df = new_vax_21,
                                      coverage=0.75,
                                      incidence=0.05
)

vax120_run_25 <- simulate_serosurvey_vax(n_sim = 100,
                                        n_survey = 1000,
                                        sens_df=new_sens_obj,
                                        spec_df=new_spec_obj,
                                        vax_df = new_vax_120,
                                        coverage=0.25,
                                        incidence=0.05
) 

vax120_run_50 <- simulate_serosurvey_vax(n_sim = 100,
                                        n_survey = 1000,
                                        sens_df=new_sens_obj,
                                        spec_df=new_spec_obj,
                                        vax_df = new_vax_120,
                                        coverage=0.5,
                                        incidence=0.05
)

vax120_run_75 <- simulate_serosurvey_vax(n_sim = 100,
                                        n_survey = 1000,
                                        sens_df=new_sens_obj,
                                        spec_df=new_spec_obj,
                                        vax_df = new_vax_120,
                                        coverage=0.75,
                                        incidence=0.05
)


bind_rows(
        bind_rows(vax21_run_25,.id="simulation") %>% mutate(campaign_day="21 days prior"),
        bind_rows(vax21_run_50,.id="simulation")%>% mutate(campaign_day="21 days prior"),
        bind_rows(vax21_run_75,.id="simulation")%>% mutate(campaign_day="21 days prior"),
        bind_rows(vax120_run_25,.id="simulation")%>% mutate(campaign_day="120 days prior"),
        bind_rows(vax120_run_50,.id="simulation")%>% mutate(campaign_day="120 days prior"),
        bind_rows(vax120_run_75,.id="simulation")%>% mutate(campaign_day="120 days prior")
)%>%
        group_by(simulation,incidence,coverage,campaign_day) %>%
        summarize(
                crude_seropos=mean(seropos),
                true_inc=mean(R)
        ) %>%
        mutate(adjusted_seropos=(crude_seropos + new_spec_estim-1)/(new_sens_estim + new_spec_estim - 1)) %>%
        mutate(adjusted_seropos=ifelse(adjusted_seropos<0,0,adjusted_seropos)) %>%
        mutate(adjusted_seropos=ifelse(adjusted_seropos>1,1,adjusted_seropos)) %>%
        gather(type,percent,-c(simulation,incidence,coverage,campaign_day)) %>%
        ggplot(aes(y=percent,x=type,col=type))+
        ggbeeswarm::geom_beeswarm()+
        geom_hline(aes(yintercept=incidence),lty=2)+
        facet_grid(campaign_day~coverage)+
        theme(axis.text.x = element_text(angle=45,hjust=1))+
        ylab("Seroprevalence")

##### New stan fitting of simulations

#specificity (non vax)

spec_dist <- loocv_spec_list$`Reduced IgG Panel`$smicpic_fitdata$`200-day`$`Outside Window`$fit%>% 
        spread_draws(`(Intercept)`,
                     b[term,group]) %>%
        mutate(logit_theta_j=`(Intercept)` + b
        ) %>%
        mutate(theta_j = 1/(1+exp(-logit_theta_j)))  %>%
        group_by(.draw) %>%
        summarize(spec=1-mean(theta_j))
        
spec_beta <- MASS::fitdistr(spec_dist$spec,
        densfun = "beta", list(shape1=1,shape2=1))

# sensitivity

sens_beta_draws <-data.frame()

for(t in 1:200){
        
        cat(t,"\n")
        
        time <- log(t) -mean(log(sens_fitdata$new_day))
        
        tmp_dist <- sens_spread %>%
                mutate(logit_theta_j=`(Intercept)` + b+
                               day1*time +
                               day2*time^2+
                               day3*time^3
                       
                ) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                group_by(.draw) %>%
                summarize(theta=mean(theta_j),n()) %>%
                mutate(days_ago=t)
        
        
        
        sens_beta_draws <-bind_rows(
                        sens_beta_draws,
                        tmp_dist
        )
}



sens_beta <- sens_beta_draws %>% group_by(.draw)%>%
                        summarize(sens=mean(theta)) %>%
                        pull(sens) %>%
                MASS::fitdistr(
                            densfun = "beta", list(shape1=1,shape2=1))

### specificity for vaccinees


vax_beta_draws <- data.frame()
vax_beta <- data.frame()

for(t in c(21,45,90,120,180)){
        cat(t,"\n")
        
        time <- log(t) -mean(log(vax_fitdata$new_day))
        
        tmp_df <- vax_spread %>%
                mutate(logit_theta_j=`(Intercept)` + b+
                               day1*time +
                               day2*time^2+
                               day3*time^3
                       
                ) %>%
                mutate(theta_j = 1/(1+exp(-logit_theta_j)))%>%
                group_by(.draw) %>%
                summarize(theta=1-mean(theta_j),n()) %>% #average across individuals here
                mutate(days_ago=t)
        
        vax_beta_draws <- bind_rows(vax_beta_draws,tmp_df)
        
        tmp_beta <- tmp_df$theta %>%
                MASS::fitdistr(
                        densfun = "beta", list(shape1=1,shape2=1))
        
        vax_beta <- bind_rows(vax_beta,
                              tmp_beta$estimate %>% as.matrix()%>%
                                      t() %>% data.frame() %>%
                                      mutate(days_ago=t)
                              )
        
        
}


vax_beta_draws %>%
        ggplot(aes(x=theta,col=factor(days_ago)))+
                geom_density()



#run it

cat_counts<- vax21_run_25$sim1 %>% 
        ungroup()%>%
        count(seropos,V) %>%
        mutate(category=glue::glue("S{seropos}Q{V}")) %>%
        select(n,category) %>%
        spread(category,n)

stan_obj <- list(
          S0Q0 = cat_counts$S0Q0+cat_counts$S0Q1,
          S1Q0 = cat_counts$S1Q0+cat_counts$S1Q1,
          S0Q1 = 0,
          S1Q1 = 0,
        
         alpha_p1 = sens_beta$estimate[1],
         alpha_p2 = sens_beta$estimate[2],
         beta_p1 = spec_beta$estimate[1],
         beta_p2 = spec_beta$estimate[2],
         epsilon_p1 = vax_beta %>%
                 filter(days_ago==21) %>%
                 pull(shape1),
         epsilon_p2 =vax_beta %>%
                 filter(days_ago==21) %>%
                 pull(shape2)
)

library(rstan)

second_model <- stan_model("code/stan/vax-model-simplified.stan")
fit <- fit_second_model <- sampling(second_model, data = stan_obj)

shinystan::launch_shinystan(fit)

