
library(tidyverse)

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





mean(vax_run$sim1$seropos)
(mean(vax_run$sim1$seropos) + spec_estim-1)/(sens_estim + spec_estim - 1)
mean(vax_run$sim1$R) 





