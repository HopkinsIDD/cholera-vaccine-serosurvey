
path <- "source/final_code/vax_strategy/generated_rds/simulations/"

#bring in the estimates
serologyonly_params <- read_rds("source/final_code/vax_strategy/generated_rds/sim-params/serologyonly/serologyonly_params.rds")
questionairre_params <- read_rds("source/final_code/vax_strategy/generated_rds/sim-params/questionairre_params.rds")

#questionairre params base
base_quest_params <- list(
        title= "IgG Reduced Panel",
        
        sens_novax_mu = questionairre_params$estimates[[1]]$`mean(sens)`,
        sens_vax_mu =  questionairre_params$estimates[[1]]$`mean(sens)`,
        spec_novax_mu =  questionairre_params$estimates[[2]]$spec_mu,
        
        sens_novax_sd = questionairre_params$estimates[[1]]$`sd(sens)`,
        sens_vax_sd =  questionairre_params$estimates[[1]]$`sd(sens)`,
        spec_novax_sd =  questionairre_params$estimates[[2]]$spec_sd
        
)

tmp_sero_params %>%
        lapply(FUN = function(x){
                
                tmp <- x
                tmp[[1]] <- x[[1]][order(names(x[[1]]))]
                tmp[[2]] <- x[[2]][order(names(x[[1]])),order(names(x[[1]]))]
                
                return(x)

        })


quest_stanmodel <- stan_model("source/final_code/vax_strategy/stan/vax-model-simplified2.stan")
sero_stanmodel <- stan_model("source/final_code/vax_strategy/stan/vax-model-multi-class.stan")

options(mc.cores=1)

sim_df <- data.frame()
for(d in c("day21","day44","day90")){

        #questionairre params ready
        tmp_quest_params  <- base_quest_params
        tmp_quest_specvax_df <- questionairre_params$estimates[[3]] %>%
                filter(new_day==str_remove(d,"day"))
        
        tmp_quest_params$spec_vax_mu <- tmp_quest_specvax_df$spec_mu
        tmp_quest_params$spec_vax_sd <-tmp_quest_specvax_df$spec_sd
        
        #get serology params ready
        tmp_sero_params <- serologyonly_params$estimates[[1]][[d]]%>%
                lapply(FUN = function(x){
                        
                        #order the names
                        #Neither
                        #RecentlyInfected
                        #RecentlyVaccinated
                        tmp_order  <- order(names(x[[1]]))
                        
                        tmp <-  list()
                        tmp[[1]] <- x[[1]][tmp_order]
                        tmp[[2]] <- x[[2]][tmp_order,tmp_order]
                        
                        return(tmp)
                        
                })

        
        
        for(inc in c(0.05)){
                for(cov in c(0.25,0.5,0.75)){

                        
                        cat("day: ",d," incidence: ",inc," ",
                            " coverage: ",cov,"\n")
                        
                        tmp_sim <-  simVaxSurvey_new(
                                unvaxxed_inc= inc,
                                cov= cov,
                                sero_params = tmp_sero_params,
                                quest_params = tmp_quest_params,
                                n_sim=10)
                        
                        count_data <- tmp_sim$ind_data %>%
                                group_by(incidence,coverage,simulation)  %>%
                                summarize(
                                          sample_inc=mean(infected_recently),
                                          S1Q1=sum(seropos*questionnaire),
                                          S0Q1=sum((1-seropos)*questionnaire),
                                          S1Q0=sum(seropos*(1-questionnaire)),
                                          S0Q0=sum((1-seropos)*(1-questionnaire)),
                                          Z1= sum(sero_class3=="Neither"),
                                          Z2= sum(sero_class3=="RecentlyInfected"),
                                          Z3= sum(sero_class3=="RecentlyVaccinated")
                                ) %>%
                                mutate(day=str_remove(d,"day"))
                        
                        #FIT SEROINCIDENCE MODELS
                        lambda_data1 <- lambda_data2<- lambda_data3 <- data.frame()
                        for(i in 1:nrow(count_data)){
                                
                                #strategy 1: ignore the problem
                                input_list1 <- list(S0Q0=count_data$S0Q0[i] + count_data$S0Q1[i],
                                                   S1Q0=count_data$S1Q0[i]+count_data$S1Q1[i],
                                                   S0Q1=0,
                                                   S1Q1=0,
                                                   
                                                   alpha_p1 = tmp_quest_params$sens_novax_mu,
                                                   alpha_p2 = tmp_quest_params$sens_novax_sd,
                                                   
                                                   beta_p1 = tmp_quest_params$spec_novax_mu,
                                                   beta_p2 = tmp_quest_params$spec_novax_sd,
                                                   
                                                   epsilon_p1 = tmp_quest_params$spec_vax_mu,
                                                   epsilon_p2 = tmp_quest_params$spec_vax_sd
                                                   
                                )
                                
                                tmp_fit1 <- sampling(quest_stanmodel,input_list1,refresh=0)
                                
                                tmp_lambda1 <- tmp_fit1 %>% spread_draws(lambda) %>%
                                        median_qi()%>%
                                        mutate(method="Ignore")
                                
                                lambda_data1 <- bind_rows(lambda_data1,tmp_lambda1) 
                                
                                #strategy 2: questionairre
                                input_list2 <- list(S0Q0=count_data$S0Q0[i],
                                                   S1Q0=count_data$S1Q0[i],
                                                   S0Q1=count_data$S0Q1[i],
                                                   S1Q1=count_data$S1Q1[i],
                                                   
                                                   alpha_p1 = tmp_quest_params$sens_novax_mu,
                                                   alpha_p2 = tmp_quest_params$sens_novax_sd,
                                                   
                                                   beta_p1 = tmp_quest_params$spec_novax_mu,
                                                   beta_p2 = tmp_quest_params$spec_novax_sd,
                                                   
                                                   epsilon_p1 = tmp_quest_params$spec_vax_mu,
                                                   epsilon_p2 = tmp_quest_params$spec_vax_sd
                                                   
                                )
                                
                                tmp_fit2 <- sampling(quest_stanmodel,input_list2,refresh=0)
                                
                                tmp_lambda2 <- tmp_fit2 %>% spread_draws(lambda,gamma) %>%
                                        median_qi() %>%
                                        mutate(method="Questionairre")
                                
                                lambda_data2 <- bind_rows(lambda_data2,tmp_lambda2)
                                
                                
                                #strategy 3: multiclassification
                                
                                input_list3 <- list(Z1=count_data$Z1[i],
                                                    Z2=count_data$Z2[i],
                                                    Z3=count_data$Z3[i],

                                                    alpha1_mu = tmp_sero_params$Neither[[1]],
                                                    alpha1_vcv = tmp_sero_params$Neither[[2]],
                                                    
                                                    alpha2_mu = tmp_sero_params$RecentlyInfected[[1]],
                                                    alpha2_vcv = tmp_sero_params$RecentlyInfected[[2]],
                                                    
                                                    alpha3_mu = tmp_sero_params$RecentlyVaccinated[[1]],
                                                    alpha3_vcv = tmp_sero_params$RecentlyVaccinated[[2]]
                                                    
                                )
                                
                                tmp_fit3 <- sampling(sero_stanmodel,input_list3,refresh=0)
                                
                                tmp_lambda3 <- tmp_fit3 %>% spread_draws(lambda,gamma) %>%
                                        median_qi() %>%
                                        mutate(method="Serological")
                                
                                lambda_data3 <- bind_rows(lambda_data3,tmp_lambda3)
                                
                        }
                        
                        #bring all the data together
                        out_data <- 
                                bind_rows(
                                        bind_cols(count_data,lambda_data1),
                                        bind_cols(count_data,lambda_data2),
                                        bind_cols(count_data,lambda_data3)
                                        
                                )

                        sim_df <- bind_rows(
                                sim_df,out_data
                        )
                }
        }
}


write_rds(sim_df,paste0(path,"new_sim_df.rds"))


sim_df %>% ggplot(aes(y=lambda,col=method,x=factor(as.numeric(day)))) +
        geom_boxplot()+
        ylim(0,0.5)+
        facet_wrap(.~coverage)+
        geom_hline(yintercept=0.05,lty=2)

sim_df %>% ggplot(aes(y=gamma,col=method,x=factor(as.numeric(day)))) +
        geom_boxplot()+
        facet_wrap(.~coverage)+
        ylim(0,1)




sim_df %>% 
        mutate(lambda.lower=ifelse(is.na(lambda.lower),.lower,lambda.lower)) %>%
        mutate(lambda.upper=ifelse(is.na(lambda.upper),.upper,lambda.upper)) %>%
        
        ggplot(aes(y=lambda,col=method,x=factor(as.numeric(day)),
                      group=paste0(day,method,simulation)
)) +
        geom_point(position=position_dodge(width = 0.5))+
        geom_linerange(aes(ymin=lambda.lower,ymax=lambda.upper),position=position_dodge(width = 0.5),
                       alpha=0.3
        )+
        facet_wrap(.~coverage)+
        theme_cowplot()+
        geom_hline(yintercept=0.05,lty=2)






sim_df %>% ggplot(aes(y=gamma,col=method,x=factor(as.numeric(day)),
                      group=paste0(day,method,simulation)
)) +
        geom_point(position=position_dodge(width = 0.5))+
        geom_linerange(aes(ymin=gamma.lower,ymax=gamma.upper),position=position_dodge(width = 0.5),
                       alpha=0.3
        )+
        facet_wrap(.~coverage)+
        ylim(0,1)+
        theme_cowplot()

serologyonly_params$estimates[[2]] %>%
        filter(day_cat!="day180") %>%
        ggplot(aes(x=category,y=prob,fill=day_cat))+
        geom_col(position=position_dodge())+
        facet_wrap(~truth)+
        theme(axis.text.x = element_text(angle=45,hjust=1))

# 
# sim_count <- sim_df %>% group_by(day,incidence,coverage,simulation)  %>%
#         summarize(S1Q1=sum(seropos*questionnaire),
#                   S0Q1=sum((1-seropos)*questionnaire),
#                   S1Q0=sum(seropos*(1-questionnaire)),
#                   S0Q0=sum((1-seropos)*(1-questionnaire)),
#                   Z1= sum(sero_class3=="Neither"),
#                   Z2= sum(sero_class3=="RecentlyInfected"),
#                   Z3= sum(sero_class3=="RecentlyVaccinated")
#         )
#         




#### old simulation stuff


# tvfpr_df <- read_rds(
#         paste0("source/final_code/vax_strategy/generated_rds/fpr-estimates/","tvfpr_df.rds")
# )
# 
# newspecdf <- tvfpr_df %>% filter(end_window==200) %>%
#         filter(seropos_type=="spec95_seropos") %>%
#         filter(variables=="Reduced Panel IgG") %>%
#         mutate(newtime=time-14) %>%
#         filter(newtime>0,newtime<=200) %>%
#         mutate(spec=1-median) %>%
#         mutate(spec=ifelse(spec>0.95,0.95,spec)) %>%
#         filter(newtime %in% c(1,7*(1:28),200)) %>%
#         select(newtime,spec)
# 
# 
# options(mc.cores = 1)
# final_df <- data.frame()
# 
# for(inc in c(0.01,0.05,0.1)){
#         for(cov in c(0.25,0.5,0.75)){
#                 
#                 final_df <- bind_rows(
#                         final_df,
#                         vaxSimSeroInc(inc,cov,newspecdf,sims=10)
#                 )
#                 
#                 
#         }
# }



# write_rds(final_df,paste0(path,"seroinc-estimates.rds"))


#--------------------------


# ##### Simulate the serosurveys -----------
# #look at different incidences
# #look at different coverage
# # look at different strategies
# # assume perfect vaccination status gathering
# 
# 
# reduced_param <- list(
#         title= "IgG Reduced Panel",
#         sens_novax = 0.5,
#         sens_vax = 0.5,
#         spec_vax = 0.7,
#         spec_novax = 0.9
#         
# )
# 
# ctxbonly_param <- list(
#         title= "CtxB IgG only",
#         sens_novax = 0.4,
#         sens_vax = 0.4,
#         spec_vax = 0.9,
#         spec_novax = 0.9
#         
# )
# 
# #this method of simulation will not work without perfect 
# # ascertainment of vaccination status
# mixed_param <- list(
#         title= "Mixed strategy",
#         sens_novax = 0.5,
#         sens_vax = 0.4,
#         spec_vax = 0.9,
#         spec_novax = 0.9
#         
# )
# 
# 
# param_list <- list(
#         `IgG Reduced Panel`=reduced_param,
#         `CtxB IgG only`=ctxbonly_param,
#         `Mixed strategy`=mixed_param
#         
# )
# 
# sim_df <- data.frame()
# sim_list <- list()
# 
# sims <-100
# 
# for(cov in c(0,0.25,0.5,0.75)){
#                 for(unvaxxed_inc in c(0.05,0.15)){
#                         
#                         cat("Coverage: ",cov," Unvaxxed Inc: ",unvaxxed_inc,"\n")
#                         
#                         reduced_sim <- simVaxSurvey(
#                                         n_sim=sims,
#                                         unvaxxed_inc = unvaxxed_inc,
#                                         cov = cov,
#                                         sens_quest=1,
#                                         spec_quest = 1,
#                                         effectiveness =0,
#                                         seropos_param=reduced_param
#                                      )
#                         
#                         ctxb_sim <-  simVaxSurvey(
#                                 n_sim=sims,
#                                 unvaxxed_inc = unvaxxed_inc,
#                                 cov = cov,
#                                 sens_quest=1,
#                                 spec_quest = 1,
#                                 effectiveness =0,
#                                 seropos_param=ctxbonly_param
#                         )
#                         
#                         mixed_sim <-  simVaxSurvey(
#                                 n_sim=sims,
#                                 unvaxxed_inc = unvaxxed_inc,
#                                 cov = cov,
#                                 sens_quest=1,
#                                 spec_quest = 1,
#                                 effectiveness =0,
#                                 seropos_param=mixed_param
#                         )
#                         
#                         tmp_df <- bind_rows(
#                                         reduced_sim$ind_data,
#                                         ctxb_sim$ind_data,
#                                         mixed_sim$ind_data
#                                         ) %>%
#                                 mutate(cov=cov,
#                                        unvaxxed_inc=unvaxxed_inc
#                                        )
#                         
#                         sim_df <- bind_rows(sim_df,tmp_df)
#                                               
#                         
#                         sim_list[[as.character(cov)]][[as.character(unvaxxed_inc)]][["IgG Reduced Panel"]] <- reduced_sim
#                         sim_list[[as.character(cov)]][[as.character(unvaxxed_inc)]][["CtxB IgG Only"]] <- ctxb_sim
#                         sim_list[[as.character(cov)]][[as.character(unvaxxed_inc)]][["Mixed Strategy"]] <- mixed_sim
# 
#                         
#                         
# }}
# 
# 
# write_rds(sim_df,paste0(path,"sim_df.rds"))
# write_rds(sim_list,paste0(path,"sim_list.rds"))
# 
# sim_df %>% group_by(cov,unvaxxed_inc,simulation,variables) %>%
#                 summarize(seropos=mean(seropos)) %>% 
#                 ggplot(aes(x=factor(variables),y=seropos,col=factor(unvaxxed_inc)))+
#                 geom_boxplot()+
#                 facet_grid(.~cov)+
#                 geom_hline(yintercept = c(0.05,0.15),lty=2)+
#                 scale_y_continuous("Seropositivity",limits=c(0,.4))+
#                 theme(axis.text.x = element_text(angle=45,hjust=1))
# 
# 
# 
# 
# 
# 
# 
# ##### fit the models -----------
# 
# sim_df <- read_rds(paste0(path,"sim_df.rds"))
# 
# vax_model <- stan_model("source/final_code/vax_strategy/stan/vax-model-simplified2.stan")
# 
# 
# newsim <- sim_df %>% group_by(cov,unvaxxed_inc,variables,simulation) %>%
#         summarize(S1Q1=sum(seropos*questionnaire),
#                   S0Q1=sum((1-seropos)*questionnaire),
#                   S1Q0=sum(seropos*(1-questionnaire)),
#                   S0Q0=sum((1-seropos)*(1-questionnaire)),
#                   study_seropos=mean(seropos),
#                   study_inc=mean(infected_recently)
#         ) %>%
#         mutate(n=S1Q1+S0Q1+S1Q0+S0Q0)
# 
# fitsummary_df <- data.frame()
# 
# for(coverage in c(0,0.25,0.5,0.75)){
#         for(unvaxxed_incidence in c(0.05,0.15)){
#                 for (var in unique(sim_df$variables)){
#                         
#                         tmp_df <- newsim %>% filter(cov==coverage) %>%
#                                 filter(unvaxxed_inc==unvaxxed_incidence) %>%
#                                 filter(variables==var)
#                         
#                         this_df <- data.frame()
#                         
#                         for(s in (1:100)){
#                                 
#                                 
#                                 s <- paste0("sim_",str_pad(s,width=5,pad="0",side="left"))
#                                 
#                                 cat("Coverage: ",coverage," Unvaxxed Inc: ",unvaxxed_incidence,"\n")
#                                 cat("Variables: ",var," Simulation: ",s,"\n")
#                                 
#                                 #estimate seroincidence model
#                                 tmp_fit <- tmp_df %>% filter(simulation==s) %>% 
#                                         estimateSeroInc_vaxpop_simple(model=vax_model,
#                                                                       params=param_list[[var]]
#                                         )
#                                 
#                                 #store summary statistics for key values
#                                 out_df <- tmp_fit %>%
#                                         spread_draws(lambda) %>% 
#                                         mutate(seroinc=lambda) %>%
#                                         median_qi() %>%
#                                         mutate(cov=coverage,
#                                                unvaxxed_inc=unvaxxed_incidence,
#                                                variables=var,
#                                                simulation=s)
#                                 
#                                 this_df <- bind_rows(this_df,out_df)
#                         }
#                         
#                         
#                         fitsummary_df <- bind_rows(fitsummary_df,this_df)
#                         
#                         write_rds(fitsummary_df,paste0(path,"running_fitsummary_df1.rds"))
#                         
#                         
#                 }}}
# 



#--------------------------
# for(coverage in c(0,0.25,0.5,0.75)){
#         for(unvaxxed_incidence in c(0.05,0.15)){
#                 for (var in unique(sim_df$variables)){
#                         
#         # coverage <- 0.5
#         # unvaxxed_incidence <- 0.08
#         # var <- "IgG Reduced Panel"
#                         
#                         tmp_df <- sim_df %>% filter(cov==coverage) %>%
#                                          filter(unvaxxed_inc==unvaxxed_incidence) %>%
#                                          filter(variables==var)
#                         
#                         this_df <- data.frame()
#                         
#                         for(s in (1:10)){
#                                 
#                                 
#                                 s <- paste0("sim_",str_pad(s,width=5,pad="0",side="left"))
#                                 
#                                 cat("Coverage: ",coverage," Unvaxxed Inc: ",unvaxxed_incidence,"\n")
#                                 cat("Variables: ",var," Simulation: ",s,"\n")
#                                 
#                                 #estimate seroincidence model
#                                 tmp_fit <- tmp_df %>% filter(simulation==s) %>% 
#                                                 estimateSeroInc_vaxpop(model=vax_model,
#                                                                        params=param_list[[var]]
#                                                                        )
#                                 
#                                 #store summary statistics for key values
#                                 out_df <- tmp_fit %>%
#                                         spread_draws(lambda) %>% 
#                                         mutate(seroinc=lambda) %>%
#                                         # spread_draws(lambda,gamma,nu) %>% 
#                                         # mutate(seroinc=gamma*lambda+(1-gamma)*lambda) %>%
#                                         median_qi() %>%
#                                         mutate(cov=coverage,
#                                                unvaxxed_inc=unvaxxed_incidence,
#                                                variables=var,
#                                                simulation=s)
#                                 
#                                 this_df <- bind_rows(this_df,out_df)
#                                 }
#                         
#                         
#                         fitsummary_df <- bind_rows(fitsummary_df,this_df)
#                         
#                         write_rds(fitsummary_df,paste0(path,"running_fitsummary_df.rds"))
# 
# 
#                         }}}

#--------------------------

# fitsummary_df<- read_rds(paste0(path,"running_fitsummary_df.rds"))
# 
# 
# fitsummary_df %>%
#         mutate(variables=factor(variables,levels=c("IgG Reduced Panel",
#                                                    "CtxB IgG only",
#                                                    "Mixed strategy"))) %>%
#         ggplot(aes(x=variables,y=seroinc,group=factor(simulation)))+
#         geom_point(position = position_dodge(0.5))+
#         geom_linerange(aes(ymin=seroinc.lower,ymax=seroinc.upper),position = position_dodge(0.5))+
#         facet_grid(cov~unvaxxed_inc)+
#         geom_hline(yintercept = c(0.15,0.05),lty=3)
# 
# fitsummary_df %>%
#         mutate(variables=factor(variables,levels=names(param_list))) %>%
#         mutate(lambda=unvaxxed_inc)%>%
#         ggplot(aes(x=factor(lambda),y=seroinc,group=factor(simulation),col=factor(lambda)))+
#         geom_point(position = position_dodge(0.5))+
#         geom_linerange(aes(ymin=seroinc.lower,ymax=seroinc.upper),position = position_dodge(0.5))+
#         facet_grid(cov~variables)+
#         xlab("True Seroincidence")+
#         scale_y_continuous("Estimated Seroincidence")+
#         geom_hline(yintercept = c(0.15,0.05),lty=3)+
#         theme_bw()
# 
# 
# 
# #get the unadjusted values of sero-incidence
# uandjusted_df<- sim_df %>% group_by(cov,unvaxxed_inc,simulation,variables) %>%
#         summarize(unadj_seropos=mean(seropos))
# 
# 
# hline_inc <- data.frame(unvaxxed_inc=c(0.15,0.05),inc=c(0.15,0.05)) 
# 
# 
# plotdf <- fitsummary_df %>% select(seroinc,cov,unvaxxed_inc,simulation,variables) %>%
#         left_join(uandjusted_df,by = c("cov", "unvaxxed_inc", "simulation", "variables")) %>%
#         gather(Incidence,`Estimated Seroincidence`,-c(cov,unvaxxed_inc,simulation,variables))%>%
#         mutate(Incidence=case_when(Incidence == "seroinc" ~"Adjusted",
#                                    Incidence == "unadj_seropos" ~"Crude"
#         )) %>%
#         mutate(Incidence= factor(Incidence,levels=c("Crude","Adjusted"))) %>%
#         mutate(`Testing Strategy`=factor(variables,levels=names(param_list))) %>%
#         # mutate(true_inc=paste0(round(0.75*unvaxxed_inc*100,1)," cases per 100")) %>%
#         # mutate(true_inc=factor(true_inc,levels=c("6 cases per 100","15 cases per 100"))) %>%
#         mutate(coverage=paste0(cov*100,"% Coverage")) 
# 
# 
# plotdf %>%
#         ggplot(aes(x=`Testing Strategy`,y=`Estimated Seroincidence`,col=Incidence))+
#         geom_boxplot()+
#         facet_grid(unvaxxed_inc~coverage)+
#         geom_hline(data=hline_inc, aes(yintercept=inc),lty=2)+
#         theme_bw()+
#         theme(axis.text.x = element_text(angle=45,hjust=1))










# 
# uandjusted_df<- sim_df %>% group_by(cov,unvaxxed_inc,simulation,variables) %>%
#         summarize(unadj_seropos=mean(seropos))
# 
# 
# hline_inc <- data.frame(true_inc=c("15 cases per 100","5 cases per 100"),inc=c(0.15,0.05)) %>%
#         mutate(true_inc=factor(true_inc,levels=c("5 cases per 100","15 cases per 100"))) 
#         
# 
# plotdf <- fitsummary_df %>% select(seroinc,cov,unvaxxed_inc,simulation,variables) %>%
#         left_join(uandjusted_df) %>%
#         gather(Incidence,`Estimated Seroincidence`,-c(cov,unvaxxed_inc,simulation,variables))%>%
#         mutate(Incidence=case_when(Incidence == "seroinc" ~"Adjusted",
#                                    Incidence == "unadj_seropos" ~"Crude"
#                                    )) %>%
#         mutate(Incidence= factor(Incidence,levels=c("Crude","Adjusted"))) %>%
#         mutate(`Testing Strategy`=factor(variables,levels=names(param_list))) %>%
#         mutate(true_inc=paste0(round(0.75*unvaxxed_inc*100,1)," cases per 100")) %>%
#         mutate(true_inc=factor(true_inc,levels=c("6 cases per 100","15 cases per 100"))) %>%
#         mutate(coverage=paste0(cov*100,"% Coverage")) 
#         
#         
# plotdf %>%
#         ggplot(aes(x=`Testing Strategy`,y=`Estimated Seroincidence`,col=Incidence))+
#         geom_boxplot()+
#         facet_grid(true_inc~coverage)+
#         geom_hline(data=hline_inc, aes(yintercept=inc),lty=2)+
#         theme_bw()+
#         theme(axis.text.x = element_text(angle=45,hjust=1))
#         


# sim_df %>% filter(cov==0.5,
#                   unvaxxed_inc==0.2,
#                   simulation==1) %>%
#                 group_by(variables,questionnaire) %>%
#                 summarize(n=n(),seropos=sum(seropos)) %>%
#                 ungroup() %>%
#                 mutate(sens=c(0.4,0.4,0.5,0.5,0.5,0.4)) %>%
#                 mutate(spec=c(0.9,0.9,0.9,0.7,0.9,0.9)) %>%
#                 group_by(variables,questionnaire) %>%
#                 mutate(inc=adjSeropos(sens=sens,spec=spec,n_seropos=seropos,n_total=n)) %>%
#                 group_by(variables) %>%
#                 summarize(sum(inc)/sum(n))
        
                

#OLDCODE

# tvfpr_df %>% filter(end_window==200) %>%
#         filter(age_group=="Adult")%>%
#         filter(variables=="Reduced Panel IgG") %>%
#         filter(seropos_type %in% c("spec90_seropos","spec95_seropos")) %>%
#         group_by(seropos_type) %>%
#         filter(max(median)==median) 
# 
# #### Function
# #get the recent infection status        
# # getInfectionStatus <-function(x, incidence, ve){
# #         
# #         if(x==1) out <- rbinom(1,1,incidence*ve)
# #         if(x==0) out <- rbinom(1,1,incidence)
# #         
# #         return(out)
# #         
# # }
# 
# 
# #get the questionairre status
# getQuestionairre<- function(x, sens, spec){
#         
#         if(x==1) out <- rbinom(1,1,sens)
#         if(x==0) out <- rbinom(1,1,1-spec)
#         
#         return(out)
# }
# 
# getSeroPos<- function(v,spec,sens){
#         
#         if(v==1) out <- rbinom(1,1,sens)
#         if(v==0) out <- rbinom(1,1,1-spec)
#         
#         return(out)
# }
# 
# 
# 
# adjSeropos <- function(sens,spec,n_seropos,n_total){
#         max((n_seropos+n_total*(spec-1))/(sens+spec-1),0)
# }
# 
# n_pop <- 2000
# 
# #hold these constant (or maybe these parameters should have some variability)
# effectiveness =0.5 #vaccine effectiveness
# spec_quest = 0.90 # specificity of questionnaire
# 
# # sens_novax <- 0.7 # sensitivity of RF among unvaxxed
# # sens_vax <- sens_novax # sensitivity of RF among unvaxxed
# # spec_vax <- 0.75 # specificity of RF among vaxxed
# # spec_novax <- 0.95 # specificity of RF among unvaxxed
# 
# calc_df <- data.frame()
# 
# for(cov in c(0.25,0.5,0.75)){
#         for(sens_quest in c(0,0.5,0.9,1)){
#                 for(true_inc in c(0.05,0.15)){
#                         
#                         
#                         cat(cov," ", sens_quest," " ,true_inc,"\n")
#                         
#                         #incidence among unvaxxed
#                         unvaxxed_inc<- true_inc/(cov*effectiveness+(1-cov)) 
#                         
#                         #first see how many people are vaxxed
#                         unvaxxed  <- rep(0,n_pop-round(n_pop*cov))
#                         vaxxed  <- rep(1,round(n_pop*cov)) 
#                         
#                         #get infection status
#                         unvaxxed_infected <-rep(1,round(length(unvaxxed)*unvaxxed_inc))
#                         unvaxxed_uninfected <- rep(0,length(unvaxxed)-length(unvaxxed_infected))
#                         vaxxed_infected <- rep(1,round(length(vaxxed)*unvaxxed_inc*effectiveness))
#                         vaxxed_uninfected <- rep(0,length(vaxxed)-length(vaxxed_infected))
#                         
#                         sim_df <- data.frame()
#                         #simulate questionairre status
#                         for(i in 1:100){
#                                 
#                                 
#                                 tmp_df  <- data.frame(
#                                                 vax_status= c(unvaxxed,vaxxed),
#                                                 infected_recently =c(
#                                                         unvaxxed_infected,
#                                                         unvaxxed_uninfected,
#                                                         vaxxed_infected,
#                                                         vaxxed_uninfected )
#                                         ) %>%   
#                                         #get the simulation number
#                                         mutate(simulation=i) %>%
#                                         
#                                         #get questionnaire status
#                                         mutate(questionnaire=sapply(vax_status, 
#                                                                     FUN=getQuestionairre,
#                                                                     sens = sens_quest,
#                                                                     spec = spec_quest
#                                         ))
#                                 
#                                 
#                                 #get seropositivity for different tests
#                                 
#                                 # FPR for 200 reduced panel
#                                 spec_vax <- rnorm(1,mean=0.7,sd=0.034)
#                                 
#                                 # specificity
#                                 spec_novax <- rnorm(1,mean=0.9,sd=0.020)
#                                 
#                                 #sensitivity
#                                 sens_novax <- rnorm(1,mean=0.5,sd=0.052)
#                                 
#                                 #for now just say senstivity drops a little with the new test
#                                 # but specificity increases among non-vaccinated
#                                 sens_novax2 <- rnorm(1,mean=0.4,sd=0.052)
#                                 
#                                 
#                                 
#                                 sero_tmp <- bind_rows(
#                                         tmp_df %>%
#                                                 mutate(seropos= ifelse(vax_status ==0,
#                                                         sapply(infected_recently,
#                                                                FUN=getSeroPos,
#                                                                sens=sens_novax,
#                                                                spec= spec_novax),
#                                                           sapply(infected_recently,
#                                                                  FUN=getSeroPos,
#                                                                  sens=sens_novax,
#                                                                  spec= spec_vax))) %>%
#                                                  mutate(variables = "Reduced Panel"),
#                                         tmp_df %>%
#                                                 mutate(seropos= ifelse(vax_status ==0,
#                                                                        sapply(infected_recently,
#                                                                               FUN=getSeroPos,
#                                                                               sens=sens_novax2,
#                                                                               spec= spec_novax),
#                                                                        sapply(infected_recently,
#                                                                               FUN=getSeroPos,
#                                                                               sens=sens_novax2,
#                                                                               spec= spec_novax))) %>%
#                                                 mutate(variables = "CtxB Panel")
#                                         
#                                 )
#                                 
#                                 sim_df <- bind_rows(sim_df,sero_tmp)
#                         }
#                         
#                         # store conditions
#                         calc_df <- bind_rows(calc_df,
#                                              sim_df %>%
#                                                 mutate(
#                                                      cov=cov,
#                                                      sens_quest=sens_quest,
#                                                      true_inc=true_inc
#                                                 ))
#                 }
#         }
# }
# 
# senspec_df<- data.frame(
#         variables= c("Reduced Panel", "Reduced Panel","CtxB Panel","CtxB Panel"),
#         questionnaire = c(0,1,0,1),
#         sens = c(0.498,0.498, 0.498-0.1, 0.498-.1),
#         spec = c(0.922,0.691,0.922,0.922)
# )
# 
# 
#  # ref_df<- calc_df  %>% group_by(simulation,variables,vax_status) %>% 
#  #                summarize(n=n(),
#  #                          true_inc=mean(infected_recently),
#  #                          sero_inc=mean(seropos)
#  #                          ) %>% left_join(senspec_df)
# 
#  
# 
#  df1 <-calc_df %>% filter(variables=="Reduced Panel") %>%  
#         group_by(simulation, questionnaire,variables,cov,true_inc,sens_quest) %>%
#          summarize(
#                  n=n(),
#                   seropos=sum(seropos)
#          )  %>%  
#          left_join(senspec_df) %>%  
#          group_by(simulation, questionnaire,variables,cov,true_inc,sens_quest) %>%
#          mutate(adj_casesOLD = adjSeropos(sens,0.922,seropos,n),
#                 adj_cases = adjSeropos(sens,spec,seropos,n)
#                 
#                                     ) %>% 
#          group_by(simulation,cov,true_inc,sens_quest) %>%
#          summarize(n=sum(n),
#                    strategy0 = sum(seropos)/n,
#                    strategy1 = sum(adj_casesOLD)/n,
#                    strategy2 = sum(adj_cases)/n
#                    )  %>%
#          gather(strategy,value,-c(simulation,n,cov,true_inc,sens_quest)) 
#          
#  
# 
#  df2 <- calc_df %>% filter((questionnaire==0 & variables=="Reduced Panel") |
#                             (questionnaire==1 & variables=="CtxB Panel")
#                     )  %>%
#  group_by(simulation, questionnaire,variables,cov,true_inc,sens_quest) %>%
#          summarize(
#                  n=n(),
#                  seropos=sum(seropos)
#          )  %>%  
#          left_join(senspec_df) %>%  
#          group_by(simulation, questionnaire,variables,cov,true_inc,sens_quest) %>%
#          mutate(adj_cases = adjSeropos(sens,spec,seropos,n)
#          ) %>% 
#          group_by(simulation,cov,true_inc,sens_quest) %>%
#          summarize(n=sum(n),
#                    strategy3 = sum(seropos)/n,
#                    strategy4 = sum(adj_cases)/n
#          )  %>%
#          gather(strategy,value,-c(simulation,n,cov,true_inc,sens_quest))
#  
#  
#  
#  
# bind_rows(df1,df2) %>% filter(sens_quest %in% c(0.5,0.9)) %>%
#          ggplot(aes(x=factor(sens_quest),y=value))+
#          geom_boxplot(aes(col=strategy),fill=NA)+
#          facet_grid(true_inc~cov)+
#         geom_hline(yintercept = c(0.05,0.15),lty=2)
#  
# 
# 
# 
# 
# 
# 
# 
# # # calculate incidence
# # tmp_calc_df <- sim_df %>%
# #         group_by(simulation,questionnaire) %>%
# #         summarize(n=n(),
# #                   infected_recently=sum(infected_recently),
# #                   seropos=sum(seropos),
# #                   true_inc = infected_recently/n,
# #                   crude_inc=seropos/n,
# #                   adj_inc1= adjSeropos(sens_novax,spec_novax,seropos,n)/n,
# #                   adj_inc2= adjSeropos(sens_novax,spec_vax,seropos,n)/n
# #         ) %>%
# #         mutate(final_rate = ifelse(questionnaire==0,adj_inc1,adj_inc2)) %>%
# #         group_by(simulation) %>%
# #         summarize(
# #                 true_cases =sum(infected_recently),
# #                 # old_adj_cases = sum(adj_inc1*n),
# #                 new_adj_cases = sum(final_rate*n),
# #                 n=sum(n)
# #         ) %>%
# #         mutate(
# #                 calc_inc= new_adj_cases/n,
# #                 cov=cov,
# #                 sens_quest=sens_quest,
# #                 true_inc=true_inc
# #         )
# # 
# # calc_df <- bind_rows(tmp_calc_df,calc_df)
# 
# calc_df %>% ggplot(aes(y=calc_inc))+
#                 geom_boxplot()
# 
# 
# plot_df <- calc_df %>% mutate(true_inc_fac=factor(true_inc,
#                                                   labels=c("5 cases per 100 pop",
#                                                            "15 cases per 100 pop"))) %>%
#         mutate(cov_fac=factor(cov,
#                               labels=c("25% Coverage",
#                                        "50% Coverage",
#                                        "75% Coverage")))  %>%
#         mutate(calc_inc=100*calc_inc)
# 
# 
# 
# calc_df%>%
#         ggplot(aes(x=factor(sens_quest),y=calc_inc,col=factor(sens_quest)))+
#         geom_boxplot()+
#         xlab("Sensitivity of Questionnaire")+
#         ylab("Estimated Incidence (Cases per 100 pop)")+
#         # geom_hline(data = plot_df %>% distinct(true_inc,cov_fac,
#         #                                        true_inc_fac),
#         #            aes(yintercept = true_inc*100),
#         #            lty=2
#         # )+
#         scale_color_viridis_d()+
#         # facet_grid(true_inc_fac~cov_fac)+
#         theme_bw()+
#         theme(legend.position = "none",
#               panel.grid.major.x = element_blank()
#         )
