
# we want to estimate 200-day seroincidence in a population

#look at the adult pop value after two doses on these days
days_eval <- c(21,44,90,180) #days post first dose

# Strategy 1: questionairre
#apply the same model over and over again, but have different adjustments
#use the same one you used in the esmate FPR

recent_raw_df <- read_rds("source/final_code/luminex_recommendation/generated_rds/random_forest_results/tvsens-spec/raw_df.rds") %>%
        filter(variables=="Reduced Panel")%>%
        filter(end_window==200) %>%
        filter(age>=18) %>%
        filter(truth==1)%>%
        mutate(new_day=ifelse(day_actual==0,1,day_actual))%>%
        mutate(lday=log(new_day)) %>%
        mutate(lday2=lday^2,
               lday3=lday^3
        )


notrecent_raw_df <- read_rds("source/final_code/luminex_recommendation/generated_rds/random_forest_results/tvsens-spec/raw_df.rds") %>%
        filter(variables=="Reduced Panel")%>%
        filter(end_window==200) %>%
        filter(age>=18) %>%
        filter(truth==0)
        
vax_raw_df <- read_rds("source/final_code/vax_strategy/generated_rds/fpr-estimates/raw_df.rds")%>%
        filter(end_window==200) %>%
        filter(variables=="Reduced Panel IgG") %>%
        filter(age>=18)  %>%
        mutate(new_day=ifelse(day_actual==0,1,day_actual))%>%
        mutate(lday=log(new_day)) %>%
        mutate(lday2=lday^2,
               lday3=lday^3
        )


library(brms)

#estimate sensitivity among recently infected (last 200 days)
sens_fit <-  brm(spec95_seropos ~ lday  + (1|id),
                                   data = recent_raw_df, 
                                   family = bernoulli(link = "logit")
)

sens_linpred <- posterior_linpred(sens_fit,
                                  newdata = distinct(recent_raw_df,id) %>%
                                          expand_grid(new_day=1:200)%>%
                                          mutate(lday=log(new_day)) %>%
                                          mutate(lday2=lday^2,
                                                 lday3=lday^3
                                          ),
                                  ndraws=1000
                )%>% t() %>%
        as.data.frame() %>%
        bind_cols(distinct(recent_raw_df,id) %>%
                          expand_grid(new_day=1:200)) %>%
        gather(.draw,lin_pred,-c(id,new_day)) %>%
        mutate(prob=exp(lin_pred)/(1+exp(lin_pred)))

# tvsens <- sens_linpred %>% group_by(.draw,new_day) %>%
#         summarize(sens=mean(prob)) %>%
#         group_by(new_day) %>%
#         summarise(mean=mean(sens),
#                   lb=quantile(sens,0.025),
#                   ub=quantile(sens,0.975)
#         )
# 
# tvsens %>% ggplot(aes(x=new_day,y=mean))+
#         geom_line()+
#         geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.5)



sens_estim <- sens_linpred %>% group_by(.draw,new_day) %>%
        summarize(sens=mean(prob)) %>%
        group_by(.draw)%>%
        summarize(sens=mean(sens))%>%
        ungroup()%>%
        summarize(mean(sens),
                  sd(sens)
        )



#estimate specificity among not recently infected (last 200 days) but unvaxxed
unvax_spec_fit <-  brm(spec95_seropos ~ 1 + (1|id),
                          data = notrecent_raw_df, 
                          family = bernoulli(link = "logit")
)

unvax_spec_linpred <- posterior_linpred(unvax_spec_fit,
                                  newdata = distinct(notrecent_raw_df,id),
                                  ndraws=1000)%>% t() %>%
        as.data.frame() %>%
        bind_cols(distinct(notrecent_raw_df,id)) %>%
        gather(.draw,lin_pred,-c(id)) %>%
        mutate(prob=exp(lin_pred)/(1+exp(lin_pred)))

unvax_spec_estim <- unvax_spec_linpred %>% group_by(.draw) %>%
        summarize(spec=1-mean(prob)) %>%
        ungroup()%>%
        summarize(spec_mu=mean(spec),
                  spec_sd=sd(spec)
                  )


#estimate specificity among not recently infected (last 200 days)  vaxxed
vax_spec_fit <-  brm(spec95_seropos ~ lday +lday2 + lday3 + (1|id),
                     data = vax_raw_df, 
                     family = bernoulli(link = "logit"))



vax_spec_linpred <- posterior_linpred(vax_spec_fit,
                                  newdata = distinct(vax_raw_df,id) %>%
                                          expand_grid(new_day=days_eval)%>%
                                          mutate(lday=log(new_day)) %>%
                                          mutate(lday2=lday^2,
                                                 lday3=lday^3
                                          ),
                                  ndraws=1000)%>% t() %>%
        as.data.frame() %>%
        bind_cols(distinct(vax_raw_df,id) %>%
                          expand_grid(new_day=days_eval)%>%
                          mutate(lday=log(new_day)) %>%
                          mutate(lday2=lday^2,
                                 lday3=lday^3
                          )) %>%
        gather(.draw,lin_pred,-c(id,new_day)) %>%
        mutate(prob=exp(lin_pred)/(1+exp(lin_pred)))


vax_spec_estim <- vax_spec_linpred %>% group_by(.draw,new_day) %>%
        summarize(spec=1-mean(prob)) %>%
        group_by(new_day)%>%
        summarize(spec_mu=mean(spec),
                  spec_sd=sd(spec)
                  )


strategy1 <- list(fits=list(sens_fit,unvax_spec_fit,vax_spec_fit),
                  estimates= list(sens_estim,unvax_spec_estim,vax_spec_estim)         
                          )

# write_rds(strategy1,"source/final_code/vax_strategy/generated_rds/sim-params/questionairre_params.rds")



# Strategy 2: serology only

# for each day, we need a different dataset/model

base_multiclass <- analysisData$netMFI_data$long_data %>%
        filter(test_type=="Luminex")%>% 
        filter(age>=18) %>%
        #limit vaccinees to only baseline and "recent samples"
        filter(!(status=="Vaccinee" & (!day %in% c(0,days_eval)))) %>%
        mutate(vaxinf_class =case_when(
                status=="Case" & day !=2 & day<=200 ~ "Recently Infected",
                status=="Vaccinee" & day %in% days_eval ~ "Recently Vaccinated",
                TRUE ~ "Neither"
        )) %>%
        left_join(antigen_df) %>%
        select(sample,id,vaxinf_class,isotype,antigen_pretty,full_test,value,day,day_actual,age,status) %>% 
        select(-isotype,-antigen_pretty) %>%
        spread(full_test,value) %>%
        mutate(vaxinf_class=str_remove_all(vaxinf_class," |,")) %>%
        mutate(vaxinf_class=factor(vaxinf_class))
        


# 6 marker formula



library(caret)

data_list <- list(
        day21=base_multiclass %>%
                filter(!(vaxinf_class=="RecentlyVaccinated" & day !=21)) ,
        day44=base_multiclass %>%
                filter(!(vaxinf_class=="RecentlyVaccinated" & day !=44)) ,
        day90=base_multiclass %>%
                filter(!(vaxinf_class=="RecentlyVaccinated" & day !=90)) ,
        day180=base_multiclass %>%
                filter(!(vaxinf_class=="RecentlyVaccinated" & day !=180))
)


included15 <- analysisData$netMFI_data$long_data %>%
                filter(test_type=="Luminex")%>%
                distinct(full_test)%>%
                filter(str_detect(full_test,"CtxB|TcpA|O139|Ogawa|Inaba")) %>%
                unlist()
                

loocv_list <- lapply(data_list, FUN=function(x){
        
        #6 marker original
        # form <- make_formula("vaxinf_class",
        #                           c("Luminex_IgG_OgawaOSPBSA",
        #                             "Luminex_IgA_OgawaOSPBSA",
        #                             "Luminex_IgM_OgawaOSPBSA",
        #                             "Luminex_IgA_CtxB",
        #                             "Luminex_IgG_CtxB",
        #                             "Luminex_IgG_TcpA"
        #                           ))
        form <- make_formula("vaxinf_class", included15)
        
                list(data=x,
                     loocv =  runCaretLOOCV(dat=x,form=form)
                     )}
                  )


write_rds(loocv_list,"source/final_code/vax_strategy/generated_rds/sim-params/serologyonly/loocv_list.rds")


options(mc.cores = 4)

multinomial_list <- lapply(loocv_list, FUN=function(x){
        
        preds <- x$loocv$pred %>%
                left_join(mutate(x$data,rowIndex=1:n()))

        m1 <- brm(pred ~ 1 + (1|id),
                  data = preds %>%
                          filter(obs=="Neither"), 
                  family = categorical(link = "logit"))
        
        m2 <- brm(pred ~ 1 + (1|id),
                  data = preds %>%
                          filter(obs=="RecentlyVaccinated"), 
                  family = categorical(link = "logit")
        )
        
        m3 <- brm(pred ~ log(day_actual)  + (1|id),
                  data = preds %>%
                          filter(obs=="RecentlyInfected"), 
                  family = categorical(link = "logit"))
        
        multinomial <- list(model_Neither=m1,
                           model_RecentlyVaccinated=m2,
                           model_RecentlyInfected=m3)
        
        list(
                data=x$data,
                loocv=x$loocv,
                multinomial=multinomial
        )
})


write_rds(multinomial_list,"source/final_code/vax_strategy/generated_rds/sim-params/serologyonly/multinomial_list.rds")


estimate_list <- multinomial_list %>%
        lapply(FUN=function(x){
                
                list(Neither=estimateAlphaConstant(x$multinomial$model_Neither),
                     RecentlyVaccinated=estimateAlphaConstant(x$multinomial$model_RecentlyVaccinated),
                     RecentlyInfected=estimateAlphaTV(x$multinomial$model_RecentlyInfected)
                     )
                
        })


write_rds(estimate_list,"source/final_code/vax_strategy/generated_rds/sim-params/serologyonly/estimate_list.rds")


##### visualize the estimates
# 
# mvtnorm::rmvnorm(1000,
#                  mean=estimate_list$day21$Neither[[1]],
#                  sigma = estimate_list$day21$Neither[[2]],
# ) %>% summary()
#         # as.data.frame() %>%
#         # mutate(.draw=1:n()) %>%
#         # gather(category,prob,-.draw)%>%
#         # ggplot(aes(x=prob))+
#         #         geom_histogram()+
#         #         facet_wrap(.~category)
        


estimate_df<- estimate_list %>%
        lapply(FUN=function(y){
                tmp <- lapply(y,FUN=function(x)
                        x[[1]] 
                ) %>%
                        bind_rows(.id="truth")
                
        }) %>%
        bind_rows(.id="day_cat") %>%
        mutate(day_vax=str_remove(day_cat,"day")) %>%
        mutate(day_vax=as.numeric(day_vax)) %>%
        gather(category,prob,-c(day_cat,truth,day_vax)) 

write_rds(estimate_df,"source/final_code/vax_strategy/generated_rds/sim-params/serologyonly/estimate_df.rds")


estimate_df %>%  filter(day_vax!=180) %>%
        ggplot(aes(x=factor(day_vax),y=prob,fill=category))+
                geom_col(position = position_dodge2())+
        facet_wrap(.~factor(truth))+
        ylim(0,1)


estimate_df %>% filter(day_vax!=180) %>%
        ggplot(aes(x=factor(truth),fill=prob,y=category))+
        geom_tile()+
        facet_wrap(.~day_vax)+
        scale_fill_distiller(palette = "Oranges",direction=1)



strategy2 <- list(fits=multinomial_list,
                  estimates= list(estimate_list,estimate_df)      
)


write_rds(strategy2,"source/final_code/vax_strategy/generated_rds/sim-params/serologyonly/serologyonly_params.rds")



### compare

# original <- read_rds("source/final_code/vax_strategy/generated_rds/sim-params/serologyonly_params.rds")
# new <- read_rds("source/final_code/vax_strategy/generated_rds/sim-params/serologyonly/serologyonly_params.rds")
# 
# p_6 <- original$estimates[[2]] %>% filter(day_vax!=180) %>%
#         ggplot(aes(x=factor(day_vax),y=prob,fill=category))+
#         geom_col(position = position_dodge2())+
#         facet_wrap(.~factor(truth))+
#         ylim(0,1)+
#         theme_bw()+
#         theme(legend.position = "bottom")
# 
# 
# p_15 <- new$estimates[[2]] %>% filter(day_vax!=180) %>%
#         ggplot(aes(x=factor(day_vax),y=prob,fill=category))+
#         geom_col(position = position_dodge2())+
#         facet_wrap(.~factor(truth))+
#         ylim(0,1)+
#         theme_bw()+
#         theme(legend.position = "bottom")
# 
# plot_grid(p_6,p_15,ncol=1)

