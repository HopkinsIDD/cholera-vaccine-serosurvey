
source("source/final_code/shared/utils.R")
source("source/final_code/shared/utils_vax.R")
source("source/final_code/shared/utils_recon.R")
source("source/final_code/vax_strategy/R/load-data.R")


library(caret)


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

multi_class_data %>% distinct(status,vaxinf_class,day,id) %>%
        group_by(status,vaxinf_class,day) %>%
                count()

multi_class_data %>% distinct(status,vaxinf_class,day,id) %>%
        group_by(vaxinf_class) %>%
        summarize(n(),
                  length(unique(id))
                  )

#get data into a good form for the caret fitting
caret_data <- multi_class_data %>%
        select(sample,id,vaxinf_class,isotype,antigen_pretty,full_test,value,day,age,status) %>% 
        select(-isotype,-antigen_pretty) %>%
        spread(full_test,value) %>%
        mutate(vaxinf_class=str_remove_all(vaxinf_class," |,")) %>%
        mutate(vaxinf_class=factor(vaxinf_class))

# lets run a bunch of models and complete feature selections
# start with 21 markers including flu
# calculate the variable importance, out of bag error
# (can also do the AUC for each and compare...)


full_var <- multi_class_data %>% distinct(full_test) %>%
                filter(!str_detect(full_test,"LTB|LTh|CTHT|Flu"))

red12 <- full_var %>%
        filter(!str_detect(full_test,"VCC|Sialidase|O139"))

var_list <-list(`21 markers`=list(var=full_var),
                `12 markers`=list(var=red12))


fit_list <-   var_list %>% 
                lapply(FUN=function(x){
                        
                        formula <- make_formula("vaxinf_class",pull(x$var,full_test))
                        fit <- runCaretLOOCV(formula, dat=caret_data)
                        importance <- caret::varImp(fit,scale=FALSE)
                        
                        x$formula <- formula
                        x$fit <- fit 
                        x$importance <- importance$importance %>%
                                        mutate(marker=rownames(.))
                        return(x)
                        
                })


last_model <- fit_list$`12 markers`

for(n in 11:2){
        
        #name the model
        new_name <- paste(n,"markers")
        
        #identify least important marker
        least_important <- last_model$importance %>%
                                filter(Overall==min(Overall)) %>% pull(marker)
        
        cat(new_name," Drop: ", least_important," \n")
        
        
        #run the caret model without the least important model
        new_var <- filter(last_model$var %>%
                                  filter(!str_detect(full_test,least_important)))
        new_formula <- make_formula("factor(vaxinf_class)",pull(new_var,full_test))
        new_fit <- runCaretLOOCV(new_formula, dat=caret_data)
        new_importance <- caret::varImp(new_fit)
        
        #output the model and name for it to take in the next one
        last_model <- list(var=new_var,
                           formula=new_formula,
                           fit=new_fit,
                           importance=new_importance$importance%>%
                                   mutate(marker=rownames(.))
                           )
        fit_list[[new_name]] <- last_model
        
}

loss_list <- fit_list %>% lapply(function(x){
        
                leastimportant <-  x$importance %>%
                        filter(Overall==min(Overall)) %>% pull(marker)
                x$fit$results %>%
                        filter(logLoss==min(logLoss)) %>%
                        mutate(lose=leastimportant)
}) 




#### see if we can get away with just a few isotypes or antigen


isotype_var_list <-list(
                `12 markers`=list(var=red12),
                `IgG markers` = list(var=red12 %>%
                                             filter(str_detect(full_test,"IgG"))
                                             ),
                `IgA markers` = list(var=red12 %>%
                                             filter(str_detect(full_test,"IgA"))
                ),
                `IgM markers` = list(var=red12 %>%
                                             filter(str_detect(full_test,"IgM"))
                ),
                `IgG and IgA markers` = list(var=red12 %>%
                                             filter(str_detect(full_test,"IgG|IgA"))
                ),
                `IgG and IgM markers` = list(var=red12 %>%
                                             filter(str_detect(full_test,"IgG|IgM"))
                ))


isotype_fit_list <-   isotype_var_list %>% 
        lapply(FUN=function(x){
                
                formula <- make_formula("vaxinf_class",pull(x$var,full_test))
                fit <- runCaretLOOCV(formula, dat=caret_data)
                importance <- caret::varImp(fit,scale=FALSE)
                
                x$formula <- formula
                x$fit <- fit 
                x$importance <- importance$importance %>%
                        mutate(marker=rownames(.))
                return(x)
                
        })


write_rds(isotype_fit_list,"source/final_code/vax_strategy/generated_rds/vax-classification/isotype_fit_list.rds")


loss_df <-lapply(isotype_fit_list, FUN=function(x){
        x$fit$results %>% filter(logLoss==min(logLoss))
}
) %>% bind_rows(.id="model")


loss_df %>%
        mutate(model=factor(model,levels=model)) %>%
        ggplot(aes(x=model,y=logLoss))+
        theme_bw()+
        geom_col()+
        coord_flip()



antigen_varlist <-list(
        `12 markers`=list(var=red12),
        `CTB and Ogawa markers` = list(var=red12%>%
                                     filter(str_detect(full_test,"CtxB|Ogawa"))
                                     ),
         `Add TcpA only` = list(var=red12%>%
                                     filter(str_detect(full_test,"CtxB|Ogawa|TcpA"))
                                     ),
         `Add Inaba only` = list(var=red12%>%
                                       filter(str_detect(full_test,"CtxB|Ogawa|Inaba"))),
        `CTB and Ogawa markers (IgG and IgM)` = list(var=red12%>%
                                               filter(str_detect(full_test,"IgM|IgG")) %>%
                                               filter(str_detect(full_test,"CtxB|Ogawa"))
        ),
        `Add TcpA only (IgG and IgM)` = list(var=red12 %>%
                                       filter(str_detect(full_test,"IgM|IgG")) %>%
                                       filter(str_detect(full_test,"CtxB|Ogawa|TcpA"))
        ),
        `Add Inaba only (IgG and IgM)` = list(var=red12 %>%
                                       filter(str_detect(full_test,"IgM|IgG")) %>%
                                        filter(str_detect(full_test,"CtxB|Ogawa|Inaba")))
        
        )

antigen_fit_list <-   antigen_varlist %>% 
        lapply(FUN=function(x){
                
                formula <- make_formula("vaxinf_class",pull(x$var,full_test))
                fit <- runCaretLOOCV(formula, dat=caret_data)
                importance <- caret::varImp(fit,scale=FALSE)
                
                x$formula <- formula
                x$fit <- fit 
                x$importance <- importance$importance %>%
                        mutate(marker=rownames(.))
                return(x)
                
        })

write_rds(antigen_fit_list,"source/final_code/vax_strategy/generated_rds/vax-classification/antigen_fit_list.rds")



# red1 <- full_var %>% 
#         filter(!str_detect(full_test,"Luminex_IgM_CtxB")) %>%
#         filter(!str_detect(full_test,"VCC"))
# 
# red2 <- red1 %>% 
#         filter(!str_detect(full_test,"O139"))
# 
# red3 <- red2 %>% 
#         filter(!str_detect(full_test,"Inaba"))
# 
# red4 <- red3 %>% 
#         filter(!str_detect(full_test,"Sialidase|Flu"))
# 
# red5 <- red4 %>% 
#         filter(str_detect(full_test,"Ig[G|A]_CtxB|Ogawa|IgG_TcpA"))
# 
# red6 <- red5 %>% 
#         filter(str_detect(full_test,"IgG"))
        

# make_formula("factor(vaxinf_class)",pull(full_var,full_test))



var_list <-list(full_var,red12)

# names(var_list) <- c("full","model 1","model 2","model 3")





loss_list <-lapply(fit_list, FUN=function(x){
                x$results %>% filter(logLoss==min(logLoss))
                }
) %>% bind_rows(.id="model")




# 
# plot(fit_list$`model 2`)
# plot(loss_list$logLoss)
# plot(fit_list$`model 1` %>% varImp())




#### conduct the pca using 6 marker model
pca <- select(caret_data,matches(full_var %>%
                                         filter(str_detect(full_test,"Ogawa|Inaba|Ig[G|A]_TcpA|Ig[G|A]_CtxB"))%>%
                                         unlist() %>% paste0(collapse = "|"))) %>%
        prcomp(center = TRUE,scale. = TRUE)
var <- factoextra::get_pca_var(pca)

#get proportion of variance
summary(pca)
(pca$sdev^2) /sum((pca$sdev^2))

#supplemental figure for PCA
var$contrib[,1:2] %>% data.frame() %>% 
        mutate(marker=rownames(.)) %>%
        gather(dimension,proportion,-marker)  %>%
        ggplot(aes(x=marker,y=proportion))+
        geom_col()+
        facet_wrap(.~dimension)+
        theme(axis.text.x=element_text(angle=45,hjust=1))


pca_pred <- predict(pca) %>% as.data.frame() %>%
                bind_cols(caret_data)




### write out the files
write_rds(caret_data,"source/final_code/vax_strategy/generated_rds/vax-classification/caret_data.rds")

write_rds(fit_list,"source/final_code/vax_strategy/generated_rds/vax-classification/fit_list.rds")
write_rds(loss_list,"source/final_code/vax_strategy/generated_rds/vax-classification/loss_list.rds")

write_rds(pca,"source/final_code/vax_strategy/generated_rds/vax-classification/pca.rds")
write_rds(var,"source/final_code/vax_strategy/generated_rds/vax-classification/var.rds")
write_rds(pca_pred,"source/final_code/vax_strategy/generated_rds/vax-classification/pca_pred.rds")



