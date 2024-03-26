rm(list=ls())

source("code/R/packages.R")
source("code/R/utils.R")


final_df <- read_rds("data/generated_data/analysis_data/final_df.rds")
final_wide <-read_rds("data/generated_data/analysis_data/final_wide.rds")



# compare_data_wide <- final_df %>%
#         # only individuals from bandladesh
#         filter(cohort %in% c("BGD Vaccinee","SMIC/PIC")) %>%
#         # recent
#         filter(day>3,day<50) %>%
#         #only cases and vaccinees
#         filter(status %in% c("Case","Vaccinee"))%>%
#         select(sample,id,isotype,antigen_pretty,marker,RAU_value,day,age,status) %>% 
#         select(-isotype,-antigen_pretty) %>%
#         # mutate(RAU_value=log10(1/RAU_value)) %>%
#         spread(marker,RAU_value)
# mds_obj <- bind_rows(
#         compare_data_wide %>% 
#                 generate_MDS() %>% mutate(variables="All"),
#         compare_data_wide %>% 
#                 generate_MDS(regex_string="IgG_Ogawa|IgG_Inaba|IgG_CtxB") %>% mutate(variables="Original"),
#         compare_data_wide %>% 
#                 generate_MDS(regex_string="CtxB") %>%  mutate(variables="CtxB")
# )


compare_data_wide <- final_df %>%
        # only individuals from bandladesh
        filter(cohort %in% c("BGD Vaccinee","SMIC/PIC", "HTI Vaccinee")) %>%
        #only cases and vaccinees
        filter(status %in% c("Case","Vaccinee"))%>%
        select(sample,id,isotype,antigen_pretty,marker,RAU_value,day,age,status,cohort) %>% 
        select(-isotype,-antigen_pretty) %>%
        # mutate(RAU_value=log10(1/RAU_value)) %>%
        spread(marker,RAU_value)


mds_obj <- bind_rows(
        compare_data_wide %>% filter(day<=45) %>%
                generate_MDS() %>% mutate(time_window="<=45 days"),
        compare_data_wide %>% filter(day>45 & day<=120)%>%
                generate_MDS() %>% mutate(time_window="46-120 days"),
        compare_data_wide %>% filter(day>120 & day<=200)%>%
                generate_MDS() %>% mutate(time_window="121-200 days")
)%>%
        mutate(time_window=factor(time_window,
                                  levels=c("<=45 days",
                                           "46-120 days",
                                           "121-200 days"
                                  )
        ))


write_rds(mds_obj,"data/generated_data/analysis_objects/mds_obj.rds")


#####

mds_obj %>%
        group_by(time_window,status) %>%
        summarise(avg_D1=mean(D1),
                  avg_D2=mean(D2)
                  ) %>%
        gather(dimension,value,-c(time_window,status)) %>%
        spread(status,value) %>%
        mutate(distance=(Case-Vaccinee)^2) %>%
        group_by(time_window) %>%
        summarise(euclidean=sqrt(sum(distance)))



# Number of permutations
n_permutations <- 1000  # You can adjust the number of permutations as needed

# Create a vector to store the permuted distances
permuted_distances <- data.frame()
tmp_obj_1 <- mds_obj %>% filter(time_window=="<=45 days")
tmp_obj_2 <- mds_obj %>% filter(time_window=="46-120 days")
tmp_obj_3 <- mds_obj %>% filter(time_window=="121-200 days")

for (i in 1:n_permutations) {

              tmp_euc <- bind_rows(
                        tmp_obj_1 %>% mutate(status=sample(tmp_obj_1$status)),
                        tmp_obj_2 %>% mutate(status=sample(tmp_obj_2$status)),
                        tmp_obj_3 %>% mutate(status=sample(tmp_obj_3$status))
                ) %>%
                group_by(time_window,status) %>%
                summarise(avg_D1=mean(D1),
                          avg_D2=mean(D2)
                ) %>%
                gather(dimension,value,-c(time_window,status)) %>%
                spread(status,value) %>%
                mutate(distance=(Case-Vaccinee)^2) %>%
                group_by(time_window) %>%
                summarise(euclidean=sqrt(sum(distance)))
        
        
        permuted_distances <- bind_rows(permuted_distances,tmp_euc)
        
}


permuted_distances %>%
        group_by(time_window) %>%
        summarize(quantile(euclidean,0.975))




# #try the pca analysis
# pca_data <- final_wide %>%
#         filter(cohort %in%   c("BGD Vaccinee","HTI Vaccinee","SMIC/PIC")) %>%
#         filter(status %in% c("Case","Vaccinee"))%>%
#         filter(day<=45) %>%
#         filter(day>2)%>%
#         select(status,cohort,matches("Ogawa|Inaba|TcpA|CtxB|IgG_O139"))
# # select(vaxinf_class,matches("Ogawa|Inaba|TcpA|IgG_O139"))
# 
# standardized_data <- scale(pca_data[,-c(1,2)])
# 
# pca_result <- prcomp(standardized_data, center = TRUE, scale = TRUE)
# 
# 
# 
# 
# 
# 
# 
# loading_data <- as.data.frame(pca_result$rotation)
# loading_data$Variable <- rownames(loading_data)
# 
# ggplot(loading_data, aes(x = PC1, y = PC2, label = Variable)) +
#         geom_hline(yintercept = 0, color = "gray", linetype = 2) +
#         geom_vline(xintercept = 0, color = "gray", linetype = 2) +
#         geom_text(hjust = 0, vjust = 0) +
#         theme_minimal()
# 
# score_data <- as.data.frame(pca_result$x)
# score_data$Group <- pca_data$status  # Replace with your actual grouping variable
# 
# ggplot(score_data, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
#         geom_point() +
#         theme_minimal()
# 
# grouped_means <- score_data %>%
#         group_by(Group) %>%
#         summarize(Mean_PC1 = mean(PC1), Mean_PC2 = mean(PC2))
# 
# group1_mean <- grouped_means %>% filter(Group == "Vaccinee")
# group2_mean <- grouped_means %>% filter(Group == "Case")
# 
# euclidean_distance <- sqrt((group1_mean$Mean_PC1 - group2_mean$Mean_PC1)^2 + (group1_mean$Mean_PC2 - group2_mean$Mean_PC2)^2)
# 
# # Number of permutations
# n_permutations <- 1000  # You can adjust the number of permutations as needed
# 
# # Create a vector to store the permuted distances
# permuted_distances <- numeric(n_permutations)
# 
# for (i in 1:n_permutations) {
#         # Permute the group labels
#         score_data$Group <- sample(score_data$Group)
#         
#         # Recalculate the distance with permuted labels
#         grouped_means_permuted <- score_data %>%
#                 group_by(Group) %>%
#                 summarize(Mean_PC1 = mean(PC1), Mean_PC2 = mean(PC2))
#         
#         group1_mean_permuted <- grouped_means_permuted %>% filter(Group == "Vaccinee")
#         group2_mean_permuted <- grouped_means_permuted %>% filter(Group == "Case")
#         
#         permuted_distances[i] <- sqrt((group1_mean_permuted$Mean_PC1 - group2_mean_permuted$Mean_PC1)^2 +
#                                               (group1_mean_permuted$Mean_PC2 - group2_mean_permuted$Mean_PC2)^2)
# }
# 
# # Calculate the p-value
# observed_distance <- euclidean_distance  # Use the observed distance calculated earlier
# p_value <- sum(permuted_distances >= observed_distance) / n_permutations

# scree_data <- data.frame(PC = 1:length(pca_result$sdev), Variance = pca_result$sdev^2/sum(pca_result$sdev^2))
# 
# ggplot(scree_data, aes(x = PC, y = Variance)) +
#         geom_point(shape = 21, fill = "blue", color = "black") +
#         geom_line(color = "red") +
#         labs(x = "Principal Component", y = "Proportion of Variance Explained") +
#         theme_minimal()


# biplot_data <- as.data.frame(pca_result$x)
# biplot_data$Group <- rownames(biplot_data)
# 
# ggplot(biplot_data, aes(PC1, PC2, label = Group)) +
#         geom_hline(yintercept = 0, color = "gray", linetype = 2) +
#         geom_vline(xintercept = 0, color = "gray", linetype = 2) +
#         geom_text(hjust = 0, vjust = 0) +
#         theme_minimal()
        