#load functions
pacman::p_load(
        tidyverse,
        readxl
        )

source("code/R/qc_files/utils_qc.R")

#load the matching file
matched_layout <- readxl::read_excel("data/raw_data/layout_lookup/layout_register 2.xlsx")

#make sure all the datafiles and and layout files are available
data_files_dir <-"data/raw_data/summaryFI"
data_files <- setdiff(
        list.files(data_files_dir),
        list.dirs(data_files_dir,full.names = FALSE)
)

layout_files_dir <- "data/raw_data/layout_lookup/layouts/new_layouts"
layout_files <- setdiff(
        list.files(layout_files_dir),
        list.dirs(layout_files_dir,full.names = FALSE)
)

batches <- matched_layout %>%
                filter(file %in% layout_files) %>%
                filter(data_file %in% data_files) %>%
                select(batch,data_file,layout_file=file) 

#match directory path
match_dir_path <- "data/generated_data/quality_control/"


#create folders with all of the merged data
for(i in 1:nrow(batches)){
        
        this_batch <- batches[i,]
        cat(this_batch$batch,"\n")
        
        #create the directory to store the output
        new_dir <- paste0("data/generated_data/quality_control/",this_batch$batch)
        dir.create(new_dir)
        
        #extract the data
        this_data_path <- paste(data_files_dir,this_batch$data_file,sep = "/")
        this_data_extract <- extractExponentcsv(this_data_path)
        
        #extract the loyout
        this_layout_path <- paste(layout_files_dir,this_batch$layout_file,sep = "/")
        this_layout_extract <- extractLayout(this_layout_path)
        
        #merge the data and the layout
        if (this_layout_extract$summary$version=="1.0"){
                this_merge <- left_join(this_data_extract$out_df,
                                        this_layout_extract$out_df,
                                        by = c("location", "batch")
                ) 
        }
        if (this_layout_extract$summary$version!="1.0"){
                this_merge <- left_join(this_data_extract$out_df,
                                        this_layout_extract$out_df,
                                        by = c("location", "batch", "BeadID")
                ) 
        }
        
        this_merge <- this_merge %>%
                #make all VCC the same
                mutate(antigen_exponent=ifelse(str_detect(antigen_exponent,"VCC"),"VCC",antigen_exponent))%>%
                #create a variable for the type of measurement
                mutate(no_control=str_remove(sample,"Control\\: ")) %>%
                mutate(type=case_when(
                        str_detect(no_control,"\\:") ~ "Standard",
                        str_detect(sample,"Control\\: ")&(!str_detect(no_control,"\\:")) ~ "Control",
                        str_detect(sample,"Blank|Sialidase|Pool|PBS") ~ "Control",
                        TRUE  ~ "Sample"
                )) %>%
                select(-no_control)
        
        #run the checks
        #do the number of rows match up between the sets
        if(nrow(this_data_extract$out_df)!=nrow(this_layout_extract$out_df) ) {
                if(nrow(this_data_extract$out_df)!=nrow(this_layout_extract$out_df)*length(unique(this_layout_extract$summary$antigen$antigen_name))){
                        stop(paste0(this_data_extract$batch,
                                    ":Layout rows (",nrow(this_layout_extract$out_df),
                                    ") doesnt match .csv rows (",
                                    nrow(this_data_extract$out_df),")"))
                }}
        
        #does every sample have sample id, batch id bead count isotype
        if(select(this_merge,batch,Count,sample,isotype) %>% is.na() %>% any()){
                this_merge %>% filter(is.na(batch)|is.na(Count)|is.na(sample)|is.na(isotype)) %>%
                        select(batch,location,sample_exponent,Count,sample,isotype)
                stop("Something went wrong with the merge for the above data")
        }
        
        
        #save the extracts and merge
        this_out <- list(
                data_extract= this_data_extract,
                layout_extract= this_layout_extract,
                extract_merge= this_merge
        )
        
        write_rds(this_out,
                  glue::glue("{new_dir}/{this_batch$batch}_merge.rds"))
        
}


#calculate the averages and fit standard curves for each batch
standards_curves <- list()
control_data <- data.frame() 
sample_data <- data.frame()

for(i in 1:nrow(batches)){

        cat(batches$batch[i],"\n")

        this_batch<- batches$batch[i]
        this_merge <- glue::glue("{match_dir_path}{this_batch}/{this_batch}_merge.rds") %>%
                read_rds()
        multiple_beadlots <- this_merge$layout_extract$summary$antigen %>%
                                distinct(antigen_name,lot_number) %>%
                                count(antigen_name) %>%
                                pull(n) %>% max()
        
        if(multiple_beadlots>1) warning("multiple bead lots, qc not done")
        if(matched_layout$type[i]=="experiment") warning("experimental data, qc not done")
        
        if(matched_layout$type[i]!="experiment" & multiple_beadlots==1){
                # first identify points with low quality
                count_threshold <- 30
                low_beads_location <- filter(this_merge$extract_merge, Count<count_threshold) %>%
                        filter(!str_detect(sample,"Blank")) 
                
                cv_threshold <- 0.2
                high_cv_location <- filter(this_merge$extract_merge, Count>=count_threshold) %>%
                        #work on the log scale
                        mutate(MFI=log10(MFI))%>%
                        #calculate the average MFI and CV
                        group_by(batch,sample,isotype,subclass,antigen_exponent) %>%
                        mutate( replicates=n(),  
                                avgMFI=mean(MFI),
                                sdMFI=sd(MFI)
                        ) %>%
                        mutate(cvMFI=sdMFI/avgMFI,
                               diff_avg=MFI-avgMFI
                        )  %>%
                        filter(cvMFI> cv_threshold) 
                
                furthest_points<- high_cv_location %>%
                        #identify the MFI that is farthest from the average
                        slice_max(diff_avg,n=1,
                                  with_ties = FALSE
                        ) %>%
                        ungroup()%>%
                        #dont include the blank in here to allow net mfi calculation
                        filter(!str_detect(sample,"Blank")) %>%
                        distinct(location)
                
                
                # remove these locations
                these_averageMFI <- this_merge$extract_merge %>%
                        #work on the log scale
                        mutate(MFI=log10(MFI))%>%
                        #remove the locations with low beads and high cv
                        anti_join(distinct(low_beads_location,location),by="location") %>%
                        anti_join(furthest_points,by="location")%>%
                        group_by(date,batch,type,sample,isotype,subclass,antigen_exponent) %>%
                        summarize( replicates=n(),  
                                   avgMFI=mean(MFI),
                                   sdMFI=sd(MFI)
                        ) %>%
                        mutate(cvMFI=sdMFI/avgMFI) %>%
                        #if any with a high cv still remain, remove it
                        filter(cvMFI<=cv_threshold ) %>%
                        #only include those with more than 1 replicate
                        filter(replicates>1) %>%
                        ungroup()
                
                
                # save them in each folder
                fits_obj <- list()
                control_df <- data.frame()
                sample_df <- data.frame()
                standard_df <- data.frame()
                
                #identify the isotypes in data
                these_isotypes <- unique(these_averageMFI$isotype)
                for(iso in these_isotypes){
                        
                        #identify the blanks
                        this_blank_df <- these_averageMFI %>%
                                filter(isotype==iso) %>%
                                filter(type=="Control") %>%
                                filter(str_detect(sample,"Blank")) %>%
                                select(sample,antigen_exponent,blank=avgMFI)
                        
                        #calculate the netmfi
                        this_netmfi <- these_averageMFI %>%
                                filter(isotype==iso) %>%
                                anti_join(this_blank_df,by = c("sample", "antigen_exponent")) %>%
                                left_join(select(this_blank_df,-sample),by="antigen_exponent") %>%
                                mutate(netmfi_value=avgMFI-blank) %>%
                                mutate(netmfi_value=ifelse(netmfi_value<=0,0,netmfi_value))
                        
                        
                        #identify the samples from the standard curve
                        this_standard_df <- this_netmfi %>%
                                filter(type=="Standard") %>%
                                mutate(dilution=str_remove(sample,"\\,"))%>%
                                mutate(dilution=str_extract(dilution,"1:[0-9]{1,7}")) %>%
                                mutate(dilution=str_remove(dilution,"1:"))%>%
                                mutate(dilution=as.numeric(dilution)) 
                        
                        #on plate 22, remove the 160 dilution
                        if(this_batch=="20220202_p0022"){
                                this_standard_df <- this_standard_df %>%
                                        filter(dilution!=160)
                        }
                        
                        #store the standard data.frame
                        standard_df <- bind_rows(standard_df,this_standard_df)
                        
                        #estimate the standard curves
                        antigens <- unique(this_standard_df$antigen_exponent)
                        for(anti in antigens){
                                #fit the standard curve
                                this_fit <- this_standard_df %>%
                                        filter(antigen_exponent==anti) %>%
                                        fitLL_model()
                                
                                
                                if(this_fit$model_type=="Net MFI, 3P-log-logistic"){
                                        this_rau_df <- this_netmfi %>%
                                                filter(antigen_exponent==anti) %>%
                                                filter(type %in% c("Control","Sample"))%>%
                                                mutate(RAU_value=sapply(avgMFI, FUN=get3PLL_RAU,fit=this_fit))
                                        
                                }
                                
                                if(this_fit$model_type=="MFI, 5P-log-logistic"){
                                        this_rau_df <- this_netmfi %>%
                                                filter(antigen_exponent==anti) %>%
                                                filter(type %in% c("Control","Sample"))%>%
                                                mutate(RAU_value=sapply(avgMFI, FUN=get5PLL_RAU,fit=this_fit))
                                }
                                # fit to standard curve
                                fits_obj[[iso]][[anti]] <- this_fit
                                # store the RAU of samples
                                if(this_fit$model_type!="Error"){
                                        control_df <- bind_rows(control_df,filter(this_rau_df,type=="Control"))
                                        sample_df <- bind_rows(sample_df,filter(this_rau_df,type=="Sample"))
                                }
                                
                        }
                        
                }
                
                #save individual objects in folder
                write_rds(fits_obj, glue::glue('{match_dir_path}{this_batch}/{this_batch}_fits.rds'))
                write_rds(standard_df, glue::glue('{match_dir_path}{this_batch}/{this_batch}_standards.rds'))
                write_rds(control_df, glue::glue('{match_dir_path}{this_batch}/{this_batch}_RAU_controls.rds'))
                write_rds(sample_df, glue::glue('{match_dir_path}{this_batch}/{this_batch}_RAU_sample.rds'))
                
                #identify the R^2 values for the plate
                r_squareds <- bind_lapply(your_list=fits_obj, your_id="isotype",
                                          your_function=function(x) bind_lapply(x,your_id="antigen",
                                                                                your_function = function(y) data.frame(r_squared=y$R.squared)
                                          ))
                #only calculate those with all curves with high R^2
                if(nrow(r_squareds)>0){
                        if(all(pull(r_squareds,r_squared)>0.9)){
                                #save curve objects in one place
                                standards_curves[[this_batch]] <- fits_obj
                                write_rds(standards_curves, glue::glue('{match_dir_path}standards_curves.rds'))
                                
                                #save controls data in one place
                                if(nrow(control_df)>0) control_data <- bind_rows(control_data,control_df)
                                write_rds(control_data, glue::glue('{match_dir_path}control_data.rds'))
                                
                                #save samples data in one place
                                if(nrow(sample_df)>0) sample_data <- bind_rows(sample_data,sample_df)
                                write_rds(sample_data, glue::glue('{match_dir_path}sample_data.rds'))
                        }
                }
        }
}

#predicted curves to plot against
predicted_curves <- bind_lapply(your_list=standards_curves, your_id="batch",
            your_function=function(x) bind_lapply(x,your_id="isotype",
                        your_function = function(y) bind_lapply(y, your_id="antigen",
                            your_function = pred_curve
                    ))) %>%
        mutate(subclass=case_when(
                str_detect(batch,"IGA_1") ~ "1",
                str_detect(batch,"IGA_2") ~ "2",
                str_detect(batch,"IGG_1") ~ "1",
                str_detect(batch,"IGG_2") ~ "2",
                str_detect(batch,"IGG_3") ~ "3",
                str_detect(batch,"IGG_4") ~ "4",
                TRUE ~ "total"
        ) )

write_rds(predicted_curves,glue::glue('{match_dir_path}predicted_curves.rds'))
                        




#run qc report
for(i in 1:nrow(batches)){

        cat(batches$batch[i],"\n")
        
        this_batch<- batches$batch[i]
        this_folder <-here::here(glue::glue("{match_dir_path}{this_batch}/"))
        this_merge <- glue::glue("{this_folder}{this_batch}_merge.rds") %>%
                read_rds()
        
                ##do this for every isotype
                #identify the isotypes in data
                these_isotypes <- unique(this_merge$extract_merge$isotype)
                for(this_iso in these_isotypes){
                        
                        #run the QA/QC analysis
                        rmarkdown::render(here::here('code/Rmd/qc_report.Rmd'),
                                          knit_root_dir = here::here(),
                                          params = list(
                                                  batch_id = this_batch,
                                                  merge_obj = this_merge,
                                                  iso=this_iso,
                                                  curves=read_rds(glue::glue(here::here('{match_dir_path}predicted_curves.rds'))),
                                                  controls=read_rds(glue::glue(here::here("{match_dir_path}control_data.rds")))
                                          ),
                                          output_file = here::here(glue::glue('{match_dir_path}{this_batch}/{this_batch}_{this_iso}_qc_report.html')
                                          ),
                                          quiet = TRUE
                                          )
                        
                        
                }       
}



# go through all the .html files to find the ones that are worth keeping 
i=1

#note in an excel file which data are worth keeping
this_batch<- batches$batch[i]
this_folder <-here::here(glue::glue("{match_dir_path}{this_batch}/"))
this_merge <- glue::glue("{this_folder}{this_batch}_merge.rds") %>%
        read_rds()

##do this for every isotype
#identify the isotypes in data
these_isotypes <- unique(this_merge$extract_merge$isotype)
for(this_iso in these_isotypes){
        browseURL(here::here(glue::glue('{match_dir_path}{this_batch}/{this_batch}_{this_iso}_qc_report.html')))
}
i <- i+1


#clean up the sample data
accept_reject_df <- read_excel("data/raw_data/accept_reject_plates.xlsx")
accepted_plates <-accept_reject_df %>%
        filter(accept_reject=="Accept")%>%
        pull(batch_id)
sample_data <- read_rds('data/generated_data/quality_control/sample_data.rds')


#create a clean version of the data sets
clean_sample_data <- sample_data %>%
        filter(batch %in% accepted_plates) %>%
        mutate(sample=str_replace(sample," ","_")) %>%
        mutate(sample=toupper(sample))%>%
        mutate(sample=str_replace(sample,"CHX_RB","RB"))%>%
        mutate(id=str_remove(sample,"_D[0-9]{1,4}")) %>% 
        mutate(studycode=str_extract(id,"RB|E01JH|IMS|R1|R2|P|S"))%>%
        mutate(day=str_extract(sample,"_D[0-9]{1,4}")) %>% 
        mutate(day=str_remove(day,"_D")) %>%
        mutate(day=as.numeric(day)) %>%
        group_by(isotype,subclass,id,day,antigen_exponent) %>%
        filter(date==max(date))%>%
        ungroup()

write_rds(clean_sample_data,'data/generated_data/quality_control/clean_sample_data.rds')


#get a full dataset of all the MFIs
full_mfi <- data.frame()

for(i in 1:nrow(batches)){

        cat(batches$batch[i],"\n")
        
        this_batch<- batches$batch[i]
        this_merge <- glue::glue("{match_dir_path}{this_batch}/{this_batch}_merge.rds") %>%
                read_rds()
        
        this_mfi <- this_merge$extract_merge
        
        full_mfi <- bind_rows(full_mfi,this_mfi)

}

write_rds(full_mfi,"data/generated_data/quality_control/complete_raw_data.rds")



