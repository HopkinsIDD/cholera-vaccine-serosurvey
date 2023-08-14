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
                select(batch,data_file,layout_file=file) %>%
                filter(batch %in% c(
                         "20220728_p0040"))

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
        
        #save this_merge
        write_rds(this_merge,glue::glue("{new_dir}/{this_batch$batch}_merge.rds"))

}

        
        

        #for each isotype, conduct QA/QC
                #bead count
                #coerricient of variation
                #separate troublesome wells
                #calculate the netMFI
                #create dilution object and fit model
                #calculate RAU for samples and controls
                #create html report

        #output object in folder
