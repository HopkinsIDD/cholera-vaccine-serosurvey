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


#create a for loop for each plate
i =1 
this_batch <- batches[i,]

        #create the directory to store the output
        dir.create(paste0("data/generated_data/quality_control/",this_batch$batch[i]))

        #extract the data
        this_data_path <- paste(data_files_dir,this_batch$data_file,sep = "/")
        this_data_extract <- extractExponentcsv(this_data_path)
        
        #extract the loyout
        this_layout_path <- paste(layout_files_dir,this_batch$layout_file,sep = "/")
        #this_layout_extract <- extract_layout(this_layout_path)
        
        #merge the data and the layout
        this_merge <- left_join(this_data_extract,
                                this_layout_extract
                                )
        #write_rds(this_merge,)
        

        #for each isotype, conduct QA/QC
                #bead count
                #coerricient of variation
                #separate troublesome wells
                #calculate the netMFI
                #create dilution object and fit model
                #calculate RAU for samples and controls
                #create html report

        #output object in folder
