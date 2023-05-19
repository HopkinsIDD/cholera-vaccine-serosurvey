

# find the layout and file name
matched_layout <- readxl::read_excel("data/raw_data/layout_lookup/layout_register 2.xlsx") %>% 
        filter(batch==batch_id)
path_csv <- paste0("data/raw_data/summaryFI/",matched_layout  %>% pull(data_file))
extract <- extractData(path_csv)


layout_sheets <- excel_sheets(paste0("data/raw_data/layout_lookup/layouts/new_layouts/",
                                     matched_layout$file))
version <- "1.0"
if("version" %in% layout_sheets){
        version <- read_excel(paste0("data/raw_data/layout_lookup/layouts/new_layouts/",
                                     matched_layout$file),sheet="version",
                              range=c("A1:A2")
        )        
}


# run extraction code based on that version
source(here::here(paste0("code/R/qc_files/extract_v",unlist(version),".R")))