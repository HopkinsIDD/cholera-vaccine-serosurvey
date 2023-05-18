## runs all plates for fullplate_qc

batch_ids <- c("20221012_p01_boxes_1_and_2",
               "20221013_p02_boxes_3_and_4",
               "20221013_p03_boxes_3_and_4",
               "20221014_p04_boxes_7_and_8",
               "20221014_p05_boxes_9_and_10")
               
for (i in seq_along(batch_ids)){
rmarkdown::render(here::here('code/fullplate_qc.Rmd'),
                  params = list(
                  batch_id = batch_ids[i]),
                  output_file = paste0(set_paths()$data,'/qc_reports/plate_report_',batch_ids[i], 
                                      '.html', sep=''))
}
