

#run main figures as a word document
rmarkdown::render("code/Rmd/main_fig.Rmd",
                  knit_root_dir="../..",
                  output_dir="figures",
                  output_format="word_document"
                  )

#run main figures as a html document
rmarkdown::render("code/Rmd/main_fig.Rmd",
                  knit_root_dir="../..",
                  output_dir="figures",
                  output_format="html_document"
)

# #run main text statistics
# rmarkdown::render("code/Rmd/main_text.Rmd",
#                   knit_root_dir="../..",
#                   output_dir="figures",
#                   output_format="html_document"
# )
# 
# #run supplement
# rmarkdown::render("code/Rmd/supp_fig.Rmd",
#                   knit_root_dir="../..",
#                   output_dir="figures",
#                   output_format="word_document"
# )