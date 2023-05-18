set_paths <- function(info = Sys.info()) {
  # you can add your user and node if you want to use other than the generic paths
  
  if (info["user"] == "andrewazman") {
    # paths for Andrew's laptop
    paths = list(
      dropbox = "/Users/andrewazman/Dropbox/Haiti Serosurvey Cholera Typhoid"
    )
  } else if (info["user"] == "hanmeng") {
    paths = list(
      dropbox = "/Users/hanmeng/Dropbox/Haiti Serosurvey Cholera Typhoid"
    )
  } else if (info["user"] == "forrestjones") {
          paths = list(
                  dropbox = "/Users/forrestjones/Dropbox/Haiti Serosurvey Cholera Typhoid"
          )
  }   else {
    # generic paths, if the machine is not above
    paths = list(
      dropbox = "~/Dropbox/Haiti Serosurvey Cholera Typhoid"
    )
    warning("User and computer not recognized, please edit R/set_paths.R.")
  }
  
  paths$data <- file.path(paths$dropbox)
  #paths$cleaning <- file.path(paths$dropbox, "data-cleaning")
  #paths$data_clean <- file.path(paths$dropbox, "data-clean")
  
  return(paths)
}
