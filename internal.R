source("util.R")

tryCatch({
  config <- yaml.load_file("local_config.yaml")
  registerDoParallel(cores=config$cores - 1)
}, error=function(err) {
  registerDoSEQ()
})

sec16a.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/internal/section16a.dta"))  

ident.merged.data <- read.dta(paste0(config$data_path, "/mobarak_seasonal_migration/internal/identity_1a-16a-16b_merged.dta"))  
