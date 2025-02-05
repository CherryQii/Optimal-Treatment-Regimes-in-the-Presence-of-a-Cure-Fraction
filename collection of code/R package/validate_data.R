# 文件名：R/validate_data.R
# data should be a data frame with columns T, Status, and A
# covriates should be a vector of column names, such as c("X1", "X2", "X3")
validate_data <- function(data, covariates = NULL) {
  required_columns <- c("T", "Status", "A")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("The input data is missing the following required columns:", 
               paste(missing_columns, collapse = ", ")))
  }
  
  if (!all(data$Status %in% c(0, 1))) {
    stop("The `Status` column must only contain values 0 (censored) and 1 (event).")
  }
  
  if (!all(data$A %in% c(0, 1))) {
    stop("The `A` column must only contain values 0 (control) and 1 (treatment).")
  }
  
  if (any(data$T <= 0)) {
    stop("The `T` column must contain positive values only.")
  }
  
  if (!is.null(covariates)) {
    missing_covariates <- setdiff(covariates, colnames(data))
    if (length(missing_covariates) > 0) {
      stop(paste("The input data is missing the following covariate columns:", 
                 paste(missing_covariates, collapse = ", ")))
    }
  }
  return(TRUE)
}
