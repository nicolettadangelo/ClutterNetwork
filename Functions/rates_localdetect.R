rates.localdetect <- function(x){
  if(is.multitype(x$X)){
      calculate_statistics(table(x$id, x$X$data$pattern))
  } else {
    stop("No rates to compute, as no true classification is known")
  }
}


