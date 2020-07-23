


asreml4pin <- function(model1){
  
  x <- model1$vparameters
  xr <- which(names(x) == "units!R")
  vsum <- paste(paste0("V", (1:length(x))[-xr]), collapse = "+")
  vsum <- paste0("(", vsum, ")")
  
  restab <- NULL
  
  for(i in 1:length(x)){
    
    xadd <- eval(parse(text = paste0("vpredict(model1, h", i, " ~ V", i, "/", vsum, ")")))
    xadd$Effect <- names(x)[i]
    
    restab <- rbind(restab, xadd)
    rm(xadd)
    
  }
  
  restab <- subset(restab, Effect != "units!R")
  restab
}
