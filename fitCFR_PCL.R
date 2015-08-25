accuracyErrorFcn <- function(pars, fcn, fix, obs) {
  if (!paramBounds(c(pars,fix))) {
    return(1000000)
  }
  preds <- fcn(free=pars,fixed=fix,return_dist = FALSE)
  obs <- obs %>% group_by(class,order) %>%
    summarise(acc = sum(score)/4) %>%
    ungroup() %>%
    complete(class,order,fill=list(acc=0))
  err <- sum((preds$acc-obs$acc^2))
  return(err)
}

RT_ErrorFcn <- function(fcn, pars, fix, obs) {
  if (!paramBounds(c(pars,fix))) {
    return(1000000)
  }
  preds <- fcn(free=pars,fixed=fix,return_dist = TRUE)
  err <- FAM:::LL(obs=obs[,c("class","order","RTrounded")],pred = preds$distribution)  
  return(err)
}

fitCFR_PCL<- function(model,inpar = FALSE,...,debugLevel = 0) {
  library(optimx)
  library(dplyr)
  library(reshape2)
  library(whoppeR)
  library(tidyr)
  
  checkRequiredParams(c(model$free,model$fixed), model$fn)
  
  if (inpar) {
    cluster <- tryDoParallelCluster(detectCores()-1)
  }
  
  if (inpar && cluster$success) {
    logfile <- tempfile("parlog", fileext = ".txt")
    writeLines(paste('[',Sys.time(),']',"INIT parlog"),
               con=logfile,sep='\n')
    clusterExport(cluster$handle,"logfile",envir = environment())
    model$results <- foreach(j =unique(model$data$subject), .verbose=T,
                             .packages=c("optimx","PCL","whoppeR","dplyr")) %dopar% {
      sink(logfile, append=TRUE)
      cat(paste("Fitting subject", j,"\n"))
      sink()
      data <- model$data[model$data$subject == j,]
      freePars <- names(model$free)
      interim<- optimx(par=model$free[c("ER","LR")],
                       fn = accuracyErrorFcn,
                       method = model$method,
                       control = list(maxit=1),
                       fcn = model$fn,
                       fix=c(model$fixed,
                             model$free[!freePars %in% c("ER","LR")]),
                       obs=data)
      newStart <- c(unlist(interim[,c("ER","LR")]), 
                    model$free[!freePars %in% c("ER","LR")])
      fit <- optimx(par = newStart,
                    fn = RT_ErrorFcn,
                    method = m$method,
                    control = list(maxit=1),
                    fcn = model$fn,
                    fix = model$fixed,
                    obs = data)      
      best <- as.vector(coef(fit))
      names(best) <- colnames(coef(fit))
      preds <- model$fn(free = best,fixed = model$fixed)
      end <- list(fit=fit, preds = preds)
    }
    
    stopCluster(cluster$handle)
    cat(paste('[',Sys.time(),']',"Finished Fitting, Goodbye"),
        file=logfile,sep='\n')
    
  } else {
    
    model$results <- vector(mode = 'list', length = length(unique(model$data$subject)))
    for (j in unique(model$data$subject)) {
      message(paste("Fitting subject", j))
      data <- model$data[model$data$subject == j,]
      freePars <- names(model$free)
      interim<- optimx(par=model$free[c("ER","LR")],
                       fn = accuracyErrorFcn,
                       method = model$method,
                       control = list(maxit=1),
                       fcn = model$fn,
                       fix=c(model$fixed,
                             model$free[!freePars %in% c("ER","LR")]),
                       obs=data)
      newStart <- c(unlist(interim[,c("ER","LR")]), 
                    model$free[!freePars %in% c("ER","LR")])
      fit <- optimx(par = newStart,
                    fn = RT_ErrorFcn,
                    method = model$method,
                    control = list(maxit=1),
                    fcn = model$fn,
                    fix = model$fixed,
                    obs = data)      
      best <- as.vector(coef(fit))
      names(best) <- colnames(coef(fit))
      preds <- model$fn(free = best,fixed = model$fixed)
      model$results[[j]] <- list(fit=fit, preds = preds )
    }
  }
  
  return(model)
}
