accuracyErrorFcn <- function(pars, fcn, fix, obs) {
  if (!paramBounds(c(pars,fix)) | anyNA(c(pars,fix))) {
    return(1000000)
  }
  preds <- fcn(free=pars,fixed=fix,return_dist = FALSE)
  obs <- obs %>% group_by(class,order) %>%
    summarise(acc = sum(score)/4) %>%
    ungroup() %>%
    complete(class,order,fill=list(acc=0))
  err <- sum((preds$acc-obs$acc)^2)
  return(err)
}

RT_ErrorFcn <- function(pars,fcn, fix, obs) {
  if (!paramBounds(c(pars,fix)) | anyNA(c(pars,fix))) {
    return(1000000)
  }
  preds <- fcn(free=pars,fixed=fix,return_dist = TRUE)
  likelihoods <- inner_join(obs[,c("class","order","RTrounded")],
                            preds$distribution,
                            by = c("class", "order", "RTrounded"))
#   likelihoods$RTdist[likelihoods$RTdist == 0] <- .Machine$double.xmin
  err <- -sum(log(likelihoods$RTdist))
  return(err)
}

fitCFR_PCL<- function(model,inpar = FALSE,...,debugLevel = 0) {
  library(optimx)
  library(dplyr)
  library(reshape2)
  library(whoppeR)
  library(tidyr)

  if (inpar) {
    cluster <- tryDoParallelCluster(detectCores()-1)
  }
  
  routine <- function() {
    data <- model$data[model$data$subject == j,]
    freePars <- names(model$free[[j]])
    interim<- optimx(par=model$free[[j]][c("ER","LR","alpha")],
                     fn = accuracyErrorFcn,
                     method = model$method,
                     control = list(maxit=1000,
                                    parscale = c(1,1,model$free[[j]]['alpha'])),
                     fcn = model$fn,
                     fix=c(model$fixed,
                           model$free[[j]][!freePars %in% c("ER","LR","alpha")]),
                     obs=data)
    newStart <- c(unlist(interim[,c("ER","LR","alpha")]), 
                  model$free[[j]][!freePars %in% c("ER","LR","alpha")])
    fit <- optimx(par = newStart,
                  fn = RT_ErrorFcn,
                  method = model$method,
                  control = list(maxit=1000,
                                 parscale = c(1,1,1,1,model$free[[j]]['alpha'],
                                              1,1,model$free[[j]]['Tmax'])),
                  fcn = model$fn,
                  fix = model$fixed,
                  obs = data)      
    best <- as.vector(coef(fit))
    names(best) <- colnames(coef(fit))
    preds <- model$fn(free = best,fixed = model$fixed)
    results <- list(fit=fit, preds = preds)
    return(results)
  }
  
  if (inpar && cluster$success) {
    logfile <- tempfile("parlog", fileext = ".txt")
    writeLines(paste('[',Sys.time(),']',"INIT parlog"),
               con=logfile,sep='\n')
    clusterExport(cluster$handle,c("accuracyErrorFcn", "RT_ErrorFcn"))
    model$results <- foreach(j =unique(model$data$subject), .verbose=T,
                             .packages=c("optimx","PCL","whoppeR",
                                         "dplyr","FAM","reshape2","tidyr")) %dopar% {
      sink(logfile, append=TRUE)
      cat(paste("Fitting subject", j,"\n"))
      results <- routine()
      sink()
      list(data.frame(subject = j, fit=results$fit),
           preds = data.frame(subject = j, results$preds$preds), 
           distribution = data.frame(subject = j,
                                     results$preds$distributions))
    }
    
    stopCluster(cluster$handle)
    cat(paste('[',Sys.time(),']',"Finished Fitting, Goodbye"),
        file=logfile,sep='\n')
    
  } else {
    
    model$results <- vector(mode = 'list',
                            length = length(unique(model$data$subject)))
    for (j in unique(model$data$subject)) {
      message(paste("Fitting subject", j))
      results <-  routine()
      model$results[[j]] <- list(data.frame(subject = j, fit=results$fit),
                                 preds = data.frame(subject = j, results$preds$preds), 
                                 distribution = data.frame(subject = j,
                                                           results$preds$distribution))
    }
  }
  
  return(model)
}
