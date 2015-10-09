fitSCKT_PCR<- function(model,inpar = FALSE,...,debugLevel = 0) {
  library(optimx)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(whoppeR)
  library(FAM)
  
  errorFcn <- function(pars, fcn, obs, N, fix) {
    if (!paramBounds(c(pars,fix))) {
      return(10000000)
    }
    
    preds <- fcn(free=pars,fixed=fix)
    preds <- unlist(preds[names(obs)])
    err <- binomialLL(obs=obs[1:16],pred= preds[1:16], N=4 ) + #N= N[1:16]
      multinomialLL(obs = obs[17:48],pred = preds[17:48],
                    N = 4)
    return(err)
  }
  if (inpar) {
    cluster <- tryDoParallelCluster(detectCores()-1)
  }

  if (inpar && cluster$success) {
    logfile <- tempfile("parlog", fileext = ".txt")
    writeLines(paste('[',Sys.time(),']',"INIT parlog"),
               con=logfile,sep='\n')
    clusterExport(cluster$handle,"logfile",envir = environment())
    model$results <- foreach(j =unique(model$IVdata$subject), .verbose=T,
                             .packages=c("optimx","PCR","whoppeR")) %dopar% {
      sink(logfile, append=TRUE)
      cat(paste("Fitting subject", j,"\n"))
      sink()
      IVdata <- model$IVdata[model$IVdata$subject == j,]
      IVobs <- IVdata$final_score[is.na(IVdata$prac_score)]
      IVn <-  IVdata$n[is.na(IVdata$prac_score)]
      names(IVn) <- names(IVobs)
      names(IVobs) <- IVdata$name[is.na(IVdata$prac_score)]
      
      # Get joint data and counts
      jointData <- model$jointData[model$jointData$subject == j,]
      jointObs <- jointData$cell_percent[jointData$practice %in% 'T']
      names(jointObs) <- jointData$name[jointData$practice %in% 'T']
      jointN <- jointData$maxCount[jointData$practice %in% 'T']
      names(jointN) <- names(jointObs)
      fit <- optimx(par=model$free,fn = errorFcn,
                    method = model$method,  
                    fcn = model$fn,
                    obs = c(IVobs, jointObs), N = c(IVn, jointN),
                    fix = model$fixed,
                    control = list(maxit = 1000))
      best <- as.vector(coef(fit))
      names(best) <- colnames(coef(fit))
      preds <- model$fn(free = best,fixed = model$fixed)
      end <- list(fit=fit, preds = preds)
    }
    
    stopCluster(cluster$handle)
    cat(paste('[',Sys.time(),']',"Finished Fitting, Goodbye"),
        file=logfile,sep='\n')
    
  } else {
    
    model$results <- vector(mode = 'list', length = length(unique(model$IVdata$subject)))
    for (j in unique(model$IVdata$subject)) {
      message(paste("Fitting subject", j))
      # Get IV data and counts
      IVdata <- model$IVdata[model$IVdata$subject == j,]
      IVobs <- IVdata$final_score[is.na(IVdata$prac_score)]
      IVn <-  IVdata$n[is.na(IVdata$prac_score)]
      names(IVn) <- names(IVobs)
      names(IVobs) <- IVdata$name[is.na(IVdata$prac_score)]
      
      # Get joint data and counts
      jointData <- model$jointData[model$jointData$subject == j,]
      jointObs <- jointData$cell_percent[jointData$practice %in% 'T']
      names(jointObs) <- jointData$name[jointData$practice %in% 'T']
      jointN <- jointData$maxCount[jointData$practice %in% 'T']
      names(jointN) <- names(jointObs)
      fit <- optimx(par=model$free,fn = errorFcn,
                    method = model$method,  
                    fcn = model$fn,
                    obs = c(IVobs, jointObs), N = c(IVn, jointN),
                    fix = model$fixed,
                    control = list(maxit = 1000))
      best <- as.vector(coef(fit))
      names(best) <- colnames(coef(fit))
      preds <- model$fn(free = best,fixed = model$fixed)
      model$results[[as.numeric(j)]] <- list(fit=fit, preds = preds )
    }
  }
  
  return(model)
}
