fitCFR_PCL<- function(model,inpar = FALSE,...,debugLevel = 0) {
  library(optimx)
  library(dplyr)
  library(reshape2)
  if (inpar) {
    packs <- requireNamespace(c('foreach', 'doParallel'), quietly = TRUE)
    if (all(packs)) {
      library(foreach)
      library(doParallel)
      cl <- makeCluster(detectCores()-1)
      registerDoParallel(cl)
      #     library(doRNG)
      #     registerDoRNG(456)
    } else {
      warning("Packages foreach and doParallel not found to do parallel processing, falling back")
      inpar = FALSE
    }
  }
  
  checkRequiredParams(c(model$par,model$fixed), model$fn)
  
  if (inpar) {
    logfile <- tempfile("parlog", fileext = ".txt")
    writeLines(paste('[',Sys.time(),']',"INIT parlog"),
               con=logfile,sep='\n')
    clusterExport(cl,"logfile",envir = environment())
    model$results <- foreach(j =unique(model$data$subject), .verbose=T,
                             .packages=c("optimx","PCL","whoppeR")) %dopar% {
      sink(logfile, append=TRUE)
      cat(paste("Fitting subject", j,"\n"))
      sink()
      m  <- model
      m$data <- model$data[model$data$subject == j,]
      interim<- optimx(par=m$par[c("ER","LR")],fn=m$fn, method = m$method,
                       itnmax = m$itnmax, data = m$data,
                       fixed=c(m$fixed, m$par[!names(m$par) %in% c("ER","LR")]),
                       fitAcc = TRUE,  fitRT= FALSE, fitting= m$fitting)
      newStart <- c(unlist(interim[,c("ER","LR")]), 
                    m$par[!(names(m$par) %in% c("ER","LR"))])
      fit <- optimx(par=newStart, fn=m$fn, method = m$method,
                    itnmax = m$itnmax, fixed=m$fixed, data = m$data,
                    fitting= m$fitting)      
      best <- as.vector(coef(fit))
      names(best) <- colnames(coef(fit))
      preds <- m$fn(free=best,fixed=m$fixed,data=m$data,
                    fitting=FALSE)[c("preds","dist")]
      end <- list(fit=fit, preds = preds)
      
      }
    
    stopCluster(cl)
    cat(paste('[',Sys.time(),']',"Finished Fitting, Goodbye"),
        con=logfile,sep='\n')
    
  } else {
    
    model$results <- vector(mode = 'list', length = length(unique(model$IVdata$subject)))
    for (j in 1:2){ #unique(model$data$subject)) {
      message(paste("Fitting subject", j))
      m  <- model
      m$data <- model$data[model$data$subject == j,]
      interim<- optimx(par=m$par[c("ER","LR")],fn=m$fn, method = m$method,
                       itnmax = m$itnmax, data = m$data,
                       fixed=c(m$fixed, m$par[!names(m$par) %in% c("ER","LR")]),
                       fitAcc = TRUE,  fitRT= FALSE, fitting= m$fitting)
      newStart <- c(unlist(interim[,c("ER","LR")]), 
                    m$par[!(names(m$par) %in% c("ER","LR"))])
      fit <- optimx(par=newStart, fn=m$fn, method = m$method,
             itnmax = m$itnmax, fixed=m$fixed, data = m$data,
             fitting= m$fitting)      
      best <- as.vector(coef(fit))
      names(best) <- colnames(coef(fit))
      preds <- m$fn(free=best,fixed=m$fixed,data=m$data,
                    fitting=FALSE)[c("preds","dist")]
      model$results[[j]] <- list(fit=fit, preds = preds )
    }
  }


  return(model)
}
