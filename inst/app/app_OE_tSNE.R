## run raw tSNE
observeEvent(input$runtSNE, {
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
    envv$InfoBox_sub = "Checking paths for tSNE"
    
    if(envv$Origin == "folder"){
      head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
      base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
      
      head.path <- gsub("/", "", head.path)
      
    }
    if(envv$Origin == "prime"){
      head.path <- "sda_results"
      base.path <- ""
      #TODO
    }
    
    SDAres <- envv$SDAres
    
    # length(which(duplicated(SDAres$scores[,])))
    
    if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_AllComps_pp50.rds"))){
      
      # print(nrow(SDAres$scores[,]))
      tsnepp <- round((nrow(SDAres$scores[,]) - 1)/3.2)
      tsnepp <- ifelse(tsnepp < 50, tsnepp, 50)
      print("perplexity:")
      print(tsnepp)
      
      
      
      envv$InfoBox_sub = "Starting tSNE with cell scores - all comps.. wait"
      tsne_CS_all <- Rtsne::Rtsne(SDAres$scores[,], verbose=TRUE, pca=FALSE, 
                                  perplexity = tsnepp, 
                                  max_iter=1000, 
                                  num_threads = 8, 
                                  check_duplicates = F)
      
      
      envv$InfoBox_sub = "Saving tSNE with cell scores - all comps.. wait"
      saveRDS(tsne_CS_all, file=paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_AllComps_pp50.rds"))
      
      envv$InfoBox_sub = paste0("tSNE complete; pp=", tsnepp)
      if(length(which(duplicated(SDAres$scores[,])))>0) envv$InfoBox_sub = paste0("tSNE complete : dups=", length(which(duplicated(SDAres$scores[,]))))
      
      
    } else {
      tsne_CS_all <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_AllComps_pp50.rds"))
    }
    envv$InfoBox_sub = "raw tSNE with cell scores complete"
    envv$tsne_CS_all <- tsne_CS_all
    
  }
  
  
  
})

##Qc tsne
observeEvent(input$runtSNEQCfilt, {
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    SDAres <- envv$SDAres
    
    suffix <- paste(setdiff(1:as.numeric(SDAres$command_arguments$num_comps), envv$QC_components), collapse = "_")
    
    envv$InfoBox_sub = "Checking paths for tSNE"
    
    if(envv$Origin == "folder"){
      head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
      base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
      
      head.path <- gsub("/", "", head.path)
      
    }
    if(envv$Origin == "prime"){
      head.path <- "sda_results"
      base.path <- ""
      #TODO
    }
    
    tsnepp <- round((nrow(SDAres$scores[,]) - 1)/3.2)
    tsnepp <- ifelse(tsnepp < 50, tsnepp, 50)
    
    tSNE_n.iter = 1000
    # tsnepp
    
    if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tsnepp,"", suffix, ".rds"))){
      
      
      
      
      
      
      envv$InfoBox_sub = "Starting tSNE with cell scores - qc comps.. wait"
      tsne_CS_qc <- Rtsne::Rtsne(SDAres$scores[,envv$QC_components], verbose=TRUE, pca=FALSE, 
                                 perplexity = tsnepp, 
                                 max_iter=tSNE_n.iter, num_threads = 8, check_duplicates = F)
      
      envv$InfoBox_sub = "Saving tSNE with cell scores - qc comps.. wait"
      saveRDS(tsne_CS_qc, file=paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tsnepp,"", suffix, ".rds"))
      
      
      
    } else {
      tsne_CS_qc <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tsnepp,"", suffix, ".rds"))
    }
    envv$InfoBox_sub = "qc tSNE with cell scores complete"
    envv$tsne_CS_qc <- tsne_CS_qc
    
  }
  
  
  
  
})