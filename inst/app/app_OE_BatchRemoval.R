observeEvent(input$nextSDA_br, {
  SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
  
  if(which(SDAorder==envv$QC_compIter) < length(SDAorder)){
    envv$QC_compIter = SDAorder[which(SDAorder==envv$QC_compIter) + 1]
  }
  
  
  
  
})

observeEvent(input$prevSDA_br, {
  
  SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
  
  if(which(SDAorder==envv$QC_compIter) < length(SDAorder)){
    envv$QC_compIter = SDAorder[which(SDAorder==envv$QC_compIter) - 1]
  }
  
  
})

observeEvent(input$SDAScoresChi_clus, {
  
  if(is.null(envv$SDAScoresChi_clusBTN)) {
    envv$SDAScoresChi_clusBTN = "ON"
  } else if( envv$SDAScoresChi_clusBTN == "OFF"){
    envv$SDAScoresChi_clusBTN = "ON"
  } else if(envv$SDAScoresChi_clusBTN == "ON"){
    envv$SDAScoresChi_clusBTN = "OFF"
  }
  
  
})



observeEvent(input$CompBatchCheckBoxSelect, {
  
  envv$Remove_comps <- input$CompBatchCheckBoxSelect
  
  
})

observeEvent(input$run_tSNE_CS_batch, {
  
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    # envv$Remove_comps = c(13, 40)
    SDAres <- envv$SDAres
    if(length(envv$Remove_comps) > 50) {
      # suffix <- paste(envv$Remove_comps, collapse = "")
      
      suffix <- paste0("LargeSetOfComps_", length(envv$Remove_comps))
      
    } else if(length(envv$Remove_comps) > 30 & length(envv$Remove_comps) <= 50 ) {
      
      suffix <- paste0(findIntRuns(as.numeric(unlist(strsplit(envv$Remove_comps, ",")))), collapse="")
      
    } else if (length(envv$Remove_comps) <= 30) {
      suffix <- paste(envv$Remove_comps, collapse = "_")
    }
    
    tSNE_n.iter <- as.numeric(input$tSNEiter) # tSNE_n.iter = 1000
    tSNE_pp <- as.numeric(input$tSNEpp) # tSNE_pp = 50
    
    choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
    selected <- setdiff(choice, envv$Remove_comps)
    
    
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
    
    
    if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tSNE_pp,"", suffix, ".rds"))){
      
      envv$InfoBox_sub = "Starting tSNE with cell scores - batch-removal.. wait"
      tsne_CS_batch <- Rtsne::Rtsne(SDAres$scores[,selected], verbose=TRUE, pca=FALSE, 
                                    perplexity = tSNE_pp, 
                                    max_iter=tSNE_n.iter, num_threads = 8, check_duplicates = F)
      
      envv$InfoBox_sub = "Saving tSNE with cell scores - batch-removal.. wait"
      saveRDS(tsne_CS_batch, file=paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tSNE_pp,"", suffix, ".rds"))
      
      
      
    } else {
      tsne_CS_batch <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tSNE_pp,"", suffix, ".rds"))
    }
    envv$InfoBox_sub = "batch-removal tSNE with cell scores complete"
    envv$tsne_CS_batch <- tsne_CS_batch
    
  }

  
})

observeEvent(input$save_batch_selection, {
  
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
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
    
    
    # choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
    selected <- envv$Remove_comps
    # print(head(selected))
    saveRDS(selected, file=paste0(envv$path2SDA_dyn, "/", head.path,"_BatchSelectedComps", ".rds"))
    
    
  }
})

observeEvent(input$load_batch_selection, {
  
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
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
    
    
    if(file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_BatchSelectedComps", ".rds"))){
      selected <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_BatchSelectedComps", ".rds"))
    } else {
      choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) 
      selected <- setdiff(choice, envv$QC_components)
      
    }
    
    
    # print(head(selected))
    envv$Remove_comps <- selected
  }
})

observeEvent(input$reset_batch_selection, {
  choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) 
  selected <- setdiff(choice, envv$QC_components)
  
  envv$Remove_comps <- selected
  
})
observeEvent(input$select_all_selection, {
  envv$Remove_comps <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
  
})