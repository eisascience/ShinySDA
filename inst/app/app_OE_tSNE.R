## run raw tSNE
observeEvent(input$runtSNE, {
  
  envv$InfoBox_sub = "Starting tSNE with cell scores - all comps.. wait"
  
  envv = ShinySDA:::Run_tSNE_full_evv(envv, input = input)
  
  
})

##Qc tsne
observeEvent(input$runtSNEQCfilt, {
  
  envv$InfoBox_sub = "Starting tSNE with cell scores - qc comps.. wait"
  
  envv = ShinySDA:::Run_tSNE_QC_envv(envv, input = input)
  
  
  
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



## Gene Expr SDA tSNE ------

output$GeneExprSDAtSNE <- renderPlot({
  
  GeneSet <- input$GeneSet
  
  
  if(length(grep(",", GeneSet)) == 0){
    
    if(length(grep('"', GeneSet)) + length(grep("'", GeneSet))>0) {
      GeneSet <- unlist(strsplit(gsub("'", '', gsub('"', '', GeneSet)), " "))
    } else {
      GeneSet <- unlist(strsplit(GeneSet, " "))
    }
    
    
  } else {
    GeneSet <- (unlist(strsplit(gsub(" ", "", gsub("'", '', gsub('"', '', GeneSet))), ",")))
    
  }
  
  
  
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    tempDFX <- as.data.frame(envv$tsne_CS_batch$Y)
    rownames(tempDFX)  <- rownames(envv$SDAres$scores)
    colnames(tempDFX) <- c("tSNE1_batch", "tSNE2_batch")
    
    tempDFX$GeneExpr <- rep(0, nrow(tempDFX))
    
    
    SDAres <- envv$SDAres
    
    GeneSet <- GeneSet[GeneSet %in% colnames(SDAres$loadings[[1]])]
    
    if(length(GeneSet)>1){
      
      GeneExpr <- SDAres$scores %*% SDAres$loadings[[1]][,as.character(GeneSet)]
      GeneExpr <- as.data.frame(rowSums(GeneExpr))
      
      TitleX = paste0("Sum-Expr of :", paste(GeneSet, collapse = "_") )
      
      LoadOrdVal <- round(SDAres$loadings[[1]][,as.character(GeneSet[1])][order(abs(SDAres$loadings[[1]][,as.character(GeneSet[1])]), decreasing = T)], 3)
      
      
      
    } else if(length(GeneSet)==1){
      GeneExpr <- SDAres$scores %*% SDAres$loadings[[1]][,as.character(GeneSet)]
      TitleX = paste0("Expr of :", GeneSet )
      LoadOrdVal <- round(SDAres$loadings[[1]][,as.character(GeneSet)][order(abs(SDAres$loadings[[1]][,as.character(GeneSet)]), decreasing = T)], 3)
      
    } else if(!length(GeneSet)>=1)  {
      GeneExpr <- SDAres$scores %*% rep(0, nrow(SDAres$loadings[[1]]))
      TitleX = "No genes in input"
      LoadOrdVal = paste("g",1:20)
    }
    
    
    
    tempDFX[rownames(GeneExpr), ]$GeneExpr <- GeneExpr[,1]
    
    
    # tempDFX <- (envv$tSNEGEx_br)
    # print(head(tempDFX))
    # TitleX <- envv$tSNEGEx_tit
    
    
    ggplot(tempDFX, aes(tSNE1_batch, tSNE2_batch,  color=cut(asinh(GeneExpr^3),
                                                             breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
      geom_point(size = 1) + theme_bw() +
      scale_color_manual("Expr", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) +
      guides(colour = guide_legend(override.aes = list(size = 2, alpha=1))) +
      theme(legend.position = "bottom", aspect.ratio=1) +
      simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE)  +
      labs(title = paste0("SDA-Batch-removed DGE\n", TitleX), 
           subtitle = paste("Found in comps: \n",
                            paste(names(LoadOrdVal)[1:5], collapse = ", "), 
                            "\n",
                            paste(LoadOrdVal[1:5], collapse = ", "), 
                            "\n",
                            paste(names(LoadOrdVal)[6:10], collapse = ", "), 
                            "\n",
                            paste(LoadOrdVal[6:10], collapse = ", "), 
                            "\n"), 
           caption = "Caption here") + 
      ylab("asinh(GeneExpr^3)")
    # ggtitle(paste0("SDA-Batch-removed DGE\n", TitleX))+
    # ylab("asinh(GeneExpr^3)")
    
    
    
    
  }
  
  
  
  
})
