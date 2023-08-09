
## UMAP DGE---------

output$DGE_SDA_UMAP <- renderPlot({
  
  if(is.null(envv$UMAP_CS_batch)){
    plot(x=0, y=0, main="UMAP CS Batch not found")
    
  } else {
    
    UMAPDF <- as.data.frame(envv$UMAP_CS_batch)
    rownames(UMAPDF)  <- rownames(envv$SDAres$scores)
    colnames(UMAPDF) <- c("UMAP1_batch", "UMAP2_batch")
    
    
    
    
    if(is.null(envv$MetaDF)){
      UMAPDF$SumScore <- rowSums(abs(envv$SDAres$scores))
      UMAPDF$SumScore <- UMAPDF$SumScore/mean(UMAPDF$SumScore)
      
      ggplot(UMAPDF, aes(UMAP1_batch, UMAP2_batch, color=(SumScore))) +
        geom_point(size = 1) + theme_bw() +
        scale_color_distiller(palette = "Spectral")  +
        ggtitle("UMAP SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n ")+
        theme(legend.position = "bottom", aspect.ratio=1)
      
    } else {
      MetaDF <- envv$MetaDF
      if("library_size" %in% colnames(MetaDF)){
        UMAPDF$library_size <- MetaDF[rownames(UMAPDF), ]$library_size
        
        ggplot(UMAPDF, aes(UMAP1_batch, UMAP2_batch, color=log10(library_size))) +
          geom_point(size = 1) + theme_bw() +
          scale_color_distiller(palette = "Spectral")  +
          ggtitle("UMAP SDA batch removed\n log10 library size \n ")+
          theme(legend.position = "bottom", aspect.ratio=1)
      } else {
        
        UMAPDF$SumScore <- rowSums(abs(envv$SDAres$scores))
        UMAPDF$SumScore <- UMAPDF$SumScore/mean(UMAPDF$SumScore)
        
        ggplot(UMAPDF, aes(UMAP1_batch, UMAP2_batch, color=(SumScore))) +
          geom_point(size = 1) + theme_bw() +
          scale_color_distiller(palette = "Spectral")  +
          ggtitle("UMAP SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n ")+
          theme(legend.position = "bottom", aspect.ratio=1)
        
      }
    }
  }
  
})





output$SDAumap_br2 <- renderPlot({
  
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    zN = envv$QC_compIter
    
    SDAres <- envv$SDAres
    # UMAPDF <- as.data.frame(envv$tsne_CS_batch$Y)
    # rownames(UMAPDF)  <- rownames(envv$SDAres$scores)
    # colnames(UMAPDF) <- c("tSNE1_batch", "tSNE2_batch")
    UMAPDF <- as.data.frame(envv$UMAP_CS_batch)
    rownames(UMAPDF)  <- rownames(envv$SDAres$scores)
    colnames(UMAPDF) <- c("UMAP1_batch", "UMAP2_batch")
    
    
    if(zN %in% envv$Remove_comps) RemoveTag = "removed" else RemoveTag = "kept"
    
    
    UMAPDF$SDAComp <- cut(asinh(SDAres$scores[,paste0("SDAV", zN, sep="")]^3), 
                           breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf))
    
    UMAPDF$SDAComp <- factor(UMAPDF$SDAComp, 
                              levels = c("(-Inf,-1]", "(-1,-0.5]", "(-0.5,0]",  "(0,0.5]",   "(0.5,1]",   "(1, Inf]" ) )
    # print(UMAPDF$SDAComp)
    # print(factor(UMAPDF$SDAComp))
    # print(factor(UMAPDF$SDAComp, 
    #              levels = c("(-Inf,-1]", "(-1,-0.5]", "(-0.5,0]",  "(0,0.5]",   "(0.5,1]",   "(1, Inf]" ) ))
    
    ggplot(UMAPDF, aes(UMAP1_batch, UMAP2_batch,  color=UMAPDF$SDAComp)) +
      geom_point(size = 1) + theme_bw() +
      scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
      guides(colour = guide_legend(override.aes = list(size = 2, alpha=1))) +
      theme(legend.position = "bottom", aspect.ratio=1) + 
      simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE) +
      ggtitle(paste0("SDAV", zN, " :: ", RemoveTag))+
      ylab("asinh(SDAscore^3)")
    
    
    
    
    
  }
})



output$GeneExprSDAUMAP <- renderPlot({
  
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
    
    UMAPDF <- as.data.frame(envv$UMAP_CS_batch)
    rownames(UMAPDF)  <- rownames(envv$SDAres$scores)
    colnames(UMAPDF) <- c("UMAP1_batch", "UMAP2_batch")
    
    UMAPDF$GeneExpr <- rep(0, nrow(UMAPDF))
    
    
    SDAres <- envv$SDAres
    
    GeneSet <- GeneSet[GeneSet %in% colnames(SDAres$loadings[[1]])]
    
    if(length(GeneSet)>1){
      
      GeneExpr <- SDAres$scores %*% SDAres$loadings[[1]][,as.character(GeneSet)]
      GeneExpr <- as.data.frame(rowSums(GeneExpr))
      
      TitleX = paste0("Sum-Expr of :", paste(GeneSet, collapse = "_") )
      
      if(is.null(rownames(SDAres$loadings[[1]]))){
        rownames(SDAres$loadings[[1]]) = paste0("C", 1:nrow(SDAres$loadings[[1]]))
      }
      
    
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
    
    
    
    UMAPDF[rownames(GeneExpr), ]$GeneExpr <- GeneExpr[,1]
    
    
    # UMAPDF <- (envv$UMAPGEx_br)
    # print(head(UMAPDF))
    # TitleX <- envv$UMAPGEx_tit
    
    
    ggplot(UMAPDF, aes(UMAP1_batch, UMAP2_batch,  color=cut(asinh(GeneExpr^3),
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
