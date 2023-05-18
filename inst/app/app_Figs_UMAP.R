
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
