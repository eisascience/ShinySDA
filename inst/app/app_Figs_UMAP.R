
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
