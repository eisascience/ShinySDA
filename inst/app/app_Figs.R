
# Main Figs -------

## tSNE CS all---------

output$tSNE_CS_all <- renderPlot({
  if(is.null(envv$tsne_CS_all)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    tsneDF <- as.data.frame(envv$tsne_CS_all$Y)
    rownames(tsneDF)  <- rownames(envv$SDAres$scores)
    colnames(tsneDF) <- c("tSNE1_all", "tSNE2_all")
    tsneDF$SumScore <- rowSums(abs(envv$SDAres$scores))
    tsneDF$SumScore <- tsneDF$SumScore/mean(tsneDF$SumScore)
    
    
    ggplot(tsneDF, aes(tSNE1_all, tSNE2_all, color=(SumScore))) +
      geom_point(size = 1) + theme_bw() + 
      scale_color_distiller(palette = "Spectral") +
      theme(legend.position = "bottom") +
      ggtitle("tSNE SDA qc Components\n Sum absolute-cell-scores normalized by its mean \n ")
    
  }
})

## tSNE CS QC---------

output$tSNE_CS_qc <- renderPlot({
  if(is.null(envv$tsne_CS_qc)){
    plot(x=0, y=0, main="Load an SDA")
    
  } else {
    
    tsneDF <- as.data.frame(envv$tsne_CS_qc$Y)
    rownames(tsneDF)  <- rownames(envv$SDAres$scores)
    colnames(tsneDF) <- c("tSNE1_qc", "tSNE2_qc")
    
    # print(head(tsneDF))
    
    
    
    if(is.null(envv$MetaDF)){
      print("MetaDF is NULL")
      tsneDF$SumScore <- rowSums(abs(envv$SDAres$scores))
      tsneDF$SumScore <- tsneDF$SumScore/mean(tsneDF$SumScore)
      
      ggplot(tsneDF, aes(tSNE1_qc, tSNE2_qc, color=(SumScore))) +
        geom_point(size = 1) + theme_bw() + 
        scale_color_distiller(palette = "Spectral") +
        theme(legend.position = "bottom") +
        ggtitle("tSNE SDA qc Components\n  Sum absolute-cell-scores normalized by its mean \n ")
      
    } else {
      MetaDF <- envv$MetaDF
      tsneDF$library_size <- MetaDF[rownames(tsneDF), ]$library_size
      
      ggplot(tsneDF, aes(tSNE1_qc, tSNE2_qc, color=log10(library_size))) +
        geom_point(size = 1) + theme_bw() + 
        scale_color_distiller(palette = "Spectral") +
        theme(legend.position = "bottom") +
        ggtitle("tSNE SDA qc Components\n log10 library size \n ")
    }
    
    
    
  }
})


## SDAScoresACross---------

output$SDAScoresAcross <- renderPlot({
  
  # ColFac_DONR.ID <- CDID()
  
  if(is.null(envv$MetaDF)){
    print("No Comp")
  } else {
    SDAScores <- envv$SDAres$scores
    # ComponentN <- as.numeric(envv$QC_compIter)
    MetaDF <- envv$MetaDF
    MetaDF <- MetaDF[rownames(SDAScores),]
    
    
    # colnames(MetaDF)
    
    
    
    ShinySDA:::Plot_CS_AcrossMeta(SDAScores = SDAScores, col_vector = col_vector,
                       MetaDF = MetaDF, MetaSelect= input$Metaselect6)
    
    
    
    
  }
  
})

output$SDAScoresAcross_Legend <- renderPlot({
  
  # ColFac_DONR.ID <- CDID()
  
  if(is.null(envv$MetaDF)){
    print("No Comp")
  } else {
    SDAScores <- envv$SDAres$scores
    # ComponentN <- as.numeric(envv$QC_compIter)
    MetaDF <- envv$MetaDF
    MetaDF <- MetaDF[rownames(SDAScores),]
    
    
    # colnames(MetaDF)
    
    
    
    ShinySDA:::Plot_CS_AcrossMeta_Legend(SDAScores = SDAScores, col_vector = col_vector,
                                  MetaDF = MetaDF, MetaSelect= input$Metaselect6)
    
    
    
    
  }
  
})



## SDAScoresBoxplot---------

output$SDAScoresBoxplot <- renderPlot({
  
  # ColFac_DONR.ID <- CDID()
  
  if(is.null(envv$MetaDF)){
    print("No Comp")
  } else {
    SDAScores <- envv$SDAres$scores
    # ComponentN <- as.numeric(envv$QC_compIter)
    MetaDF <- envv$MetaDF
    MetaDF <- MetaDF[rownames(SDAScores),]
    
    
    # colnames(MetaDF)
    
    
    
    ShinySDA::Plot_CS_BoxplotMeta(SDAScores = SDAScores, col_vector = col_vector,
                        MetaDF = MetaDF, MetaSelect= input$Metaselect5)
    
    
    
    
  }
  
})


# output$SDAScoresAcross2 <- renderPlot({
#   
#   # ColFac_DONR.ID <- CDID()
#   
#   if(is.null(envv$MetaDF)){
#     print("No Comp")
#   } else {
#     SDAScores <- envv$SDAres$scores
#     ComponentN <- as.numeric(envv$QC_compIter)
#     MetaDF <- envv$MetaDF
#     MetaDF <- MetaDF[rownames(SDAScores),]
#     
#     
#     # colnames(MetaDF)
#     
#     
#     
#     
#     tempDF <- data.frame(cell_index = 1:nrow(SDAScores), 
#                          score = asinh((SDAScores[, paste0("SDAV", ComponentN)])^3), 
#                          experiment = MetaDF$Population, 
#                          ColFac = MetaDF$SubjectId)
#     
#     tempDF <- tempDF[order(tempDF$score),]
#     tempDF$cell_index <- 1:nrow(tempDF)
#     
#     # tempDF$cell_index_cut <- cut(tempDF$cell_index, quantile(tempDF$cell_index))
#     
#     # print(levels(tempDF$cell_index_cut))
#     
#     ggplot(tempDF, 
#            aes(cell_index, score, colour = ColFac)) + 
#       geom_jitter(size = 1, width=0, height = 3, alpha = .6) +
#       # geom_boxplot(aes(x= factor(cut(cell_index, quantile(tempDF$cell_index))), y=score, 
#       # colour = ColFac), outlier.colour = "red", outlier.shape = 8) +
#       # geom_point(size = 0.5, stroke = 0) + 
#       xlab("Cell Index") + ylab("asinh(Score^3)") + 
#       #scale_color_brewer(palette = "Paired") + 
#       theme_bw() + 
#       theme(legend.position = "none") + 
#       guides(colour = guide_legend(ncol = 4, override.aes = list(size = 2, alpha=1))) +
#       scale_colour_manual(values =(col_vector),
#                           guide = guide_legend(nrow=2)) +
#       # guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
#       ggtitle(paste0("SDAV", ComponentN)) + 
#       geom_smooth(method = "lm", formula = y ~ x, size = 2, colour="red")+ 
#       geom_smooth(method = "loess", formula = y ~ x, size = 2, colour="dodgerblue")  +
#       facet_wrap(~experiment, ncol = 2)
#     
#     
#     
#   }
#   
# })



output$GL_Cor_HMplot <- renderPlot({
  
  ShinySDA::Plot_CorSDA_Loadings(SDAres = envv$SDAres)
  
})


# Batch removed Figs ----

## tSNE batch removed 1---------

output$SDAtsne_br1 <- renderPlot({
  
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    print("QC_compIter: ")
    # print(envv$QC_compIter)
    zN = envv$QC_compIter
    
    SDAres <- envv$SDAres
    tempDFX <- as.data.frame(envv$tsne_CS_qc$Y)
    colnames(tempDFX) <- c("tSNE1_qc", "tSNE2_qc")
    
    print(head(SDAres$scores))
    
    tempDFX$SDAComp <- SDAres$scores[,paste0("SDAV", zN, sep="")]
    
    # print(head(tempDFX))
    
    ggplot(tempDFX, aes(tSNE1_qc, tSNE2_qc,  color=cut(asinh(SDAComp^3), breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf)))) +
      geom_point(size = 1) + theme_bw() +
      scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
      guides(colour = guide_legend(override.aes = list(size = 2, alpha=1))) +
      theme(legend.position = "bottom", aspect.ratio=1) + 
      simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE) + 
      ggtitle(paste0("SDAV", zN, sep=""))+ylab("asinh(SDAscore^3)")
    
  }
})

## tSNE DGE---------

output$DGE_SDA_tSNE <- renderPlot({
  
  if(is.null(envv$tsne_CS_batch)){
    plot(x=0, y=0, main="tsne CS Batch not found")
    
  } else {
    
    tsneDF <- as.data.frame(envv$tsne_CS_batch$Y)
    rownames(tsneDF)  <- rownames(envv$SDAres$scores)
    colnames(tsneDF) <- c("tSNE1_batch", "tSNE2_batch")
    
    
    
    
    if(is.null(envv$MetaDF)){
      tsneDF$SumScore <- rowSums(abs(envv$SDAres$scores))
      tsneDF$SumScore <- tsneDF$SumScore/mean(tsneDF$SumScore)
      
      ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=(SumScore))) +
        geom_point(size = 1) + theme_bw() +
        scale_color_distiller(palette = "Spectral")  +
        ggtitle("tSNE SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n ")+
        theme(legend.position = "bottom", aspect.ratio=1)
      
    } else {
      MetaDF <- envv$MetaDF
      tsneDF$library_size <- MetaDF[rownames(tsneDF), ]$library_size
      
      ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=log10(library_size))) +
        geom_point(size = 1) + theme_bw() +
        scale_color_distiller(palette = "Spectral")  +
        ggtitle("tSNE SDA batch removed\n log10 library size \n ")+
        theme(legend.position = "bottom", aspect.ratio=1)
    }
  }
  
})


## tSNE CS batch 1---------

output$tSNE_CS_batch1 <- renderPlot({
  
  if(is.null(envv$tsne_CS_batch)){
    plot(x=0, y=0, main="tsne CS Batch not found")
    
  } else {
    
    tsneDF <- as.data.frame(envv$tsne_CS_batch$Y)
    rownames(tsneDF)  <- rownames(envv$SDAres$scores)
    colnames(tsneDF) <- c("tSNE1_batch", "tSNE2_batch")
    
    
    
    
    if(is.null(envv$MetaDF)){
      tsneDF$SumScore <- rowSums(abs(envv$SDAres$scores))
      tsneDF$SumScore <- tsneDF$SumScore/mean(tsneDF$SumScore)
      
      ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=(SumScore))) +
        geom_point(size = 1) + theme_bw() +
        scale_color_distiller(palette = "Spectral")  +
        ggtitle("tSNE SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n No Meta loaded")+
        theme(legend.position = "bottom", aspect.ratio=1)
      
    } else {
      MetaDF <- envv$MetaDF
      
      tsneDF$Meta <- MetaDF[rownames(tsneDF), input$Metaselect1]
      
      ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=Meta)) +
        geom_point(size = 1, alpha=.4)+ theme_bw() +
        theme(legend.position = "bottom", aspect.ratio=1) +
        ggtitle(paste0("tSNE - batch removed cell scores\n", input$Metaselect1)) +
        scale_color_manual(values = col_vector
                           #c(rep(colorRampPalette(brewer.pal(12,"Paired"))(30),2),"black","grey")
        ) + 
        guides(colour = guide_legend(override.aes = list(size = 2, alpha=1), ncol=5))
      
      
    }
  }
  
})

## tSNE CS batch 2---------

output$tSNE_CS_batch2 <- renderPlot({
  
  if(is.null(envv$tsne_CS_batch)){
    plot(x=0, y=0, main="tsne CS Batch not found")
    
  } else {
    
    tsneDF <- as.data.frame(envv$tsne_CS_batch$Y)
    rownames(tsneDF)  <- rownames(envv$SDAres$scores)
    colnames(tsneDF) <- c("tSNE1_batch", "tSNE2_batch")
    
    
    
    
    if(is.null(envv$MetaDF)){
      tsneDF$SumScore <- rowSums(abs(envv$SDAres$scores))
      tsneDF$SumScore <- tsneDF$SumScore/mean(tsneDF$SumScore)
      
      ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=(SumScore))) +
        geom_point(size = 1) + theme_bw() +
        scale_color_distiller(palette = "Spectral")  +
        ggtitle("tSNE SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n No Meta loaded")+
        theme(legend.position = "bottom", aspect.ratio=1)
      
    } else {
      MetaDF <- envv$MetaDF
      
      tsneDF$Meta <- MetaDF[rownames(tsneDF), input$Metaselect2]
      
      ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=Meta)) +
        geom_point(size = 1, alpha=.4)+ theme_bw() +
        theme(legend.position = "bottom", aspect.ratio=1) +
        ggtitle(paste0("tSNE - batch removed cell scores\n", input$Metaselect2)) +
        scale_color_manual(values = col_vector
                           #c(rep(colorRampPalette(brewer.pal(12,"Paired"))(30),2),"black","grey")
        ) + 
        guides(colour = guide_legend(override.aes = list(size = 2, alpha=1), ncol=5))
      
      
    }
  }
  
})

## tSNE CS batch 3---------


output$tSNE_CS_batch3 <- renderPlot({
  
  if(is.null(envv$tsne_CS_batch)){
    plot(x=0, y=0, main="tsne CS Batch not found")
    
  } else {
    
    tsneDF <- as.data.frame(envv$tsne_CS_batch$Y)
    rownames(tsneDF)  <- rownames(envv$SDAres$scores)
    colnames(tsneDF) <- c("tSNE1_batch", "tSNE2_batch")
    
    
    
    
    if(is.null(envv$MetaDF)){
      tsneDF$SumScore <- rowSums(abs(envv$SDAres$scores))
      tsneDF$SumScore <- tsneDF$SumScore/mean(tsneDF$SumScore)
      
      ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=(SumScore))) +
        geom_point(size = 1) + theme_bw() +
        scale_color_distiller(palette = "Spectral")  +
        ggtitle("tSNE SDA batch removed\n  Sum absolute-cell-scores normalized by its mean \n No Meta loaded")+
        theme(legend.position = "bottom", aspect.ratio=1)
      
    } else {
      MetaDF <- envv$MetaDF
      
      tsneDF$Meta <- MetaDF[rownames(tsneDF), input$Metaselect3]
      
      ggplot(tsneDF, aes(tSNE1_batch, tSNE2_batch, color=Meta)) +
        geom_point(size = 1, alpha=.4)+ theme_bw() +
        theme(legend.position = "bottom", aspect.ratio=1) +
        ggtitle(paste0("tSNE - batch removed cell scores\n", input$Metaselect3)) +
        scale_color_manual(values = col_vector
                           #c(rep(colorRampPalette(brewer.pal(12,"Paired"))(30),2),"black","grey")
        ) + 
        guides(colour = guide_legend(override.aes = list(size = 2, alpha=1), ncol=5))
      
      
    }
  }
  
})


# batch removed Figs -----

## tSNE batch removed 2---------

output$SDAtsne_br2 <- renderPlot({
  
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    zN = envv$QC_compIter
    
    SDAres <- envv$SDAres
    tempDFX <- as.data.frame(envv$tsne_CS_batch$Y)
    rownames(tempDFX)  <- rownames(envv$SDAres$scores)
    colnames(tempDFX) <- c("tSNE1_batch", "tSNE2_batch")
    
    if(zN %in% envv$Remove_comps) RemoveTag = "removed" else RemoveTag = "kept"
    
    
    tempDFX$SDAComp <- cut(asinh(SDAres$scores[,paste0("SDAV", zN, sep="")]^3), 
                           breaks = c(-Inf, -1, -.5, 0, .5, 1, Inf))
    
    tempDFX$SDAComp <- factor(tempDFX$SDAComp, 
                              levels = c("(-Inf,-1]", "(-1,-0.5]", "(-0.5,0]",  "(0,0.5]",   "(0.5,1]",   "(1, Inf]" ) )
    # print(tempDFX$SDAComp)
    # print(factor(tempDFX$SDAComp))
    # print(factor(tempDFX$SDAComp, 
    #              levels = c("(-Inf,-1]", "(-1,-0.5]", "(-0.5,0]",  "(0,0.5]",   "(0.5,1]",   "(1, Inf]" ) ))
    
    ggplot(tempDFX, aes(tSNE1_batch, tSNE2_batch,  color=tempDFX$SDAComp)) +
      geom_point(size = 1) + theme_bw() +
      scale_color_manual("CS", values = rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue")) ) + 
      guides(colour = guide_legend(override.aes = list(size = 2, alpha=1))) +
      theme(legend.position = "bottom", aspect.ratio=1) + 
      simplify2 + coord_cartesian(xlim = NULL, ylim = NULL, expand = FALSE) +
      ggtitle(paste0("SDAV", zN, " :: ", RemoveTag))+
      ylab("asinh(SDAscore^3)")
    
    
    
    
    
  }
})

## tSNE br2Tab---------

output$SDAtsne_br2Tab <- renderPlot({
  
  if(is.null(envv$SDAres)){
    plot(x=0, y=0, main="Load an SDA")
  } else {
    
    zN = envv$QC_compIter
    
    SDAres <- envv$SDAres
    tempDFX <- as.data.frame(envv$tsne_CS_batch$Y)
    rownames(tempDFX)  <- rownames(envv$SDAres$scores)
    colnames(tempDFX) <- c("tSNE1_batch", "tSNE2_batch")
    
    if(zN %in% envv$Remove_comps) RemoveTag = "removed" else RemoveTag = "kept"
    
    
    tempDFX$SDAComp <- SDAres$scores[,paste0("SDAV", zN, sep="")]
    
    
    
    if(!is.null(envv$MetaDF)){
      MetaDF <- envv$MetaDF
      tempDFX$Meta <- MetaDF[rownames(tempDFX), input$Metaselect3]
      tempDFX <- table(cut(asinh(tempDFX$SDAComp^3), 
                           c(-Inf, -1, -.5, 0, .5, 1, Inf)), tempDFX$Meta)
      # print(tempDFX)
      # print(rownames(tempDFX))
      
      
      tempDFX <-  tempDFX[rowSums(tempDFX)!=0, ]
      ppg2 <- ggplot(reshape2::melt(tempDFX)) +
        geom_bar(aes(x=as.character(Var2), y=value, fill=factor(Var1, levels=c("(-Inf,-1]", "(-1,-0.5]", "(-0.5,0]",  "(0,0.5]",   "(0.5,1]",   "(1, Inf]" ) )), 
                 stat="identity", width = 0.7, position="fill") +
        theme_bw()  + scale_fill_manual(values=rev(c("red", "orange", "yellow", "lightblue", "dodgerblue", "blue"))) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                           labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
        ggtitle(paste0("Relative Contribution\n","SDA", zN, " :: ", RemoveTag)) + ylab("Relative % cells")
      
      print(ppg2)
      
      
    } else {
      plot(x=0, y=0, main="No Meta")
    }
    
    
    
  }
})


## GO figs -----
### GO pos ---------

output$GOpos <- renderPlot({
  
  
  
  if(is.null(envv$GO_data)){
    plot(x=0, y=0, main="Load an SDA the GO data")
    
  } else {
    GO_data  <- envv$GO_data
    zN = envv$QC_compIter
    go_volcano_plot(x=GO_data, 
                    component = paste("V", zN, "P", sep=""))+ 
      theme_bw()+ theme(aspect.ratio = 1)
    
  }
  
})

### GO neg ---------

output$GOneg <- renderPlot({
  
  if(is.null(envv$GO_data)){
    plot(x=0, y=0, main="Load an SDA the GO data")
    
  } else {
    GO_data  <- envv$GO_data
    zN = envv$QC_compIter
    go_volcano_plot(x=GO_data, component = paste("V", zN, "N", sep=""))+ 
      theme_bw()+ theme(aspect.ratio = 1)
    
  }
  
})


# Enrich & Explore figs -----


## Gene Enrichmend SDA pos ------
output$GenesEnrichSDAPos <- renderPlot({
  
  SDAres <- envv$SDAres
  SDA_TopNpos <- envv$SDA_TopNpos
  
  envv$TopN
  
  # N = total number of genes (usually not entire genome, since many have unk func)
  N=length(colnames(SDAres$loadings[[1]]))
  # k = number of genes submitted, top N
  k = envv$TopN
  
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
  
  GeneSetNot <- GeneSet[!GeneSet %in% colnames(SDAres$loadings[[1]][,])]
  
  print("length of your genes:")
  print(length(GeneSet))
  GeneSet <- GeneSet[GeneSet %in% colnames(SDAres$loadings[[1]][,])]
  print("length of your genes in this dataset:")
  print(length(GeneSet))
  
  
  
  
  # print("length of your genes in this dataset:")
  # print(length(GeneSet))
  
  plotEnrich(GeneSetsDF=SDA_TopNpos, 
             GeneVec = GeneSet, 
             plotTitle= paste0("Gene-set enrichment\n SDA top ", k, " pos loadings\nGene universe size: ", N, "\n Hypergeometric test: * adj.p < 0.01 \n Genes not found: ",
                               paste0(GeneSetNot, collapse = ", ")),
             xLab = "SDA Comps",
             N=N,
             k=k)
  
  
})

## Gene Enrichmend SDA neg ------

output$GenesEnrichSDANeg <- renderPlot({
  
  SDAres <- envv$SDAres
  SDA_TopNneg <- envv$SDA_TopNneg
  
  # envv$TopN
  
  # N = total number of genes (usually not entire genome, since many have unk func)
  N=length(colnames(SDAres$loadings[[1]]))
  # k = number of genes submitted, top N
  k = envv$TopN
  
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
  
  
  GeneSetNot <- GeneSet[!GeneSet %in% colnames(SDAres$loadings[[1]][,])]
  
  print("length of your genes:")
  print(length(GeneSet))
  GeneSet <- GeneSet[GeneSet %in% colnames(SDAres$loadings[[1]][,])]
  print("length of your genes in this dataset:")
  print(length(GeneSet))
  
  
  
  
  # print("length of your genes in this dataset:")
  # print(length(GeneSet))
  
  plotEnrich(GeneSetsDF=SDA_TopNneg, 
             GeneVec = GeneSet, 
             plotTitle= paste0("Gene-set enrichment\n SDA top ", k, " neg loadings\nGene universe size: ", N, "\n Hypergeometric test: * adj.p < 0.01 \n Genes not found: ",
                               paste0(GeneSetNot, collapse = ", ")),
             xLab = "SDA Comps",
             N=N,
             k=k)
  
  
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