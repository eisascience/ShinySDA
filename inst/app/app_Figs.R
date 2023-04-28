
# Main Figs -------




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



## GO figs -----
### GO pos ---------

output$GOpos <- renderPlot({
  
  
  
  if(is.null(envv$GO_data)){
    plot(x=0, y=0, main="Load an SDA the GO data")
    
  } else {
    GO_data  <- envv$GO_data
    zN = envv$QC_compIter
    if(any(grepl("Pos", names(envv$GO_data)))){
      ShinySDA::go_volcano_plot(x=envv$GO_data, component = paste(zN, "-Pos", sep=""))+ 
        theme_bw()+ theme(aspect.ratio = 1)
    } else {
      if(any(grepl("V", names(envv$GO_data)))){
        ShinySDA::go_volcano_plot(x=envv$GO_data, component = paste("V", zN, "P", sep=""))+ 
          theme_bw()+ theme(aspect.ratio = 1)
      }
    }
    
  }
  
})

### GO neg ---------

output$GOneg <- renderPlot({
  
  if(is.null(envv$GO_data)){
    plot(x=0, y=0, main="Load an SDA the GO data")
    
  } else {
    GO_data  <- envv$GO_data
    zN = envv$QC_compIter
    
    if(any(grepl("Neg", names(envv$GO_data)))){
      ShinySDA::go_volcano_plot(x=envv$GO_data, component = paste(zN, "-Neg", sep=""))+ 
        theme_bw()+ theme(aspect.ratio = 1)
    } else {
      if(any(grepl("V", names(envv$GO_data)))){
        ShinySDA::go_volcano_plot(x=envv$GO_data, component = paste("V", zN, "N", sep=""))+ 
          theme_bw()+ theme(aspect.ratio = 1)
      }
    }
    
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


