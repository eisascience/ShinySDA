
output$SDAScoresChiPos <- renderPlot({
  
  # ColFac_DONR.ID <- CDID()
  
  if(is.null(envv$MetaDF)){
    print("No Comp")
  } else {
    
  
    SDAScores <- envv$SDAres$scores
    ComponentN <- as.numeric(envv$QC_compIter)
    MetaDF <- envv$MetaDF
    MetaDF <- MetaDF[rownames(SDAScores),]
    
    
    #TODO: seperate tab to have this chisqr with drop down to select from envv$PossibleMetaVec
    
    ShinySDA::plot_SDA_ChiSqrHM(MetaDF = MetaDF, MetaSelect=input$Metaselect4, 
                      SDAScores=SDAScores, direction="pos", clustStat = T)
      

    
  }
  
})


output$SDAScoresChiNeg <- renderPlot({
  
  # ColFac_DONR.ID <- CDID()
  
  if(is.null(envv$MetaDF)){
    print("No Comp")
  } else {
    SDAScores <- envv$SDAres$scores
    ComponentN <- as.numeric(envv$QC_compIter)
    MetaDF <- envv$MetaDF
    MetaDF <- MetaDF[rownames(SDAScores),]
    
    
    
    ShinySDA::plot_SDA_ChiSqrHM(MetaDF = MetaDF, MetaSelect=input$Metaselect4, 
                      SDAScores=SDAScores, direction="neg", clustStat = T)
    
    
    
    
    
  }
  
})