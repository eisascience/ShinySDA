


output$table_SDAScoresChiPos <- DT::renderDataTable({

  if(!is.null(envv$ChiSqrMeta_tabDF)){
    envv$ChiSqrMeta_tabDF = as.data.frame(envv$ChiSqrMeta_tabDF)
    # print(envv$ChiSqrMeta_tabDF)
    print(nrow(envv$ChiSqrMeta_tabDF))
    # print(length(ifelse(envv$FailingKertosis, "fail", "pass")))
    envv$ChiSqrMeta_tabDF$Kertosis = ifelse(envv$SDAres$FailingKertosis, "fail", "pass")
    envv$ChiSqrMeta_tabDF$ScoreThreshold = ifelse(envv$SDAres$FailScoreThreshold, "fail", "pass")
    envv$ChiSqrMeta_tabDF$MaxScore = ifelse(envv$SDAres$FailingMaxScore, "fail", "pass")
    envv$ChiSqrMeta_tabDF$FinalCall = ifelse(envv$SDAres$FailingFilters, "fail", "pass")
    
      
  } else {
    envv$ChiSqrMeta_tabDF = data.frame(A1=c("", "", "","", "", ""), A2=c("", "", "","", "", ""))
  }
  
  
  
  if(envv$Origin == "folder"){
    head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
    base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
    
    head.path <- gsub("/", "", head.path)
    
  }
  if(envv$Origin == "prime"){
    head.path <- "sda_results"
    base.path <- ""
    #TODO simplify?
  }

  DT::datatable(envv$ChiSqrMeta_tabDF, 
                extensions = "Buttons", 
            options = list(pageLength = nrow(envv$ChiSqrMeta_tabDF), #
                           order = list(list(0, 'asc')), #
                           lengthChange = F,
                           autoWidth = T, #
                           fixedHeader = T, #
                           paging = T,
                           scrollX=T, 
                           searching = F,
                           ordering = T,
                           columnDefs = list(list(className = 'dt-center', targets = "_all")),#
                           css = c("table.dataTable { background-color: white; }", "table.dataTable { color: black; }"),  #
                           dom = 'Bfrtip',
                           buttons = c('copy', 'csv', 'excel', 'pdf'),
                           pageLength=5, 
                           lengthMenu=c(3,5,10) )) %>% DT::formatStyle(columns = c(1:6), fontSize = '10px')
   
  # DT::datatable(envv$ChiSqrMeta_tabDF,
  #               options = list(pageLength = 30, #searching = F, lengthChange = F,
  #               extensions = c("Buttons"),
  #               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
  #               columnDefs = list(list(className = 'dt-center', targets = "_all")),
  #               # language = list(
  #               #   info = 'Showing _START_ to _END_ of _TOTAL_ entries',
  #               #   lengthMenu = 'Show _MENU_ entries',
  #               #   paginate = list(
  #               #     next = '>',
  #               #     previous = '<'
  #               #   )
  #               # ),
  #               css = c("table.dataTable { background-color: white; }", "table.dataTable { color: black; }"),
  #               dom = 'tB',
  #               autoWidth = TRUE,
  #               order = list(list(0, 'asc')),
  #               fixedHeader = TRUE#,
  #               # download = paste0(envv$path2SDA_dyn, "/", head.path,"_ChiSqrMetaDataAnalysis.csv")
  #               ) ) %>% DT::formatStyle(columns = c(1,2,3), fontSize = '10px')
  
})



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
    
    ShinySDA:::plot_SDA_ChiSqrHM(MetaDF = MetaDF, MetaSelect=input$Metaselect4, 
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