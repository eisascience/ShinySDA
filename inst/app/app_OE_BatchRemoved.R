observeEvent(input$nextSDA_br2, {
  
  
  SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
  
  
  
  if(envv$QC_compIter + 1 %in% SDAorder){
    # envv$QC_compIter = SDAorder[which(SDAorder==envv$QC_compIter) + 1]
    # updateTextInput(session, "SDAVn", value = envv$QC_compIter)
    updateTextInput(session, "SDAVn", value = SDAorder[which(SDAorder==envv$QC_compIter) + 1])
    
  }
  
  
  
  
})

observeEvent(input$prevSDA_br2, {
  
  # choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
  # SDAorder <- setdiff(choice, envv$Remove_comps)
  
  SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
  
  
  if(envv$QC_compIter - 1 %in% SDAorder){
    #envv$QC_compIter = SDAorder[which(SDAorder==envv$QC_compIter) - 1]
    # updateTextInput(session, "SDAVn", value = envv$QC_compIter)
    updateTextInput(session, "SDAVn", value = SDAorder[which(SDAorder==envv$QC_compIter) - 1])
  }
  
  
  
  
})


observeEvent(input$SDAVn, {
  
  
  
  if(!is.null(envv$QC_compIter)){
    
    SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
    
    if(is.null(input$SDAVn) | (!(as.numeric(input$SDAVn) %in% SDAorder))) {
      
      
      print("Its null")
      updateTextInput(session, "SDAVn", value = envv$QC_compIter)
    } else {
      envv$QC_compIter = as.numeric(input$SDAVn)
    }
  }
  
  
  
  # # SDAorder <- 1:as.numeric(envv$SDAres$command_arguments$num_comps)
  # if(!is.null(envv$QC_compIter)){
  #   
  #   if(!is.null(input$SDAVn)) {
  #     
  #     
  #   } else {
  #     # input$SDAVn <- envv$QC_compIter
  #   }
  # } else {
  #     
  #   }
  
  # choice <-  1:as.numeric(envv$SDAres$command_arguments$num_comps) #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
  # SDAorder <- setdiff(choice, envv$Remove_comps)
  # 
  # input$SDAVn <- as.character(envv$QC_compIter)
  # 
  # if(as.numeric(input$SDAVn) %in% SDAorder) {
  # 
  #   envv$QC_compIter = as.numeric(input$SDAVn)
  # 
  # } else {
  
  
  
  
  
})


observeEvent(input$C2Cpos, {
  
  
  Out1 <- print_gene_list(results=envv$SDAres, as.numeric(envv$QC_compIter), PosOnly = T) %>%
    #group_by(package) %>%
    #tally() %>%
    #arrange(desc(n), tolower(package)) %>%
    #mutate(percentage = n / nrow(pkgData()) * 100) %>%
    #select("Package name" = package, "% of downloads" = percentage) %>%
    as.data.frame() %>%
    head(as.numeric(input$NoOfGenes)) 
  Out1 <- Out1$Gene.Name
  
  # print(Out1)
  clipr::write_clip(Out1)
  
})


observeEvent(input$C2Cneg, {
  
  
  Out2 <- print_gene_list(results=envv$SDAres, as.numeric(envv$QC_compIter), NegOnly = T) %>%
    #group_by(package) %>%
    #tally() %>%
    #arrange(desc(n), tolower(package)) %>%
    #mutate(percentage = n / nrow(pkgData()) * 100) %>%
    #select("Package name" = package, "% of downloads" = percentage) %>%
    as.data.frame() %>%
    head(as.numeric(input$NoOfGenes)) 
  Out2 <- Out2$Gene.Name
  
  # print(Out1)
  clipr::write_clip(Out2)
  
})