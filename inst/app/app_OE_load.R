
observeEvent(input$Download_SDA_primeseq, {
  
  if(input$SDA_OFId != "" & nchar(input$SDA_OFId)==6){
    envv$Origin = "prime"
    envv$InfoBox_sub = "Downloading from Prime-seq"

  
    
    labkey.setDefaults(apiKey=input$apiKey)
    Rdiscvr::SetLabKeyDefaults(baseUrl = input$baseURL, 
                               defaultFolder = input$defaultFolder)
    

    envv$sda_path = ShinySDA::download_SDA(outf = input$SDA_PS_save, outFileID=input$SDA_OFId)
    envv$sda_dims_path = ShinySDA::download_SDA_dimnames(outf = input$SDA_PS_save, outFileID=input$SDA_OFId)
    
    envv$path2SDA_dyn <- paste0(input$SDA_PS_save, "sda.", input$SDA_OFId)
    
    if(!dir.exists(envv$path2SDA_dyn)) {
      dir.create(envv$path2SDA_dyn)
      dir.create(paste0(envv$path2SDA_dyn, "/", "sda_results"))
    }
    
    envv$InfoBox_sub = "Downloading Completed"
    
  } else {
    envv$Origin = "unk"
    envv$InfoBox_sub = "OutFile ID error, needs to be a 6-digit number"
    envv$sda_path = NULL
    envv$sda_dims_path = NULL
  }
  
  
})


observeEvent(input$Load_SDA_primeseq, {

  envv$InfoBox_sub = "Loading from Prime-seq download"
  
  if(!is.null(envv$sda_path)){
    print(envv$sda_path)
    SDAres = readRDS(envv$sda_path)
    
    if(!is.null(envv$sda_dims_path)){
      NamesDimsDF = readRDS(envv$sda_dims_path)
      rownames(SDAres$scores) <- NamesDimsDF[[1]]
      colnames(SDAres$loadings[[1]]) <- NamesDimsDF[[2]]
      
    }
    colnames(SDAres$scores) = gsub("V", "SDAV", colnames(SDAres$scores))
    rownames(SDAres$loadings[[1]]) = gsub("V", "SDAV", rownames(SDAres$loadings[[1]]))
    
    if(is.null(envv$MetaDF)){
      
      
      #TODO: grab metaDF somehow
      
      if(!file.exists(paste0(envv$path2SDA_dyn, "/", input$SDA_OFId, "_MetaDF.rds"))){
        envv$InfoBox_sub = "Downloading Metadata"
        # metaDF = get_MetaDF(cellbarcodes = rownames(SDAres$scores))
        MetaDF = ShinySDA::DownloadMetadataForSdaResults(input$SDA_OFId)
        saveRDS(MetaDF, paste0(envv$path2SDA_dyn, "/", input$SDA_OFId, "_MetaDF.rds"))
      } else {
        MetaDF = readRDS(paste0(envv$path2SDA_dyn, "/", input$SDA_OFId, "_MetaDF.rds"))
      }
      
      rownames(MetaDF) = MetaDF$cellbarcode
      
      envv$InfoBox_sub = "Metadata Loaded"
      
      
      PossibleMetaVec = c("SampleDate", "SubjectId", "ExpID",
        "SingleR_Labels","SingleR_Labels_Fine",
        "hpca.label", "hpca.label.fine",
        "blueprint.label", "blueprint.label.fine",
        "predicted_labels", "majority_voting",
        "Phase", "BarcodePrefix", "DatasetId",
        "Population", "WorkbookId",
        "RNA_snn_res.0.2", "RNA_snn_res.0.6", "RNA_snn_res.1.2")

      PossibleMetaVec = PossibleMetaVec[PossibleMetaVec %in% colnames(MetaDF)]

      envv$PossibleMetaVec = PossibleMetaVec

      updateSelectInput(session, 
                        "Metaselect1",
                        choices = PossibleMetaVec, 
                        selected = PossibleMetaVec[1])
      
      updateSelectInput(session, 
                        "Metaselect2",
                        choices = PossibleMetaVec, 
                        selected = PossibleMetaVec[1])
      
      
      updateSelectInput(session, 
                        "Metaselect3",
                        choices = PossibleMetaVec, 
                        selected = PossibleMetaVec[1])
      
      updateSelectInput(session, 
                        "Metaselect4",
                        choices = PossibleMetaVec, 
                        selected = PossibleMetaVec[1])
      
      #Cell scores Boxplot
      updateSelectInput(session, 
                        "Metaselect5",
                        choices = PossibleMetaVec, 
                        selected = PossibleMetaVec[1])
      
      #Cell scores Across
      updateSelectInput(session, 
                        "Metaselect6",
                        choices = PossibleMetaVec, 
                        selected = PossibleMetaVec[1])
      
      
      
      #$SubjectId, $ExpID, $EXP.ID, $SampleDate, $SingleR_Labels, $BarcodePrefix
      
      for(featX in PossibleMetaVec){
        #TODO: seprate out cleaning of meta in a fx
        
        if(featX == "SubjectId") { MetaDF[,featX] = paste0("Rh", MetaDF[,featX]) }
        if(featX %in% c("RNA_snn_res.0.2", "RNA_snn_res.0.6", "RNA_snn_res.1.2")) { MetaDF[,featX] = paste0("Cl", MetaDF[,featX]) }
        
        MetaDF[,featX] = factor( MetaDF[,featX])
      }
      
      # print(head(MetaDF))
      
      MetaDF$library_size = MetaDF$nCount_RNA
      MetaDF$ExpID = MetaDF$Population
      MetaDF$EXP.ID = MetaDF$RNA_snn_res.0.2 #TODO this is a patch EXP.ID needs to be removed and consquent fx updated to use envv$PossibleMetaVec
      
      envv$MetaDF = MetaDF
    }
    
    
    
    # print(head(MetaDF))
    
    envv$InfoBox_sub = "Load from Prime-seq completed, starting pre-processing"
    
    #adds some stats
    SDAres      <- ShinySDA::AddCompStats(SDAres)
    
    print("CompStats added")
    
    envv$SDAres <- SDAres
    
    envv = ShinySDA::preprocess_SDA(SDAres = SDAres, QuantThr = 0.95, envv = envv,
                   TopN = 150, MetaDF = NULL)
    
    
    
    envv$InfoBox_sub = "Load and pre-processing complete on Prime-seq sourced sda"
    
  }
  
 
  
  
})



observeEvent(input$loadSDA, {
  
  envv$Origin = "folder"
  envv$InfoBox_sub = "Loading from SDA folder (traditional)"
  
  
  # print(head(paste0(input$SDAroot, "/", input$folder.name)))
  
  envv$path2SDA_dyn <- paste0(input$SDAroot, "/", input$folder.name)
  
  if(file.exists(envv$path2SDA_dyn)) {
    
    print(envv$path2SDA_dyn)
    
    envv$InfoBox_sub = "Loading SDA"
    
    
    if(sum(grepl("_dimnames", list.files(dirname(envv$path2SDA_dyn), recursive = T)))==0){
      DimNamesPath <- paste0(dirname(dirname(envv$path2SDA_dyn)), "/")
    } else {
      DimNamesPath <- paste0(dirname(envv$path2SDA_dyn), "/")
      log
    }
    print(DimNamesPath)
    
    
    SDAres <- SDAtools::load_results(
      results_folder = envv$path2SDA_dyn,
      data_path =  DimNamesPath)
    
    print("SDA results loaded")
    
    
    #update the names
    colnames(SDAres$scores) <- paste("SDA", 1:ncol(SDAres$scores), sep="")
    rownames(SDAres$loadings[[1]]) <- paste("SDA", 1:ncol(SDAres$scores), sep="")
    
    if(file.exists(paste0(envv$path2SDA_dyn, "_dimnames.rds"))){
      print(head(NamesDimsDF))
      NamesDimsDF = readRDS(paste0(envv$path2SDA_dyn, "_dimnames.rds"))
      rownames(SDAres$scores) <- NamesDimsDF[[1]]
      colnames(SDAres$loadings[[1]]) <- NamesDimsDF[[2]]
    }
    
    # print(SDAres$loadings[[1]][1:10,1:10])
    

    envv$InfoBox_sub = "Load from folder complete, initiating pre-processing"
    
    #adds some stats
    SDAres      <- AddCompStats(SDAres)
    
    print("CompStats added")
    
    envv$SDAres <- SDAres
    
    
    envv = preprocess_SDA(SDAres = SDAres, QuantThr = 0.95, envv = envv,
                   TopN = 150, MetaDF = NULL)
    
    
    
    envv$InfoBox_sub = "Load and pre-processing complete on Prime-seq sourced sda"
    
  } else { 
    # updateTextInput(session, "loadSDAmsg", value = "File not found")
    
    
  }
})


