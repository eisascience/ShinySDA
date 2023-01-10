## get GO
observeEvent(input$getSDAGo, {
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
    envv$InfoBox_sub = "Checking paths for tSNE"
    
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
    
   
    NComps <- as.numeric(envv$SDAres$command_arguments$num_comps)
    

    print("starting GO")
    
    SDAres <- envv$SDAres
    
    envv$InfoBox_sub <- paste0(NComps, " comps, ~30 sec per comp")
    
    
    print(paste0(envv$path2SDA_dyn, "/", head.path,"_SDA_GO_Comps",input$species, ".rds"))
    
    if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,"_SDA_GO_Comps",input$species, ".rds"))){
      
      library(AnnotationHub) # source("https://bioconductor.org/biocLite.R"); biocLite("AnnotationHub")
      library(clusterProfiler) # source("https://bioconductor.org/biocLite.R"); biocLite("clusterProfiler")
      
      hub <- AnnotationHub()
      
      if(input$species == "mouse") qrhub="org.MM.eg"
      if(input$species == "human") qrhub="org.Hs.eg.db"
      if(input$species == "rhesus") qrhub="org.Mmu.eg.db"
      
      RefGenome.names <- query(hub, qrhub)#org.Hs.eg.db  org.MM.eg
      
      
      # if(input$species == "mouse") qrhub.id="AH84123"
      # if(input$species == "human") qrhub.id="org.Hs.eg.db"
      # if(input$species == "rhesus") qrhub.id="org.Mmu.eg.db"
      
      print(RefGenome.names$ah_id)
      RefGenome <- hub[[RefGenome.names$ah_id]] 
      
      
      envv$GOAnn <- list(RefGenome = RefGenome, RefGenome.names = RefGenome.names)
      
      GO_data <- list()
      
      for (i in 1:NComps){
        print(i)
        envv$InfoBox_sub <- paste0("Getting GO: %", round(i/NComps,3)*100)
        
        print("...negatives")
        GO_data[[paste0("V",i,"N")]] <- GO_enrichment(results =SDAres, i, side="N", geneNumber = 100, threshold=0.05, OrgDb =RefGenome)
        print("...positives")
        GO_data[[paste0("V",i,"P")]] <- GO_enrichment(results =SDAres, i, side="P", geneNumber = 100, threshold=0.05, OrgDb =RefGenome)
      }
      
      saveRDS(GO_data, paste0(envv$path2SDA_dyn, "/", head.path,"_SDA_GO_Comps",input$species, ".rds"))
      
    } else{
      print("Loaded previously saved GO annotations.")
      
      GO_data <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_SDA_GO_Comps",input$species, ".rds"))
    }
    
    envv$InfoBox_sub <- "GO data loaded"
    envv$GO_data  <- GO_data
    
  }
  
  
  
  
})