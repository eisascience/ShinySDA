
## get biomart
observeEvent(input$getGeneAnn, {
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
    
    envv$InfoBox_sub = "Checking paths"
    
    library(biomaRt)
    
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
    
    # print("head path:")
    # print(head.path)
    # print("base path:")
    # print(base.path)
    # print("path2SDA_dyn:")
    # print(envv$path2SDA_dyn)
    
    
    SDAres <- envv$SDAres
    
    # stringr::str_split("../../../../Conrad/R/Utah/sda_results/Testis_Run2", "sda_results/")
    
    # basename("../../../../Conrad/R/Utah/sda_results/Testis_Run2")
    
    print(paste0("looking for files: ", envv$path2SDA_dyn, "/", head.path, "_SDAtools_GeneLoc_", input$species, ".rds"))
    print(paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_", input$species,".chromosome.lengths.rds"))
    print(paste0(envv$path2SDA_dyn, "/", head.path,  "_biomaRt_gene_loc_", 
                 input$species, ".rds"))
    
    if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path,  "_biomaRt_gene_loc_",
                           input$species, ".rds"))){
      
      envv$InfoBox_sub = "Downloading gene locations"

      if(input$species == "mouse") ens_ds="mmusculus_gene_ensembl"
      if(input$species == "human") ens_ds="hsapiens_gene_ensembl"
      if(input$species == "rhesus") ens_ds="mmulatta_gene_ensembl"

      print("getting gene locations")
      gene_locations <- ShinySDA::get.location(gene.symbols=colnames(SDAres$loadings[[1]]),
                                     data_set = ens_ds,
                                     gene_name = "external_gene_name")

      print("saving gene locations")
      saveRDS(gene_locations, paste0(envv$path2SDA_dyn, "/", head.path,  "_biomaRt_gene_loc_",
                                     input$species, ".rds"))

    } else {
      envv$InfoBox_sub = "Loading prev. downloaded gene locations"
      print("loading local biomaRt_gene_loc")
      gene_locations <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,  "_biomaRt_gene_loc_",
                                       input$species, ".rds"))
    }



    if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_GeneLoc_", input$species, ".rds"))){

      if(input$species == "mouse") ens_ds="mmusculus_gene_ensembl"
      if(input$species == "human") ens_ds="hsapiens_gene_ensembl"
      if(input$species == "rhesus") ens_ds="mmulatta_gene_ensembl"
      
      envv$InfoBox_sub = "Downloading Chr Lengths"

      print("loading SDAtools Chr Lengths")
      GeneLoc         <- SDAtools::load_gene_locations(path = base.path,
                                                       genes = colnames(SDAres$loadings[[1]]),
                                                       organism = ens_ds,
                                                       name="human")
      print("loading SDAtool chrom lengths")
      chromosome.lengths <- SDAtools::load_chromosome_lengths(organism = ens_ds)

      print("saving SDAtools chrom")
      saveRDS(GeneLoc, paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_GeneLoc_", input$species, ".rds"))
      saveRDS(chromosome.lengths, paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_", input$species,".chromosome.lengths.rds"))

    } else {
      envv$InfoBox_sub = "Loading prev. downloaded Chr Lengths"
      GeneLoc                <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_GeneLoc_", input$species, ".rds"))
      chromosome.lengths     <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_", input$species,".chromosome.lengths.rds"))
    }

    envv$chromosome.lengths <- chromosome.lengths
    envv$GeneLoc <- GeneLoc
    envv$gene_locations <- gene_locations
    
    envv$InfoBox_sub = "Chromosome location and lengths fully loaded"
    
  }
  
  
  
  
})