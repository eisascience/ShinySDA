
# Envv update fx ---------

#' This function to update the environment list when loading from local
#'
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @return the updated environment list
Run_LoadSDA_local_evv <- function(envv, input){

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
  envv$SDAres      <- AddCompStats(SDAres)
  
  # print("CompStats added")
  
  # envv$SDAres <- SDAres
  
  
  envv = preprocess_SDA(SDAres = envv$SDAres, QuantThr = 0.95, envv = envv,
                        TopN = 150, MetaDF = NULL)
  
  envv$Ncomp = ncol(SDAres$scores)
  envv$Ncell = nrow(SDAres$scores)
  envv$Ngene = ncol(SDAres$loadings[[1]]) 
  
  
  envv$InfoBox_sub = paste0("Load and pre-processing complete\n",
                            "Ncomp: ", envv$Ncomp, " Ncells: ", envv$Ncell, " Ngenes: ", envv$Ngene)
  
  
} else { 
  # updateTextInput(session, "loadSDAmsg", value = "File not found")
  envv$InfoBox_sub = "File not found loading error"
  
  
  
}
return(envv)
}

#' This function to update the environment list when doanloading SDA obj from prime-seq
#'
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @return the updated environment list
Run_DowLoadSDA_evv <- function(envv, input){
  
  if(input$SDA_OFId != "" & nchar(input$SDA_OFId)==6){
    envv$Origin = "prime"
    # envv$InfoBox_sub = "Downloading from Prime-seq"
    
    
    
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
    
    # envv$InfoBox_sub = "Downloading Completed"
    
  } else {
    envv$Origin = "unk"
    envv$InfoBox_sub = "OutFile ID error, needs to be a 6-digit number"
    envv$sda_path = NULL
    envv$sda_dims_path = NULL
  }
  
  return(envv)
}

#' This function to update the environment list when grabbing the meta data DF
#'
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @param session session  associated with ShinySDA
#' @return the updated environment list
Run_MetaDF_evv <- function(envv, input, session){
  
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

updateCheckboxGroupInput(session, "chkbx_MetaSave", choices = PossibleMetaVec, selected =  PossibleMetaVec[1])



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

return(envv)

}


#' This function to update the environment list when loading SDA from primeseq
#'
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @param session session  associated with ShinySDA
#' @return the updated environment list
Run_LoadSDA_evv <- function(envv, input, session){


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
    
    envv = ShinySDA:::Run_MetaDF_evv(envv, input, session)
  }
  
  
  
  # print(head(MetaDF))
  
  # envv$InfoBox_sub = "Load from Prime-seq completed, starting pre-processing"
  #adds some stats
  envv$SDAres      <- ShinySDA::AddCompStats(SDAres)
  
  # print("CompStats added")
  
  # envv$SDAres <- SDAres

  envv = ShinySDA::preprocess_SDA(SDAres = envv$SDAres, 
                                  QuantThr = 0.95, 
                                  envv = envv,
                                  TopN = 150, MetaDF = NULL)
  

  envv$Ncomp = ncol(SDAres$scores)
  envv$Ncell = nrow(SDAres$scores)
  envv$Ngene = ncol(SDAres$loadings[[1]])
 
  
  envv$InfoBox_sub = paste0("Load and pre-processing complete\n",
                            "Ncomp: ", envv$Ncomp, " Ncells: ", envv$Ncell, " Ngenes: ", envv$Ngene)
  
}
  
  return(envv)
  
}

#' This function to update the environment list when getting gene annotations
#'
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @return the updated environment list
Run_GeneAnn_evv <- function(envv, input){
  if(is.null(envv$SDAres)){ 
    # envv$InfoBox_sub = "Load SDA"
  } else {
    
    
    # envv$InfoBox_sub = "Checking paths"
    
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
      
      # envv$InfoBox_sub = "Downloading gene locations"
      
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
      # envv$InfoBox_sub = "Loading prev. downloaded gene locations"
      print("loading local biomaRt_gene_loc")
      gene_locations <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,  "_biomaRt_gene_loc_",
                                       input$species, ".rds"))
    }
    
    
    
    if(!file.exists(paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_GeneLoc_", input$species, ".rds"))){
      
      if(input$species == "mouse") ens_ds="mmusculus_gene_ensembl"
      if(input$species == "human") ens_ds="hsapiens_gene_ensembl"
      if(input$species == "rhesus") ens_ds="mmulatta_gene_ensembl"
      
      # envv$InfoBox_sub = "Downloading Chr Lengths"
      
      print("loading SDAtools Chr Lengths")
      GeneLoc         <- SDAtools::load_gene_locations(path = envv$path2SDA_dyn,
                                                       genes = colnames(SDAres$loadings[[1]]),
                                                       organism = ens_ds,
                                                       name="human")
      print("loading SDAtool chrom lengths")
      chromosome.lengths <- SDAtools::load_chromosome_lengths(organism = ens_ds)
      
      print("saving SDAtools chrom")
      saveRDS(GeneLoc, paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_GeneLoc_", input$species, ".rds"))
      saveRDS(chromosome.lengths, paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_", input$species,".chromosome.lengths.rds"))
      
    } else {
      # envv$InfoBox_sub = "Loading prev. downloaded Chr Lengths"
      GeneLoc                <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_GeneLoc_", input$species, ".rds"))
      chromosome.lengths     <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path, "_SDAtools_", input$species,".chromosome.lengths.rds"))
    }
    
    envv$chromosome.lengths <- chromosome.lengths
    envv$GeneLoc <- GeneLoc
    envv$gene_locations <- gene_locations
    
    envv$InfoBox_sub = "Chromosome location and lengths fully loaded"
    
  }
  
  return(envv)
}


#' This function to update the environment list when running GO enrichmet analysis
#'
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @return the updated environment list
Run_GO_evv <- function(envv, input){
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
    # envv$InfoBox_sub = "Checking paths for tSNE"
    
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
    # envv$InfoBox_sub <- "Stating GO"
    
    SDAres <- envv$SDAres
    
    
    
    
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
      
      
      if(!is.null(envv$SDAres$goEnrichment)){
        
        GO_data = envv$SDAres$goEnrichment
        # envv$InfoBox_sub <- paste0("GO was found in the SDA obj")
        
        
      } else{
        # envv$InfoBox_sub <- paste0(NComps, " comps, ~30 sec per comp")
        
        GO_data <- list()
        
        for (i in 1:NComps){
          print(i)
          # envv$InfoBox_sub <- paste0("Getting GO: %", round(i/NComps,3)*100)
          
          print("...negatives")
          GO_data[[paste0("V",i,"N")]] <- GO_enrichment(results =SDAres, i, side="N", geneNumber = 100, threshold=0.05, OrgDb =RefGenome)
          print("...positives")
          GO_data[[paste0("V",i,"P")]] <- GO_enrichment(results =SDAres, i, side="P", geneNumber = 100, threshold=0.05, OrgDb =RefGenome)
        }
      }
      
      saveRDS(GO_data, paste0(envv$path2SDA_dyn, "/", head.path,"_SDA_GO_Comps",input$species, ".rds"))
      
    } else{
      print("Loaded previously saved GO annotations.")
      
      GO_data <- readRDS(paste0(envv$path2SDA_dyn, "/", head.path,"_SDA_GO_Comps",input$species, ".rds"))
    }
    
    envv$InfoBox_sub <- "GO data loaded"
    envv$GO_data  <- GO_data
    
  }
  
  
  return(envv)
}



#' This function to update the environment list when running the fulls tSNE
#'
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @return the updated environment list
Run_tSNE_full_evv <- function(envv, input){
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
    
    SDAres <- envv$SDAres
    
    tsnepp <- round((nrow(SDAres$scores[,]) - 1)/3.2)
    tsnepp <- ifelse(tsnepp < 50, tsnepp, 50)
    # print("perplexity:")
    # print(tsnepp)
    
   
    
    envv$tsne_CS_all = ShinySDA:::.run_tSN_scores(SDAres = SDAres,
                                                  save_path=paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_AllComps_pp50.rds"), 
                                                  tsnepp=tsnepp, tSNE_n.iter=500, 
                                                  num_threads = 8, check_duplicates = F, 
                                                  selectComps = NULL )
    
    envv$InfoBox_sub = "raw tSNE with cell scores complete"
    
    
    
  }
  return(envv)
  
}

#' This function to update the environment list when runing the QC tSNE
#'
#' @param envv environment list associated with ShinySDA
#' @param input environment list associated with innput paramters in ShinySDA
#' @return the updated environment list
Run_tSNE_QC_envv <- function(envv, input){
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    SDAres <- envv$SDAres
    
    suffix <- paste(setdiff(1:as.numeric(SDAres$command_arguments$num_comps), envv$QC_components), collapse = "_")
    
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
    
    tsnepp <- round((nrow(SDAres$scores[,]) - 1)/3.2)
    tsnepp <- ifelse(tsnepp < 30, tsnepp, 30)
    
    tSNE_n.iter = 500
    # tsnepp
    
    
    
    # envv$InfoBox_sub = "Starting tSNE with cell scores - qc comps.. wait"
    
    envv$tsne_CS_qc = ShinySDA:::.run_tSN_scores(SDAres = SDAres,
                                                 save_path=paste0(envv$path2SDA_dyn, "/", head.path,"_tSNE_CellScore_QCfil_",tSNE_n.iter,"_pp",tsnepp,"", suffix, ".rds"), 
                                                 tsnepp=tsnepp, tSNE_n.iter=500, 
                                                 num_threads = 8, check_duplicates = F, 
                                                 selectComps = envv$QC_components )
    
    # envv$InfoBox_sub = "Saving tSNE with cell scores - qc comps.. wait"

    envv$InfoBox_sub = "qc tSNE with cell scores complete"
    
  }
  return(envv)
}




#Mcl Dev ---------


#' This function to run tsne 
#'
#' @param SDAres SDA score matrix where cols are the comps and rows are cells
#' @param save_path path to save the tsne rds obj
#' @param tsnepp tSNE preplexity
#' @param tSNE_n.iter N iter for tSNE
#' @param num_threads number of threads
#' @param check_duplicates check for dup rows
#' @param selectComps select SDA compos to run the tSNE
#' @return A list of two vector text of enrichments, one for pos one for neg scored cells
.run_tSN_scores <- function(SDAres, save_path, tsnepp, tSNE_n.iter, num_threads = 8, check_duplicates = F, selectComps = NULL){
  
  if(is.null(selectComps)) selectComps = colnames(SDAres$scores)
  
  if(!file.exists(save_path)){
    
    
    tsne_res <- Rtsne::Rtsne(SDAres$scores[,selectComps], 
                             verbose=TRUE, pca=FALSE, 
                             perplexity = tsnepp, 
                             max_iter=tSNE_n.iter, num_threads = num_threads, check_duplicates = check_duplicates)
    
    saveRDS(tsne_res, file=save_path)
    
    
    
  } else {
    tsne_res <- readRDS(save_path)
  }
  return(tsne_res)
}


.GetCompsDF <- function(MetaDF, MetaSelect, SDAScores, direction){
  CompsDF <- as.data.frame(lapply(levels(factor(MetaDF[,MetaSelect,drop = TRUE])), function(CondX){
    apply(SDAScores[rownames(MetaDF)[which(MetaDF[,MetaSelect,drop = TRUE] == CondX)], ,drop = FALSE], 2,
          function(x){
            if (direction == 'neg'){
              round(sum(x<0)/nrow(SDAScores)*100, 4)
            } else if (direction == 'pos'){
              round(sum(x>0)/nrow(SDAScores)*100, 4)
            } else {
              stop(paste0('Direction must be pos or neg: ', direction))
            }
          })
  }))
  
  colnames(CompsDF) <- levels(factor(MetaDF[,MetaSelect]))
  CompsDF[CompsDF==0] = 0.0001 #0 values produce NAs in Chisqr this is a minor adjustment
  
  CompsDF <- as.data.frame(CompsDF[naturalsort::naturalsort(rownames(CompsDF)),])
  CompsDF$SDA <- factor(rownames(CompsDF), levels=rownames(CompsDF))
  
  return(CompsDF)
}

#' Plot SDA Results Using Chi-Square Heatmap
#'
#' This function plots SDA results using a chi-square heatmap.
#'
#' @param MetaDF A data frame containing metadata for the samples.
#' @param MetaSelect A character string representing the metadata column to use for plotting.
#' @param SDAScores A data frame containing the SDA scores.
#' @return A list of two vector text of enrichments, one for pos one for neg scored cells
#' @export
SDA_ChiSqrHMresDF <- function(MetaDF, MetaSelect, SDAScores){
  
  direction = "pos"
  
  CompsDF = .GetCompsDF(MetaDF, MetaSelect, SDAScores, direction=direction)
  
  
  ChiT <- chisq.test(CompsDF[,1:(ncol(CompsDF)-1)])
  
  ChiTres <- ChiT$residuals
  ChiTres[which(is.na(ChiTres))] = 0

  ChiTres.pos = ChiTres
  ChiTres.pos[ChiTres.pos<1] = 0 
  nonzero_positions.pos <- which(ChiTres.pos != 0, arr.ind = TRUE)
  table.pos <- ChiTres.pos
  table.pos[nonzero_positions.pos] <- colnames(ChiTres.pos)[nonzero_positions.pos[,2]]
  
  
  OutVec.pos = apply(table.pos, 1, function(x) {
    y = paste(x[x != "0"], collapse = ", ") 
    ifelse(y=="", y, paste0("E4 ", y))
  })
  
  ChiTres.pos = ChiTres
  ChiTres.pos[ChiTres.pos>-1] = 0 
  nonzero_positions.pos <- which(ChiTres.pos != 0, arr.ind = TRUE)
  table.pos <- ChiTres.pos
  table.pos[nonzero_positions.pos] <- colnames(ChiTres.pos)[nonzero_positions.pos[,2]]
  
  OutVec.pos = paste0(OutVec.pos, " | ",  apply(table.pos, 1, function(x) {
    y = paste(x[x != "0"], collapse = ", ") 
    ifelse(y=="", y, paste0("D4 ", y))
  })) 
  OutVec.pos[OutVec.pos ==" | "] = ""
  
  
  
  
  
  
  direction = "neg"
  
  CompsDF = .GetCompsDF(MetaDF, MetaSelect, SDAScores, direction=direction)
  
  
  ChiT <- chisq.test(CompsDF[,1:(ncol(CompsDF)-1)])
  
  ChiTres <- ChiT$residuals
  ChiTres[which(is.na(ChiTres))] = 0
  
  ChiTres.neg = ChiTres
  ChiTres.neg[ChiTres.neg<1] = 0 
  nonzero_negitions.neg <- which(ChiTres.neg != 0, arr.ind = TRUE)
  table.neg <- ChiTres.neg
  table.neg[nonzero_negitions.neg] <- colnames(ChiTres.neg)[nonzero_negitions.neg[,2]]
  
  
  OutVec.neg = apply(table.neg, 1, function(x) {
    y = paste(x[x != "0"], collapse = ", ") 
    ifelse(y=="", y, paste0("E4 ", y))
  })
  
  ChiTres.neg = ChiTres
  ChiTres.neg[ChiTres.neg>-1] = 0 
  nonzero_negitions.neg <- which(ChiTres.neg != 0, arr.ind = TRUE)
  table.neg <- ChiTres.neg
  table.neg[nonzero_negitions.neg] <- colnames(ChiTres.neg)[nonzero_negitions.neg[,2]]
  
  OutVec.neg = paste0(OutVec.neg, " | ",  apply(table.neg, 1, function(x) {
    y = paste(x[x != "0"], collapse = ", ") 
    ifelse(y=="", y, paste0("D4 ", y))
  })) 
  OutVec.neg[OutVec.neg ==" | "] = ""
  
  
  
  
return(list(Pos=OutVec.pos,
            Neg=OutVec.neg 
            ))
}






# Plots ----------

#' Plot_CS_BoxplotMeta
#'
#' This function plots Box plots of the cell scores relative to a selected meta feature
#'
#' @param SDAScores SDA score matrix
#' @param MetaDF Metadata Dataframe
#' @param MetaSelect A factored feature of Metadata
#' @param KeepCols a numeric vector of features to keep, default NULL = all
#' @param scaleMean0 Boolean default T to scale each compoent to mean 0
#' @return A pheatmap cor heatmap
#' @export
Plot_CS_BoxplotMeta <- function(SDAScores = NULL, MetaDF, 
                                MetaSelect= "Population", KeepCols=NULL,
                                scaleMean0 = T, col_vector) {
  
  if(!is.null(KeepCols)) SDAScores = SDAScores[,KeepCols]
  
  if(scaleMean0) SDAScores = scale(SDAScores)

  
  ggls = lapply(1:ncol(SDAScores), function(ComponentN){
    ggplot(data.table(cell_index = 1:nrow(SDAScores), 
                      score = asinh((SDAScores[, paste0("SDAV", ComponentN)])^3), 
                      experiment = MetaDF$BarcodePrefix, 
                      ColFac = MetaDF[,MetaSelect]), 
           aes(ColFac, score)) + 
      geom_boxplot(outlier.size = 0) + 
      xlab("") + 
      ylab("asinh(Score^3)") + 
      #scale_color_brewer(palette = "Paired") + 
      theme_bw() + 
      theme(legend.position = "bottom", 
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      guides(colour = guide_legend(ncol = 4, override.aes = list(size = 2, alpha=1))) +
      scale_colour_manual(values =(col_vector),
                          guide = guide_legend(nrow=2)) +
      # guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
      ggtitle(paste0("SDAV", ComponentN))
  })
  
  cowplot::plot_grid(plotlist = ggls, ncol=5)
}


#' Plot_CS_AcrossMeta
#'
#' This function plots scatter pf the cell scores relative to a selected meta feature
#'
#' @param SDAScores SDA score matrix
#' @param MetaDF Metadata Dataframe
#' @param MetaSelect A factored feature of Metadata
#' @param KeepCols a numeric vector of features to keep, default NULL = all
#' @param scaleMean0 Boolean default T to scale each compnoent to mean 0
#' @return A pheatmap cor heatmap
#' @export
Plot_CS_AcrossMeta <- function(SDAScores = NULL, MetaDF, 
                                MetaSelect= "Population", KeepCols=NULL,
                                scaleMean0 = T, col_vector) {
  
  if(!is.null(KeepCols)) SDAScores = SDAScores[,KeepCols]
  
  if(scaleMean0) SDAScores = scale(SDAScores)
  

  ggls = lapply(1:ncol(SDAScores), function(ComponentN){
    
    plottingDF = data.table(
      score = asinh((SDAScores[, paste0("SDAV", ComponentN)])^3), 
      experiment = MetaDF$BarcodePrefix, 
      ColFac = MetaDF[,MetaSelect])
    
    # plottingDF = plottingDF[naturalsort::naturalorder(plottingDF$ColFac),] #slows things
    plottingDF = plottingDF[order(plottingDF$ColFac),]
    
    plottingDF$cell_index = 1:nrow(plottingDF)
    
    ggplot(plottingDF,
    aes(cell_index, score, colour = ColFac)) + 
      geom_point(size = 1, stroke = 0) + 
      xlab("Cell Index") + ylab("asinh(Score^3)") + 
      #scale_color_brewer(palette = "Paired") + 
      theme_bw() + 
      theme(legend.position = "none") + 
      # guides(colour = guide_legend(ncol = 4, override.aes = list(size = 2, alpha=1))) +
      scale_colour_manual(values =(col_vector),
                          guide = guide_legend(nrow=2)) +
      # guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
      ggtitle(paste0("SDAV", ComponentN)) + ylim(range(-10.5, 10.5))
  })
  
  cowplot::plot_grid(plotlist = ggls, ncol=5)
}



#' Plot_CS_AcrossMeta_Legend
#'
#' This function grabs the legend of Plot_CS_AcrossMeta to plot
#'
#' @param SDAScores SDA score matrix
#' @param MetaDF Metadata Dataframe
#' @param MetaSelect A factored feature of Metadata
#' @param KeepCols a numeric vector of features to keep, default NULL = all
#' @param scaleMean0 Boolean default T to scale each compnoent to mean 0
#' @return A pheatmap cor heatmap
#' @export
Plot_CS_AcrossMeta_Legend <- function(SDAScores = NULL, MetaDF, 
                               MetaSelect= "Population", KeepCols=NULL,
                               scaleMean0 = T, col_vector) {
  
  if(!is.null(KeepCols)) SDAScores = SDAScores[,KeepCols]
  
  if(scaleMean0) SDAScores = scale(SDAScores)
  
  ggls = lapply(1:1, function(ComponentN){
    ggplot(data.table(cell_index = 1:nrow(SDAScores), 
                      score = asinh((SDAScores[, paste0("SDAV", ComponentN)])^3), 
                      experiment = MetaDF$BarcodePrefix, 
                      ColFac = MetaDF[,MetaSelect]),
           aes(cell_index, score, colour = ColFac)) + 
      geom_point(size = 1, stroke = 0) + 
      # xlab("Cell Index") + ylab("asinh(Score^3)") + 
      #scale_color_brewer(palette = "Paired") + 
      # theme_bw() + 
      theme(legend.position = "bottom") + 
      guides(colour = guide_legend(ncol = 5, override.aes = list(size = 3, alpha=1))) +
      scale_colour_manual(values =(col_vector),
                          guide = guide_legend(nrow=2)) #+
      # guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) +
      # ggtitle(paste0("SDAV", ComponentN)) + ylim(range(-10, 10))
  })
  plot(cowplot::get_legend(ggls[[1]]))
  # cowplot::plot_grid(plotlist = ggls, ncol=5)
}



#' Plot_CorSDA_Loadings
#'
#' This function plots Correlatio Heatmap of the gene loadings
#'
#' @param SDAres SDA results
#' @return A pheatmap cor heatmap
#' @export
Plot_CorSDA_Loadings <- function(SDAres = NULL) {
  rownames(SDAres$loadings[[1]]) <- paste0('SDAV', 1:nrow(SDAres$loadings[[1]]))
  
  pheatmap::pheatmap(cor(t(SDAres$loadings[[1]][,])))
}





#' Plot SDA Results Using Chi-Square Heatmap
#'
#' This function plots SDA results using a chi-square heatmap.
#'
#' @param MetaDF A data frame containing metadata for the samples.
#' @param MetaSelect A character string representing the metadata column to use for plotting.
#' @param SDAScores A data frame containing the SDA scores.
#' @param direction A character string indicating the direction of the analysis (either "pos" or "neg").
#' @param clustStat A logical value indicating whether to include clustering statistics in the plot.
#' @return A pheatmap plot using a chi-square
#' @export
plot_SDA_ChiSqrHM <- function(MetaDF, MetaSelect, SDAScores, direction, clustStat = T){
  

  CompsDF <- as.data.frame(lapply(levels(factor(MetaDF[,MetaSelect,drop = TRUE])), function(CondX){
    apply(SDAScores[rownames(MetaDF)[which(MetaDF[,MetaSelect,drop = TRUE] == CondX)], ,drop = FALSE], 2,
          function(x){
            if (direction == 'neg'){
              round(sum(x<0)/nrow(SDAScores)*100, 4)
            } else if (direction == 'pos'){
              round(sum(x>0)/nrow(SDAScores)*100, 4)
            } else {
              stop(paste0('Direction must be pos or neg: ', direction))
            }
          })
  }))
  
  colnames(CompsDF) <- levels(factor(MetaDF[,MetaSelect]))
  CompsDF[CompsDF==0] = 0.0001 #0 values produce NAs in Chisqr this is a minor adjustment
  
  CompsDF <- as.data.frame(CompsDF[naturalsort::naturalsort(rownames(CompsDF)),])
  CompsDF$SDA <- factor(rownames(CompsDF), levels=rownames(CompsDF))
  
  
  
  ChiT <- chisq.test(CompsDF[,1:(ncol(CompsDF)-1)])
  
  ChiTres <- ChiT$residuals
  ChiTres[which(is.na(ChiTres))] = 0
  
  ChiResSD <- round(apply(ChiTres, 1, sd),2)
  ChiResSD[which(is.na(ChiResSD))] <- 0
  ChiResSD[ChiResSD < 0.2] <- ""
  
 
  pheatmap::pheatmap((t(ChiTres)),
                     cluster_cols = clustStat, cluster_rows = clustStat,
                     color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10), 
                     labels_col = paste0(rownames(CompsDF), " sd_", ChiResSD)
  )
  
  
}






# get_MetaDF <- function(cellbarcodes) { 
#   loupeFileIds <- unique(sapply(cellbarcodes, function(x){    
#     return(unlist(strsplit(x, split = '_'))[1]) 
#   })) 
#   
#   # Make a pseudo-seurat object: 
#   counts <- matrix(1:length(loupeFileIds), nrow = 1, dimnames = list('Gene1', paste0(loupeFileIds, '_CB'))) 
#   meta <- data.frame(BarcodePrefix = loupeFileIds, row.names = paste0(loupeFileIds, '_CB')) 
#   seuratObj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta) 
#   seuratObj <- Rdiscvr::QueryAndApplyCdnaMetadata(seuratObj) 
#   return(seuratObj@meta.data)
# }


#' go_volcano_plot
#'
#' @return a volcano-like plot
#' @export
go_volcano_plot <- function(x=GO_data, component="V5N", extraTitle=""){
  if(extraTitle=="") extraTitle = paste("Component : ", component, sep="")
  
  #print(
  ggplot(data.table(x[[component]]), aes(GeneOdds/BgOdds, -log(pvalue), size = Count)) +
    geom_point(aes(colour=p.adjust<0.05)) +
    scale_size_area() +
    geom_label_repel(data = data.table(x[[component]])[order(p.adjust)][1:30][p.adjust<0.7], aes(label = Description, size = 0.25), size = 3, force=2) + 
    ggtitle(paste("",extraTitle, sep="\n") ) +
    xlab("Odds Ratio") +
    scale_x_log10(limits=c(1,NA), breaks=c(1,2,3,4,5,6,7,8))
  #)
}


#' plotEnrich
#'
#' @return plot of enrichment
#' @export
plotEnrich <- function(GeneSetsDF, GeneVec, plotTitle="", xLab="", N = NULL, k = NULL, ReturnPval=T, BiPlot = F){
  
  
  GeneVec_overlap <- apply(GeneSetsDF, 2, function(x){
    length(which(x %in% GeneVec))
  })
  table(GeneVec_overlap)
  
  
  
  m = length(GeneVec)
  n = N - m
  
  ## Random expectation
  marked.proportion <- m / N; marked.proportion
  exp.x <- k * marked.proportion; exp.x
  
  x = GeneVec_overlap
  
  ## Fold enrichment, as computed by David
  fold.enrichment <-  (x / k ) / (m / N)
  
  # barplot(fold.enrichment, las=2)
  
  
  
  # ggplot(data=data.frame(x=names(fold.enrichment),
  #            y=fold.enrichment), aes(x=x, y=y)) + theme_bw()  +
  #   geom_bar(stat="identity") +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("Fold-enrichment \n Sertoli cell-only syndrome vs control DE genes")
  
  p.value <-  phyper(q=x-1, m=m, n=n, k=k, lower.tail=FALSE)
  # p.value <-  p.adjust(p.value, "BH")
  p.value <-  p.adjust(p.value, "fdr")
  
  fold.enrichment <- fold.enrichment[gtools::mixedsort(names(fold.enrichment))]
  
  p.value <-  p.value[names(fold.enrichment)]
  
  tDF <- data.frame(x=(fold.enrichment),
                    y=p.value)
  rownames(tDF) <- names(fold.enrichment)
  # tDF[tDF$y <0.01,]
  
  # print(head(p.value))
  # print(head(names(fold.enrichment)))
  # print(names(p.value)[p.value<0.01])
  # print(paste(names(p.value)[p.value<0.01]))
  # print(rownames(tDF[tDF$y <0.01,]))
  
  tempSigComps <- rownames(tDF[tDF$y <0.01,])
  
  tempLen <- length(tempSigComps)
  
  tempKeepLen <- tempLen - length(grep("Removed", tempSigComps))
  
  if(tempLen> 16){
    plotTitle <- paste0(plotTitle, "\nSig Comps ", tempKeepLen, "/", tempLen, ":", 
                        paste(rownames(tDF[tDF$y <0.01,])[1:7], collapse = ", "), "\n", 
                        paste(rownames(tDF[tDF$y <0.01,])[8:15], collapse = ", "), "\n", 
                        paste(rownames(tDF[tDF$y <0.01,])[16:tempLen], collapse = ", "))
    
  } else  if(tempLen> 8){
    plotTitle <- paste0(plotTitle, "\nSig Comps: ", tempKeepLen, "/", tempLen, ":",
                        paste(rownames(tDF[tDF$y <0.01,])[1:7], collapse = ", "), "\n", 
                        paste(rownames(tDF[tDF$y <0.01,])[8:tempLen], collapse = ", "))
    
  } else {
    plotTitle <- paste0(plotTitle, "\nSig Comps: ", tempKeepLen, "/", tempLen, ":",
                        paste(rownames(tDF[tDF$y <0.01,]), collapse = ", "))
    
  }
  
  if(BiPlot) print(ggplot(data=tDF, aes(x=x, y=y, label = ifelse(y < 0.01, "*", ""))) + theme_bw()  +
                     geom_point() + geom_line() +
                     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                     geom_text(vjust = 0) +
                     ggtitle(plotTitle) +
                     xlab("Fold-enrichment") +
                     ylab("1 - P-value = 1 - P(X>=x) = P(X<x)"))
  
  tDF <- data.frame(x=factor(names(fold.enrichment), levels = names(fold.enrichment)),
                    y=fold.enrichment,
                    p=p.value)
  
  rownames(tDF) <- names(fold.enrichment)
  
  
  
  print(ggplot(data=tDF, aes(x=x, y=y, label = ifelse(p < 0.01, "*", ""))) + theme_bw()  +
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          geom_text(vjust = 0) + 
          ggtitle(plotTitle) + 
          xlab(xLab) + 
          ylab("Fold-enrichment"))
  
  
  if(ReturnPval) return(1-p.value)
  
}





# Web/Prime Access ----------


#' DownloadMetadataForSdaResults
#'
#' @param sdaResultOutputId A numeric value often 6 digits long
#' @return A dataframe of metadata
#' @export
DownloadMetadataForSdaResults <- function(sdaResultOutputId) {
  fn <- tempfile()
  
  rows <- labkey.selectRows(
    baseUrl=Rdiscvr:::.getBaseUrl(),
    folderPath=Rdiscvr:::.getLabKeyDefaultFolder(),
    schemaName="sequenceanalysis",
    queryName="outputfiles",
    viewName="",
    colSort="-rowid",
    colSelect="rowid,workbook/workbookid,dataid/webdavurlrelative",
    colFilter=makeFilter(c("rowid", "EQUAL", sdaResultOutputId)),
    containerFilter=NULL,
    colNameOpt="rname"
  )
  
  if (nrow(rows) > 1) {
    stop(paste0('More than one matching file found, this should not occur.  RowId: ', sdaResultOutputId))
  } else if (nrow(rows) == 0) {
    stop(paste0('File not found. RowId: ', sdaResultOutputId))
  }
  
  wb <- rows[['workbook_workbookid']]
  if (is.na(wb) || is.null(wb)){
    wb <- ''
  }
  
  remotePath <- rows[['dataid_webdavurlrelative']]
  remotePath <- dirname(remotePath)
  
  files <- labkey.webdav.listDir(
    baseUrl=Rdiscvr:::.getBaseUrl(),
    folderPath=paste0(Rdiscvr:::.getLabKeyDefaultFolder(), wb),
    remoteFilePath = remotePath
  )
  
  fileUrl <- NULL
  for (f in files$files) {
    if (endsWith(f$text, '.seurat.meta.txt')) {
      fileUrl <- f$text
    }
  }
  
  if (is.null(fileUrl)) {
    stop('Unable to find seurat metadata for output file')
  }
  
  labkey.webdav.get(
    baseUrl=Rdiscvr:::.getBaseUrl(),
    folderPath=paste0(Rdiscvr:::.getLabKeyDefaultFolder(), wb),
    remoteFilePath = paste0(remotePath, '/', fileUrl),
    localFilePath = fn
  )
  
  df <- read.table(fn, sep = ',', header = TRUE)
  unlink(fn)
  
  return(df)
}


#' Download SDA Results
#'
#' This function downloads SDA results from a remote location and saves them to a file.
#'
#' @param outf A character string representing the file path to save the SDA results.
#' @param outFileID A numeric often 6 digits representing the file identifier for the SDA results.
#' @return the path where the file is saved is returned.
#' @export
download_SDA <- function(outf = NULL, outFileID = NULL){
  outf = sub("/$", "", as.character(outf))
  outf = paste0(outf, "/sda", ".", outFileID, ".rds")
  Rdiscvr::DownloadOutputFile(outFileID, outFile = outf, overwrite = F)
  print("sda downloaded")
  return(outf)
}

#' Download SDA Dimnames 
#'
#' This function downloads SDA Dimnames file from a remote location and saves them to a file.
#'
#' @param outf A character string representing the file path to save the SDA results.
#' @param outFileID A numeric often 6 digits representing the file identifier for the SDA results.
#' @return the path where the file is saved is returned.
#' @export
download_SDA_dimnames <- function(outf = NULL, outFileID = NULL){
  outf = sub("/$", "", as.character(outf))
  outf = paste0(outf, "/sda", ".", outFileID, ".rawData_dimnames.rds")
  Rdiscvr::DownloadOutputFile(outFileID, outFile = outf, 
                              overwrite = F, pathTranslator = function(x){
                                x <- gsub(x, pattern = '/sda', replacement = '/sdaOutput')
                                x <- gsub(x, pattern = '.rds', replacement = '')
                                x <- paste0(x, '/rawData_dimnames.rds')
                                
                                return(x)
                              })  
  print("sda dimnames downloaded")
  return(outf)
}





# Pipeline Fxs ----------

#' Preprocess SDA Results
#'
#' This function preprocesses SDA results by filtering and summarizing the data.
#'
#' @param SDAres A data frame containing the SDA results.
#' @param QuantThr A numeric value between 0 and 1 representing the quantile threshold for filtering the data.
#' @param envv An environment or list containing the data used in the SDA analysis.
#' @param TopN An integer value representing the number of top features to retain.
#' @param MetaDF A data frame containing metadata for the samples.
#' @return A list containing the filtered and summarized SDA results.
#' @export
preprocess_SDA <- function(SDAres = NULL, 
                           QuantThr = 0.95,
                           envv = NULL,
                           TopN = 150,
                           MetaDF = NULL
                           ){
  
  envv$TopN <- TopN
  
  
  # some stats
  envv$MaxScore.thr <- quantile(SDAres$component_statistics$max_score, c(QuantThr))
  envv$QC_components  <- SDAres$component_statistics[SDAres$component_statistics$max_score<envv$MaxScore.thr,]
  envv$QC_components  <- envv$QC_components[order(envv$QC_components$Component), ]$Component
  
  envv$QC_compIter <- min(as.numeric(envv$QC_components))
  
  envv$NComps <- as.numeric(SDAres$command_arguments$num_comps)
  
  
  
  #top loaded genes
  SDA_TopNpos <- (as.data.frame(lapply(1:envv$NComps, function(xN){
    as.data.frame(print_gene_list(results = SDAres, i=xN, PosOnly = T, NegOnly = F, TopN = TopN))[1:TopN,1]
  })))
  colnames(SDA_TopNpos) <- paste0("SDAV", 1:envv$NComps)
  
  SDA_TopNneg <- (as.data.frame(lapply(1:envv$NComps, function(xN){
    as.data.frame(print_gene_list(results = SDAres, i=xN, PosOnly = F, NegOnly = T, TopN = TopN))[1:TopN,1]
  })))
  colnames(SDA_TopNneg) <- paste0("SDAV", 1:envv$NComps)
  
  envv$SDA_TopNpos <- SDA_TopNpos
  envv$SDA_TopNneg <- SDA_TopNneg
  
  
  if(!is.null(MetaDF)){
    envv$InfoBox_sub = "SDA & Meta Loaded"
  } else {
    envv$InfoBox_sub = "SDA Meta not found"
  }
  
  
  return(envv)
  
}




#' Add Component Statistics to SDA object
#'
#' This function adds component statistics to SDA results.
#'
#' @param SDAres A data object of SDA results.
#' @param sdThrsh A numeric value representing the standard deviation threshold for filtering the data.
#' @param maxscoreThrsh An integer value representing the maximum score threshold for filtering the data.
#' @param maxloadThrsh A numeric value representing the maximum load threshold for filtering the data.
#' @param redoCalc A logical value indicating whether to recalculate the component statistics.
#' @return A data object of SDA results with added component statistics.
#' @export
AddCompStats <- function(SDAres, sdThrsh = 0.04, maxscoreThrsh=20, maxloadThrsh = 1, redoCalc=T){
  if (redoCalc)   SDAres$component_statistics <- NULL
  if (is.null(SDAres$component_statistics) ) {
    
    SDAres$component_statistics <- as.data.frame(data.table(
      Component = 1:SDAres$n$components, 
      Component_name = dimnames(SDAres$scores)[[2]], 
      max_score = apply(abs(SDAres$scores),  2, max),
      max_loading = apply(abs(SDAres$loadings[[1]]), 1, max),
      mean_score = apply(abs(SDAres$scores),  2, mean),
      mean_loading = apply(abs(SDAres$loadings[[1]]), 1, mean),
      sd_score = apply(abs(SDAres$scores),  2, sd),
      sd_loading = apply(abs(SDAres$loadings[[1]]), 1, sd),
      ssqrd_score = apply(SDAres$scores,  2, function(x) sum(x^2)),
      ssqrd_loading = apply(SDAres$loadings[[1]], 1, function(x) sum(x^2))
    )[order(-Component)])
    
  }
  
  SDAres$component_statistics$Component_namev2    <- SDAres$component_statistics$Component_name
  SDAres$component_statistics$Component_name_plot <- SDAres$component_statistics$Component_name
  
  SDAres$component_statistics[which(SDAres$component_statistics$sd_loading<sdThrsh),]$Component_namev2         <- rep("", length(which(SDAres$component_statistics$sd_loading<sdThrsh)))
  
  # SDAres$component_statistics[which(SDAres$component_statistics$max_score<maxscoreThrsh & SDAres$component_statistics$max_loading<maxloadThrsh),]$Component_name_plot <- rep("", length(which(SDAres$component_statistics$max_score<maxscoreThrsh & SDAres$component_statistics$max_loading<maxloadThrsh)))
  
  SDAres$component_statistics <- data.table(SDAres$component_statistics)
  return(SDAres)
  
}

#' convert mouse to human genes
#'
#'
#' @param x A vector of mouse genes
#' @return A vector of human genes
#' @export
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  #humanx <- unique(genesV2[, 2])
  humanx <- data.frame(genesV2)
  #print(head(humanx))
  return(humanx)
}

#' print_loadings_scores
#'
#' @return plots of gene loadings 
#' @export
print_loadings_scores <- function(SDAResult = NA, ComponentN=NA, ColFac = NA, Prefix="SDAV", GeneLoc=NA){
  library(ggthemes)
  library(scales)
  if(!is.factor(ColFac)) ColFac <- factor(ColFac)
  
  SDAScores    <- SDAResult$scores
  SDALoadings <- SDAResult$loadings[[1]]
  
  
  #cowplot::plot_grid(
  g1 <- ggplot(data.table(cell_index = 1:nrow(SDAScores), 
                          score = SDAScores[, paste0(Prefix, ComponentN)], experiment = gsub("_.*", 
                                                                                             "", gsub("[A-Z]+\\.", "", rownames(SDAScores))), ColFac = ColFac), 
               aes(cell_index, score, colour = ColFac)) + 
    geom_point(size = 1, stroke = 0) + 
    xlab("Cell Index") + ylab("Score") + 
    #scale_color_brewer(palette = "Paired") + 
    
    
    theme_bw() + 
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(override.aes = list(size = 2, alpha=1))) +
    scale_colour_manual(values = colorRampPalette(solarized_pal()(8))(length(levels(ColFac))),
                        guide = guide_legend(nrow=2)) +
    #guides(color = guide_legend(ncol = 2, override.aes = list(size = 2))) + 
    ggtitle(paste0(Prefix, ComponentN)) 
  #,
  
  g2 <- genome_loadings(SDALoadings[ComponentN,], 
                        label_both = T, 
                        max.items = 10, 
                        gene_locations =   GeneLoc)
  #, ncol = 1)
  print(g1)
  print(g2)
}

#' get.location
#'
#' @return tabulated gene map
#' @export
get.location <- function(gene.symbols, data_set, gene_name){
  
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", dataset = data_set)
  
  mapTab <- getBM(attributes = c(gene_name,'chromosome_name','start_position'),
                  filters = gene_name, values = gene.symbols, mart = ensembl, uniqueRows=TRUE)
  
  mapTab <- as.data.table(mapTab)
  
  setnames(mapTab, gene_name,"gene_symbol")
  
  # Remove duplicate genes!!
  # first which genes are duplicated
  duplicates <- mapTab[duplicated(mapTab, by="gene_symbol")]$gene_symbol
  # order by chr so patch versions to go bottom, then choose first unique by name
  unduplicated <- unique(mapTab[gene_symbol %in% duplicates][order(chromosome_name)], by="gene_symbol")
  
  # remove duplicates and replace with unique version
  mapTab <- mapTab[!gene_symbol %in% duplicates]
  mapTab <- rbind(mapTab, unduplicated)
  
  # change all patch chr names to NA
  mapTab[!chromosome_name %in% c(1:22,"X","Y","MT")]$chromosome_name <- NA # mice actually 19
  mapTab[is.na(chromosome_name)]$start_position <- NA
  return(mapTab)
}

#' GO_enrichment
#'
#' @return GO results
#' @export
GO_enrichment <- function(results = NULL, component, geneNumber = 100, threshold=0.01, side="N", OrgDb = NULL){
  # results = SDAres; component = 1; geneNumber = 100; threshold=0.01; side="N"; OrgDb = RefGenome
  require(data.table)
  if(side=="N"){
    top_genes <- data.table(as.matrix(results$loadings[[1]][component, ]), keep.rownames = TRUE)[order(V1)][1:geneNumber]$rn
  }else{
    top_genes <- data.table(as.matrix(results$loadings[[1]][component, ]), keep.rownames = TRUE)[order(-V1)][1:geneNumber]$rn
  }
  
  gene_universe <- data.table(as.matrix(results$loadings[[1]][component,]), keep.rownames = TRUE)$rn
  print(head(top_genes))
  ego <- enrichGO(gene = top_genes,
                  universe = gene_universe,
                  OrgDb = OrgDb,
                  keyType = 'SYMBOL',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)
  
  ego@result$Enrichment <- frac_to_numeric(ego@result$GeneRatio)/frac_to_numeric(ego@result$BgRatio)
  
  ego@result$GeneOdds <- unlist(lapply(strsplit(ego@result$GeneRatio, "/", fixed = TRUE), function(x){ x<-as.numeric(x) ;x[1] / (x[2]-x[1])}))
  ego@result$BgOdds <- unlist(lapply(strsplit(ego@result$BgRatio, "/", fixed = TRUE), function(x){ x<-as.numeric(x) ;x[1] / (x[2]-x[1])}))
  
  return(ego@result)
  
}

#' frac_to_numeric
#'
#' @return  converted frac
#' @export
frac_to_numeric <- function(x) sapply(x, function(x) eval(parse(text=x)))



#' print_gene_list
#'
#' @return ordered gene list
#' @export
print_gene_list <- function(results = results, i, PrintRes = F, PosOnly = F, NegOnly = F, AbsLoad = T, TopN = 150) {
  # results = SDAres; i=1; PosOnly = T; NegOnly = F; TopN = TopN
  
  if(AbsLoad)  tmp <- data.table(as.matrix(results$loadings[[1]][i,]), keep.rownames = TRUE)[order(-abs(V1))][1:TopN]
  
  if(PosOnly)  tmp <- data.table(as.matrix(results$loadings[[1]][i,]), keep.rownames = TRUE)[order(-(V1))][1:TopN]
  
  if(NegOnly)  tmp <- data.table(as.matrix(results$loadings[[1]][i,]), keep.rownames = TRUE)[order((V1))][1:TopN]
  
  
  setnames(tmp, c("Gene.Name","Loading"))
  setkey(tmp, Gene.Name)
  
  
  # Display Result
  if(PrintRes) print(tmp[order(-abs(Loading))]) else return(tmp[order(-abs(Loading))])
}


#' findIntRuns
#'
#' @return numeric
#' @export
findIntRuns <- function(run){
  rundiff <- c(1, diff(run))
  difflist <- split(run, cumsum(rundiff!=1))
  unlist(lapply(difflist, function(x){
    if(length(x) %in% 1:2) as.character(x) else paste0(x[1], "-", x[length(x)])
  }), use.names=FALSE)
}

