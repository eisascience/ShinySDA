observeEvent(input$pos_score_save, {
  
})

observeEvent(input$neg_score_save, {
  
})

observeEvent(input$upd_score_tabl, {
  
  SDAScores <- envv$SDAres$scores
  MetaDF <- envv$MetaDF
  
  if(length(input$chkbx_MetaSave)>1 ){
    Pos.ls = list()
    Neg.ls = list()
    for(x in input$chkbx_MetaSave){
      StrLs = ShinySDA:::SDA_ChiSqrHMresDF(MetaDF, MetaSelect=x, SDAScores)
      Pos.ls[[x]] = StrLs$Pos
      Neg.ls[[x]] = StrLs$Neg
    }
    # names(Pos.ls) = paste0("Pos_",names(Pos.ls))
    # names(Neg.ls) = paste0("Neg_",names(Neg.ls))
    
    # renderDataTable(mtcars$cyl*10)
   
    # DT::replaceData(session, "upd_score_tabl",  data.frame(cbind(as.data.frame(Pos.ls), as.data.frame(Neg.ls))),
    #             resetPaging = TRUE)
    
    paste_if_not_blank <- function(df){
      # Create an empty vector to store the output
      output_vector <- character(nrow(df))
      
      # Apply the function to each row of the dataframe
      output_vector <- apply(df, 1, function(row){
        # Paste the values of the row if not blank
        return(paste(row[row != ""], collapse = " | "))
      })
      
      return(output_vector)
    }
    

    # envv$ChiSqrMeta_tabDF = data.frame(cbind(as.data.frame(Pos.ls), as.data.frame(Neg.ls)))
    envv$ChiSqrMeta_tabDF = data.frame(Pos=paste_if_not_blank(as.data.frame(Pos.ls)), 
               Neg=paste_if_not_blank(as.data.frame(Neg.ls)))
    
    # print(data.frame(Pos=paste_if_not_blank(as.data.frame(Pos.ls)), 
    #                  Neg=paste_if_not_blank(as.data.frame(Neg.ls))))
    # print(head(data.frame(cbind(as.data.frame(Pos.ls), as.data.frame(Neg.ls)))))
    
  }
  

})




observeEvent(input$SaveAsSerObj, {
  
  
  if(is.null(envv$SDAres)){ 
    envv$InfoBox_sub = "Load SDA"
  } else {
    
    head.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][2]
    base.path <- stringr::str_split(envv$path2SDA_dyn, "sda_results/")[[1]][1]
    
    head.path <- gsub("/", "", head.path)
    
    
    SDAres <- envv$SDAres
    
    
    
    # print(names(envv))
    library(Seurat)
    
    ## create an empty Seurat object
    Mat1 <- abs(Matrix::rsparsematrix(20, nrow(SDAres$scores), density = .8))
    
    
    colnames(Mat1) <- rownames(SDAres$scores)
    rownames(Mat1) <- paste0("Gene", 1:nrow(Mat1))
    
    
    
    
    Ser1 <- CreateSeuratObject(Mat1)
    Ser1 <- NormalizeData(Ser1) #needed for FindVariableFeatures
    Ser1 <- FindVariableFeatures(Ser1) #needed for ScaleData
    Ser1 <- ScaleData(Ser1, features = rownames(x = Ser1)) #needed for RunPCA
    Ser1 <- RunPCA(Ser1, npcs = 5, verbose = T) # creates @ reduction object
    
    
    
    reduction.data <- CreateDimReducObject(
      embeddings = (SDAres$scores),
      loadings = t(SDAres$loadings[[1]]),
      assay = "RNA",
      stdev = apply(SDAres$scores, 2, sd),
      key = "SDA_",
      misc = list(command.args = SDAres$command_arguments,
                  n.cells = SDAres$n,
                  pip.mat = SDAres$pips[[1]],
                  pip.frac = SDAres$pip_fraction,
                  iterations = seq_len(ncol(SDAres$free_energy)) * as.numeric(SDAres$command$free_freq),
                  free.energu = SDAres$free_energy[1, ])
    )
    
    # Ser0@reductions[["SDA"]] <- reduction.data
    
    Ser1@reductions[["SDA"]] <- reduction.data
    
    rownames(envv$tsne_CS_batch$Y) <- rownames((SDAres$scores))
    
    tsne.reduction_CS_BR <- CreateDimReducObject(
      embeddings = (envv$tsne_CS_batch$Y),
      key = "tSNECSBR_",
      assay = "RNA"
    )
    
    Ser1@reductions[["tSNECSBR"]] <- tsne.reduction_CS_BR
  }
  
  Ser1@misc$SDA_processing_results <- list(GOAnn = envv$GOAnn,
                                           GO_data = envv$GO_data,
                                           MaxScore.thr = envv$MaxScore.thr,
                                           path2SDA_dyn = envv$path2SDA_dyn,
                                           SDA_TopNneg = envv$SDA_TopNneg,
                                           TopN = envv$TopN,
                                           SDA_TopNpos = envv$SDA_TopNpos,
                                           Remove_comps = envv$Remove_comps,
                                           QC_components = envv$QC_components,
                                           tsne_CS_raw = envv$tsne_CS_all,
                                           MetaDF = envv$MetaDF,
                                           chromosome.lengths = envv$chromosome.lengths,
                                           gene_locations = envv$gene_locations)
  
  
  print("Saving...")
  saveRDS(Ser1, file=paste0(envv$path2SDA_dyn, "/", head.path,"_FinalSerObj", ".rds"))
  print("Save complete...")
  
  
  #   1] "MaxScore.thr"  "y"             "GOAnn"         "path2SDA_dyn"  "SDAres"        "SDA_TopNneg"   "TopN"         
  # [8] "tsne_CS_all"   "tsne_CS_qc"    "Remove_comps"  "tsne_CS_batch" "QC_components" "QC_compIter"   "SDA_TopNpos"  
  # [15] "InfoBox_sub"   "MetaDF"        "GO_data"  
})