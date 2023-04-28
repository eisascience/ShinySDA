# general OE


observeEvent(input$runAllProc, {
  
  envv$InfoBox_sub = "Getting gene annotations"
  envv = ShinySDA:::Run_GeneAnn_evv(envv, input = input)
  envv$InfoBox_sub = "Starting tSNE with cell scores - all comps.. wait"
  envv = ShinySDA:::Run_tSNE_full_evv(envv, input = input)
  # envv$InfoBox_sub = "Starting tSNE with cell scores - qc comps.. wait"
  # envv = ShinySDA:::Run_tSNE_QC_envv(envv, input = input)
  envv$InfoBox_sub <- "Stating GO"
  envv = ShinySDA:::Run_GO_evv(envv, input = input)
  
})


# load in OE
# observeEvent(input$LoadFolderInput, {
# 
#   envv$Origin = "folder"
#   envv$InfoBox_sub = "load from input"
#   # updateTabItems(session, "page", selected = "PrimeSeqInput")
#  
# })
# 
# observeEvent(input$LoadPrimeSeqInput, {
#   
#   envv$Origin = "prime"
#   envv$InfoBox_sub = "load from input"
#   # updateTabItems(session, "page", selected = "FolderInput")
#   
#   
# })