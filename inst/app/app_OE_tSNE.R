## run raw tSNE
observeEvent(input$runtSNE, {
  
  envv$InfoBox_sub = "Starting tSNE with cell scores - all comps.. wait"
  
  envv = ShinySDA:::Run_tSNE_full_evv(envv, input = input)
  
  
})

##Qc tsne
observeEvent(input$runtSNEQCfilt, {
  
  envv$InfoBox_sub = "Starting tSNE with cell scores - qc comps.. wait"
  
  envv = ShinySDA:::Run_tSNE_QC_envv(envv, input = input)
  
  
  
})