
observeEvent(input$Download_SDA_primeseq, {
  
  envv$InfoBox_sub = "Downloading from Prime-seq"
  envv$Origin = "prime"
  envv = ShinySDA:::Run_DowLoadSDA_evv(envv, input)
  
  
})

observeEvent(input$Down_N_Load_SDA, {
  
  envv$InfoBox_sub = "Downloading & Loading from Prime-seq"
  envv$Origin = "prime"
  envv = ShinySDA:::Run_DowLoadSDA_evv(envv, input)
  envv = ShinySDA:::Run_LoadSDA_evv(envv, input, session)
  
  
})


observeEvent(input$Load_SDA_primeseq, {

  envv$InfoBox_sub = "Loading from Prime-seq download"
  envv$Origin = "prime"
  envv = ShinySDA:::Run_LoadSDA_evv(envv, input, session)
  
  
 
  
  
})



observeEvent(input$loadSDA, {
  
  envv$Origin = "folder"
  envv$InfoBox_sub = "Loading from SDA folder (traditional)"
  
  envv = ShinySDA:::Run_LoadSDA_local_evv(envv, input)
  
  
})


