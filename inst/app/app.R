
# options(repos=structure(BiocManager::repositories()))
# detach("package:Biobase", unload=TRUE)

library(shiny)
library(shinyWidgets)
library(shinydashboard)
# library(shinyjs)
library(shinyFiles)

library(ggforce)
library(ggplot2)
library(ggrepel)
library(viridis)
library(RColorBrewer)

library(grid)
library(gridExtra) 

library(data.table)
library(dplyr)

library(Seurat)

library(SDAtools)

library(Rlabkey)
library(Rdiscvr)

# library(roxygen2)

library(rclipboard)


library("BiocParallel")
register(MulticoreParam(4))

library(ShinySDA)

library(ggthemes)
library(scales)

library(biomaRt)

library(AnnotationHub) # source("https://bioconductor.org/biocLite.R"); biocLite("AnnotationHub")

library(clusterProfiler) # source("https://bioconductor.org/biocLite.R"); biocLite("clusterProfiler")

library(Matrix)


# if (Sys.getenv("SCRATCH_DIR") != "") {
#   cachedir <- paste0(Sys.getenv("SCRATCH_DIR"), "ShinySDA-cache")
# } else {
#   cachedir <- paste0(Sys.getenv("TMPDIR"), "ShinySDA-cache")
# }

if (Sys.getenv("SCRATCH_DIR") != "") {
  init.path = paste0(Sys.getenv("SCRATCH_DIR"), "/data/ShinySDA")
  serv = "monkeydo"
}  else {
  init.path = getwd()
  if(grepl("Maggie", init.path)) init.path = "/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/sda/primeseq/"
  serv = "local"
}


# source(system.file('app/fxs.R', package = 'ShinySDA', mustWork = TRUE), local = TRUE)

# print(Sys.getenv("SCRATCH_DIR"))

ui <- dashboardPage(skin="red",
                    dashboardHeader(title = "ShinySDA"),
                    #https://rstudio.github.io/shinydashboard/appearance.html#icons
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem("Main", tabName = "MainDash", icon = icon("dashboard"), selected = T),
                        menuItem("Prime-seq Input", tabName = "PrimeSeqInput", icon = icon("dashboard")),
                        menuItem("Folder Input", tabName = "FolderInput", icon = icon("dashboard")),
                        menuItem("PreProcessing", tabName = "PreProcess", icon = icon("dashboard")),
                        menuItem("Run QC", tabName = "QCplots", icon = icon("wrench")),
                        menuItem("Cell-Score ChiSqr Heatmap w/ Meta", tabName = "ChiSqrHM", icon = icon("wrench")),
                        menuItem("Cell-Score ChiSqr Export w/ Meta", tabName = "ChiSqrTab", icon = icon("wrench")),
                        menuItem("Cell-score across w/ Meta", tabName = "CSplots", icon = icon("wrench")),
                        menuItem("Cell-score boxplots w/ Meta", tabName = "CSBoxPlots", icon = icon("wrench")),
                        menuItem("Cell-score paired scatter w/ Meta", tabName = "CSScatter", icon = icon("wrench")),
                        # menuItem("Cell-score tSNE", tabName = "CStSNEPlots", icon = icon("wrench")),
                        menuItem("Gene-loading Cor HM", tabName = "GL_cor_HM", icon = icon("wrench")),
                        menuItem("Batch removal", tabName = "BatchRemove", icon = icon("toolbox")),
                        menuItem("DGE Batch-Removed", tabName = "DGEsda", icon = icon("autoprefixer")),
                        menuItem("Gene Explorer", tabName = "GeneExplorer", icon = icon("dna")),
                        menuItem("Save Out", tabName = "SaveOut", icon = icon("save")),
                        menuItem("@eisamahyari", icon = icon("heart"), 
                                 href = "https://eisascience.github.io")
                      )
                    ),
                    
                    dashboardBody(
                      # useShinyjs(),
                      tags$head(
                        tags$style(HTML("
                                        .content-wrapper {
                                        background-color: black !important;
                                        }
                                        .main-sidebar {
                                        background-color: black !important;
                                        }
                                        .multicol .shiny-options-group{
                                        -webkit-column-count: 5; /* Chrome, Safari, Opera */
                                        -moz-column-count: 5;    /* Firefox */
                                        column-count: 5;
                                        -moz-column-fill: balanced;
                                        -column-fill: balanced;
                                        }
                                        .checkbox{
                                        margin-top: 0px !important;
                                        -webkit-margin-after: 0px !important; 
                                        }
                                        "))),
                      tabItems(
                        
                        
                        # Main---------------
                        tabItem(tabName = "MainDash",
                                h2("Main Dashboard"),
                                fluidRow(
                                  valueBoxOutput("InfoBox_Main", width = 6),
                                  
                                  box(textInput("apiKey", "Prime-seq API Kit", 
                                                value =""),
                                      textInput("baseURL", "baseURL", 
                                                value ="https://prime-seq.ohsu.edu"),
                                      textInput("defaultFolder", "defPSfold", 
                                                value ="Labs/Bimber"),
                                      width = 5, background = "olive"
                                      
                                  ))),
                        
                        # Prime-seq input---------
                        tabItem(tabName = "PrimeSeqInput",
                                h2("Download from Prime-seq: Use the OutputFileId for an SDA object"),
                                fluidRow(
                                  
                                  valueBoxOutput("InfoBox_Prime", width = 6),
                                  
                                    box(textInput("SDA_OFId", "OutputFileId to an SDA Obj", 
                                                  value ="506842"),
                                        textInput("SDA_PS_save", "path to download and SDA objects", 
                                                  value =init.path),
                                        column(
                                          width = 7,
                                          style = "float: left;",
                                          actionButton("Download_SDA_primeseq", "Download SDA Obj"),
                                          actionButton("Load_SDA_primeseq", "Load SDA Obj"),
                                          actionButton("Down_N_Load_SDA", "** Download & Load **"),
                                          
                                          
                                        ),
                                        # column(
                                        #   width = 3,
                                        #   style = "float: right;",
                                        #   actionButton("DnL_SDA_primeseq", "Do Both", align="right"),
                                        # ),
                                    width = 10, background = "olive"
  
                                ))),
       
                        # Folder input------
                        tabItem(tabName = "FolderInput",
                                h2("Folder input: Make sure _MetaDF.rds includes $SubjectId, $ExpID, $EXP.ID, $SampleDate, $SingleR_Labels, $BarcodePrefix"),
                                fluidRow(
                                  valueBoxOutput("InfoBox_Folder", width = 6),
                                  
                                  box(textInput("SDAroot", "Path to SDA folders. Expects dimnames in one dir up.", 
                                                value =init.path),
                                      uiOutput("select.folder"),
                                      actionButton("loadSDA", "0. Load SDA"),
                                      width = 10, background = "teal"
                                      
                                  ))),
                        
                        # Pre-Process---------
                        tabItem(tabName = "PreProcess",
                                h2("Pre-Process the SDA object"),
                                fluidRow(
                                  
                                  valueBoxOutput("InfoBox_PP", width = 6),
                                  
                                  box(title = "SpeciesSelect", status = "primary", 
                                      solidHeader = TRUE,collapsible = TRUE,
                                      selectInput(inputId = 'species',
                                                  label = 'Species',
                                                  choices =  c('human', 'mouse', 'rhesus'),
                                                  selected = , multiple = F),
                                      width = 3, background = "red"
                                  ),
                                  
                                  box(
                                    # textInput("SDAroot", "Path to SDA folders. Expects dimnames in one dir up.",
                                    #             value =init.path),
                                  #     uiOutput("select.folder"),
                                      # actionButton("loadSDA", "0. Load SDA"),
                                      actionButton("getGeneAnn", "1. Get Gene Annotations"),
                                      actionButton("getSDAGo", "2. Get SDA GO Enrichments"),
                                      actionButton("runtSNE", "3. Run tSNE (cs-all)"),
                                      actionButton("runtSNEQCfilt", "4.5 Run tSNE (cs-qc)"),
                                      actionButton("runAllProc", "** Run all **"),
                                      
                                      # textInput("loadSDAmsg", "File Status", "not loaded"),
                                      width = 10
                                      #fileInput("SDAin", "Browse")
                                  ),
                                  
                                  box(
                                    title = "QC_MaxScore_filt", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqcMaxScorefilt"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "tSNE CS All", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("tSNE_CS_all"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "tSNE_CS_QC", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("tSNE_CS_qc"), 
                                    width = 5, background = "black"
                                  ))),
                        
                        # QC plots-------
                        tabItem(tabName = "QCplots",
                                h2("Full QC plots content"),
                                fluidRow(
                                  
                                  box(
                                    title = "Kurtosis", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("KurtosisPlot"), 
                                    width = 10, background = "black"
                                  ),
                                  box(
                                    title = "QC1", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc1"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "QC2", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc2"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "QC3", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc3"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "QC4", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc4"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "QC5", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc5"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "QC6", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAqc6"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "Convergence", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("convergence"), 
                                    width = 5, background = "black"
                                  ) ,
                                  box(
                                    title = "Gene Loading Hist.", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("loadhist"), 
                                    width = 5, background = "black"
                                  ) ,
                                  box(
                                    title = "Score Dist.", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("scoredist"), 
                                    width = 5, background = "black"
                                  ) ,
                                  #The PIP is the mean of the posterior. You can think of it as a measure of how likely it is that this variable is in the true model.
                                  box(
                                    title = "Post. inclus. prob. dist", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("pipdist"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "Slack-Slab Prior", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("slackslabprior"), 
                                    width = 10, background = "black"
                                  )
                                  
                                  
                     
                                  
                                )
                        ),
                        
                        # Cell-score Boxplot-------
                        tabItem(tabName = "CSBoxPlots",
                                h2("Cell Score Box Plots, scaled to mean = 0"),
                                fluidRow(
                                  box(
                                    title = "Metadata Selection", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE, background = "black",
                                    selectInput("Metaselect5", "Meta select:",
                                                c("Population" = "Population",
                                                  "SampleDate" = "SampleDate",
                                                  "SubjectId" = "SubjectId",
                                                  "Phase" = "Phase"), selected = "Population"),
                                  ),
                                  box(
                                    title = "Scores order by Meta", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAScoresBoxplot", height = 1500, width = 900),
                                    width = 10,
                                    background = "black"
                                  )
                                )
                        ),
                        
                        # Cell-score paired-------
                        tabItem(tabName = "CSScatter",
                                h2("Cell Score Box Plots, scaled to mean = 0"),
                                fluidRow(
                                  
                                )
                        ),

                        # # Cell-score Feature plot-------
                        # tabItem(tabName = "CStSNEPlots",
                        #         h2("Full Cell Score Feature Plots"),
                        #         fluidRow(
                        # 
                        #           # box(
                        #           #   title = "Scores order by Meta", status = "primary", solidHeader = TRUE,
                        #           #   collapsible = TRUE,
                        #           #   plotOutput("SDAScoresBoxplot", height = 1500, width = 900),
                        #           #   width = 10,
                        #           #   background = "black"
                        #           # )
                        #         )
                        # ),

                        # Cell Score Across plots-------
                        tabItem(tabName = "CSplots",
                                h2("Full Cell Score Feature Plots"),
                                fluidRow(
                                  box(
                                    title = "Metadata Selection", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE, background = "black",
                                    selectInput("Metaselect6", "Meta select:",
                                                c("Population" = "Population",
                                                  "SampleDate" = "SampleDate",
                                                  "SubjectId" = "SubjectId",
                                                  "Phase" = "Phase"), selected = "Population"),
                                  ),
                                  box(
                                    title = "Scores order by Meta", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAScoresAcross_Legend"),
                                    width = 5,
                                    background = "black"
                                  ),
                                  box(
                                    title = "Scores order by Meta", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAScoresAcross", height = 1500, width = 900),
                                    width = 10,
                                    background = "black"
                                  )
                                )
                        ),
                        
                        # ChiSqrTab ------
                        tabItem(tabName = "ChiSqrTab",
                                h2("ChiSqr Heatmap Enrichment of Cell-scores per Meta"),
                                fluidRow(
                                 
                                  box(
                                    title = "Select to save", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE, background = "black",
                                    checkboxGroupInput("chkbx_MetaSave", "Select items:", 
                                                       c("Item 1", "Item 2", "Item 3")),
                                    # actionButton("pos_score_save", "Save Pos Score"),
                                    # actionButton("neg_score_save", "Save Neg Score"),
                                    actionButton("upd_score_tabl", "Update Tables"),
                                    
                                  ),
                                  box(
                                    title = "ChiSqrRes Scores Pos cellscores", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    DT::dataTableOutput("table_SDAScoresChiPos"),
                                    width = 10, background = "light-blue"
                                  ),
                                  
                                )
                        ),



                        # ChiSqrHM --------
                        tabItem(tabName = "ChiSqrHM",
                                h2("ChiSqr Heatmap Enrichment of Cell-scores per Meta"),
                                fluidRow(
                                  box(
                                    title = "Select to plot", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE, background = "black",
                                    selectInput("Metaselect4", "Meta select:",
                                                c("Population" = "Population",
                                                  "SampleDate" = "SampleDate",
                                                  "SubjectId" = "SubjectId",
                                                  "Phase" = "Phase"), selected = "Population"),
                                  ),
                                  
                                  box(
                                    title = "ChiSqrRes Scores Pos cellscores", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAScoresChiPos", height = 500),
                                    width = 10, background = "black"
                                  ),
                                  box(
                                    title = "ChiSqrRes Scores Neg cellscores", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("SDAScoresChiNeg", height = 500),
                                    width = 10, background = "black"
                                  ),
                                )
                        ),

                        # Gene-loadig Cor HM-------
                        tabItem(tabName = "GL_cor_HM",
                                h2("Gene Loadings Correlation Heatmap"),
                                fluidRow(
                                  box(
                                    title = "Cor. HM", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GL_Cor_HMplot", height = 500),
                                    width = 10, background = "black"
                                  ),
                                )
                        ),
                        
                        # Batch removal content----------
                        tabItem(tabName = "BatchRemove",
                                h2("Batch Removal"),
                                fluidRow(
                                  box(
                                    title = "SDA projected on tSNE", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    actionButton("prevSDA_br", "Prev comp"),
                                    actionButton("nextSDA_br", "Next comp"),
                                    plotOutput("SDAtsne_br1"), 
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "Batch Removal Selection", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE, 
                                    column(10,
                                           uiOutput("CompBatchCheckBoxSelect")
                                    ),
                                    textInput("CompSelecTextIn", "Paste fail/pass", ""),
                                    actionButton("save_batch_selection", "Save Selection"),
                                    actionButton("load_batch_selection", "Load last Selection"),
                                    actionButton("reset_batch_selection", "Reset last Selection"),
                                    actionButton("select_all_selection", "Select all"),
                                    width=5, background = "black"),
                                  box(
                                    title = "Run tSNE Batch Removed", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE, 
                                    
                                    selectInput("tSNEiter", "tSNE n-iter:",
                                                c("Fast-500" = "500",
                                                  "Fast2-1000" = "1000",
                                                  "Med1-2000" = "2000",
                                                  "Med2-5000" = "5000",
                                                  "Robust-10000" = "10000",
                                                  "OverKill-20000"="20000"), selected = "500"),
                                    selectInput("tSNEpp", "tSNE prplx:",
                                                c("1" = "1",
                                                  "10" = "10",
                                                  "50-default" = "50",
                                                  "100-med" = "100",
                                                  "200-max" = "200",
                                                  "300-insane!"  = "300"), selected = "50"),
                                    actionButton("run_tSNE_CS_batch", "Run tSNE (batch-removed)"),
                                    actionButton("SDAScoresChi_clus", "Show/Hide Pairwise Clustering"),
                                    width = 5, background = "black",
                                  ),
                                  
                                  
                                  
                                  box(
                                    title = "tSNE Batch removed", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("DGE_SDA_tSNE"), 
                                    width = 10, background = "black"
                                  )
                                )
                        ),
                        
                        # DGE-SDA-BatchRemove tab content--------
                        tabItem(tabName = "DGEsda",
                                h2("DGE_SDA Batched Removed DGE"),
                                fluidRow(
                                  box(
                                    title = "tSNE Meta Exploration", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    
                                    plotOutput("tSNE_CS_batch1"),
                                    selectInput("Metaselect1", "Meta select:",
                                                c("Population" = "Population",
                                                  "SampleDate" = "SampleDate",
                                                  "SubjectId" = "SubjectId",
                                                  "Phase" = "Phase"), selected = "Population"),
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "tSNE Meta Exploration", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    
                                    plotOutput("tSNE_CS_batch2"),
                                    selectInput("Metaselect2", "Meta select:",
                                                c("Population" = "Population",
                                                  "SampleDate" = "SampleDate",
                                                  "SubjectId" = "SubjectId",
                                                  "Phase" = "Phase"), selected = "Population"),
                                    width = 5, background = "black"
                                  ),
                                  box(
                                    title = "tSNE Meta Exploration", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    
                                    plotOutput("tSNE_CS_batch3"),
                                    selectInput("Metaselect3", "Meta select:",
                                                c("Population" = "Population",
                                                  "SampleDate" = "SampleDate",
                                                  "SubjectId" = "SubjectId",
                                                  "Phase" = "Phase"), selected = "Population"),
                                    width = 5, background = "black"
                                  ),
                                  box(title = "tSNE SDA projection", status = "primary", solidHeader = TRUE,
                                      collapsible = TRUE,
                                      actionButton("prevSDA_br2", "Prev comp"),
                                      actionButton("nextSDA_br2", "Next comp"),
                                      actionButton("C2Cpos", "Copy2ClipPosGenes"),
                                      actionButton("C2Cneg", "Copy2ClipNegGenes"),
                                      textInput("NoOfGenes", "No. of Genes to output:", "20"),
                                      textInput("SDAVn", "SDA comp. No:"),
                                      plotOutput("SDAtsne_br2"),
                                      width=5, background = "black"
                                  ),
                                  box(title = "SDA score tabulation", status = "primary", solidHeader = TRUE,
                                      collapsible = TRUE,
                                      plotOutput("SDAtsne_br2Tab"),
                                      width=10, background = "black"
                                  )
                                ),
                                box(
                                  title = "Pos. Loadings GO", status = "primary", solidHeader = TRUE,
                                  collapsible = TRUE,
                                  plotOutput("GOpos"), #plotlyOutput
                                  width = 5, background = "black"
                                ),
                                box(
                                  title = "Neg. Loadings GO", status = "primary", solidHeader = TRUE,
                                  collapsible = TRUE,
                                  plotOutput("GOneg"),
                                  width = 5, background = "black"
                                ),
                                
                                box(title = "Pos. Top Genes", status = "info", solidHeader = TRUE, width = 4,
                                    tableOutput("packageTablePos")
                                ),
                                box(title = "Neg. Top Genes", status = "info", solidHeader = TRUE, width = 4,
                                    tableOutput("packageTableNeg")
                                )
                        ),
                        
                        # Batch-Removed Exploration--------
                        tabItem(tabName = "GeneExplorer",
                                h2("DGE_SDA Batched Removed DGE Explorer"),
                                fluidRow(
                                  box(
                                    title = "Inputs", status = "warning", solidHeader = TRUE,
                                    "Multiple formatting of gene sets accepted", 
                                    br(), "List can be seperated by comma e.g. from ", 
                                    br(), "   or spaces e.g. from Excel", 
                                    br(), "Also, single or double quotes or not",
                                    #sliderInput("ComponentN", "Slider input:", 1, 150, 1),
                                    textInput("GeneSet", "A set of genes", "'CD19', 'CD20', 'MS4A1', 'IGHM', 'IGHA2'"),
                                    width = 5
                                  ),
                                  
                                  box(
                                    title = "BatchRemoved-DGE Expr", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GeneExprSDAtSNE"),
                                    width = 5
                                  ),
                                  
                                  
                                  box(
                                    title = "Positive Loadings", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GenesEnrichSDAPos"),
                                    width = 10
                                  ),
                                  
                                  box(
                                    title = "Negative Loadings", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    plotOutput("GenesEnrichSDANeg"),
                                    width = 10
                                  )
                                  
                                )
                        ),
                        
                        # Save out----------
                        tabItem(tabName = "SaveOut",
                                h2("Save the results for downstream analysis"),
                                fluidRow(
                                  box(
                                    title = "Save as Seurat Object", status = "primary", solidHeader = TRUE,
                                    collapsible = TRUE,
                                    actionButton("SaveAsSerObj", "Save as Seurat Obj"),
                                    width = 5, background = "black"
                                  )))
                        
                      ) #end of tabItems
                ) #end of body
        ) #end UI

# Server --------

server <- function(input, output, session) {
  

  ## loading local files
  shinyFileChoose(input,'file', session=session,roots=c(wd='.'))
  
  
  ## environment defaults
  envv=reactiveValues(y=NULL)
  envv$InfoBox_sub = "Load in SDA, use either the Prime-seq or Folder input tabs"
  
  ## controls how to handle pipe depending on data
  envv$Origin = "unk"
  envv$serv = serv
  
  
  ## Feature Enrichment; used in returning enrichment profile from chi sqr analysis
  # envv$FeatEnrich$pos = list()
  # envv$FeatEnrich$neg = list()

  # Use the try() function to run the function and handle the error
  result <- try(source("apik.R",local = TRUE), silent = TRUE)
  
  # Check if an error occurred
  if (inherits(result, "try-error")) {
    # either make apik.R and put your api key in it as envv$apik = "" or type it directly
  } else {
    # The function ran successfully. Use the result here.
    source("apik.R",local = TRUE)
    updateTextInput(session, "apiKey", value = isolate(envv$apik))
  }
  
  
  ### SDA local folder
  output$select.folder <-
    renderUI(expr = selectInput(inputId = 'folder.name',
                                label = 'Folder Name',
                                choices = list.dirs(path = input$SDAroot,
                                                    full.names = FALSE,
                                                    recursive = FALSE)))
  
  
  
  source("app_Theme.R",local = TRUE)
  
  source("app_InfoBox.R",local = TRUE)
  
  # source("fxs.R",local = TRUE)
  
  
  ## ObserveEvents--------------------------------------
  source("app_OE.R",local = TRUE)
  source("app_OE_load.R",local = TRUE)
  source("app_OE_biomart.R",local = TRUE)
  source("app_OE_GO.R",local = TRUE)
  source("app_OE_tSNE.R",local = TRUE)
  
  
  ## Plots--------------------------------------
  source("app_Figs.R",local = TRUE)
  
  ## QC tabs--------------------------------------
  source("app_Figs_ChiSqr.R",local = TRUE)
  
  source("app_Figs_QC.R",local = TRUE)
  
  ## Batch removal tab--------------------------------------
  
  source("app_OE_BatchRemoval.R",local = TRUE)
  

  output$CompBatchCheckBoxSelect <- renderUI({
    
    split_text = ""
    
    MaxCompN = as.numeric(envv$SDAres$command_arguments$num_comps)
    
    choice <-  1:MaxCompN #paste0("SDA", 1:as.numeric(envv$SDAres$command_arguments$num_comps)) # envv$QC_components
    
    split_text <- unlist(strsplit(input$CompSelecTextIn, " "))
    
    if(all(split_text %in% c("pass", "fail")) & length(split_text)==MaxCompN){
      # print(split_text)
      # print(length(split_text))
      #this feels reversed but you are picking comps to remove so pass is not 2 be removed
      logical_text <- as.logical(ifelse(split_text == "pass", F, T))
      selected <- choice[logical_text]
      # envv$Remove_comps = selected
    } else {
      if(is.null(envv$Remove_comps)){
        selected <- setdiff(choice, envv$QC_components) #setdiff(choice, paste0("SDA",envv$QC_components))
        
      } else {
        selected <- envv$Remove_comps
        
      }
    }
    
    print(selected)
    
    tags$div(align = 'left',
             class = 'multicol',
             checkboxGroupInput("CompBatchCheckBoxSelect", "components",
                                choices=choice, selected = selected))
    
  })
  

  
  
  
  
  
  ## Batch Removed DGE --------------------------------------
  
  source("app_OE_BatchRemoved.R",local = TRUE)
  
  
  output$packageTablePos <- renderTable({
    
    print_gene_list(results=envv$SDAres, as.numeric(envv$QC_compIter), PosOnly = T) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes))
  }, digits = 1)
  
  output$packageTableNeg <- renderTable({
    print_gene_list(results=envv$SDAres, as.numeric(envv$QC_compIter), NegOnly = T) %>%
      as.data.frame() %>%
      head(as.numeric(input$NoOfGenes))
  }, digits = 1)
  
  
  ## Enrich N Explore-----------
  
  
  
  
  ## SAve out ------------
  
  source("app_OE_SaveOut.R",local = TRUE)
  
  
  
  
  
}

shinyApp(ui, server)

