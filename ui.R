# Load packages ----
library(dplyr)
library(DT)
library(ggplot2)
library(Seurat)
library(shiny)
library(shinythemes)
library(reticulate)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    # shiny theme
    theme = shinythemes::shinytheme("simplex"),
    
    navbarPage("CITE Analysis",
               # Page1: Data Preprocessing 
               tabPanel("Preprocessing",
                    navlistPanel(
                      id="navlistSet",
                        
                        # Step1: Upload RNA Data ----
                        tabPanel("Upload RNA Data",
                                 fileInput("data.rna", h4("Choose CSV File"),
                                           multiple = FALSE,
                                           accept = c("text/csv",
                                                      "text/comma-separated-values,text/plain",
                                                      ".csv")),
                                 fluidRow(
                                     column(3,checkboxInput("header", "Header", TRUE)),
                                     column(3,radioButtons("sep", "Separator",
                                                           choices = c(Comma = ",",
                                                                       Semicolon = ";",
                                                                       Tab = "\t"),
                                                           selected = ",")),
                                     column(3,radioButtons("quote", "Quote",
                                                           choices = c(None = "",
                                                                       "Double Quote" = '"',
                                                                       "Single Quote" = "'"),
                                                           selected = '"')),
                                     column(3, radioButtons("disp", "Display",
                                                            choices = c(Head = "head",
                                                                        All = "all"),
                                                            selected = "head"))
                                 ),
                                 

                                 # Horizontal line ----
                                 tags$hr(),
                                 
                                 h4("Observations"),
                                 fluidRow(
                                   column(12,DT::dataTableOutput("rna.contents"))
                                 ),
                        ),
                        
                        # Step2: QC Filter(genes) ----
                        tabPanel("QC Filter (genes)",
                                 h4("Initialize the RNA Seurat object "),
                                 fluidRow(
                                   column(4,numericInput("min.cells", label = h5("Minimum number of cells per gene"), value = 3)),
                                   column(4,numericInput("min.features", label = h5("Minimum number of genes per cell"), value = 200)),
                                 ),
                                 actionButton("create", label = "Create"),
                                 
                                 # Horizontal line ----
                                 tags$hr(),
                                 
                                 h4("Filter Options"),
                                 fluidRow(
                                   column(8,wellPanel(
                                     fluidRow(
                                       column(6,textInput("regex", label = h5("Regex expression input"), placeholder = "Eg. ^MT-")),
                                       column(6,textInput("qcLabel",label = h5("QC Label"),placeholder = "Eg. percent.mt")),
                                     ),
                                     fluidRow(
                                       column(6,actionButton("testRegex", label = "Test Regex",class ="btn-info")),
                                       column(6,actionButton("addFilter", label = "Add Filter",class ="btn-success"))
                                     ),
                                     br(),
                                     helpText('Genes that match regex:'),
                                     helpText('-------------------------------------------'),
                                     htmlOutput("genesMatch")
                                   ),),
                                   column(4,verbatimTextOutput("filterExp",placeholder = FALSE))
                                 ),
                                 actionButton("plotqc",label = "Plot"),
                        ),
                        
                        # Step3: VlnPlot(Filter Cells) ----
                        tabPanel("VlnPlot (Filter Cells)",
                                 h4("Vln Plot (Filter Cells)"),
                                 fluidRow(
                                   column(6,plotOutput("nFeaturePlot")),
                                   column(6,plotOutput("MitoPlot")),
                                 ),
                                 fluidRow(
                                   column(6,wellPanel(
                                     sliderInput("nFeaRange", "Threshold:",
                                                 min = 1, max = 3000,
                                                 value = c(min,max)),
                                     htmlOutput("genesbynFeaNum"),
                                     br(),
                                     downloadButton('downloadnFea', label = 'Download Plot')
                                   )),
                                   column(6,wellPanel(
                                     sliderInput("mitoRange", "Threshold:",
                                                 min = 1, max = 3000,
                                                 value = c(min,max)),
                                     htmlOutput("genesbyMitoNum"),
                                     br(),
                                     downloadButton('downloadnMito', label = 'Download Plot')
                                   ))
                                 ),
                                 div(style="text-align:center",actionButton("doFilier",label = "Filter Cells(within thesholds)",class ="btn-success")),
                                 
                                 # Horizontal line ----
                                 tags$hr(),
                                 
                                 h4("Feature Scatter Plot"),
                                 fluidRow(
                                   column(6,plotOutput("ReaPlot1")),
                                   column(6,plotOutput("ReaPlot2")),
                                 ),
                                 fluidRow(
                                   column(6,downloadButton('downReaPlot1', label = 'Download Plot')),
                                   column(6,downloadButton('downReaPlot2', label = 'Download Plot')),
                                 ),
                                 br(),br()
                        ),
                        
                      # Step4: Norm/Dect HVGs ----
                        tabPanel("Norm/Dect HVGs",
                                 h4("Normalizing the data"),
                                 fluidRow(
                                   column(4,selectInput("normMed", label = h5("Normalization method"), 
                                                        choices = list("LogNormalize", "CLR", "RC"), 
                                                        selected = ),),
                                   column(4,numericInput("scaleFactor", label = h5("Scale factor"), value = 10000))
                                 ),
                                 actionButton("normalize",label = "Normalize"),
                                 
                                 # Horizontal line ----
                                 tags$hr(),
                                 
                                 h4("Identification of highly variable features (feature selection)"),
                                 wellPanel(
                                   fluidRow(
                                     column(4,selectInput("dectMed",label = h5("Selection method"),
                                                          choices = list("vst" = "vst","mean.var.plot (mvp)" = "mvp","dispersion (disp)" = "disp"),
                                                          selected = )),
                                     column(3,numericInput("featureNum", label = h5("Features number"), value = 1000)),
                                     column(3,numericInput("topFeature", label = h5("Top features number"), value = 10)),
                                     column(2,actionButton("detecte",label = "Detecte",class = "btn-success"))
                                   ),
                                   tags$style(type='text/css', "#detecte { width:100%; margin-top: 35px;}"),
                                 ),
                                 br(),
                                 fluidRow(
                                   column(6,plotOutput("HvgPlot1")),
                                   column(6,plotOutput("HvgPlot2")),
                                 ),
                                 fluidRow(
                                   column(6,downloadButton('downloadHvg1', label = 'Download Plot')),
                                   column(6,downloadButton('downloadHvg2', label = 'Download Plot'))
                                 ),
                                 br(),br()
                        ),
                        
                        # Protein Prediction(optional) ---
                        tabPanel('Protein Prediction (Optional)'),
                        
                        #step5: PCA Reduction ----
                        tabPanel('PCA Reduction',
                                 h4("Perform linear dimensional reduction"),
                                 fluidRow(
                                   column(4,numericInput("pcsNum", label = h5("Number of PCs to show"), value = 4),
                                          numericInput("pcsGenesNum", label = h5("Number of Genes to show"), value = 5),
                                          actionButton("computePCA",label = "Compute")
                                          ),
                                   column(8,br(),br(),wellPanel(
                                     htmlOutput("pcaPrintTitle"),
                                     "-------------------------------------",
                                     htmlOutput("pcaPrint")
                                   ))
                                 ),
                                 
                                 # Horizontal line ----
                                 tags$hr(),
                                 
                                 h4("Determine the dimensionality of the dataset"),
                                 h5("1) PC Elbow Plot"),
                                 plotOutput("pcElbowPlot"),
                                 h5("2) PC JackStraw Plot (Optional)"),
                                 br(),br()
                                 ),
                        
                        # Step6: Upload ADT data ----
                        tabPanel("Upload Protein Data",
                                 fileInput("data.adt", h4("Choose CSV File"),
                                           multiple = FALSE,
                                           accept = c("text/csv",
                                                      "text/comma-separated-values,text/plain",
                                                      ".csv")),
                                 fluidRow(
                                   column(3,checkboxInput("adt.header", "Header", TRUE)),
                                   column(3,radioButtons("adt.sep", "Separator",
                                                         choices = c(Comma = ",",
                                                                     Semicolon = ";",
                                                                     Tab = "\t"),
                                                         selected = ",")),
                                   column(3,radioButtons("adt.quote", "Quote",
                                                         choices = c(None = "",
                                                                     "Double Quote" = '"',
                                                                     "Single Quote" = "'"),
                                                         selected = '"')),
                                   column(3, radioButtons("adt.disp", "Display",
                                                          choices = c(Head = "head",
                                                                      All = "all"),
                                                          selected = "head"))
                                 ),
                                 
                                 # Horizontal line ----
                                 tags$hr(),
                                 
                                 h4("Observations"),
                                 fluidRow(
                                   column(12,DT::dataTableOutput("adt.contents"))
                                 ),
                                 
                                 # Horizontal line ----
                                 tags$hr(),
                                 
                                 h4("Normalizing the data"),
                                 fluidRow(
                                   column(4,selectInput("adt.normMed", label = h5("Normalization method"), 
                                                        choices = list("LogNormalize", "CLR", "RC"), 
                                                        selected = "CLR"),),
                                   column(4,numericInput("adt.scaleFactor", label = h5("Scale factor"), value = 10000))
                                 ),
                                 br(),
                                 actionButton("adt.normalize",label = "Normalize"),
                                 br(),br()
                        )
                    ),
               ),
               
               # Page 2: Protein Prediction
               tabPanel("Protein Prediction",
                        sidebarLayout(
                          sidebarPanel(
                            fileInput("rna.predict", h4("Choose RNA File (Optional)"),
                                      multiple = FALSE,
                                      accept = c("text/csv",
                                                 "text/comma-separated-values,text/plain",
                                                 ".csv")),
                            fluidRow(
                              column(6,checkboxInput("pred.header", "Header", TRUE)),
                              column(6,radioButtons("pred.sep", "Separator",
                                                    choices = c(Comma = ",",
                                                                Semicolon = ";",
                                                                Tab = "\t"),
                                                    selected = ",")),
                            ),
                            fluidRow(
                              column(6,radioButtons("pred.quote", "Quote",
                                                    choices = c(None = "",
                                                                "Double Quote" = '"',
                                                                "Single Quote" = "'"),
                                                    selected = '"')),
                              column(6, radioButtons("pred.disp", "Display",
                                                     choices = c(Head = "head",
                                                                 All = "all"),
                                                     selected = "head"))
                            ),
                            
                            # Horizontal line ----
                            tags$hr(),
                            
                            actionButton("predictADT",label = "Compute",class ="btn-success")
                          ),
                          mainPanel(
                            h4("Prediction Output"),
                            fluidRow(
                              column(12,DT::dataTableOutput("adt.predict"))
                            ),
                            downloadButton("downloadPred", "Download Data")
                          )
                        )),
               
               # Page3: Downstream Analysis
               tabPanel("Downstream Analysis",
                    sidebarLayout(
                        sidebarPanel(
                            h4("FindClusters paramters"),
                            fluidRow(
                              column(6,numericInput("pcsToUseNum", label = h5("Choose PCs number to use"), value = 10)),
                              column(6,numericInput("resolution", label = h5("Scale factor"), value = 0.8))
                            ),
                            actionButton("doCluster",label = "Cluster Cells"),
                            h4("Find biomarkers"),
                            fluidRow(
                              column(6,numericInput("minDiffPct", label = h5("Minimum percentage difference"), value = 0.3)),
                              column(6,numericInput("TopGenestoShow", label = h5("Top genes to show per cluster"), value = 5))
                            ),
                            actionButton("findMarkers",label = "Find AllMarkers"),
                            h4("Vlz protein levels"),
                            fluidRow(
                              column(6,textInput("adt&rnaToShow1", label = h5("ADT Input"), value = "adt_CD3,adt_CD11c,adt_CD8,adt_CD16")),
                              column(6,textInput("adt&rnaToShow2", label = h5("Corresponding RNA Input"), value = "CD3E,ITGAX,CD8A,FCGR3A")),
                            ),
                            actionButton("doADTLev",label = "Plot ADT levels")
                        ),
                        mainPanel(
                            tabsetPanel(
                              # Analysis1: Clustering ----
                              tabPanel("Clustering",
                                       radioButtons("cluterPlotChoice", label = h4("Choose reduction method to proceed with:"),
                                                    choices = list("tsne" = "tsne", "umap" = "umap"), 
                                                    selected = "umap"),
                                       plotOutput("cluterPlot"),
                                       br(),
                                       downloadButton('downloadCluster', label = 'Download Plot')),
                              
                              # Analysis2: Finding markers ----
                              tabPanel("Finding markers",
                                       br(),
                                       h4("Finding cluster biomarkers"),
                                       DT::dataTableOutput("rnaMarkers.contents"),
                                       downloadButton('downloadallMarkers', label = 'Download File'),
                                       
                                       # Horizontal line ----
                                       tags$hr(),
                                       
                                       h4("Visualizing Marker Expression:"),
                                       fluidRow(
                                         column(6,
                                                wellPanel(
                                                  fluidRow(
                                                    column(8,selectInput('markersToPlot', label = h5('Genes to plot:'), choices = NULL, multiple=TRUE, selectize=TRUE)),
                                                    column(4,actionButton("doMarkerVln",label = "Plot",class = "btn-success"))
                                                  ),
                                                  tags$style(type='text/css', "#doMarkerVln { width:100%; margin-top: 35px;}"),
                                                ),
                                                )
                                       ),
                                       plotOutput("markerVlnPlot"),
                                       
                                       # Horizontal line ----
                                       tags$hr(),
                                       
                                       h4("Assigning cell type identity to clusters: (Optional)"),
                                       fluidRow(
                                         column(6,
                                                wellPanel(
                                                  fluidRow(
                                                    column(8,textInput("cellTypeInput", label = h5("Cell type input"), placeholder = "Eg. 'Naive CD4 T', 'NK'")),
                                                    column(4,actionButton("cellTypeIdent",label = "Plot"))
                                                  ),
                                                  tags$style(type='text/css', "#cellTypeIdent { width:100%; margin-top: 35px;}"),
                                                ))
                                       ),
                                       plotOutput("typeIdentPlot"),
                                       br()
                                       ),
                              
                              # Analysis3: Visualize protein levels on RNA clusters ----
                              tabPanel("Vlz protein levels",
                                       h4("Visualize protein levels on RNA clusters"),
                                       h5("1) Feature Plot"),
                                       plotOutput("adtLevelFeaPlot"),
                                       h5("2) Ridge Plot"),
                                       plotOutput("adtLevelRidgePlot"),
                                       br(),br())
                            )
                        ),      
                    ),
                )
               )
))
