shinyServer(function(input, output, session) {
  options(shiny.maxRequestSize = 1000*1024^2)
  
  # Generate a summary of the RNA dataset ---- 
  output$summary.rna <- renderPrint({
    df <- read.csv(input$data.rna$datapath)
    summary(df)
  })
  
  # Generate a table of the RNA dataset ----
  output$rna.contents <- DT::renderDataTable({
  
    req(input$data.rna)
    
    tryCatch(
      {
        df.rna <<- read.csv(input$data.rna$datapath,
                       header = input$header,
                       sep = input$sep,
                       quote = input$quote,row.names = 1)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    if(ncol(df.rna) >20){
      df.rna.disp = df.rna[,0:20]
    }
    else{
      df.rna.disp = df.rna
    }
    
    if(input$disp == "head") {
      return(DT::datatable(head(df.rna.disp),options = list(scrollX = TRUE)))
    }
    else {
      return(DT::datatable(df.rna.disp,options = list(scrollX = TRUE)))
    }
    
  })
  
  # Generate a Seurat object of the RNA dataset ----
    cbmc <- eventReactive(input$create,{
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...', {
                     observe({ print("create Seurat begin") })
                     df.rna <<- CollapseSpeciesExpressionMatrix(df.rna)
                     cbmcTmp <- CreateSeuratObject(counts = df.rna,
                                                   min.cells = input$min.cells,
                                                   min.features = input$min.features)
                     shiny::setProgress(value = 0.7)
                     cbmcTmp[["percent.mt"]] <- PercentageFeatureSet(cbmcTmp, pattern = "^MT-")
                   })
     
     cbmcTmp
  })
  
  observe({ print(cbmc()) })
  
  # Test Regex, return genes that match regex ----
  testreg <- eventReactive(input$testRegex,{
    observe({ print("filter begin") })
    rownames(cbmc())[grep(input$regex,rownames(cbmc()))]
  })
  observe({ print(testreg()) })
  output$genesMatch <-renderText({
    genes <- c("<h5>num of genes: ",length(testreg()),"</h5>",testreg())
    # paste("Num of genes:" , length(testreg()) ,'\n',testreg())
    HTML(genes)
  })
  
  # Add QC filters ,save filter----
  filterItem <- reactiveValues()
  observe({
    if(input$addFilter > 0){
      # print(input$addFilter)
      filterItem$items <- c(isolate(filterItem$items), isolate(input$qcLabel))
      filterItem$items <- c(isolate(filterItem$items), isolate(input$regex))
      
      nFeat <<- cbmc()@meta.data$nFeature_RNA
      # cbmc()[["percent.mt"]] <- PercentageFeatureSet(cbmc(), pattern = "^MT-")
      nMito <<- cbmc()@meta.data$percent.mt
    }
  })
  
  # Return filter name to screen ----
  output$filterExp <-renderText({
    if(input$addFilter < 1)
      return()
    filterItem$items[1]
  })
  
  # Plot Vln Plot(Filter Cells) ----
  output$nFeaturePlot <- renderPlot({
    input$plotqc
    
    nFeaDown  <<- VlnPlot(cbmc(), features = "nFeature_RNA") + geom_hline(yintercept=input$nFeaRange[1], linetype="dashed", color = "green",size=1.1) + geom_hline(yintercept=input$nFeaRange[2], linetype="dashed", color = "blue",size=1.1)
    nFeaDown
  })
  
  output$MitoPlot <- renderPlot({
    input$plotqc
    MitoDown <<- VlnPlot(cbmc(), features = filterItem$items[1]) + geom_hline(yintercept=input$mitoRange[1], linetype="dashed", color = "green",size=1.1) + geom_hline(yintercept=input$mitoRange[2], linetype="dashed", color = "blue",size=1.1)
    MitoDown
  })
  
  # Change the min and max of a slider input ----
  observe({
    if(input$plotqc > 0){
      updateSliderInput(session, "nFeaRange",min = floor(min(nFeat))-1, 
                        max = ceiling(max(nFeat))+1,value = c(min,max))
      updateSliderInput(session, "mitoRange",min = floor(min(nMito))-1, 
                        max = ceiling(max(nMito))+1,value = c(min,max))
    }
  })
  
  # Write genes numbers of a slider input ----
  output$genesbynFeaNum <- renderUI({
    if(input$plotqc < 1){
      return()
    }
    
    filter1 <- nFeat[nFeat >input$nFeaRange[1] & nFeat <input$nFeaRange[2]]
    HTML(paste0("Number of cells within gene detection thresholds ", "<b>",length(filter1), "</b>"," out of ", "<b>",length(nFeat),"</b>"))
  })
  
  output$genesbyMitoNum <- renderUI({
    if(input$plotqc < 1){
      return()
    }
    
    filter2 <- nMito[nMito >input$mitoRange[1] & nMito <input$mitoRange[2]]
    HTML(paste0("Number of cells within gene detection thresholds ", "<b>",length(filter2), "</b>"," out of ", "<b>",length(nMito),"</b>"))
  })
  
  # Download QC plot ----
  output$downloadnFea <- downloadHandler(
    filename ='nFeature_RNA.png',
    content = function(file) {
      ggsave(file, plot = nFeaDown, device = "png")
    }
  )
  output$downloadnMito <- downloadHandler(
    filename ='percent.mt.png',
    content = function(file) {
      ggsave(file, plot = MitoDown, device = "png")
    }
  )
  output$downReaPlot1 <- downloadHandler(
    filename ='Feature Scatter.png',
    content = function(file) {
      ggsave(file, plot = corPlot1, device = "png")
    }
  )
  output$downReaPlot2 <- downloadHandler(
    filename ='Feature Scatter.png',
    content = function(file) {
      ggsave(file, plot = corPlot2, device = "png")
    }
  )
  
  # Plot Feature Scatter Plot ----
  output$ReaPlot1 <- renderPlot({
    input$plotqc
    corPlot1  <<- FeatureScatter(cbmc(), feature1 = "nCount_RNA", feature2 = "percent.mt")
    corPlot1
  })
  output$ReaPlot2 <- renderPlot({
    input$plotqc
    corPlot2  <<- FeatureScatter(cbmc(), feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    corPlot2
  })
  
  # Click button, do cell Filter (within thesholds) ----
  observeEvent(input$doFilier, {
    tmp1 <<- input$nFeaRange[1]
    tmp2 <<- input$nFeaRange[2]
    tmp3 <<- input$mitoRange[1]
    tmp4 <<- input$mitoRange[2]
    cbmcQc <<- subset(cbmc(), subset = nFeature_RNA > tmp1 & nFeature_RNA < tmp2 & percent.mt > tmp3 & percent.mt < tmp4)
    # updateTabsetPanel(session, "navlistSet",selected = "Norm/Dect HVGs")
    
    print(cbmcQc)
  })
  
  # Click button, do normalization ----
  observeEvent(input$normalize,{
    cbmcNor <<- NormalizeData(cbmcQc, normalization.method = input$normMed, scale.factor = input$scaleFactor)
  })
  
  # Click button, do find HVGs ----
  # observeEvent(input$detecte,{
  #   cbmcFind <<- FindVariableFeatures(cbmcNor, selection.method = input$dectMed, nfeatures = input$featureNum)
  #   top10 <- head(VariableFeatures(cbmcFind), input$topFeature)
  #   
  #   hvgplot1 <<- VariableFeaturePlot(cbmcFind)
  #   hvgplot2 <<- LabelPoints(plot = hvgplot1, points = top10, repel = TRUE)
  # })
  
  # Plot HVGs Plot ----
  output$HvgPlot1 <- renderPlot({
    if(input$detecte < 1){
      return()
    }
    cbmcFind <<- FindVariableFeatures(cbmcNor, selection.method = input$dectMed, nfeatures = input$featureNum)
    top10 <<- head(VariableFeatures(cbmcFind), input$topFeature)
    hvgplot1 <<- VariableFeaturePlot(cbmcFind)
    hvgplot1
  })
  output$HvgPlot2 <- renderPlot({
    if(input$detecte < 1){
      return()
    }
    hvgplot2 <<- LabelPoints(plot = hvgplot1, points = top10, repel = TRUE)
    hvgplot2
  })
  
  # Download HVGs plot ----
  output$downloadHvg1 <- downloadHandler(
    filename ='HVGs.png',
    content = function(file) {
      ggsave(file, plot = hvgplot1, device = "png")
    }
  )
  output$downloadHvg2 <- downloadHandler(
    filename ='Top HVGs.png',
    content = function(file) {
      ggsave(file, plot = hvgplot2, device = "png")
    }
  )
  
  # Generate a table of the ADT dataset ----
  output$adt.contents <- DT::renderDataTable({
    
    req(input$data.adt)
    
    tryCatch(
      {
        df.adt <<- read.csv(input$data.adt$datapath,
                            header = input$adt.header,
                            sep = input$adt.sep,
                            quote = input$adt.quote,row.names = 1)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    if(ncol(df.adt) >20){
      df.adt.disp = df.adt[,0:20]
    }
    else{
      df.adt.disp = df.adt
    }
    
    if(input$adt.disp == "head") {
      return(DT::datatable(head(df.adt.disp),options = list(scrollX = TRUE)))
    }
    else {
      return(DT::datatable(df.adt.disp,options = list(scrollX = TRUE)))
    }
    
  })
  
  # Click button, do ADT normalization ----
  observeEvent(input$adt.normalize,{
    curcellnames <- colnames(as.matrix(cbmcFind[["RNA"]]@data))
    print(length(curcellnames))
    df.adt <- subset(as.matrix(df.adt), select=c(curcellnames))
    
    cbmcFind[["ADT"]] <- CreateAssayObject(counts = df.adt)
    cbmcFind <<- NormalizeData(cbmcFind, assay = "ADT",normalization.method = input$adt.normMed, scale.factor = input$adt.scaleFactor)
    cbmcFind <<- ScaleData(cbmcFind, assay = "ADT")
  })
  
  # # Click button, compute PCA ----
  # observeEvent(input$computePCA,{
  #   # standard scaling (no regression)
  #   cbmcFind <<- ScaleData(cbmcFind)
  #   
  #   # Run PCA
  #   cbmcFind <<- RunPCA(cbmcFind, verbose = FALSE)
  # })
  
  # PrintPCA ----
  output$pcaPrint <- renderText({
    if(input$computePCA < 1)
      return()
    
    # standard scaling (no regression)
    cbmcFind <<- ScaleData(cbmcFind)
    # Run PCA
    cbmcFind <<- RunPCA(cbmcFind, verbose = FALSE)
    
    # print(cbmcFind[["pca"]], dims = 1:input$pcsNum, nfeatures = input$pcsGenesNum)
    printStr <- capture.output(print(cbmcFind[["pca"]],dims = 1:input$pcsNum, nfeatures = input$pcsGenesNum))
    printStr <- paste(printStr, collapse = "<br>")
    HTML(printStr)
  })
  output$pcaPrintTitle <- renderText({
    paste("<b>PCA Print Output:</b>")
  })
  
  
  # PC Elbow Plot ----
  output$pcElbowPlot <- renderPlot({
    if(input$computePCA < 1)
      return()
    
    ElbowPlot(cbmcFind, ndims = 50)
  })
  
  # Protein Prediction ----
  # ()
  
  # Click button, predict ADT by Python file
  predictADT <- eventReactive(input$predictADT,{
    source_python('E://Python Virtualenv/CITE-analysis/adt_prediction.py')
    predict.ADT <- rf_multi()
  })
  
  output$adt.predict <- DT::renderDataTable({
    if(ncol(predictADT()) >20) {
      return(DT::datatable(predictADT()[,0:20],options = list(scrollX = TRUE)))
    }
    else {
      return(DT::datatable(predictADT(),options = list(scrollX = TRUE)))
    }
  })
  
  # Downloadable csv of predict ADT ----
  output$downloadPred <- downloadHandler(
    filename = function() {
      paste("ADT Prediction", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(predictADT(), file, row.names = FALSE)
    }
  )
  
  # Cluster cells ----
  observeEvent(input$doCluster,{
    req(cbmcFind)
    cbmcFind <<- FindNeighbors(cbmcFind, dims = 1:input$pcsToUseNum)
    cbmcFind <<- FindClusters(cbmcFind, resolution = input$resolution)
    cbmcFind <<- RunTSNE(cbmcFind, dims = 1:input$pcsToUseNum, method = "FIt-SNE")
    cbmcFind <<- RunUMAP(cbmcFind, dims = 1:input$pcsToUseNum)
  })

  # Plot UMAP ----
  output$cluterPlot <- renderPlot({
    if(input$doCluster < 1)
      return()
    
    clusterplot <<- DimPlot(cbmcFind, reduction = input$cluterPlotChoice)
    clusterplot
  })
  
  # Download Cluster plot ----
  output$downloadCluster <- downloadHandler(
    filename ='cluster.png',
    content = function(file) {
      ggsave(file, plot = clusterplot, device = "png")
    }
  )
  
  # Click button, find all markers
  allmarkers <- eventReactive(input$findMarkers,{
    
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', {
                   cbmc.rna.markers <- FindAllMarkers(cbmcFind, max.cells.per.ident = 100, min.diff.pct = input$minDiffPct, only.pos = TRUE)
                   shiny::setProgress(value = 0.7)
                   cbmc.rna.markers %>% group_by(cluster) %>% top_n(n = input$TopGenestoShow, wt = avg_logFC)
                 })
  })
  
  output$rnaMarkers.contents <- DT::renderDataTable({
    if(ncol(allmarkers()) >20) {
      return(DT::datatable(allmarkers()[,0:20],options = list(scrollX = TRUE)))
    }
    else {
      return(DT::datatable(allmarkers(),options = list(scrollX = TRUE)))
    }
  })
  
  # Download RNA all markers CSV file ----
  output$downloadallMarkers <- downloadHandler(
    filename ='rna.allmarkers.csv',
    content = function(file) {
      write.csv(allmarkers(), file)
    }
  )
  
  # unpdate select choice by found marker genes names ---- 
  observe({
    updateSelectizeInput(session,'markersToPlot',
                         choices = allmarkers()['gene'] , selected=NULL)
  })
  
  # Click button, plot markers Vln Plot ---
  output$markerVlnPlot <- renderPlot({
    if(input$doMarkerVln < 1)
      return()
    
    VlnPlot(cbmcFind, features = input$markersToPlot)
  })
  
  # Click button, assign cell types to cluster ----
  output$cellTypeIdent <- renderPlot({
    
  })
  
  # click button, plot protein levels on ADT clusters ----
  output$adtLevelFeaPlot <- renderPlot({
    if(input$doADTLev < 1)
      return()
    
    featureStr <- paste(input$`adt&rnaToShow1`,input$`adt&rnaToShow2`, sep = ',')
    featureVec <- unlist(strsplit(featureStr,split=','))
    FeaturePlot(cbmcFind, features = featureVec, min.cutoff = "q05", max.cutoff = "q95", ncol = 4)
  })
  
  output$adtLevelRidgePlot <- renderPlot({
    if(input$doADTLev < 1)
      return()
    
    featureVec <- unlist(strsplit(input$'adt&rnaToShow1',split=','))
    RidgePlot(cbmcFind, features = featureVec, ncol = 2)
  })

})
