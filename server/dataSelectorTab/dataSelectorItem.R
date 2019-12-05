dataSelectorTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    dataSelectorMessages <- allReactiveMsgs$dataSelectorMessages
    
    updateCurrentSource <- eventReactive(input$dataSource,{
        currentMetadata$source <- input$dataSource
        allReactiveVars <- clearReactiveVars(allReactiveVars)
    })
    
    updateCurrentDataset <- eventReactive(input$dataSource,{
        currentMetadata$dataset <- input$dataDataset
        allReactiveVars <- clearReactiveVars(allReactiveVars)
    })

    updateCurrentMetadata <- eventReactive(input$dataDataset,{
        s <- currentMetadata$source
        d <- input$dataDataset
        currentMetadata$dataset <- d
        currInd <- paste0("FROM metadata WHERE source='",s,
            "' AND dataset='",d,"'")
        summInd <- paste0("FROM datasets WHERE dataset='",d,"'")
        
        currentMetadata$class <- as.character(dbGetQuery(metadata,
            paste0("SELECT DISTINCT class ",currInd))$class)
        currentMetadata$metadata <- 
            dbGetQuery(metadata,paste0("SELECT * ",currInd))
        currentMetadata$genome <- as.character(dbGetQuery(metadata,
            paste0("SELECT DISTINCT genome ",currInd))$genome)[1]
        currentMetadata$short_summary <- as.character(dbGetQuery(metadata,
            paste0("SELECT DISTINCT short_summary ",summInd))$short_summary)[1]
        currentMetadata$title <- as.character(dbGetQuery(metadata,
            paste0("SELECT DISTINCT title ",summInd))$title)[1]
        currentMetadata$link <- as.character(dbGetQuery(metadata,
            paste0("SELECT DISTINCT link ",summInd))$link)[1]
        
        if (!is.na(currentMetadata$genome)
            && is.null(loadedGenomes[[currentMetadata$genome]]$dbGene)) {
            load(file.path("genome",currentMetadata$genome,"gene.rda"))
            loadedGenomes[[currentMetadata$genome]]$dbGene <<- gene
            gn <- names(gene)
            names(gn) <- as.character(elementMetadata(gene)$gene_name)
            loadedGenomes[[currentMetadata$genome]]$geneNames <<- gn
            dataSelectorMessages <- updateMessages(
                dataSelectorMessages,
                type="SUCCESS",
                msg=paste(getTime("SUCCESS"),"Genome ",
                    currentMetadata$genome," genes loaded!", sep="")
            )
        }
        if (!isEmpty(s) && !isEmpty(d))
            dataSelectorMessages <- updateMessages(
                dataSelectorMessages,
                type="INFO",
                msg=paste(getTime("INFO"),"Current source: ",s,
                    ", Current dataset: ",d,", Classes: ", 
                    paste(currentMetadata$class,collapse=", "),
                    ", Total samples: ",nrow(currentMetadata$metadata),
                    sep="")
            )
    })
    
    createDataset <- eventReactive(input$createDataset,{
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        c <- currentMetadata$class
        m <- currentMetadata$metadata
        currentMetadata$final <- do.call("rbind",lapply(c,function(x) {
            if (input[[paste("sampleSelectType_",x,sep="")]]=="all")
                return(m[m$class==x,])
            else if (input[[paste("sampleSelectType_",x,sep="")]]=="none") {
                m <- match(x,currentMetadata$class)
                currentMetadata$class <- currentMetadata$class[-m]
                return(NULL)
            }
            else {
                sel <- 
                    input[[paste("classTable_",x,"_rows_selected",sep="")]]
                if (length(sel)==0)
                    return(m[m$class==x,])
                else {
                    mm <- m[m$class==x,]
                    return(mm[sel,])
                }
            }
        }))
        
        msgpart <- paste(lapply(c,function(x,m) {
            paste("Number of samples in ",x,": ",nrow(m[m$class==x,]))
        },currentMetadata$final),collapse=", ")
        dataSelectorMessages <- updateMessages(
            dataSelectorMessages,
            type="SUCCESS",
            msg=paste(getTime("SUCCESS"),"Dataset created! Source: ",
                currentMetadata$source,", Dataset: ",
                currentMetadata$dataset,", Classes: ",
                paste(currentMetadata$class,collapse=", "),", ",msgpart,
                ", Total samples: ",nrow(currentMetadata$final),sep="")
        )
    })
    
    clearDataset <- eventReactive(input$clearDataset,{
        currentMetadata$final <- NULL
        allReactiveVars <- clearReactiveVars(allReactiveVars)
        dataSelectorMessages <- updateMessages(
            dataSelectorMessages,
            type="WARNING",
            msg=paste(getTime("WARNING"),"Dataset cleared!")
        )
    })
    
    loadDataset <- eventReactive(input$loadDataset,{
        s <- input$dataSource
        d <- input$dataDataset
        c <- currentMetadata$class
        if (!isEmpty(s) && !isEmpty(d)) {
            if (is.null(dataFiles[[s]][[d]]))
                dataSelectorMessages <- updateMessages(
                    dataSelectorMessages,
                    type="WARNING",
                    msg=paste(getTime("WARNING"),
                        "Dataset ",d," from ",s," does not exist! Please wait ",
                        "until the loading of existing dataset names from ",s,
                        " is complete. If the problem persists, contact the ",
                        "administrator.",sep="")
                )
            else {
                if (is.null(loadedData[[s]][[d]])) {
                    dataSelectorMessages <- updateMessages(
                        dataSelectorMessages,
                        type="INFO",
                        msg=paste(getTime("INFO"),
                            "Loading dataset ",d," from ",s,sep="")
                    )                    
                    load(dataFiles[[s]][[d]])
                    # ^$#$#@@#$^%$$!!!!
                    loadedData[[s]][[d]] <<- b2c.out
                    dataSelectorMessages <- updateMessages(
                        dataSelectorMessages,
                        type="SUCCESS",
                        msg=paste(getTime("SUCCESS"),
                            "Loaded ",d," from ",s,"!",sep="")
                    )
                }
                else
                    dataSelectorMessages <- updateMessages(
                        dataSelectorMessages,
                        type="INFO",
                        msg=paste(getTime("INFO"),
                            "Dataset ",d," from ",s," already loaded.",sep="")
                    )
            }
        }
    })

    return(list(
        updateCurrentSource=updateCurrentSource,
        updateCurrentDataset=updateCurrentDataset,
        updateCurrentMetadata=updateCurrentMetadata,
        createDataset=createDataset,
        clearDataset=clearDataset,
        loadDataset=loadDataset
    ))
}

dataSelectorTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    dataSelectorMessages <- allReactiveMsgs$dataSelectorMessages
    
    #classDataTables <- reactive({
    #   s <- input$dataSource
    #   d <- input$dataDataset
    #   c <- currentMetadata$class
    #   lapply(c,function(x,s,d) {
    #       tab <- metadata[which(as.character(metadata$source)==s 
    #           & as.character(metadata$dataset)==d
    #           & as.character(metadata$class)==x),
    #           c("sample_id","alt_id","norm_factor","library_strategy")]
    #       output[[paste("classTable",x,sep="_")]] <- 
    #           DT::renderDataTable(
    #               tab,
    #               class="display compact",
    #               rownames=FALSE,
    #               options=list(
    #                   searchHighlight=TRUE,
    #                   pageLength=10,
    #                   lengthMenu=c(10,20,50,100)
    #               )
    #           )
    #   },s,d)
    #})
    
    handleSampleSelection <- reactive({
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        c <- currentMetadata$class

        lapply(c,function(x) {
            observeEvent(input[[paste("clearSelection_",x,sep="")]],{
                proxy <- dataTableProxy(paste("classTable_",x,sep=""))
                selectRows(proxy,NULL)
            })
        })
        lapply(c,function(x) {
            relevant_rows <- as.numeric(dbGetQuery(metadata,
                paste0("SELECT COUNT(*) FROM metadata WHERE source='",s,
                    "' AND dataset='",d,"' AND class='",x,"'")))
            N <- 1:relevant_rows
            observeEvent(input[[paste("invertSelection_",x,sep="")]],{
                sel <- input[[paste("classTable_",x,
                    "_rows_selected",sep="")]]
                if (length(sel)>0) {
                    N <- N[-sel]
                    proxy <- dataTableProxy(paste("classTable_",x,sep=""))
                    selectRows(proxy,N)
                }
            })
        })
    })
    
    #updateDataSelectorMessages <- reactive({
    #   output$dataSelectorMessages <- renderUI({
    #       lapply(dataSelectorMessages$messages,function(x) {
    #           switch(x$type,
    #               INFO = {
    #                   div(class="info-box",x$msg)
    #               },
    #               SUCCESS = {
    #                   div(class="success-box",x$msg)
    #               },
    #               WARNING = {
    #                   div(class="warn-box",x$msg)
    #               },
    #               ERROR = {
    #                   div(class="error-box",x$msg)
    #               }
    #           )
    #       })
    #   }) 
    #})
    
    return(list(
        #classDataTables=classDataTables,
        handleSampleSelection=handleSampleSelection
    ))
}

dataSelectorTabPanelRenderUI <- function(output,session,allReactiveVars,
    allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    dataSelectorMessages <- allReactiveMsgs$dataSelectorMessages
    
    output$currentDatasetTable <- DT::renderDataTable(
        if (is.null(currentMetadata$final))
            data.frame(
                sample_id=character(0),
                alt_id=character(0),
                class=character(0),
                norm_factor=numeric(0),
                quality=integer(0)
            )
        else    
            as.data.frame(currentMetadata$final[,
                c("sample_id","alt_id","class","norm_factor")]),
        class="display compact",
        rownames=FALSE,
        options=list(
            searchHighlight=TRUE,
            pageLength=10,
            lengthMenu=c(10,20,50,100)
        )
    )
        
    output$dataSource <- renderUI({
        selectInput(
            inputId="dataSource",
            label="Select data source",
            choices=sources
        )
    })
    
    # Init the dataset
    output$dataDataset <- renderUI({
        selectInput(
            inputId="dataDataset",
            label="Select dataset",
            choices=as.character(dbGetQuery(metadata,
                paste0("SELECT DISTINCT(dataset) FROM metadata WHERE source='",
                    input$dataSource,"'"))$dataset)
        )
    })
    
    # Fill genome
    output$dataGenome <- renderUI({
        disabled(textInput(
            inputId="dataGenome",
            label="Genome",
            value=currentMetadata$genome
        ))
    })
    
    # Fill short summary
    output$dataShortSummary <- renderUI({
        s <- input$dataSource
        d <- input$dataDataset
        url <- a(currentMetadata$title ,href=currentMetadata$link)
        if (!isEmpty(s) && !isEmpty(d)) {
            list(
                h4(tagList("Title: ", url)),
                helpText(currentMetadata$short_summary)
                )
        }
    })

    output$dataSelectHint <- renderUI({
        s <- input$dataSource
        d <- input$dataDataset
        if (!isEmpty(s) && !isEmpty(d)) {
            if (!is.null(loadedData[[s]][[d]]))
                list(
                    h4("Select samples for selected dataset"),
                    helpText(paste("Select 'Custom samples' for ",
                        "customized sample selections or 'All samples' to ",
                        "select all samples from each class."))
                )
        }
    })
    
    output$dataCustomSamples <- renderUI({
        s <- input$dataSource
        d <- input$dataDataset
        c <- as.character(dbGetQuery(metadata,
            paste0("SELECT DISTINCT(class) FROM metadata WHERE source='",s,
                "' AND dataset='",d,"'"))$class)

        if (!isEmpty(s) && !isEmpty(d)) {
            if (!is.null(loadedData[[s]][[d]])) {
                #classDataTables()
                s <- input$dataSource
                d <- input$dataDataset
                c <- currentMetadata$class
                lapply(c,function(x,s,d) {
                    tab <- dbGetQuery(metadata,
                        paste0("SELECT sample_id,alt_id,norm_factor,",
                            "library_strategy FROM metadata WHERE dataset='",d,
                            "' AND class='",x,"' AND source='",s,"'"))
                    
                    output[[paste("classTable",x,sep="_")]] <- 
                        DT::renderDataTable(
                            tab,
                            class="display compact",
                            rownames=FALSE,
                            options=list(
                                searchHighlight=TRUE,
                                pageLength=10,
                                lengthMenu=c(10,20,50,100)
                            )
                        )
                },s,d)
                lapply(c,function(x,s,d) {
                    list(
                        radioButtons(
                            inputId=paste("sampleSelectType",x,sep="_"),
                            label=paste("Select",x,"samples"),
                            inline=TRUE,
                            choices=list(
                                "All samples"="all",
                                "Custom samples"="custom",
                                "No samples"="none"
                            )
                        ),
                        conditionalPanel(
                            condition=paste("input['sampleSelectType_",x,
                                "']=='custom'",sep=""),
                            div(
                                class="small table-container",
                                DT::dataTableOutput(paste("classTable",x,
                                    sep="_")),
                                div(
                                    class="pull-left",
                                    style="display:block; margin:5px;",
                                    actionButton(
                                        paste("clearSelection_",x,sep=""),
                                        paste("Clear selection"),
                                        class="btn-xs"
                                    )
                                ),
                                div(
                                    class="pull-left",
                                    style="display:block; margin:5px;",
                                    actionButton(
                                        paste("invertSelection_",x,sep=""),
                                        paste("Invert selection"),
                                        class="btn-xs"
                                    )
                                )
                            )
                        ),
                        hr()
                    )
                },s,d)
            } 
        }
    })
    
    output$dataSelectorMessages <- renderUI({
        lapply(dataSelectorMessages$messages,function(x) {
            switch(x$type,
                INFO = {
                    div(class="info-box",x$msg)
                },
                SUCCESS = {
                    div(class="success-box",x$msg)
                },
                WARNING = {
                    div(class="warn-box",x$msg)
                },
                ERROR = {
                    div(class="error-box",x$msg)
                }
            )
        })
    })
}

dataSelectorTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    dataSelectorMessages <- allReactiveMsgs$dataSelectorMessages
    
    dataSelectorTabPanelReactiveEvents <- 
        dataSelectorTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
       
    updateCurrentSource <- 
        dataSelectorTabPanelReactiveEvents$updateCurrentSource
    updateCurrentDataset <- 
        dataSelectorTabPanelReactiveEvents$updateCurrentDataset
    updateCurrentMetadata <- 
        dataSelectorTabPanelReactiveEvents$updateCurrentMetadata
    createDataset <- dataSelectorTabPanelReactiveEvents$createDataset
    clearDataset <- dataSelectorTabPanelReactiveEvents$clearDataset
    loadDataset <- dataSelectorTabPanelReactiveEvents$loadDataset
    
    dataSelectorTabPanelReactiveExprs <- 
        dataSelectorTabPanelReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
    
    classDataTables <- dataSelectorTabPanelReactiveExprs$classDataTables
    handleSampleSelection <- 
        dataSelectorTabPanelReactiveExprs$handleSampleSelection
    
    dataSelectorTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
    
    observe({
        #tryCatch({
        #   shinyjs::disable("loadDataset")
            loadDataset()
        #},error=function(e) {
        #   s <- currentMetadata$source
        #   d <- currentMetadata$dataset
        #   dataSelectorMessages <- updateMessages(
        #       dataSelectorMessages,
        #       type="ERROR",
        #       msg=paste(getTime("ERROR"),
        #           "Error while loading dataset ",d," from ",s,". Please ",
        #           "try loading again and if it fails contact the ",
        #           "administrator with the following error: ",e,sep="")
        #   )
        #},
        #finally={
        #   shinyjs::enable("loadDataset")
        #})
    })
    
    observe({
        updateCurrentSource()
        updateCurrentMetadata()
        handleSampleSelection()
        createDataset()
        clearDataset()
    })
    
    observe({
        s <- currentMetadata$source
        d <- currentMetadata$dataset
        if (is.null(loadedData[[s]][[d]])) {
            shinyjs::disable("createDataset")
            shinyjs::disable("showFastqc")

        }else {
            shinyjs::enable("createDataset")
            shinyjs::enable("showFastqc")   
        }         
    })

    observeEvent(input$showFastqc, {
        d <- currentMetadata$dataset
        #addResourcePath(d,paste0(getwd(),"/data/",d))
        reportPath <- paste0(appConfig$urls$base,"/data/",d,"/",
            paste0(d,"_multiqc_report.html"))
        showModal(modalDialog(
            title=paste("Fastq sample quality control for dataset",d),
            #includeHTML("www/GSE59017_multiqc_report.html")
            tags$iframe(
                #src=paste0(d,"/",d,"_multiqc_report.html"),
                src=reportPath,
                width="100%",height=640,
                #sandbox=paste("allow-forms allow-popups allow-pointer-lock",
                #   "allow-same-origin allow-scripts allow-top-navigation"),
                frameborder=0
            ),
            size="l"
        ))
    })
}
