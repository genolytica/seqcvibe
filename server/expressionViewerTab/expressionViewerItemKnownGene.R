expressionExplorerTabPanelEventReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    return(NULL)
}

expressionExplorerTabPanelReactive <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentTables <- allReactiveVars$currentTables
    
    createCurrentTables <- reactive({
        for (c in currentMetadata$class)
            currentTables[[c]] <- NULL
    })
    
    handleExpressionSampleSelection <- reactive({
        c <- currentMetadata$class
        lapply(c,function(x) {
            observeEvent(input[[paste("clearCountSelection_",x,sep="")]],{
                proxy <- dataTableProxy(paste("countTable_",x,sep=""))
                selectRows(proxy,NULL)
            })
        })
        lapply(c,function(x) {
            observeEvent(input[[paste("invertCountSelection_",x,sep="")]],{
                N <- input[[paste("countTable_",x,"_rows_all",
                    sep="")]]
                sel <- input[[paste("countTable_",x,"_rows_selected",
                    sep="")]]
                if (length(sel)>0) {
                    N <- N[-sel]
                    proxy <- dataTableProxy(paste("countTable_",x,sep=""))
                    selectRows(proxy,N)
                }
            })
        })
    })
    
    handleExpressionDownloadSelection <- reactive({
        c <- currentMetadata$class
        lapply(c,function(x) {
            output[[paste("exportCountSelection_",x,sep="")]] <- 
                downloadHandler(
                    filename=function() {
                        tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                        paste(x,"_gene_expression_",tt,".txt", sep='')
                    },
                    content=function(con) {
                        sel <- input[[paste("countTable_",x,
                            "_rows_selected",sep="")]]
                        if (length(sel)>0)
                            write.table(currentTables[[x]][sel,,drop=FALSE],
                                file=con,sep="\t",quote=FALSE,col.names=NA)
                    }
                )
            output[[paste("exportCountAll_",x,sep="")]] <- downloadHandler(
                filename=function() {
                    tt <- format(Sys.time(),format="%Y%m%d%H%M%S")
                    paste(x,"_gene_expression_",tt,".txt", sep='')
                },
                content=function(con) {
                    write.table(currentTables[[x]],file=con,sep="\t",
                            quote=FALSE,col.names=NA)
                }
            )
        })
    })
    
    return(list(
        createCurrentTables=createCurrentTables,
        #countDataTables=countDataTables,
        handleExpressionSampleSelection=handleExpressionSampleSelection,
        handleExpressionDownloadSelection=handleExpressionDownloadSelection
    ))
}

expressionExplorerTabPanelRenderUI <- function(output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata
    currentTables <- allReactiveVars$currentTables
    
    output$rnaExpressionTables <- renderUI({
        M <- currentMetadata$final
        if (is.null(M))
            div(
                style="display:inline-block; margin:5px;",
                h4(paste("Please create a dataset first using the",
                    "'Data selector' tab."))
            )
        else {
            s <- unique(as.character(M$source))
            d <- unique(as.character(M$dataset))
            c <- unique(as.character(M$class))
            #countDataTables()
            dbGene <- loadedGenomes[[currentMetadata$genome]]$dbGene
            lapply(c,function(x,M,D) {
                samples <- as.character(M[which(as.character(M$class)==x),
                    "sample_id"])
                switch(input$rnaExpressionGeneList,
                    all = {
                        switch(input$rnaExpressionMeasureRadio,
                            raw = {
                                tab <- D$counts[,samples,drop=FALSE]
                                average <- round(apply(tab,1,
                                    input$rnaExpressionAverageRadio))
                                deviation <- round(apply(tab,1,
                                    input$rnaExpressionDeviationRadio))
                            },
                            norm = {
                                tab <- D$norm[,samples,drop=FALSE]
                                average <- round(apply(tab,1,
                                    input$rnaExpressionAverageRadio))
                                deviation <- round(apply(tab,1,
                                    input$rnaExpressionDeviationRadio))
                            },
                            rpkm = {
                                tab <- round(edgeR::rpkm(
                                    D$counts[,samples,drop=FALSE],
                                    gene.length=D$length,
                                    lib.size=unlist(D$libsize[samples]),
                                    ),digits=6)
                                average <- round(apply(tab,1,
                                    input$rnaExpressionAverageRadio),digits=6)
                                deviation <- round(apply(tab,1,
                                    input$rnaExpressionDeviationRadio),digits=6)
                            },
                            rpgm = {
                                tab <- round(D$norm[,samples,
                                    drop=FALSE]/D$length,digits=6)
                                average <- round(apply(tab,1,
                                    input$rnaExpressionAverageRadio),digits=6)
                                deviation <- round(apply(tab,1,
                                    input$rnaExpressionDeviationRadio),digits=6)
                            }
                        )
                    },
                    select = {
                        g <- input$selectExpressionGeneName
                        if (isEmpty(g))
                            tab <- data.frame(
                                name=character(0),
                                average=numeric(0),
                                deviation=numeric(0),
                                sample=character(0)
                            )
                        else {
                            switch(input$rnaExpressionMeasureRadio,
                                raw = {
                                    tab <- D$counts[g,samples,drop=FALSE]
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio))
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio))
                                },
                                norm = {
                                    tab <- D$norm[g,samples,drop=FALSE]
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio))
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio))
                                },
                                rpkm = {
                                    tab <- round(edgeR::rpkm(
                                        D$counts[g,samples,drop=FALSE],
                                        gene.length=D$length[g],
                                        lib.size=unlist(D$libsize[samples])
                                        ),digits=6)
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio),
                                            digits=6)
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio),
                                            digits=6)
                                },
                                rpgm = {
                                    tab <- round(D$norm[g,
                                        samples,drop=FALSE]/D$length[g],
                                            digits=6)
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio),
                                            digits=6)
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio),
                                            digits=6)
                                }
                            )
                        }
                    },
                    custom = {
                        g <- input$rnaExpressionCustomList
                        if (isEmpty(g))
                            tab <- data.frame(
                                name=character(0),
                                average=numeric(0),
                                deviation=numeric(0),
                                sample=character(0)
                            )
                        else {
                            g <- strsplit(g,split="\n")[[1]]
                            m <- match(g,as.character(dbGene$gene_name))
                            na <- which(is.na(m))
                            if (length(na)>0)
                                m <- m[-na]
                            gene_id <- names(dbGene)[m]
                            switch(input$rnaExpressionMeasureRadio,
                                raw = {
                                    tab <- D$counts[gene_id,samples,drop=FALSE]
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio))
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio))
                                },
                                norm = {
                                    tab <- D$norm[gene_id,samples,drop=FALSE]
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio))
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio))
                                },
                                rpkm = {
                                    tab <- round(edgeR::rpkm(
                                        D$counts[gene_id,samples,drop=FALSE],
                                        gene.length=D$length[gene_id],
                                        lib.size=unlist(D$libsize[samples])
                                        ),digits=6)
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio),digits=6)
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio),digits=6)
                                },
                                rpgm = {
                                    tab <- round(D$norm[gene_id,
                                        samples,drop=FALSE]/D$length[gene_id],
                                        digits=6)
                                    average <- round(apply(tab,1,
                                        input$rnaExpressionAverageRadio),digits=6)
                                    deviation <- round(apply(tab,1,
                                        input$rnaExpressionDeviationRadio),digits=6)
                                }
                            )
                        }
                    }
                )
                if (input$rnaExpressionGeneList=="all" || !isEmpty(g)) {
                    switch(input$rnaExpressionScaleRadio,
                        natural = {
                            tab <- tab
                            average <- average
                            deviation <- deviation
                        },
                        log2 = {
                            tab <- round(log2(tab+1),digits=6)
                            average <- round(apply(tab,1,
                                input$rnaExpressionAverageRadio),digits=6)
                            deviation <- round(apply(tab,1,
                                input$rnaExpressionDeviationRadio),digits=6)
                        }
                    )
                    
                    # Use alternative ids if present
                    if (!is.null(M$alt_id))
                        colnames(tab) <- as.character(M[which(
                            as.character(M$class)==x),"alt_id"])
                    # Attach some annotation
                    meta <- data.frame(
                        #chromosome=as.character(seqnames(dbGene[rownames(tab)])),
                        #start=start(dbGene[rownames(tab)]),
                        #end=end(dbGene[rownames(tab)]),
                        name=as.character(dbGene[rownames(tab)]$gene_name),
                        average=average,
                        deviation=deviation
                    )
                    currentTables[[x]] <- cbind(meta,tab)
                }
                else 
                    meta <- NULL
                output[[paste("countTable",x,sep="_")]] <- 
                    DT::renderDataTable(
                        cbind(meta,tab),
                        class="display compact",
                        rownames=FALSE,
                        options=list(
                            searchHighlight=TRUE,
                            pageLength=10,
                            lengthMenu=c(10,20,50,100)
                        )
                    )
            },M,loadedData[[s]][[d]])
            
            lapply(c,function(x,s,d) {
                list(
                    div(
                        class="small table-container",
                        h4(paste("Gene expression for ",x)),
                        DT::dataTableOutput(paste("countTable",x,sep="_")),
                        div(
                            style="display:inline-block; margin:5px;",
                            actionButton(
                                paste("clearCountSelection_",x,sep=""),
                                paste("Clear selection"),
                                class="btn-xs"
                            )
                        ),
                        div(
                            style="display:inline-block; margin:5px;",
                            actionButton(
                                paste("invertCountSelection_",x,sep=""),
                                paste("Invert selection"),
                                class="btn-xs"
                            )
                        ),
                        div(
                            style="display:inline-block; margin:5px;",
                            downloadButton(
                                paste("exportCountSelection_",x,sep=""),
                                paste("Export selection"),
                                class="btn-xs"
                            )
                        ),
                        div(
                            style="display:inline-block; margin:5px;",
                            downloadButton(
                                paste("exportCountAll_",x,sep=""),
                                paste("Export all"),
                                class="btn-xs"
                            )
                        )
                    ),
                    hr()
                )
            },s,d)
        }
    })
}

expressionExplorerTabPanelObserve <- function(input,output,session,
    allReactiveVars,allReactiveMsgs) {
    currentMetadata <- allReactiveVars$currentMetadata  
        
    expressionExplorerTabPanelReactiveEvents <- 
        expressionExplorerTabPanelEventReactive(input,output,session,
            allReactiveVars,allReactiveMsgs)
            
    expressionExplorerTabPanelReactiveExprs <- 
        expressionExplorerTabPanelReactive(input,output,session,allReactiveVars,
            allReactiveMsgs)
            
    createCurrentTables <- 
        expressionExplorerTabPanelReactiveExprs$createCurrentTables
    handleExpressionSampleSelection <- 
        expressionExplorerTabPanelReactiveExprs$handleExpressionSampleSelection
    handleExpressionDownloadSelection <- 
        expressionExplorerTabPanelReactiveExprs$handleExpressionDownloadSelection
    
    expressionExplorerTabPanelRenderUI(output,session,allReactiveVars,
        allReactiveMsgs)
    
    observe({
        if (is.null(currentMetadata$final)) {
            updateSelectizeInput(session,"selectExpressionGeneName",
                choices=NULL,
                server=TRUE
            )
            shinyjs::reset("rnaExpressionCustomList")
        }
        else {
            geneNames <- loadedGenomes[[currentMetadata$genome]]$geneNames
            g <- isolate({input$selectExpressionGeneName})
            i <- grep(paste0("^",g),geneNames,perl=TRUE)
            if (length(i)>0) {
                updateSelectizeInput(session,"selectExpressionGeneName",
                    choices=geneNames[i],
                    #selected=input$selectExpressionGeneName,
                    selected=g,
                    server=TRUE
                )
            }
        }
    })
    
    observe({
        createCurrentTables()
        handleExpressionSampleSelection()
        handleExpressionDownloadSelection()
    })
}
